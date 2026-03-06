use rfd::FileDialog;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};
use tauri::Manager;

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
struct InspectRequest {
    rscript_path: String,
    input_path: String,
    metadata_input_path: Option<String>,
    assay: Option<String>,
    gene_query: Option<String>,
    gene_limit: Option<u32>,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
struct BuildRequest {
    rscript_path: String,
    config: Value,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
struct PickPathRequest {
    mode: String,
    title: Option<String>,
    default_name: Option<String>,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct InspectResponse {
    summary: Value,
    stdout: String,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct BuildResponse {
    config_path: String,
    output_path: Option<String>,
    stdout: String,
    stderr: String,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct SuggestedPaths {
    rscript_path: Option<String>,
    backend_version: Option<String>,
    backend_revision: Option<String>,
}

#[derive(Debug)]
struct CommandOutput {
    status: i32,
    stdout: String,
    stderr: String,
}

fn canonicalize_lossy(path: impl AsRef<Path>) -> String {
    let candidate = path.as_ref();
    candidate
        .canonicalize()
        .unwrap_or_else(|_| candidate.to_path_buf())
        .to_string_lossy()
        .to_string()
}

fn script_path(root: &Path, script_name: &str) -> Result<PathBuf, String> {
    let path = root.join("scripts").join(script_name);
    if !path.exists() {
        return Err(format!(
            "Could not find {} inside bundled KaroSpaceR backend: {}",
            script_name,
            path.to_string_lossy()
        ));
    }
    Ok(path)
}

fn run_rscript<I, S>(
    rscript_path: &str,
    working_dir: &Path,
    script: &Path,
    args: I,
) -> Result<CommandOutput, String>
where
    I: IntoIterator<Item = S>,
    S: AsRef<OsStr>,
{
    let output = Command::new(rscript_path)
        .current_dir(working_dir)
        .arg(script)
        .args(args)
        .output()
        .map_err(|err| format!("Failed to run Rscript: {err}"))?;

    Ok(CommandOutput {
        status: output.status.code().unwrap_or(-1),
        stdout: String::from_utf8_lossy(&output.stdout).to_string(),
        stderr: String::from_utf8_lossy(&output.stderr).to_string(),
    })
}

fn ensure_rscript_path(rscript_path: &str) -> Result<(), String> {
    if !Path::new(rscript_path).exists() {
        return Err(format!("Rscript path does not exist: {rscript_path}"));
    }
    Ok(())
}

fn suggest_rscript_path() -> Option<String> {
    let common_paths = [
        "/opt/homebrew/bin/Rscript",
        "/usr/local/bin/Rscript",
        "/Library/Frameworks/R.framework/Resources/bin/Rscript",
    ];

    for candidate in common_paths {
        if Path::new(candidate).exists() {
            return Some(candidate.to_string());
        }
    }

    std::env::var_os("PATH")
        .and_then(|raw| {
            std::env::split_paths(&raw)
                .map(|dir| dir.join("Rscript"))
                .find(|path| path.exists())
        })
        .map(canonicalize_lossy)
}

fn bundled_backend_dev_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("vendor")
        .join("KaroSpaceR")
}

fn description_field(root: &Path, field_name: &str) -> Option<String> {
    let description = fs::read_to_string(root.join("DESCRIPTION")).ok()?;
    let prefix = format!("{field_name}:");
    description.lines().find_map(|line| {
        line.strip_prefix(&prefix)
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
    })
}

fn vendored_revision(root: &Path) -> Option<String> {
    let metadata = fs::read_to_string(root.join("VENDORED_FROM.txt")).ok()?;
    metadata.lines().find_map(|line| {
        line.strip_prefix("source_commit_short=")
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
    })
}

fn resolve_bundled_backend_root(app: &tauri::AppHandle) -> Result<PathBuf, String> {
    let dev_candidate = bundled_backend_dev_root();
    if dev_candidate.exists() {
        return Ok(dev_candidate);
    }

    let resource_candidate = app
        .path()
        .resource_dir()
        .map_err(|err| format!("Could not resolve the app resource directory: {err}"))?
        .join("backend")
        .join("KaroSpaceR");

    if resource_candidate.exists() {
        return Ok(resource_candidate);
    }

    Err(format!(
        "Could not resolve the bundled KaroSpaceR backend. Looked in:\n- {}\n- {}",
        dev_candidate.to_string_lossy(),
        resource_candidate.to_string_lossy()
    ))
}

#[tauri::command]
fn guess_backend_paths(app: tauri::AppHandle) -> SuggestedPaths {
    let backend_root = resolve_bundled_backend_root(&app).ok();
    SuggestedPaths {
        rscript_path: suggest_rscript_path(),
        backend_version: backend_root
            .as_deref()
            .and_then(|root| description_field(root, "Version")),
        backend_revision: backend_root.as_deref().and_then(vendored_revision),
    }
}

#[tauri::command]
fn pick_path(request: PickPathRequest) -> Result<Option<String>, String> {
    let mut dialog = FileDialog::new();
    if let Some(title) = request.title {
        dialog = dialog.set_title(&title);
    }
    let result = match request.mode.as_str() {
        "file" => dialog.pick_file(),
        "folder" => dialog.pick_folder(),
        "save" => {
            if let Some(default_name) = request.default_name {
                dialog = dialog.set_file_name(&default_name);
            }
            dialog.save_file()
        }
        other => {
            return Err(format!("Unsupported pick mode: {other}"));
        }
    };
    Ok(result.map(canonicalize_lossy))
}

#[tauri::command]
fn inspect_dataset(
    app: tauri::AppHandle,
    request: InspectRequest,
) -> Result<InspectResponse, String> {
    ensure_rscript_path(&request.rscript_path)?;
    let backend_root = resolve_bundled_backend_root(&app)?;
    let script = script_path(&backend_root, "karospace_inspect_r.R")?;

    let mut args = vec!["--input".to_string(), request.input_path.clone()];
    if let Some(metadata_input_path) = request.metadata_input_path.as_ref() {
        if !metadata_input_path.trim().is_empty() {
            args.push("--metadata-input".to_string());
            args.push(metadata_input_path.clone());
        }
    }
    if let Some(assay) = request.assay.as_ref() {
        if !assay.trim().is_empty() {
            args.push("--assay".to_string());
            args.push(assay.clone());
        }
    }
    if let Some(gene_query) = request.gene_query.as_ref() {
        if !gene_query.trim().is_empty() {
            args.push("--gene-query".to_string());
            args.push(gene_query.clone());
        }
    }
    if let Some(gene_limit) = request.gene_limit {
        args.push("--gene-limit".to_string());
        args.push(gene_limit.to_string());
    }

    let output = run_rscript(&request.rscript_path, &backend_root, &script, args)?;
    if output.status != 0 {
        return Err(format!(
            "Inspect failed with status {}.\nstdout:\n{}\nstderr:\n{}",
            output.status, output.stdout, output.stderr
        ));
    }

    let summary: Value = serde_json::from_str(&output.stdout).map_err(|err| {
        format!(
            "Inspect command returned invalid JSON: {err}\nstdout:\n{}",
            output.stdout
        )
    })?;

    Ok(InspectResponse {
        summary,
        stdout: output.stdout,
    })
}

#[tauri::command]
fn build_viewer(app: tauri::AppHandle, request: BuildRequest) -> Result<BuildResponse, String> {
    ensure_rscript_path(&request.rscript_path)?;
    let backend_root = resolve_bundled_backend_root(&app)?;
    let script = script_path(&backend_root, "karospace_build_r.R")?;

    let output_path = request
        .config
        .get("output")
        .and_then(|value| value.as_str())
        .map(str::to_string);

    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|err| format!("System clock error: {err}"))?
        .as_millis();
    let config_path = std::env::temp_dir().join(format!("karospacer-builder-{timestamp}.json"));
    let config_json = serde_json::to_vec_pretty(&request.config)
        .map_err(|err| format!("Failed to serialize config JSON: {err}"))?;
    fs::write(&config_path, config_json)
        .map_err(|err| format!("Failed to write temporary config file: {err}"))?;

    let output = run_rscript(
        &request.rscript_path,
        &backend_root,
        &script,
        [
            "--config".to_string(),
            config_path.to_string_lossy().to_string(),
        ],
    )?;
    if output.status != 0 {
        return Err(format!(
            "Build failed with status {}.\nstdout:\n{}\nstderr:\n{}",
            output.status, output.stdout, output.stderr
        ));
    }

    Ok(BuildResponse {
        config_path: canonicalize_lossy(&config_path),
        output_path,
        stdout: output.stdout,
        stderr: output.stderr,
    })
}

#[tauri::command]
fn open_path(path: String) -> Result<(), String> {
    let status = if cfg!(target_os = "macos") {
        Command::new("open").arg(&path).status()
    } else if cfg!(target_os = "windows") {
        Command::new("cmd")
            .args(["/C", "start", "", &path])
            .status()
    } else {
        Command::new("xdg-open").arg(&path).status()
    }
    .map_err(|err| format!("Failed to launch opener for {path}: {err}"))?;

    if status.success() {
        Ok(())
    } else {
        Err(format!("Open command exited with status: {status}"))
    }
}

fn main() {
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![
            guess_backend_paths,
            pick_path,
            inspect_dataset,
            build_viewer,
            open_path
        ])
        .run(tauri::generate_context!())
        .expect("error while running KaroSpaceRBuilder");
}
