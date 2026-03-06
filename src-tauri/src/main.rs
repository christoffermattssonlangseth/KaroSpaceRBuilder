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
    rscript_path: Option<String>,
    input_path: String,
    metadata_input_path: Option<String>,
    assay: Option<String>,
    gene_query: Option<String>,
    gene_limit: Option<u32>,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
struct BuildRequest {
    rscript_path: Option<String>,
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
    runtime_mode: String,
    bundled_rscript_path: Option<String>,
    bundled_r_home: Option<String>,
    bundled_r_version: Option<String>,
    backend_version: Option<String>,
    backend_revision: Option<String>,
}

#[derive(Debug)]
struct CommandOutput {
    status: i32,
    stdout: String,
    stderr: String,
}

#[derive(Debug)]
struct ResolvedRuntime {
    rscript_path: PathBuf,
    r_home: Option<PathBuf>,
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
    runtime: &ResolvedRuntime,
    working_dir: &Path,
    script: &Path,
    args: I,
) -> Result<CommandOutput, String>
where
    I: IntoIterator<Item = S>,
    S: AsRef<OsStr>,
{
    let mut command = Command::new(&runtime.rscript_path);
    command.current_dir(working_dir).arg(script).args(args);

    if let Some(r_home) = runtime.r_home.as_ref() {
        let library_path = r_home.join("library");
        command
            .env("R_HOME", r_home)
            .env("R_LIBS", &library_path)
            .env("R_LIBS_SITE", &library_path)
            .env("R_LIBS_USER", &library_path)
            .env("R_PROFILE_USER", "")
            .env("R_ENVIRON_USER", "");
    }

    let output = command
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

fn bundled_r_dev_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("bundle-resources")
        .join("macos-arm64")
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

fn bundled_r_framework_root(root: &Path) -> PathBuf {
    root.join("R.framework")
}

fn bundled_r_home(root: &Path) -> Option<PathBuf> {
    let framework_root = bundled_r_framework_root(root);
    let direct_resources = framework_root.join("Resources");
    if direct_resources.join("bin").join("exec").join("R").exists() {
        return Some(direct_resources);
    }

    let current_resources = framework_root
        .join("Versions")
        .join("Current")
        .join("Resources");
    if current_resources
        .join("bin")
        .join("exec")
        .join("R")
        .exists()
    {
        return Some(current_resources);
    }

    let versions_dir = framework_root.join("Versions");
    let mut version_dirs = fs::read_dir(&versions_dir)
        .ok()?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| path.is_dir())
        .collect::<Vec<_>>();
    version_dirs.sort();

    version_dirs.into_iter().find_map(|version_dir| {
        let resources = version_dir.join("Resources");
        resources
            .join("bin")
            .join("exec")
            .join("R")
            .exists()
            .then_some(resources)
    })
}

fn bundled_rscript_path(root: &Path) -> PathBuf {
    root.join("bin").join("karospacer-rscript")
}

fn bundled_r_version(root: &Path) -> Option<String> {
    fs::read_link(
        bundled_r_framework_root(root)
            .join("Versions")
            .join("Current"),
    )
    .ok()
    .and_then(|path| {
        path.file_name()
            .map(|value| value.to_string_lossy().to_string())
    })
    .or_else(|| {
        let versions_dir = bundled_r_framework_root(root).join("Versions");
        let mut versions = fs::read_dir(versions_dir)
            .ok()?
            .filter_map(Result::ok)
            .map(|entry| entry.path())
            .filter(|path| path.is_dir())
            .filter_map(|path| {
                path.file_name()
                    .map(|value| value.to_string_lossy().to_string())
            })
            .collect::<Vec<_>>();
        versions.sort();
        versions.into_iter().next()
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

fn resolve_bundled_r_root(app: &tauri::AppHandle) -> Result<PathBuf, String> {
    let dev_candidate = bundled_r_dev_root();
    if bundled_rscript_path(&dev_candidate).exists() {
        return Ok(dev_candidate);
    }

    let resource_candidate = app
        .path()
        .resource_dir()
        .map_err(|err| format!("Could not resolve the app resource directory: {err}"))?
        .join("runtime")
        .join("macos-arm64");

    if bundled_rscript_path(&resource_candidate).exists() {
        return Ok(resource_candidate);
    }

    Err(format!(
        "Could not resolve a bundled Apple Silicon R runtime. Looked in:\n- {}\n- {}",
        dev_candidate.to_string_lossy(),
        resource_candidate.to_string_lossy()
    ))
}

fn resolve_runtime(
    app: &tauri::AppHandle,
    requested_rscript_path: Option<&str>,
) -> Result<ResolvedRuntime, String> {
    if let Ok(root) = resolve_bundled_r_root(app) {
        let r_home = bundled_r_home(&root).ok_or_else(|| {
            format!(
                "Found a bundled Apple Silicon R runtime, but could not resolve its R home inside {}",
                root.to_string_lossy()
            )
        })?;
        return Ok(ResolvedRuntime {
            rscript_path: bundled_rscript_path(&root),
            r_home: Some(r_home),
        });
    }

    let requested_path = requested_rscript_path
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(str::to_string)
        .or_else(suggest_rscript_path);

    let rscript_path = requested_path.ok_or_else(|| {
        "No bundled Apple Silicon R runtime was found, and no external Rscript path was provided."
            .to_string()
    })?;
    ensure_rscript_path(&rscript_path)?;

    Ok(ResolvedRuntime {
        rscript_path: PathBuf::from(rscript_path),
        r_home: None,
    })
}

#[tauri::command]
fn guess_backend_paths(app: tauri::AppHandle) -> SuggestedPaths {
    let backend_root = resolve_bundled_backend_root(&app).ok();
    let bundled_runtime = resolve_bundled_r_root(&app).ok();
    SuggestedPaths {
        rscript_path: suggest_rscript_path(),
        runtime_mode: if bundled_runtime.is_some() {
            "bundled".to_string()
        } else {
            "external".to_string()
        },
        bundled_rscript_path: bundled_runtime
            .as_deref()
            .map(bundled_rscript_path)
            .map(canonicalize_lossy),
        bundled_r_home: bundled_runtime
            .as_deref()
            .and_then(bundled_r_home)
            .map(canonicalize_lossy),
        bundled_r_version: bundled_runtime.as_deref().and_then(bundled_r_version),
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
    let runtime = resolve_runtime(&app, request.rscript_path.as_deref())?;
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

    let output = run_rscript(&runtime, &backend_root, &script, args)?;
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
    let runtime = resolve_runtime(&app, request.rscript_path.as_deref())?;
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
        &runtime,
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
