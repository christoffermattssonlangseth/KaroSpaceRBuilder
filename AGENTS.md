# AGENTS.md

## Purpose

This repository is the desktop wrapper for `KaroSpaceR`.

It must remain a thin UI/process layer. It should help a user choose files, inspect an input dataset, edit export settings, run the backend export, and open the resulting HTML. It must not reimplement payload construction, object normalization, viewer packaging, or analytics logic.

## Hard Boundaries

- Do not move export logic from `KaroSpaceR` into this repo.
- Do not add dataset normalization, Seurat parsing, or payload-building logic here.
- Do not modify the Python `KaroSpace` repository from this repo.
- Do not assume shell startup files will populate `$PATH` inside a desktop app.
- Treat `Rscript` as an explicit dependency.
- Treat the bundled `vendor/KaroSpaceR` snapshot as the runtime backend, while keeping the sibling `../KaroSpaceR` repo as the canonical upstream source.

## Backend Contract

This repo depends on the following stable contracts from `KaroSpaceR`.
At runtime it uses the vendored snapshot in `vendor/KaroSpaceR`.
When the contract changes, update `KaroSpaceR` first and then refresh the vendored snapshot here.

- Inspect:
  - `scripts/karospace_inspect_r.R`
  - machine-readable JSON output
- Build:
  - `scripts/karospace_build_r.R --config build.json`
  - config-driven standalone HTML export

The intended flow is:

1. choose `Rscript`
2. choose an `.rds` input
3. run inspect
4. load `default_config` from inspect JSON
5. let the user edit config values
6. write config JSON
7. run build
8. open the exported HTML

If the backend contract needs to change, prefer changing `KaroSpaceR` first and then adapting this repo. Do not invent a parallel builder-only contract.

## Current Scope

Current scaffolded scope:

- Tauri desktop shell
- Static frontend in `src/`
- Native file/folder pickers
- Bundled `KaroSpaceR` backend snapshot in `vendor/KaroSpaceR`
- GitHub Actions macOS release packaging/signing flow
- Rust commands for:
  - guessing runtime paths
  - choosing files/folders
  - running inspect
  - running build
  - opening the output HTML
- JSON config preview/edit flow

Not implemented yet:

- persistent settings storage
- progress streaming from long-running R jobs
- richer validation and better field-level error rendering
- bundled R runtime

## Repo Structure

- `src/`
  Desktop UI.
- `src/app.js`
  Frontend state, inspect flow, config editing, build flow.
- `src/index.html`
  Main layout.
- `src/styles.css`
  Visual design.
- `src-tauri/src/main.rs`
  Tauri backend commands and process launching.
- `src-tauri/tauri.conf.json`
  Tauri app configuration.
- `src-tauri/capabilities/default.json`
  Tauri window capability declaration.
- `.github/workflows/desktop-release.yml`
  Tagged macOS release packaging/signing pipeline.

## Working Rules For Agents

- Keep this repo a wrapper over `KaroSpaceR`, not a second backend.
- Prefer explicit file paths and machine-readable contracts over heuristics.
- If a UI field maps to a backend config field, keep the naming close to the backend config.
- If adding new export options, first verify that `KaroSpaceR` exposes them through inspect/build.
- Use the inspect JSON `default_config` as the source of truth for initial UI state.
- Keep the builder resilient to missing tools and bad paths with clear error messages.
- Keep the vendored backend snapshot traceable to an upstream `KaroSpaceR` revision.
- Prefer calling existing Rust commands over adding frontend-only assumptions.

## Tauri Notes

- This repo is using Tauri v2 concepts and config structure.
- GUI apps often do not inherit shell environment setup, so do not assume `Rscript` is discoverable automatically.
- Prefer letting the user browse to `Rscript` explicitly.
- Resolve bundled backend files from app resources in packaged builds instead of assuming a sibling checkout exists.
- Avoid broad command execution surfaces; keep Rust commands narrow and task-specific.
- If adding OS integrations, keep macOS first and avoid unnecessary platform-specific branching until needed.

## Testing

Before finishing builder changes, run what is relevant:

```bash
node --check src/app.js
python -m json.tool package.json >/dev/null
python -m json.tool src-tauri/tauri.conf.json >/dev/null
python -m json.tool src-tauri/capabilities/default.json >/dev/null
cargo fmt --check --manifest-path src-tauri/Cargo.toml
```

If Tauri dependencies are installed, also run:

```bash
npm install
npm run dev
```

Before blaming the builder, verify the backend directly against the vendored snapshot or the sibling upstream repo:

```bash
Rscript vendor/KaroSpaceR/scripts/karospace_inspect_r.R --input /path/to/object.rds
Rscript vendor/KaroSpaceR/scripts/karospace_build_r.R --config /path/to/build.json
```

## Near-Term Priorities

Priority order for future work:

1. add better inspect/build error presentation in the UI
2. persist last-used `Rscript` path
3. add progress/log streaming for long exports
4. improve config editing UX for colors, genes, and advanced options
5. validate the signed macOS GitHub release flow with real Apple secrets
6. evaluate whether bundling an R runtime is worth the complexity

## Output Expectations

When making changes, optimize for:

- correct use of the `KaroSpaceR` backend contracts
- no duplication of exporter logic
- clear, recoverable desktop errors
- explicit paths and reproducible config output
- a UI that is practical for real datasets, not just toy demos
