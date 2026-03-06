# KaroSpaceRBuilder

`KaroSpaceRBuilder` is a separate desktop wrapper for `KaroSpaceR`.

It does not contain export logic. It drives a pinned bundled snapshot of `KaroSpaceR`
through two backend contracts:

- `scripts/karospace_inspect_r.R` for machine-readable dataset inspection
- `scripts/karospace_build_r.R --config build.json` for export

## MVP assumptions

- `Rscript` is installed on the machine
- the app bundles the required `KaroSpaceR` backend files
- this app shells out to those scripts instead of bundling an R runtime

## Current flow

1. choose `Rscript`
2. choose an `.rds` input
3. run inspect
4. edit the generated build config in the form
5. build the viewer
6. open the exported HTML

## Dev commands

```bash
npm install
npm run dev
npm run build
```

## Notes

- This scaffold uses a plain static frontend and Tauri's built-in dev server.
- It prefers explicit `Rscript` selection over hidden assumptions about shell `$PATH`.
- The bundled backend snapshot lives in `vendor/KaroSpaceR`.
- The packaged app bundles `KaroSpaceR`, but it does not bundle an R runtime.

## GitHub Releases

macOS release builds are published from GitHub Actions when you push a `v*` tag,
for example:

```bash
git tag v0.1.1
git push origin v0.1.1
```

The workflow in `.github/workflows/desktop-release.yml` creates or updates the
GitHub Release for that tag and uploads signed/notarized macOS DMGs for:

- Apple Silicon: `aarch64-apple-darwin`
- Intel Mac: `x86_64-apple-darwin`

Downloaded releases include the bundled `KaroSpaceR` backend, but end users
still need:

- `Rscript`
- the required R packages for their datasets

## GitHub Signing Secrets

To enable signed and notarized macOS releases, configure these repository
secrets:

- `APPLE_CERTIFICATE`: Base64-encoded exported `Developer ID Application` certificate (`.p12`)
- `APPLE_CERTIFICATE_PASSWORD`: Password used when exporting the `.p12`
- `APPLE_ID`: Apple ID email used for notarization
- `APPLE_PASSWORD`: Apple app-specific password used for notarization
- `APPLE_TEAM_ID`: Apple Developer Team ID
- `KEYCHAIN_PASSWORD`: Temporary CI keychain password used while importing the certificate

The workflow imports the certificate into a temporary macOS keychain, derives
the `Developer ID Application` identity from that keychain, and passes the
resulting signing identity plus the notarization credentials to Tauri.
