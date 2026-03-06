# KaroSpaceRBuilder

`KaroSpaceRBuilder` is a separate desktop wrapper for `KaroSpaceR`.

It does not contain export logic. It drives a pinned bundled snapshot of `KaroSpaceR`
through two backend contracts:

- `scripts/karospace_inspect_r.R` for machine-readable dataset inspection
- `scripts/karospace_build_r.R --config build.json` for export

## Runtime modes

- The app always bundles the required `KaroSpaceR` backend files from `vendor/KaroSpaceR`.
- On Apple Silicon, you can stage a bundled R runtime into `bundle-resources/macos-arm64`
  and package a self-contained macOS app.
- If no bundled runtime is present, the app falls back to a user-provided `Rscript`.

## Current flow

1. use the bundled Apple Silicon R runtime when present, otherwise choose `Rscript`
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

To stage a local Apple Silicon R runtime into the app resources before packaging:

```bash
npm run stage:r
npm run build -- --target aarch64-apple-darwin
```

## Notes

- This scaffold uses a plain static frontend and Tauri's built-in dev server.
- It prefers a bundled Apple Silicon runtime when one is present.
- Fallback `Rscript` selection is still supported for dev mode and non-bundled builds.
- The bundled backend snapshot lives in `vendor/KaroSpaceR`.
- The staged Apple Silicon runtime lives in `bundle-resources/macos-arm64` and is intentionally gitignored.
- The runtime staging script rewrites library paths so the packaged `.app` can run without `/Library/Frameworks/R.framework`.

## Apple Silicon bundled-R prototype

The local Apple Silicon packaging path is now working end to end:

- stage `R.framework` plus the current package library with `npm run stage:r`
- build the Tauri app/DMG for `aarch64-apple-darwin`
- run inspect/build through the packaged app's embedded `karospacer-rscript`

This has been smoke-tested locally against:

- `/Users/chrislangseth/Downloads/28124747/IntegratedHeartST.rds`

and successfully produced a packaged-app export HTML.

## GitHub Releases

The current reliable release path for the bundled-R Apple Silicon app is manual
from this Mac, using the DMG built after `npm run stage:r`.

Example:

```bash
gh release create v0.2.0 \
  src-tauri/target/aarch64-apple-darwin/release/bundle/dmg/KaroSpaceRBuilder_0.2.0_aarch64.dmg \
  --title "KaroSpaceRBuilder v0.2.0" \
  --notes "Apple Silicon DMG with bundled KaroSpaceR backend and bundled R runtime."
```

The GitHub Actions workflow in `.github/workflows/desktop-release.yml` is kept
for future CI validation, but automatic tagged releases are intentionally paused
until CI-side runtime staging exists. That avoids publishing a fallback build
that silently drops the bundled R runtime.

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
