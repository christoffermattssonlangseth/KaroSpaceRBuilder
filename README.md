# KaroSpaceRBuilder

`KaroSpaceRBuilder` is a separate desktop wrapper for `KaroSpaceR`.

It does not contain export logic. It drives `KaroSpaceR` through two backend contracts:

- `scripts/karospace_inspect_r.R` for machine-readable dataset inspection
- `scripts/karospace_build_r.R --config build.json` for export

## MVP assumptions

- `Rscript` is installed on the machine
- the local `KaroSpaceR` repo exists and contains the current scripts
- this app shells out to those scripts instead of bundling an R runtime

## Current flow

1. choose `Rscript`
2. choose the local `KaroSpaceR` repo
3. choose an `.rds` input
4. run inspect
5. edit the generated build config in the form
6. build the viewer
7. open the exported HTML

## Dev commands

```bash
npm install
npm run dev
```

## Notes

- This scaffold uses a plain static frontend and Tauri's built-in dev server.
- It currently prefers explicit file paths and native file pickers over hidden assumptions about shell `$PATH`.
