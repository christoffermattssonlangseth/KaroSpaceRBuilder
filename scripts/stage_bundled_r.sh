#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FRAMEWORK_SRC="${R_FRAMEWORK_SRC:-/Library/Frameworks/R.framework}"
DEST_ROOT="${BUNDLED_R_DEST:-$ROOT_DIR/bundle-resources/macos-arm64}"
DEST_FRAMEWORK="$DEST_ROOT/R.framework"
WRAPPER_DIR="$DEST_ROOT/bin"
WRAPPER_PATH="$WRAPPER_DIR/karospacer-rscript"
INCLUDE_SEURAT="${INCLUDE_SEURAT:-1}"
EXTRA_PACKAGES="${EXTRA_PACKAGES:-}"
BUNDLE_FULL_LIBRARY="${BUNDLE_FULL_LIBRARY:-1}"

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "This staging script only supports macOS." >&2
  exit 1
fi

if [[ "$(uname -m)" != "arm64" ]]; then
  echo "This staging script is currently Apple Silicon only." >&2
  exit 1
fi

if [[ ! -d "$FRAMEWORK_SRC" ]]; then
  echo "R.framework not found at $FRAMEWORK_SRC" >&2
  exit 1
fi

CURRENT_VERSION="$(readlink "$FRAMEWORK_SRC/Versions/Current")"
if [[ -z "$CURRENT_VERSION" ]]; then
  echo "Could not resolve $FRAMEWORK_SRC/Versions/Current" >&2
  exit 1
fi

SOURCE_VERSION_ROOT="$FRAMEWORK_SRC/Versions/$CURRENT_VERSION"
SOURCE_RESOURCES="$SOURCE_VERSION_ROOT/Resources"
DEST_VERSION_ROOT="$DEST_FRAMEWORK/Versions/$CURRENT_VERSION"
DEST_RESOURCES="$DEST_VERSION_ROOT/Resources"
DEST_LIBRARY="$DEST_RESOURCES/library"
METADATA_PATH="$DEST_ROOT/R-BUNDLE-METADATA.txt"

mkdir -p "$DEST_ROOT"
rm -rf "$DEST_FRAMEWORK"
rm -rf "$WRAPPER_DIR"
mkdir -p "$DEST_FRAMEWORK/Versions"
mkdir -p "$DEST_RESOURCES"
mkdir -p "$WRAPPER_DIR"

ln -s "$CURRENT_VERSION" "$DEST_FRAMEWORK/Versions/Current"
ln -s "Versions/Current/Resources" "$DEST_FRAMEWORK/Resources"
ln -s "Versions/Current/R" "$DEST_FRAMEWORK/R"
ln -s "Versions/Current/Resources/lib" "$DEST_FRAMEWORK/Libraries"

rsync -a \
  --exclude 'Resources/library/*' \
  --exclude 'Resources/doc/***' \
  --exclude 'Resources/tests/***' \
  --exclude 'Resources/include/***' \
  --exclude 'Resources/lib/*.dSYM/***' \
  --exclude 'Resources/lib/pkgconfig/***' \
  "$SOURCE_VERSION_ROOT/" \
  "$DEST_VERSION_ROOT/"

mkdir -p "$DEST_LIBRARY"

PACKAGE_LIST=()
if [[ "$BUNDLE_FULL_LIBRARY" == "1" ]]; then
  rsync -a "$SOURCE_RESOURCES/library/" "$DEST_LIBRARY/"
  while IFS= read -r pkg; do
    [[ -n "$pkg" ]] && PACKAGE_LIST+=("$pkg")
  done < <(find "$DEST_LIBRARY" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort)
else
  while IFS= read -r pkg; do
    [[ -n "$pkg" ]] && PACKAGE_LIST+=("$pkg")
  done < <(
    INCLUDE_SEURAT="$INCLUDE_SEURAT" EXTRA_PACKAGES="$EXTRA_PACKAGES" Rscript --vanilla -e '
      seed <- c("jsonlite", "Matrix", "SeuratObject", "SummarizedExperiment", "S4Vectors", "SingleCellExperiment")
      if ("SpatialExperiment" %in% rownames(installed.packages())) {
        seed <- c(seed, "SpatialExperiment")
      }
      if (Sys.getenv("INCLUDE_SEURAT") == "1" && "Seurat" %in% rownames(installed.packages())) {
        seed <- c(seed, "Seurat")
      }
      extra <- trimws(unlist(strsplit(Sys.getenv("EXTRA_PACKAGES"), ",")))
      extra <- extra[nzchar(extra)]
      seed <- unique(c(seed, extra))
      ip <- installed.packages()
      db <- ip[, c("Package", "Depends", "Imports", "LinkingTo", "Suggests")]
      deps <- tools::package_dependencies(seed, recursive = TRUE, db = db, which = c("Depends", "Imports", "LinkingTo"))
      keep <- sort(unique(c(seed, unlist(deps, use.names = FALSE))))
      priority <- ip[, "Priority"]
      keep <- sort(unique(c(keep, rownames(ip)[!is.na(priority) & priority %in% c("base", "recommended")])))
      keep <- keep[keep %in% rownames(ip)]
      cat(paste(keep, collapse = "\n"))
    '
  )

for pkg in "${PACKAGE_LIST[@]}"; do
    pkg_src="$SOURCE_RESOURCES/library/$pkg"
    if [[ -d "$pkg_src" ]]; then
      rsync -a "$pkg_src" "$DEST_LIBRARY/"
    fi
  done
fi

rewrite_symlink() {
  local link_path="$1"
  local target_path
  local mapped_target
  local link_dir
  local relative_target

  target_path="$(readlink "$link_path")"
  link_dir="$(dirname "$link_path")"

  if [[ "$target_path" == /* ]]; then
    if [[ "$target_path" == "$FRAMEWORK_SRC/"* ]]; then
      mapped_target="$DEST_FRAMEWORK${target_path#$FRAMEWORK_SRC}"
      if [[ -e "$mapped_target" ]]; then
        relative_target="$(python3 -c 'import os, sys; print(os.path.relpath(sys.argv[1], sys.argv[2]))' "$mapped_target" "$link_dir")"
        ln -snf "$relative_target" "$link_path"
      else
        rm "$link_path"
      fi
    fi
    return
  fi

  mapped_target="$(python3 -c 'import os, sys; print(os.path.normpath(os.path.join(sys.argv[1], sys.argv[2])))' "$link_dir" "$target_path")"
  if [[ ! -e "$mapped_target" ]]; then
    rm "$link_path"
  fi
}

SYMLINKS=()
while IFS= read -r symlink; do
  [[ -n "$symlink" ]] && SYMLINKS+=("$symlink")
done < <(find "$DEST_FRAMEWORK" -type l | sort)

for symlink in "${SYMLINKS[@]}"; do
  rewrite_symlink "$symlink"
done

rewrite_binary() {
  local file="$1"
  local lib_dir="$DEST_RESOURCES/lib"
  local file_dir
  local rel_to_lib
  local dep
  local dep_name
  local new_dep

  file_dir="$(dirname "$file")"
  rel_to_lib="$(python3 -c 'import os, sys; print(os.path.relpath(sys.argv[1], sys.argv[2]))' "$lib_dir" "$file_dir")"

  while IFS= read -r dep; do
    case "$dep" in
      "$FRAMEWORK_SRC/Resources/lib/"*|"$SOURCE_RESOURCES/lib/"*)
        dep_name="$(basename "$dep")"
        if [[ "$rel_to_lib" == "." ]]; then
          new_dep="@loader_path/$dep_name"
        else
          new_dep="@loader_path/$rel_to_lib/$dep_name"
        fi
        install_name_tool -change "$dep" "$new_dep" "$file"
        ;;
    esac
  done < <(otool -L "$file" | awk 'NR > 1 { print $1 }')

  if [[ "$file_dir" == "$lib_dir" && "$(basename "$file")" == *.dylib ]]; then
    install_name_tool -id "@loader_path/$(basename "$file")" "$file"
  fi
}

MACHO_FILES=()
while IFS= read -r macho; do
  [[ -n "$macho" ]] && MACHO_FILES+=("$macho")
done < <(
  find "$DEST_FRAMEWORK" -type f ! -path '*.dSYM/*' -print0 | while IFS= read -r -d '' file; do
    if file "$file" | grep -q 'Mach-O'; then
      printf '%s\n' "$file"
    fi
  done
)

for macho in "${MACHO_FILES[@]}"; do
  rewrite_binary "$macho"
done

{
  echo "bundled_runtime=apple-silicon-r"
  echo "source_framework=$FRAMEWORK_SRC"
  echo "current_version=$CURRENT_VERSION"
  echo "include_seurat=$INCLUDE_SEURAT"
  echo "extra_packages=$EXTRA_PACKAGES"
  echo "bundle_full_library=$BUNDLE_FULL_LIBRARY"
  echo "package_count=${#PACKAGE_LIST[@]}"
} > "$METADATA_PATH"

cat > "$WRAPPER_PATH" <<'EOF'
#!/bin/sh
set -eu

BUNDLE_ROOT="$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)"
FRAMEWORK_ROOT="$BUNDLE_ROOT/R.framework"

resolve_r_home_dir() {
  if [ -d "$FRAMEWORK_ROOT/Resources" ]; then
    printf '%s\n' "$FRAMEWORK_ROOT/Resources"
    return 0
  fi

  if [ -L "$FRAMEWORK_ROOT/Versions/Current" ] && [ -d "$FRAMEWORK_ROOT/Versions/Current/Resources" ]; then
    printf '%s\n' "$FRAMEWORK_ROOT/Versions/Current/Resources"
    return 0
  fi

  for version_dir in "$FRAMEWORK_ROOT"/Versions/*; do
    if [ -d "$version_dir/Resources" ]; then
      printf '%s\n' "$version_dir/Resources"
      return 0
    fi
  done

  return 1
}

R_HOME_DIR="$(resolve_r_home_dir || true)"
R_EXEC="$R_HOME_DIR/bin/exec/R"
R_LIBRARY="$R_HOME_DIR/library"

if [ -z "$R_HOME_DIR" ] || [ ! -x "$R_EXEC" ]; then
  echo "Bundled R executable was not found at $R_EXEC" >&2
  exit 1
fi

export R_HOME="$R_HOME_DIR"
export R_SHARE_DIR="$R_HOME_DIR/share"
export R_INCLUDE_DIR="$R_HOME_DIR/include"
export R_DOC_DIR="$R_HOME_DIR/doc"
export R_LIBS="$R_LIBRARY"
export R_LIBS_SITE="$R_LIBRARY"
export R_LIBS_USER="$R_LIBRARY"
export R_PROFILE_USER=
export R_ENVIRON_USER=

case "${1:-}" in
  "" )
    exec "$R_EXEC" --vanilla --no-echo
    ;;
  -h|--help)
    exec "$R_EXEC" --help
    ;;
  --version)
    exec "$R_EXEC" --version
    ;;
  -e)
    if [ "$#" -lt 2 ]; then
      echo "Missing expression for -e" >&2
      exit 1
    fi
    expr="$2"
    shift 2
    exec "$R_EXEC" --vanilla --no-echo -e "$expr" --args "$@"
    ;;
  *)
    script_path="$1"
    shift
    exec "$R_EXEC" --vanilla --no-echo --file="$script_path" --args "$@"
    ;;
esac
EOF

chmod +x "$WRAPPER_PATH"
for macho in "${MACHO_FILES[@]}"; do
  codesign --force --sign - --timestamp=none "$macho"
done

echo "Staged bundled R runtime at $DEST_FRAMEWORK"
echo "Copied ${#PACKAGE_LIST[@]} packages into $DEST_LIBRARY"
