const invoke = window.__TAURI__?.core?.invoke;

const state = {
  backendInfo: null,
  inspect: null,
  config: null,
  lastOutputPath: null,
  busy: false
};

const elements = {
  rscriptPath: document.getElementById("rscript-path"),
  backendSummary: document.getElementById("backend-summary"),
  inputPath: document.getElementById("input-path"),
  metadataInputPath: document.getElementById("metadata-input-path"),
  inspectButton: document.getElementById("inspect-button"),
  buildButton: document.getElementById("build-button"),
  openOutputButton: document.getElementById("open-output-button"),
  guessPaths: document.getElementById("guess-paths"),
  browseRscript: document.getElementById("browse-rscript"),
  browseInput: document.getElementById("browse-input"),
  browseMetadataInput: document.getElementById("browse-metadata-input"),
  chooseOutput: document.getElementById("choose-output"),
  inspectEmpty: document.getElementById("inspect-empty"),
  inspectContent: document.getElementById("inspect-content"),
  defaultsSummary: document.getElementById("defaults-summary"),
  capabilitiesSummary: document.getElementById("capabilities-summary"),
  columnSummary: document.getElementById("column-summary"),
  configOutput: document.getElementById("config-output"),
  configTitle: document.getElementById("config-title"),
  configGroupby: document.getElementById("config-groupby"),
  configInitialColor: document.getElementById("config-initial-color"),
  configAssay: document.getElementById("config-assay"),
  configTheme: document.getElementById("config-theme"),
  configTopGenes: document.getElementById("config-top-genes"),
  configNeighborMode: document.getElementById("config-neighbor-mode"),
  configMarkerTest: document.getElementById("config-marker-test"),
  configLightweight: document.getElementById("config-lightweight"),
  additionalColors: document.getElementById("additional-colors"),
  interactionGroupby: document.getElementById("interaction-groupby"),
  configPreview: document.getElementById("config-preview"),
  logOutput: document.getElementById("log-output")
};

function log(message, mode = "append") {
  if (mode === "replace") {
    elements.logOutput.textContent = message;
    return;
  }
  elements.logOutput.textContent = `${elements.logOutput.textContent}\n${message}`.trim();
  elements.logOutput.scrollTop = elements.logOutput.scrollHeight;
}

function setBusy(isBusy) {
  state.busy = isBusy;
  elements.inspectButton.disabled = isBusy;
  elements.buildButton.disabled = isBusy || !state.config;
  elements.openOutputButton.disabled = !state.lastOutputPath || isBusy;
}

function requireInvoke() {
  if (!invoke) {
    throw new Error("Tauri bridge is not available. Run this UI inside the Tauri app.");
  }
}

function pathPayload() {
  return {
    rscriptPath: elements.rscriptPath.value.trim(),
    inputPath: elements.inputPath.value.trim(),
    metadataInputPath: elements.metadataInputPath.value.trim() || null
  };
}

function assertRequiredPaths() {
  const paths = pathPayload();
  if (!paths.rscriptPath) throw new Error("Rscript path is required.");
  if (!paths.inputPath) throw new Error("Input RDS path is required.");
  return paths;
}

function formatBackendSummary(info = null) {
  if (!info?.backendVersion) {
    return "KaroSpaceR is bundled with the app. Only Rscript is required locally.";
  }

  const revision = info.backendRevision ? ` (${info.backendRevision})` : "";
  return `Bundled KaroSpaceR ${info.backendVersion}${revision}. Only Rscript is required locally.`;
}

function setSelectOptions(select, values, allowEmpty = false) {
  select.innerHTML = "";
  if (allowEmpty) {
    const option = document.createElement("option");
    option.value = "";
    option.textContent = "<none>";
    select.append(option);
  }
  values.forEach((value) => {
    const option = document.createElement("option");
    option.value = value;
    option.textContent = value;
    select.append(option);
  });
}

function renderSummaryList(target, entries) {
  target.innerHTML = "";
  entries.forEach(([label, value]) => {
    const dt = document.createElement("dt");
    dt.textContent = label;
    const dd = document.createElement("dd");
    dd.textContent = value ?? "<none>";
    target.append(dt, dd);
  });
}

function renderChipGroup(target, values, selectedValues, onToggle) {
  target.innerHTML = "";
  values.forEach((value) => {
    const label = document.createElement("label");
    label.className = "checkbox-chip";
    const input = document.createElement("input");
    input.type = "checkbox";
    input.checked = selectedValues.includes(value);
    input.addEventListener("change", () => onToggle(value, input.checked));
    const span = document.createElement("span");
    span.textContent = value;
    label.append(input, span);
    target.append(label);
  });
}

function sortedColumnNames(type) {
  const summary = state.inspect?.columns?.summary ?? [];
  return summary
    .filter((column) => (type === "categorical" ? column.type === "categorical" : true))
    .map((column) => column.name);
}

function normalizeConfigArrays(config) {
  for (const key of ["additional_colors", "marker_genes_groupby", "interaction_markers_groupby"]) {
    if (config[key] == null) continue;
    if (!Array.isArray(config[key])) {
      config[key] = [config[key]];
    }
  }
}

function updateConfigFromForm() {
  if (!state.config) {
    return;
  }

  state.config.output = elements.configOutput.value.trim() || null;
  state.config.title = elements.configTitle.value.trim() || null;
  state.config.groupby = elements.configGroupby.value || null;
  state.config.initial_color = elements.configInitialColor.value || null;
  state.config.assay = elements.configAssay.value || null;
  state.config.theme = elements.configTheme.value || "light";
  state.config.neighbor_mode = elements.configNeighborMode.value || "spatial";
  state.config.marker_test = elements.configMarkerTest.value || "mean_diff";
  state.config.lightweight = elements.configLightweight.checked;

  const topGenes = elements.configTopGenes.value.trim();
  state.config.top_genes_n = topGenes ? Number.parseInt(topGenes, 10) : null;
  if (Number.isNaN(state.config.top_genes_n)) {
    state.config.top_genes_n = null;
  }

  normalizeConfigArrays(state.config);
  elements.configPreview.value = JSON.stringify(state.config, null, 2);
  elements.buildButton.disabled = state.busy;
}

function syncFormFromConfig() {
  if (!state.config) {
    return;
  }

  normalizeConfigArrays(state.config);
  elements.configOutput.value = state.config.output ?? "";
  elements.configTitle.value = state.config.title ?? "";
  elements.configGroupby.value = state.config.groupby ?? "";
  elements.configInitialColor.value = state.config.initial_color ?? "";
  elements.configAssay.value = state.config.assay ?? "";
  elements.configTheme.value = state.config.theme ?? "light";
  elements.configTopGenes.value = state.config.top_genes_n ?? "";
  elements.configNeighborMode.value = state.config.neighbor_mode ?? "spatial";
  elements.configMarkerTest.value = state.config.marker_test ?? "mean_diff";
  elements.configLightweight.checked = Boolean(state.config.lightweight);

  const allColors = sortedColumnNames("any");
  renderChipGroup(
    elements.additionalColors,
    allColors.filter((name) => name !== state.config.groupby && name !== state.config.initial_color),
    state.config.additional_colors ?? [],
    (value, checked) => {
      const next = new Set(state.config.additional_colors ?? []);
      if (checked) {
        next.add(value);
      } else {
        next.delete(value);
      }
      state.config.additional_colors = Array.from(next);
      updateConfigFromForm();
    }
  );

  renderChipGroup(
    elements.interactionGroupby,
    sortedColumnNames("categorical"),
    state.config.interaction_markers_groupby ?? [],
    (value, checked) => {
      const next = new Set(state.config.interaction_markers_groupby ?? []);
      if (checked) {
        next.add(value);
      } else {
        next.delete(value);
      }
      state.config.interaction_markers_groupby = Array.from(next);
      updateConfigFromForm();
    }
  );

  elements.configPreview.value = JSON.stringify(state.config, null, 2);
}

function renderInspect() {
  if (!state.inspect) {
    elements.inspectEmpty.classList.remove("hidden");
    elements.inspectContent.classList.add("hidden");
    return;
  }

  elements.inspectEmpty.classList.add("hidden");
  elements.inspectContent.classList.remove("hidden");

  renderSummaryList(elements.defaultsSummary, [
    ["Group By", state.inspect.defaults.groupby],
    ["Initial Color", state.inspect.defaults.initial_color],
    ["Additional Colors", (state.inspect.defaults.additional_colors ?? []).join(", ") || "<none>"],
    ["Assay", state.inspect.defaults.assay ?? "<none>"],
    ["Output", state.inspect.defaults.output ?? "<none>"]
  ]);

  renderSummaryList(elements.capabilitiesSummary, [
    ["Source Class", (state.inspect.source_class ?? []).join(", ")],
    ["Assays", (state.inspect.assays.available ?? []).join(", ") || "<none>"],
    ["Graphs", (state.inspect.graphs.available ?? []).join(", ") || "<none>"],
    ["Gene Count", String(state.inspect.genes.count ?? 0)],
    ["Metadata Merge", state.inspect.metadata_merge ? `${state.inspect.metadata_merge.overlap_rows}/${state.inspect.metadata_merge.primary_rows}` : "<none>"]
  ]);

  elements.columnSummary.innerHTML = "";
  (state.inspect.columns.summary ?? []).forEach((column) => {
    const chip = document.createElement("div");
    chip.className = "column-chip";
    chip.textContent = `${column.name} · ${column.type} · ${Math.round((column.coverage ?? 0) * 100)}%`;
    elements.columnSummary.append(chip);
  });

  const categoricalColumns = sortedColumnNames("categorical");
  const allColorColumns = sortedColumnNames("any");
  setSelectOptions(elements.configGroupby, categoricalColumns);
  setSelectOptions(elements.configInitialColor, allColorColumns);
  setSelectOptions(elements.configAssay, state.inspect.assays.available ?? [], true);
  syncFormFromConfig();
}

async function callBackend(command, payload) {
  requireInvoke();
  return invoke(command, payload);
}

async function pickPath(mode, title, defaultName = null) {
  const result = await callBackend("pick_path", {
    request: {
      mode,
      title,
      defaultName
    }
  });
  return result || null;
}

async function guessPaths() {
  try {
    const result = await callBackend("guess_backend_paths", {});
    state.backendInfo = result;
    if (result.rscriptPath) elements.rscriptPath.value = result.rscriptPath;
    elements.backendSummary.textContent = formatBackendSummary(result);
    log(
      result.rscriptPath
        ? `Guessed Rscript path.\n${formatBackendSummary(result)}`
        : formatBackendSummary(result),
      "replace"
    );
  } catch (error) {
    log(String(error), "replace");
  }
}

async function inspectDataset() {
  try {
    const paths = assertRequiredPaths();
    setBusy(true);
    log("Running inspect...", "replace");
    const result = await callBackend("inspect_dataset", {
      request: {
        ...paths,
        assay: null,
        geneQuery: null,
        geneLimit: 50
      }
    });
    state.inspect = result.summary;
    state.config = structuredClone(result.summary.default_config);
    normalizeConfigArrays(state.config);
    state.lastOutputPath = null;
    renderInspect();
    updateConfigFromForm();
    log(result.stdout?.trim() || "Inspect completed.", "replace");
  } catch (error) {
    log(String(error), "replace");
  } finally {
    setBusy(false);
  }
}

async function buildViewer() {
  try {
    const paths = assertRequiredPaths();
    if (!state.config) {
      throw new Error("Run inspect before building.");
    }
    updateConfigFromForm();
    setBusy(true);
    log("Building viewer...", "replace");
    const result = await callBackend("build_viewer", {
      request: {
        ...paths,
        config: state.config
      }
    });
    state.lastOutputPath = result.outputPath ?? state.config.output ?? null;
    elements.openOutputButton.disabled = !state.lastOutputPath;
    log(
      [
        "Build completed.",
        result.configPath ? `Config: ${result.configPath}` : null,
        result.outputPath ? `Output: ${result.outputPath}` : null,
        result.stdout?.trim() || null,
        result.stderr?.trim() ? `stderr:\n${result.stderr.trim()}` : null
      ].filter(Boolean).join("\n\n"),
      "replace"
    );
  } catch (error) {
    log(String(error), "replace");
  } finally {
    setBusy(false);
  }
}

async function openOutput() {
  if (!state.lastOutputPath) {
    return;
  }
  try {
    await callBackend("open_path", { path: state.lastOutputPath });
  } catch (error) {
    log(String(error), "replace");
  }
}

async function boot() {
  try {
    await guessPaths();
  } catch (error) {
    log(String(error), "replace");
  }
}

elements.guessPaths.addEventListener("click", guessPaths);
elements.inspectButton.addEventListener("click", inspectDataset);
elements.buildButton.addEventListener("click", buildViewer);
elements.openOutputButton.addEventListener("click", openOutput);

elements.browseRscript.addEventListener("click", async () => {
  const path = await pickPath("file", "Choose Rscript");
  if (path) elements.rscriptPath.value = path;
});

elements.browseInput.addEventListener("click", async () => {
  const path = await pickPath("file", "Choose Input RDS");
  if (path) elements.inputPath.value = path;
});

elements.browseMetadataInput.addEventListener("click", async () => {
  const path = await pickPath("file", "Choose Metadata Input");
  if (path) elements.metadataInputPath.value = path;
});

elements.chooseOutput.addEventListener("click", async () => {
  const path = await pickPath("save", "Choose Output HTML", "viewer.html");
  if (path) {
    elements.configOutput.value = path;
    updateConfigFromForm();
  }
});

[
  elements.configOutput,
  elements.configTitle,
  elements.configGroupby,
  elements.configInitialColor,
  elements.configAssay,
  elements.configTheme,
  elements.configTopGenes,
  elements.configNeighborMode,
  elements.configMarkerTest,
  elements.configLightweight
].forEach((element) => {
  element.addEventListener("change", updateConfigFromForm);
  element.addEventListener("input", updateConfigFromForm);
});

boot();
