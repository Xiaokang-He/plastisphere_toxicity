# Repository copy of the RDA source-data workflow.
# Repository-only adaptation: input/output paths now use data/input and results.
# The Hellinger transformation, RDA model, envfit permutations and export logic are retained.

# ============================================================
# RDA figure source data
# Prefer running AFTER the RDA and envfit sections of the original script.
# Standalone mode reads the same three CSVs and fits only the RDA analysis.
#
# Existing fitted/plot objects are reused when all required objects exist.
# Actual plotted site coordinates are taken from PP$data or plotdata.
# Original script: RDA1 site score * 0.15; RDA2 site score * 0.20.
# Arrow endpoints come from RDAE / rda$CCA$biplot, with no site multiplier.
# ============================================================

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(repo_root, "data", "input")
results_dir <- file.path(repo_root, "results")

required_packages <- c("vegan", "writexl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

export_RDA_figure_source_data <- function(
    analysis_dir = input_dir,
    output_dir = file.path(results_dir, "supp_fig24_rda_source_data"),
    object_environment = .GlobalEnv,
    prefer_existing_objects = TRUE) {

  get_if_present <- function(name) {
    if (exists(name, envir = object_environment, inherits = FALSE)) {
      get(name, envir = object_environment, inherits = FALSE)
    } else NULL
  }
  candidate_plot <- get_if_present("plotdata")
  candidate_ggplot <- get_if_present("PP")
  if (is.list(candidate_ggplot) && is.data.frame(candidate_ggplot$data)) {
    candidate_plot <- candidate_ggplot$data
  }
  reuse_ok <- isTRUE(prefer_existing_objects) &&
    inherits(get_if_present("rda"), "rda") &&
    is.list(get_if_present("rdascore")) &&
    inherits(get_if_present("envfit"), "envfit") &&
    !is.null(get_if_present("env")) &&
    !is.null(get_if_present("sampledata")) &&
    is.data.frame(candidate_plot) &&
    all(c("sample", "RDAS1", "RDAS2", "group") %in% names(candidate_plot))

  if (reuse_ok) {
    execution_mode <- "Existing RDA plot/model objects and envfit results reused."
    fit_source_note <- "Existing RDA envfit r-squared and P values; permutations are not rerun."
  } else {
    original_files <- file.path(analysis_dir, c(
      "OTU - 副本.csv", "环境因子 - 副本.csv", "分组依据 - 副本.csv"
    ))
    if (!all(file.exists(original_files))) {
      stop("Run the original RDA/envfit code first, or supply these CSV files: ",
           paste(original_files[!file.exists(original_files)], collapse = "; "))
    }
    otu_work <- read.csv(original_files[1], header = TRUE, row.names = 1)
    env_work <- read.csv(original_files[2], header = TRUE, row.names = 1)
    group_work <- read.csv(original_files[3], header = FALSE, colClasses = "character")
    hellinger_work <- vegan::decostand(t(otu_work), method = "hellinger")
    if (nrow(env_work) != nrow(hellinger_work) ||
        nrow(group_work) != nrow(hellinger_work) || ncol(group_work) < 2L) {
      stop("OTU, indicator and grouping row counts do not match the original workflow.")
    }
    model_work <- vegan::rda(hellinger_work, env_work, scale = TRUE)
    scores_work <- vegan::scores(model_work)
    if (ncol(scores_work$sites) < 2L) stop("Fewer than two RDA axes.")
    plotted_work <- data.frame(
      sample = rownames(scores_work$sites),
      RDAS1 = scores_work$sites[, 1] * 0.15,
      RDAS2 = scores_work$sites[, 2] * 0.20,
      group = group_work[[2]], stringsAsFactors = FALSE
    )
    fit_work <- vegan::envfit(model_work, env_work, permutations = 999)
    object_environment <- list2env(list(
      rda = model_work, rdascore = scores_work,
      RDAE = as.data.frame(model_work$CCA$biplot[, 1:2, drop = FALSE]),
      plotdata = plotted_work, env = env_work, sampledata = hellinger_work,
      envfit = fit_work
    ), parent = baseenv())
    execution_mode <- "Standalone RDA export from the three original CSV inputs."
    fit_source_note <- paste(
      "RDA envfit recalculated with 999 permutations in standalone mode.",
      "The supplied script has no fixed seed; permutation P values may vary across runs."
    )
  }

  # 1. Read the existing analysis and plot objects -------------
  get_object <- function(name) {
    if (!exists(name, envir = object_environment, inherits = FALSE)) {
      stop("Missing object: ", name, ". Run the original RDA code first.")
    }
    get(name, envir = object_environment, inherits = FALSE)
  }
  model <- get_object("rda")
  if (!inherits(model, "rda")) {
    stop("Object 'rda' is not a fitted vegan RDA model.")
  }
  if (exists("PP", envir = object_environment, inherits = FALSE)) {
    plot_object <- get_object("PP")
    candidate <- plot_object$data
    if (is.data.frame(candidate) &&
        all(c("sample", "RDAS1", "RDAS2", "group") %in% names(candidate))) {
      plotted <- candidate
    } else {
      plotted <- get_object("plotdata")
    }
  } else {
    plotted <- get_object("plotdata")
  }
  plotted <- as.data.frame(plotted, check.names = FALSE)
  needed <- c("sample", "RDAS1", "RDAS2", "group")
  if (!all(needed %in% names(plotted))) {
    stop("The RDA plot data must contain sample, RDAS1, RDAS2 and group.")
  }
  if (anyNA(plotted$sample) || anyDuplicated(as.character(plotted$sample)) ||
      anyNA(plotted$group)) {
    stop("Missing/duplicate plotted sample IDs or missing group labels.")
  }

  ordination_scores <- get_object("rdascore")
  sites <- as.matrix(ordination_scores$sites)
  if (ncol(sites) < 2L) stop("Two site-score axes are required.")
  sample_order <- match(as.character(plotted$sample), rownames(sites))
  if (anyNA(sample_order)) {
    stop("A plotted sample could not be matched to rdascore$sites.")
  }
  sites <- sites[sample_order, 1:2, drop = FALSE]
  arrows <- if (exists("RDAE", envir = object_environment, inherits = FALSE)) {
    as.matrix(get_object("RDAE"))[, 1:2, drop = FALSE]
  } else {
    as.matrix(model$CCA$biplot)[, 1:2, drop = FALSE]
  }
  if (any(!is.finite(sites)) || any(!is.finite(arrows)) ||
      any(!is.finite(plotted$RDAS1)) || any(!is.finite(plotted$RDAS2))) {
    stop("Nonfinite site or arrow coordinates.")
  }

  # 2. Readable group/sample labels ---------------------------
  # Unknown group labels are retained without guessing their treatment.
  original_group <- as.character(plotted$group)
  group_key <- gsub("[[:space:]_-]", "", tolower(trimws(original_group)))
  group_map <- c(
    ctr = "Control", control = "Control",
    mp1 = "Microplastic-Low", mp10 = "Microplastic-High",
    mpp1 = "Plastisphere-Low", mpp10 = "Plastisphere-High",
    microplasticlow = "Microplastic-Low",
    microplastichigh = "Microplastic-High",
    plastispherelow = "Plastisphere-Low",
    plastispherehigh = "Plastisphere-High"
  )
  experimental_group <- unname(group_map[group_key])
  experimental_group[is.na(experimental_group)] <-
    original_group[is.na(experimental_group)]
  replicate_number <- as.integer(ave(
    seq_along(experimental_group), experimental_group, FUN = seq_along
  ))

  sample_coordinates <- data.frame(
    "Sample ID" = paste0(experimental_group, "_", replicate_number),
    "Original sample ID" = as.character(plotted$sample),
    "Experimental group" = experimental_group,
    "RDA1 plotted coordinate" = as.numeric(plotted$RDAS1),
    "RDA2 plotted coordinate" = as.numeric(plotted$RDAS2),
    "RDA1 unscaled site score" = unname(sites[, 1]),
    "RDA2 unscaled site score" = unname(sites[, 2]),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  scale_matches <- c(
    RDA1 = isTRUE(all.equal(
      as.numeric(plotted$RDAS1), unname(sites[, 1]) * 0.15,
      tolerance = 1e-10
    )),
    RDA2 = isTRUE(all.equal(
      as.numeric(plotted$RDAS2), unname(sites[, 2]) * 0.20,
      tolerance = 1e-10
    ))
  )
  # The actual plot coordinates above are preserved even if plotting
  # multipliers were changed after the supplied version of the script.
  indicator_vectors <- data.frame(
    "Indicator" = rownames(arrows),
    "Arrow start RDA1" = rep(0, nrow(arrows)),
    "Arrow start RDA2" = rep(0, nrow(arrows)),
    "Arrow end RDA1" = unname(arrows[, 1]),
    "Arrow end RDA2" = unname(arrows[, 2]),
    check.names = FALSE, stringsAsFactors = FALSE
  )

  # 3. Axis percentages: exactly the denominator in the plot ---
  eigenvalues <- as.numeric(model$CCA$eig)
  axis_names <- names(model$CCA$eig)
  if (is.null(axis_names)) axis_names <- paste0("RDA", seq_along(eigenvalues))
  if (length(eigenvalues) < 2L || sum(eigenvalues) <= 0) {
    stop("Insufficient constrained RDA axes.")
  }
  constrained_percent <- 100 * eigenvalues / sum(eigenvalues)
  reference_percent <- c(16.08, 14.23) # image labels; used only for checking
  axis_variance <- data.frame(
    "Axis" = axis_names,
    "Eigenvalue" = eigenvalues,
    "Constrained variance explained (%)" = constrained_percent,
    "Total variance explained (%)" = 100 * eigenvalues / model$tot.chi,
    "Axis label percentage" = round(constrained_percent, 2),
    "Shown in figure" = seq_along(eigenvalues) <= 2L,
    check.names = FALSE
  )
  figure_correspondence <- data.frame(
    "Axis" = axis_names[1:2],
    "Calculated axis label (%)" = round(constrained_percent[1:2], 2),
    "Reference image label (%)" = reference_percent,
    "Matches reference image" = round(constrained_percent[1:2], 2) == reference_percent,
    "Expected site multiplier" = c(0.15, 0.20),
    "Plotted coordinates match multiplier" = unname(scale_matches),
    check.names = FALSE
  )

  # 4. Reuse the existing RDA envfit result --------------------
  # This is variable fitting, not the overall RDA model P value.
  fit <- get_object("envfit")
  if (!inherits(fit, "envfit")) {
    stop("Run envfit on the RDA model at the end of the original script first.")
  }
  fit_rows <- list()
  for (component in c("vectors", "factors")) {
    block <- fit[[component]]
    if (!is.null(block) && length(block$r)) {
      fit_names <- names(block$r)
      if (is.null(fit_names)) {
        stop("The existing envfit result has no variable names.")
      }
      fit_p <- block$pvals
      if (length(fit_p) != length(block$r)) {
        stop("The existing envfit r-squared and P-value lengths differ.")
      }
      if (!is.null(names(fit_p))) fit_p <- fit_p[fit_names]
      fit_rows[[component]] <- data.frame(
        "Variable" = fit_names,
        "Variable type" = if (component == "vectors") "Numeric vector" else "Factor",
        "Envfit r squared" = unname(block$r),
        "Permutation P value" = unname(fit_p),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    }
  }
  envfit_results <- if (length(fit_rows)) {
    do.call(rbind, fit_rows)
  } else {
    data.frame(
      "Variable" = character(), "Variable type" = character(),
      "Envfit r squared" = numeric(), "Permutation P value" = numeric(),
      check.names = FALSE
    )
  }
  rownames(envfit_results) <- NULL

  # 5. Existing inputs and model summary ----------------------
  env_data <- as.data.frame(get_object("env"), check.names = FALSE)
  hellinger <- as.matrix(get_object("sampledata"))
  host_indicator_data <- data.frame(
    "Original sample ID" = rownames(env_data),
    env_data, row.names = NULL, check.names = FALSE
  )
  otu_hellinger <- data.frame(
    "Original sample ID" = rownames(hellinger),
    hellinger, row.names = NULL, check.names = FALSE
  )
  r2 <- vegan::RsquareAdj(model)
  model_statistics <- data.frame(
    "Statistic" = c(
      "Number of samples", "Number of OTU features",
      "Number of explanatory variables", "Constrained model rank",
      "Total inertia", "Constrained inertia", "Unconstrained inertia",
      "R squared", "Adjusted R squared"
    ),
    "Value" = c(
      nrow(sample_coordinates), ncol(hellinger), ncol(env_data),
      model$CCA$rank, model$tot.chi, model$CCA$tot.chi,
      if (is.null(model$CA)) 0 else model$CA$tot.chi,
      r2$r.squared, r2$adj.r.squared
    ),
    check.names = FALSE
  )
  group_mapping <- data.frame(
    "Sample ID" = sample_coordinates[["Sample ID"]],
    "Original sample ID" = as.character(plotted$sample),
    "Original group label" = original_group,
    "Experimental group" = experimental_group,
    "Within-group sequence" = replicate_number,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  input_names <- c(
    "OTU - 副本.csv", "环境因子 - 副本.csv", "分组依据 - 副本.csv"
  )
  input_paths <- file.path(analysis_dir, input_names)
  input_available <- file.exists(input_paths)

  export_notes <- data.frame(
    "Item" = c(
      "Coordinate source", "Site score source", "Arrow source",
      "Expected RDA1 site multiplier", "Expected RDA2 site multiplier",
      "Arrow multiplier", "Axis percentage denominator",
      "Envfit statistics", "Overall/axis permutation tests",
      "Group assignment", "Replicate suffixes",
      "Reference figure agreement",
      "Original OTU file", "Original indicator file", "Original grouping file"
    ),
    "Description" = c(
      "Existing PP$data or RDA plotdata; plotted coordinates are exported unchanged.",
      "Existing rdascore$sites; both unscaled and plotted scores are retained.",
      "Existing RDAE, otherwise rda$CCA$biplot; not envfit unit arrows.",
      "0.15", "0.20", "1 in the supplied script.",
      "Sum of constrained RDA eigenvalues, matching the original axis-label code.",
      fit_source_note,
      "Not calculated by this exporter; the supplied code only performs envfit.",
      "The group assignment stored in the plotted data is retained.",
      "Readable suffixes follow current within-group row order; original sample IDs are retained.",
      if (all(figure_correspondence[["Matches reference image"]])) {
        "Calculated first-two-axis percentages match 16.08% and 14.23%."
      } else {
        "Calculated percentages differ from the supplied image; confirm the input/model version before submission."
      },
      if (input_available[1]) "Included as OTU_original and copied CSV." else "Not available at analysis_dir.",
      if (input_available[2]) "Copied CSV; the in-memory env data are included." else "In-memory env data are included; CSV not available.",
      if (input_available[3]) "Copied CSV; plotted sample-group mapping is included." else "Plotted sample-group mapping is included; CSV not available."
    ),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  export_notes <- rbind(
    data.frame("Item" = "Execution mode", "Description" = execution_mode,
               check.names = FALSE, stringsAsFactors = FALSE),
    export_notes
  )

  tables <- list(
    Sample_coordinates = sample_coordinates,
    Indicator_vectors = indicator_vectors,
    Axis_variance = axis_variance,
    Envfit_results = envfit_results,
    Model_statistics = model_statistics,
    Figure_correspondence = figure_correspondence,
    Sample_group_mapping = group_mapping,
    Host_indicator_data = host_indicator_data,
    OTU_Hellinger = otu_hellinger
  )
  if (input_available[1]) {
    otu_original <- read.csv(
      input_paths[1], header = TRUE, check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (!nzchar(names(otu_original)[1])) names(otu_original)[1] <- "OTU ID"
    tables$OTU_original <- otu_original
  }
  tables$Export_notes <- export_notes

  # 6. Validate and export ------------------------------------
  stopifnot(
    nrow(sample_coordinates) == nrow(plotted),
    !anyDuplicated(sample_coordinates[["Sample ID"]]),
    nrow(indicator_vectors) == nrow(arrows),
    all(is.finite(sample_coordinates[["RDA1 plotted coordinate"]])),
    all(is.finite(sample_coordinates[["RDA2 plotted coordinate"]]))
  )
  if (any(vapply(tables, nrow, integer(1)) > 1048575L)) {
    stop("An output table exceeds Excel's row limit.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create output folder.")
  xlsx_path <- file.path(output_dir, "RDA_Figure_Source_Data.xlsx")
  writexl::write_xlsx(tables, path = xlsx_path)

  for (table_name in c(
    "Sample_coordinates", "Indicator_vectors", "Axis_variance",
    "Envfit_results", "Figure_correspondence"
  )) {
    write.csv(
      tables[[table_name]], file.path(output_dir, paste0(table_name, ".csv")),
      row.names = FALSE, na = "NA", fileEncoding = "UTF-8"
    )
  }
  if (any(input_available)) {
    input_copy_dir <- file.path(output_dir, "Original_input_files")
    dir.create(input_copy_dir, recursive = TRUE, showWarnings = FALSE)
    copied <- file.copy(
      input_paths[input_available],
      file.path(input_copy_dir, input_names[input_available]),
      overwrite = TRUE
    )
    if (!all(copied)) warning("Some original input CSV files could not be copied.")
  }
  cat("\nRDA source data exported: ", xlsx_path, "\n", sep = "")
  print(figure_correspondence, row.names = FALSE)
  if (!all(figure_correspondence[["Matches reference image"]])) {
    message("Axis percentages do not match the supplied image; review Figure_correspondence.")
  }
  invisible(tables)
}

RDA_source_data <- export_RDA_figure_source_data()
