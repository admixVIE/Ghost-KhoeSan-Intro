library(demografr)
library(slendr)
init_env(uv=T)
.libPaths( c("/lisc/data/scratch/admixlab/mk_data/rlib/", .libPaths()) )  
#library(tibble, quietly = TRUE)
library(e1071, quietly = TRUE)
library(GenomicRanges, quietly = TRUE)
library(reticulate, quietly = TRUE)
library(dplyr)
library(scales)
py_require("scipy")
library(ggplot2)
library(patchwork)

'%ni%' <- Negate('%in%')
options("scipen"=100)


## load data
cv<-list(); posts<-list()
load(file="/lisc/data/scratch/admixlab/mk_data/san/simul_red/red_selection.Robject")
cv[["reduced"]]<-cv_sel
posts[["reduced"]]<-post_sel
load(file="/lisc/data/scratch/admixlab/mk_data/san/simul_wss/wss_selection.Robject")
cv[["wss"]]<-cv_sel
posts[["wss"]]<-post_sel
load(file="/lisc/data/scratch/admixlab/mk_data/san/simul_nog/nog_selection.Robject")
cv[["noghost"]]<-cv_sel
posts[["noghost"]]<-post_sel

load(file="/lisc/data/scratch/admixlab/mk_data/san/abc_inference.Robject")
load(file="/lisc/data/scratch/admixlab/mk_data/san/abc_cross_val.Robject")

######################################################################
### panel A: demographic model
######################################################################

posterior_medians <- apply(
  myabc$adj.values,
  2,
  median,
  na.rm = TRUE
)

list2env(
  as.list(posterior_medians),
  envir = .GlobalEnv
)

A <- population("Ancestor", time = 60001,    N = 18507,remove = T_AGhost-10)
Ghost <- population("Ghost", time = T_AGhost, N = Ne_Ghost, parent = A,remove = 280)
AH <- population("Human_ancestor", time = T_AGhost, N = 18507, parent = A, remove = T_KS-10)
KS <- population("Khoe_San", time = T_KS, N = Ne_KS, parent = AH,remove = 100)
RuN <- population("ARHG_RHGn", time = T_KS, N = Ne_RuN, parent = AH, remove = T_ARHG-10)
ARHG <- population("ARHG", time = T_ARHG, N = Ne_ARHG, parent = RuN, remove = T_RHG-10)
RHGn <- population("RHGn", time = T_ARHG, N = Ne_RHGn, parent = RuN,remove = 100)
wRHG <- population("wRHG", time = T_RHG, N = Ne_wRHG, parent = ARHG,remove = 100)
eRHG <- population("eRHG", time = T_RHG, N = Ne_eRHG, parent = ARHG,remove = 100)

gf<- list( gene_flow(from = KS, to = RuN, start = 9679, end = 9678, proportion = gf_KR),
           gene_flow(from = RuN, to = KS, start = 9678, end = 9677, proportion = gf_RK),
           gene_flow(from = KS, to = RHGn, start = 4490, end = 4489, proportion = gf_KN),
           gene_flow(from = RHGn, to = KS, start = 4491, end = 4490, proportion = gf_NK),
           gene_flow(from = KS, to = ARHG, start = 4408, end = 4407, proportion = gf_KA),
           gene_flow(from = ARHG, to = KS, start = 4407, end = 4406, proportion = gf_AK),
           gene_flow(from = RHGn, to = ARHG, start = 2756, end = 2755, proportion = gf_NA),
           gene_flow(from = ARHG, to = RHGn, start = 2757, end = 2756, proportion = gf_AN),
           gene_flow(from = KS, to = RHGn, start = 1383, end = 1382, proportion = gf_KN2),
           gene_flow(from = RHGn, to = KS, start = 1382, end = 1381, proportion = gf_NK2),
           gene_flow(from = wRHG, to = eRHG, start = 411, end = 410, proportion = gf_WE),
           gene_flow(from = eRHG, to = wRHG, start = 410, end = 409, proportion = gf_EW),
           gene_flow(from = eRHG, to = RHGn, start = 393, end = 392, proportion = gf_EN),
           gene_flow(from = RHGn, to = eRHG, start = 394, end = 393, proportion = gf_NE),
           gene_flow(from = wRHG, to = RHGn, start = 302, end = 301, proportion = gf_WN),
           gene_flow(from = RHGn, to = wRHG, start = 303, end = 302, proportion = gf_NW),
           gene_flow(from = Ghost, to = KS, start = 303, end = 302, proportion = gf_Ghost))

model <- compile_model(
  populations = list(A, Ghost, AH, RuN, ARHG,KS, RHGn, wRHG, eRHG),
  gene_flow = gf,
  generation_time = 1)

samples <- schedule_sampling(
  model, times = 100,
  list(KS, 25), list(RHGn, 50), list(eRHG, 11), list(wRHG, 17),
  strict = TRUE)

population_colours <- c(
  Ancestor          = "#DEDEDE", # light grey
  Human_ancestor         = "#A9A9A9", # darker grey
  Ghost      = "#45B8B0", # turquoise
  ARHG_RHGn  = "#D99AC4", # light magenta
  ARHG       = "#C4ADDC", # light violet
  Khoe_San         = "#D55E5E", # stronger red
  RHGn       = "#4C91C6", # stronger blue
  eRHG       = "#9270B8", # stronger violet
  wRHG       = "#E4933A"  # stronger orange
)

samples_plot <- samples
samples_plot$time[samples_plot$time == 0] <- 150

model_plot <- plot_model(
  model,
  log = FALSE,
  samples = samples_plot
) +
  scale_fill_manual(values = population_colours, drop = FALSE) +
  scale_colour_manual(values = population_colours, drop = FALSE) +
  scale_y_continuous(
    trans = scales::log1p_trans(),
    limits = c(90, NA),
    breaks = c(250, 500, 1000, 5000, 10000, 20000, 50000),
    labels = scales::label_comma(),
    expand = expansion(mult = c(0, 0.04))
  ) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(3, 2, 3, 2, unit = "mm"),
    legend.position = "none"
  )


panel_a <- wrap_elements(
  full = model_plot,
  clip = FALSE
)



######################################################################
### panel D: confusion matrices
######################################################################
compact_heatmaps <- function(plot_list, guide_title) {
  Map(
    function(p, i) {
      q <- p +
        guides(fill = guide_colorbar(
          title = guide_title,
          direction = "vertical",
          title.position = "top",
          barheight = unit(22, "mm"),
          barwidth = unit(3, "mm")
        )) +
        theme(
          plot.margin = margin(1, 1, 1, 1, unit = "mm"),
          axis.text.x = element_text(size = 6, lineheight = 0.8),
          axis.text.y = element_text(size = 6, lineheight = 0.8),
          axis.title = element_text(size = 7),
          plot.title = element_text(size = 8),
          plot.subtitle = element_text(size = 6.5),
          legend.title = element_text(size = 6.5),
          legend.text = element_text(size = 6),
          legend.margin = margin(0, 0, 0, 0)
        )
      
      # Use shared-looking axis titles
      if (i != 1L) q <- q + labs(y = NULL)
      if (i != 2L) q <- q + labs(x = NULL)
      
      q
    },
    plot_list,
    seq_along(plot_list)
  )
}


panel_order <- c("noghost", "wss", "reduced")

model_order <- list(
  noghost = c("Ghost", "NoGhost"),
  wss     = c("Ghost", "NoGhost", "WSS_puls", "WSS_cont"),
  reduced = c("Ghost", "NoGhost", "WSS_puls", "WSS_cont")
)

setting_labels <- c(
  noghost = "Full model",
  wss = "WSS comparison",
  reduced = "WSS reduced"
)

model_labels <- c(
  Ghost = "Ghost",
  NoGhost = "No\nghost",
  WSS_puls = "WSS\npulse",
  WSS_cont = "WSS\ncontinuous"
)

theme_pub <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    plot.margin = margin(3, 4, 3, 4)
  )

norm_name <- function(x)
  gsub("[^[:alnum:]]", "", tolower(x))

quiet_summary <- function(x) {
  ans <- NULL
  invisible(capture.output(
    ans <- suppressWarnings(summary(x))
  ))
  ans
}

pluck_first <- function(x, candidates) {
  if (!is.list(x)) return(NULL)
  
  nx <- norm_name(names(x))
  for (nm in norm_name(candidates)) {
    j <- match(nm, nx)
    if (!is.na(j)) return(x[[j]])
  }
  NULL
}

align_vector <- function(x, levels) {
  old_names <- names(x)
  x <- as.numeric(x)
  
  if (!is.null(old_names)) {
    j <- match(norm_name(levels), norm_name(old_names))
    if (all(!is.na(j))) x <- x[j]
  }
  
  if (length(x) != length(levels))
    stop("Could not align model probabilities with model names.")
  
  names(x) <- levels
  x
}

#### 
cv_key <- function(x) {
  gsub("[^[:alnum:]]", "", tolower(as.character(x)))
}

extract_mean_model_probs <- function(x, levels, tol = NULL) {
  
  if (is.null(x$model.probs))
    stop("`model.probs` was not found.")
  
  if (is.null(tol)) {
    if (length(x$model.probs) != 1L)
      stop("Multiple tolerances found; specify `tol`.")
    i <- 1L
  } else {
    i <- match(tol, x$tols)
    if (is.na(i))
      stop("Tolerance not found: ", tol)
  }
  
  probabilities <- as.matrix(x$model.probs[[i]])
  true_model <- as.character(x$true)
  
  probability_names <- colnames(probabilities)
  
  if (is.null(probability_names))
    probability_names <- x$names$models
  
  j <- match(cv_key(levels), cv_key(probability_names))
  
  if (anyNA(j)) {
    stop(
      "Missing probability columns: ",
      paste(levels[is.na(j)], collapse = ", ")
    )
  }
  
  probabilities <- probabilities[, j, drop = FALSE]
  colnames(probabilities) <- levels
  
  result <- lapply(levels, function(model) {
    
    rows <- cv_key(true_model) == cv_key(model)
    values <- probabilities[rows, , drop = FALSE]
    values <- values[complete.cases(values), , drop = FALSE]
    
    if (!nrow(values))
      stop("No valid cross-validation samples for ", model)
    
    data.frame(
      true_model = model,
      posterior_model = levels,
      mean_probability = colMeans(values),
      n = nrow(values)
    )
  })
  
  do.call(rbind, result)
}

pretty_model_labels <- function(x) {
  result <- unname(model_labels[x])
  result[is.na(result)] <- x[is.na(result)]
  result
}

make_mean_probability_plot <- function(
    x, setting, tol = NULL) {
  
  levels <- model_order[[setting]]
  
  d <- extract_mean_model_probs(
    x,
    levels = levels,
    tol = tol
  )
  
  d$true_model <- factor(
    d$true_model,
    levels = rev(levels)
  )
  
  d$posterior_model <- factor(
    d$posterior_model,
    levels = levels
  )
  
  d$light_text <- d$mean_probability > 0.55
  
  ggplot(
    d,
    aes(
      x = posterior_model,
      y = true_model,
      fill = mean_probability
    )
  ) +
    geom_tile(
      colour = "white",
      linewidth = 0.4
    ) +
    geom_text(
      aes(
        label = sprintf("%.2f", mean_probability),
        colour = light_text
      ),
      size = 2.3
    ) +
    scale_colour_manual(
      values = c(
        `FALSE` = "#222222",
        `TRUE` = "white"
      ),
      guide = "none"
    ) +
    scale_fill_gradient(
      low = "#F7FBFF",
      high = "#08519C",
      limits = c(0, 1),
      name = "Mean posterior\nprobability"
    ) +
    scale_x_discrete(labels = pretty_model_labels) +
    scale_y_discrete(labels = pretty_model_labels) +
    labs(
      title = setting_labels[[setting]],
      x = "Model receiving probability",
      y = "True model"
    ) +
    guides(
      fill = guide_colorbar(
        direction = "vertical",
        title.position = "top",
        barheight = grid::unit(20, "mm"),
        barwidth = grid::unit(3, "mm")
      )
    ) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(colour = "black", size = 6),
      axis.title = element_text(size = 7),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 8
      ),
      plot.margin = margin(1, 1, 1, 1, unit = "mm")
    )
}

cv_d <- cv[panel_order]

plots_d <- Map(
  function(setting, i) {
    p <- make_mean_probability_plot(
      cv_d[[setting]],
      setting = setting,
      tol = NULL
    )
    
    if (i != 1L) p <- p + labs(y = NULL)
    if (i != 2L) p <- p + labs(x = NULL)
    
    p
  },
  panel_order,
  seq_along(panel_order)
)

plots_d <- compact_heatmaps(
  plots_d,
  guide_title = "Row\nproportion"
)


panel_d <- (
  wrap_plots(
    plots_d,
    nrow = 1,
    guides = "collect"
  ) &
    theme(legend.position = "right")
)




######################################################################
### panel E: Bayes factors
######################################################################

panel_order <- c("noghost", "wss", "reduced")

model_order <- list(
  noghost = c("Ghost", "NoGhost"),
  wss     = c("Ghost", "NoGhost", "WSS_puls", "WSS_cont"),
  reduced = c("Ghost", "NoGhost", "WSS_puls", "WSS_cont")
)

setting_labels <- c(
  noghost = "Full model",
  wss = "WSS comparison",
  reduced = "WSS reduced"
)

model_labels <- c(
  Ghost    = "Ghost",
  NoGhost  = "No\nghost",
  WSS_puls = "WSS\npulse",
  WSS_cont = "WSS\ncontinuous"
)

theme_bf <- theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(colour = "black", size = 6),
    axis.title = element_text(colour = "black",size = 7),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
    plot.subtitle = element_text(hjust = 0.5, size = 7),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    plot.margin = margin(1, 0, 1, 0, unit = "mm")
  )


bf_key <- function(x) {
  gsub("[^[:alnum:]]", "", tolower(as.character(x)))
}

same_models <- function(x, y) {
  if (is.null(x) || length(x) != length(y))
    return(FALSE)
  
  x <- bf_key(x)
  y <- bf_key(y)
  
  !anyDuplicated(x) && setequal(x, y)
}

find_postpr_fits <- function(x) {
  fits <- list()
  
  walk <- function(y) {
    if (inherits(y, "postpr")) {
      fits[[length(fits) + 1L]] <<- y
      return(invisible(NULL))
    }
    
    if (is.list(y)) {
      for (i in seq_along(y))
        walk(y[[i]])
    }
  }
  
  walk(x)
  fits
}

read_postpr_fit <- function(fit, levels) {
  
  if (is.null(fit$pred))
    return(NULL)
  
  expected <- bf_key(levels)
  
  ## Full model set recorded in the fit
  fit_models <- NULL
  
  if (!is.null(fit$names$models))
    fit_models <- as.character(fit$names$models)
  
  if (is.null(fit_models) && !is.null(names(fit$nmodels)))
    fit_models <- names(fit$nmodels)
  
  if (!is.null(fit_models) &&
      !same_models(fit_models, levels)) {
    return(NULL)
  }
  
  ## Posterior probabilities
  pred <- fit$pred
  p_values <- as.numeric(pred)
  p_names <- names(pred)
  
  if (is.null(p_names)) {
    dn <- dimnames(pred)
    if (length(dn) == 1L)
      p_names <- dn[[1]]
  }
  
  if (is.null(p_names) &&
      !is.null(fit_models) &&
      length(fit_models) == length(p_values)) {
    p_names <- fit_models
  }
  
  if (is.null(p_names) &&
      length(p_values) == length(levels)) {
    p_names <- levels
  }
  
  if (is.null(p_names) ||
      length(p_names) != length(p_values)) {
    return(NULL)
  }
  
  found <- bf_key(p_names)
  
  if (anyDuplicated(found) ||
      !all(found %in% expected)) {
    return(NULL)
  }
  
  # Missing categories are permitted only if the fit records
  # the complete expected model set
  if (!setequal(found, expected) &&
      is.null(fit_models)) {
    return(NULL)
  }
  
  if (any(!is.finite(p_values)) ||
      any(p_values < 0)) {
    return(NULL)
  }
  
  posterior <- setNames(rep(0, length(levels)), levels)
  posterior[match(found, expected)] <- p_values
  
  if (sum(posterior) <= 0)
    return(NULL)
  
  posterior <- posterior / sum(posterior)
  
  ## Prior probabilities: use reference-table model counts
  prior <- setNames(rep(1 / length(levels), length(levels)), levels)
  
  if (!is.null(fit$nmodels)) {
    n_values <- as.numeric(fit$nmodels)
    n_names <- names(fit$nmodels)
    
    if (is.null(n_names) &&
        !is.null(fit_models) &&
        length(fit_models) == length(n_values)) {
      n_names <- fit_models
    }
    
    if (same_models(n_names, levels)) {
      n_values <- n_values[
        match(expected, bf_key(n_names))
      ]
      
      if (all(is.finite(n_values)) &&
          all(n_values > 0)) {
        prior <- setNames(n_values / sum(n_values), levels)
      }
    }
  }
  
  method <- if (is.null(fit$method)) {
    "unknown"
  } else {
    as.character(fit$method)[1]
  }
  
  list(
    posterior = posterior,
    prior = prior,
    method = method
  )
}

select_postpr_fit <- function(
    x, levels,
    prefer = c("neuralnet", "rejection")) {
  
  fits <- find_postpr_fits(x)
  
  if (!length(fits))
    stop("No objects of class 'postpr' were found.")
  
  methods <- vapply(
    fits,
    function(fit) {
      if (is.null(fit$method)) ""
      else bf_key(fit$method[1])
    },
    character(1)
  )
  
  for (wanted in bf_key(prefer)) {
    candidates <- which(methods == wanted)
    
    for (i in candidates) {
      result <- read_postpr_fit(fits[[i]], levels)
      
      if (!is.null(result))
        return(result)
    }
  }
  
  available <- unique(methods[nzchar(methods)])
  
  stop(
    "No usable neuralnet or rejection result for models: ",
    paste(levels, collapse = ", "),
    ". Methods found: ",
    paste(available, collapse = ", ")
  )
}

format_bf_labels <- function(log_bf, cap = 10) {
  lim <- log10(cap)
  lower <- format(
    signif(1 / cap, 3),
    trim = TRUE,
    scientific = FALSE
  )
  
  result <- character(length(log_bf))
  
  result[log_bf > lim]  <- paste0(">", cap)
  result[log_bf < -lim] <- paste0("<", lower)
  
  middle <- abs(log_bf) <= lim
  
  result[middle] <- vapply(
    10^log_bf[middle],
    function(x) {
      format(signif(x, 2), trim = TRUE, scientific = FALSE)
    },
    character(1)
  )
  
  result
}

make_bf_plot <- function(
    x, setting, cap = 10,
    prefer = c("neuralnet", "rejection")) {
  
  if (cap <= 1)
    stop("`cap` must be greater than one.")
  
  lev <- model_order[[setting]]
  
  selected <- select_postpr_fit(
    x,
    levels = lev,
    prefer = prefer
  )
  
  posterior <- selected$posterior
  prior <- selected$prior
  
  # Avoid infinite values for zero rejection probabilities
  posterior <- pmax(posterior, 1e-12)
  posterior <- posterior / sum(posterior)
  
  model_evidence <- log10(posterior / prior)
  
  log_bf <- outer(
    model_evidence,
    model_evidence,
    FUN = "-"
  )
  
  dimnames(log_bf) <- list(lev, lev)
  diag(log_bf) <- 0
  
  d <- as.data.frame(as.table(log_bf))
  names(d) <- c("numerator", "denominator", "logBF")
  
  lim <- log10(cap)
  
  d$fill_value <- pmax(
    -lim,
    pmin(lim, d$logBF)
  )
  
  d$label <- format_bf_labels(d$logBF, cap)
  d$light_text <- abs(d$fill_value) > 0.55 * lim
  
  d$numerator <- factor(
    d$numerator,
    levels = rev(lev)
  )
  
  d$denominator <- factor(
    d$denominator,
    levels = lev
  )
  
  method_name <- switch(
    bf_key(selected$method),
    neuralnet = "Neural network",
    rejection = "Rejection",
    selected$method
  )
  
  lower_label <- format(
    signif(1 / cap, 3),
    trim = TRUE,
    scientific = FALSE
  )
  
  ggplot(
    d,
    aes(denominator, numerator, fill = fill_value)
  ) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(
      aes(label = label, colour = light_text),
      size = 2.2
    ) +
    scale_colour_manual(
      values = c(
        `FALSE` = "#222222",
        `TRUE` = "white"
      ),
      guide = "none"
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-lim, lim),
      breaks = c(-lim, 0, lim),
      labels = c(
        paste0("≤", lower_label),
        "1",
        paste0("≥", cap)
      ),
      name = "BF (row/column)"
    ) +
    scale_x_discrete(labels = model_labels) +
    scale_y_discrete(labels = model_labels) +
    coord_equal() +
    labs(
      title = setting_labels[[setting]],
      subtitle = method_name,
      x = "Denominator model",
      y = "Numerator model"
    ) +
    guides(
      fill = guide_colorbar(
        direction = "horizontal",
        title.position = "top"
      )
    ) +
    theme_bf
}

missing_settings <- setdiff(panel_order, names(posts))

if (length(missing_settings)) {
  stop(
    "Missing settings in `posts`: ",
    paste(missing_settings, collapse = ", ")
  )
}

# Check which method will be used
methods_used_b <- setNames(
  vapply(
    panel_order,
    function(z) {
      select_postpr_fit(
        posts[[z]],
        model_order[[z]]
      )$method
    },
    character(1)
  ),
  panel_order
)


plots_e <- lapply(
  panel_order,
  function(z) {
    make_bf_plot(
      posts[[z]],
      setting = z,
      cap = 10,
      prefer = c("neuralnet", "rejection")
    )
  }
)

names(plots_e) <- panel_order

plots_e <- compact_heatmaps(
  plots_e,
  guide_title = "BF\n(row/column)"
)


panel_e <- (
  wrap_plots(
    plots_e,
    nrow = 1,
    guides = "collect"
  ) &
    theme(legend.position = "right")
)



######################################################################
### panel B: prediction error
######################################################################

extract_prediction_error <- function(x, parameter) {
  
  key <- function(z)
    gsub("[^[:alnum:]]", "", tolower(z))
  
  true <- as.matrix(x$true)
  
  parameter_names <- colnames(true)
  if (is.null(parameter_names))
    parameter_names <- x$names$parameter.names
  
  j <- match(key(parameter), key(parameter_names))
  
  if (is.na(j)) {
    stop(
      "Parameter not found. Ghost-related parameters are: ",
      paste(
        grep("ghost", parameter_names,
             ignore.case = TRUE, value = TRUE),
        collapse = ", "
      )
    )
  }
  
  estimates <- x$estim
  
  ## Obtain tolerances
  tolerance <- suppressWarnings(
    as.numeric(sub("^tol", "", names(estimates)))
  )
  
  if (length(tolerance) != length(estimates) ||
      any(!is.finite(tolerance))) {
    tolerance <- x$tols
  }
  
  if (length(tolerance) != length(estimates))
    stop("Could not match tolerances to estimate matrices.")
  
  result <- lapply(seq_along(estimates), function(i) {
    
    estimated <- as.matrix(estimates[[i]])
    estimated_names <- colnames(estimated)
    
    if (!is.null(estimated_names)) {
      jj <- match(key(parameter), key(estimated_names))
    } else {
      jj <- j
    }
    
    if (is.na(jj))
      stop("Parameter missing from ", names(estimates)[i])
    
    observed  <- true[, j]
    predicted <- estimated[, jj]
    
    keep <- is.finite(observed) & is.finite(predicted)
    
    observed  <- observed[keep]
    predicted <- predicted[keep]
    
    denominator <- sum(
      (observed - mean(observed))^2
    )
    
    error <- if (length(observed) > 1 && denominator > 0) {
      sum((predicted - observed)^2) / denominator
    } else {
      NA_real_
    }
    
    data.frame(
      parameter = parameter_names[j],
      tolerance = tolerance[i],
      prediction_error = error,
      n = length(observed)
    )
  })
  
  result <- do.call(rbind, result)
  result[order(result$tolerance), ]
}


parameters_c <- c("Ne_Ghost", "T_AGhost", "gf_Ghost")

pe <- do.call(
  rbind,
  lapply(parameters_c, function(p) {
    extract_prediction_error(cv_abc, p)
  })
)

rownames(pe) <- NULL
pe$parameter <- factor(pe$parameter, levels = parameters_c)
pe<-pe[which(pe$tolerance<0.21),]

parameter_breaks <- c("Ne_Ghost", "T_AGhost", "gf_Ghost")

parameter_labels <- c(
  Ne_Ghost = "Ne (Ghost)",
  T_AGhost = "Ghost split time",
  gf_Ghost = "Gene flow proportion"
)

parameter_shapes <- c(
  Ne_Ghost = 16,
  T_AGhost = 17,
  gf_Ghost = 15
)

panel_b <- ggplot(
  pe,
  aes(
    x = tolerance,
    y = prediction_error,
    colour = parameter,
    shape = parameter
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.35,
    colour = "grey55"
  ) +
  geom_line(linewidth = 0.7, linetype = "solid") +
  geom_point(size = 2.7) +
  scale_x_continuous(
    breaks = sort(unique(pe$tolerance)),
    labels = function(x) format(
      x,
      scientific = FALSE,
      trim = TRUE,
      drop0trailing = TRUE
    )
  ) +
  scale_colour_manual(
    values = c(
      Ne_Ghost = "#0072B2",
      T_AGhost = "#D55E00",
      gf_Ghost = "#009E73"
    ),
    breaks = parameter_breaks,
    labels = parameter_labels,
    guide = guide_legend(
      nrow = 3,
      byrow = TRUE,
      override.aes = list(
        shape = unname(parameter_shapes[parameter_breaks]),
        linetype = "solid"
      )
    )
  ) +
  scale_shape_manual(
    values = parameter_shapes,
    breaks = parameter_breaks,
    labels = parameter_labels,
    guide = "none"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.12))
  ) +
  labs(
    title = "Prediction error",
    x = "Tolerance",
    y = "Normalized prediction error",
    colour = NULL,
    shape = NULL
  ) +
  theme_pub +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.justification = "center"
  )




######################################################################
### panel C: posteriors of ghost parameters
######################################################################

posterior_parameters <- c(
  "Ne_Ghost",
  "T_AGhost",
  "gf_Ghost"
)

parameter_titles <- list(
  Ne_Ghost = "N\u2091 (Ghost)",  # Nₑ (Ghost)
  T_AGhost = "Ghost split time",
  gf_Ghost = "Gene flow proportion"
)

posterior_key <- function(x) {
  gsub("[^[:alnum:]]", "", tolower(as.character(x)))
}

parameter_colours <- c(
  Ne_Ghost = "#0072B2",
  T_AGhost = "#D55E00",
  gf_Ghost = "#009E73"
)

interval_colours <- c(
  Ne_Ghost = "#D9EAF3",
  T_AGhost = "#F9E7D9",
  gf_Ghost = "#D9F0EA"
)

median_colours <- c(
  Ne_Ghost = "#004F7C",
  T_AGhost = "#8F3F00",
  gf_Ghost = "#006B4E"
)

extract_abc_posterior <- function(
    x,
    parameters,
    use_adjusted = TRUE) {
  
  # Prefer regression-adjusted posterior samples
  samples <- NULL
  sample_type <- NULL
  
  if (use_adjusted &&
      !is.null(x$adj.values) &&
      length(x$adj.values) > 0) {
    
    samples <- x$adj.values
    sample_type <- "Adjusted posterior"
    
  } else if (!is.null(x$unadj.values) &&
             length(x$unadj.values) > 0) {
    
    samples <- x$unadj.values
    sample_type <- "Unadjusted posterior"
    
  } else {
    stop(
      "No posterior samples were found in ",
      "`x$adj.values` or `x$unadj.values`."
    )
  }
  
  samples <- as.matrix(samples)
  
  parameter_names <- colnames(samples)
  
  if (is.null(parameter_names) &&
      !is.null(x$names$parameter.names)) {
    parameter_names <- x$names$parameter.names
  }
  
  if (is.null(parameter_names)) {
    stop(
      "Could not determine the names of the posterior parameters."
    )
  }
  
  if (length(parameter_names) != ncol(samples)) {
    stop(
      "The number of parameter names does not match ",
      "the number of posterior-sample columns."
    )
  }
  
  colnames(samples) <- parameter_names
  
  j <- match(
    posterior_key(parameters),
    posterior_key(parameter_names)
  )
  
  if (anyNA(j)) {
    missing_parameters <- parameters[is.na(j)]
    
    stop(
      "Parameters not found in the ABC posterior: ",
      paste(missing_parameters, collapse = ", "),
      "\nAvailable Ghost-related parameters: ",
      paste(
        grep(
          "ghost",
          parameter_names,
          ignore.case = TRUE,
          value = TRUE
        ),
        collapse = ", "
      )
    )
  }
  
  samples <- samples[, j, drop = FALSE]
  colnames(samples) <- parameters
  
  # ABC weights
  if (!is.null(x$weights) &&
      length(x$weights) == nrow(samples)) {
    
    weights <- as.numeric(x$weights)
    
  } else {
    
    weights <- rep(1, nrow(samples))
  }
  
  if (any(!is.finite(weights)) ||
      any(weights < 0) ||
      sum(weights) <= 0) {
    warning(
      "Invalid ABC weights found; using equal weights instead."
    )
    weights <- rep(1, nrow(samples))
  }
  
  weights <- weights / sum(weights)
  
  result <- do.call(
    rbind,
    lapply(seq_along(parameters), function(k) {
      data.frame(
        parameter = parameters[k],
        value = samples[, k],
        weight = weights,
        stringsAsFactors = FALSE
      )
    })
  )
  
  # Remove missing and infinite posterior values separately
  # for each parameter.
  result <- result[
    is.finite(result$value) &
      is.finite(result$weight) &
      result$weight >= 0,
    ,
    drop = FALSE
  ]
  
  if (!nrow(result))
    stop("No finite posterior samples remain.")
  
  # Renormalize the weights within each parameter
  result <- do.call(
    rbind,
    lapply(
      split(result, result$parameter),
      function(z) {
        if (sum(z$weight) <= 0)
          z$weight <- rep(1 / nrow(z), nrow(z))
        else
          z$weight <- z$weight / sum(z$weight)
        
        z
      }
    )
  )
  
  rownames(result) <- NULL
  
  result$parameter <- factor(
    result$parameter,
    levels = parameters
  )
  
  attr(result, "sample_type") <- sample_type
  result
}


posterior_c <- extract_abc_posterior(
  myabc,
  parameters = posterior_parameters,
  use_adjusted = TRUE
)


# Extract prior samples used by demografr/abc
extract_abc_prior <- function(x, parameters) {
  prior <- attr(x, "components")$parameters
  
  if (is.null(prior))
    stop("Prior samples were not found in attr(x, 'components')$parameters.")
  
  prior <- as.matrix(prior)
  
  key <- function(z)
    gsub("[^[:alnum:]]", "", tolower(z))
  
  j <- match(key(parameters), key(colnames(prior)))
  
  if (anyNA(j))
    stop("Missing prior parameters: ",
         paste(parameters[is.na(j)], collapse = ", "))
  
  do.call(rbind, lapply(seq_along(parameters), function(i) {
    data.frame(
      parameter = parameters[i],
      value = prior[, j[i]],
      weight = 1 / nrow(prior)
    )
  }))
}

prior_c <- extract_abc_prior(myabc, posterior_parameters)

weighted_quantile <- function(x, w, probs) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  
  ord <- order(x)
  x <- x[ord]
  w <- w[ord] / sum(w)
  
  approx(
    x = c(0, cumsum(w)),
    y = c(x[1], x),
    xout = probs,
    rule = 2,
    ties = "ordered"
  )$y
}

prepare_density_data <- function(data, parameter) {
  z <- data[
    as.character(data$parameter) == parameter &
      is.finite(data$value) &
      is.finite(data$weight) &
      data$weight > 0 &
      data$value >= 0,
    ,
    drop = FALSE
  ]
  
  if (nrow(z) < 2 || length(unique(z$value)) < 2)
    stop("Insufficient non-negative values for ", parameter)
  
  z$weight <- z$weight / sum(z$weight)
  z
}

weighted_density_zero <- function(z, n = 1024, adjust = 1) {
  x <- z$value
  w <- z$weight / sum(z$weight)
  
  bw <- stats::bw.nrd0(x)
  
  if (!is.finite(bw) || bw <= 0) {
    bw <- max(abs(x), 1) * 0.001
  }
  
  bw <- bw * adjust
  
  # Reflection around zero provides boundary correction
  reflected_x <- c(x, -x)
  reflected_w <- c(w, w) / 2
  
  d <- stats::density(
    reflected_x,
    weights = reflected_w,
    bw = bw,
    n = n,
    from = 0,
    to = max(x) + 3 * bw
  )
  
  data.frame(
    value = d$x,
    density = 2 * d$y
  )
}

density_interval <- function(density_data, lower, upper) {
  x <- sort(unique(c(
    lower,
    density_data$value[
      density_data$value > lower &
        density_data$value < upper
    ],
    upper
  )))
  
  data.frame(
    value = x,
    density = approx(
      density_data$value,
      density_data$density,
      xout = x,
      rule = 2
    )$y
  )
}

make_density_subplot <- function(parameter, position) {
  prior <- prepare_density_data(prior_c, parameter)
  posterior <- prepare_density_data(posterior_c, parameter)
  
  prior_density <- weighted_density_zero(prior)
  posterior_density <- weighted_density_zero(posterior)
  
  q <- weighted_quantile(
    posterior$value,
    posterior$weight,
    probs = c(0.025, 0.5, 0.975)
  )
  
  credible_density <- density_interval(
    posterior_density,
    lower = q[1],
    upper = q[3]
  )
  
  median_height <- approx(
    posterior_density$value,
    posterior_density$density,
    xout = q[2],
    rule = 2
  )$y
  
  median_data <- data.frame(
    value = q[2],
    density = median_height
  )

  posterior_colour <- parameter_colours[[parameter]]
  interval_colour  <- interval_colours[[parameter]]
  median_colour    <- median_colours[[parameter]]
  
  ggplot() +
    geom_ribbon(
      data = credible_density,
      aes(x = value, ymin = 0, ymax = density),
      fill = interval_colour,
      colour = NA
    ) +
    geom_line(
      data = prior_density,
      aes(
        x = value,
        y = density,
        linetype = "Prior"
      ),
      colour = "grey40",
      linewidth = 0.55
    ) +
    geom_segment(
      data = median_data,
      aes(
        x = value, xend = value,
        y = 0, yend = density
      ),
      inherit.aes = FALSE,
      colour = median_colour,
      linewidth = 0.5
    ) +
    geom_line(
      data = posterior_density,
      aes(
        x = value,
        y = density,
        linetype = "Posterior"
      ),
      colour = posterior_colour,
      linewidth = 0.75
    ) +
    scale_linetype_manual(
      values = c(
        Prior = "dashed",
        Posterior = "solid"
      ),
      breaks = c("Prior", "Posterior"),
      name = NULL,
      guide = guide_legend(
        override.aes = list(
          colour = c("grey40", "grey25"),
          linewidth = c(0.55, 0.75)
        )
      )
    ) +
    scale_x_continuous(
      limits = c(0, NA),
      labels = scales::label_number(
        scale_cut = scales::cut_short_scale()
      ),
      expand = expansion(mult = c(0, 0.03))
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(
      title = parameter_titles[[parameter]],
      x = if (position == 3) "Parameter value" else NULL,
      y = if (position == 2) "Density" else NULL
    ) +
    theme_minimal(base_size = 7) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text = element_text(colour = "black", size = 6),
      axis.title = element_text(size = 6.5),
      plot.title = element_text(
        size = 7, hjust = 0.02
      ),
      plot.margin = margin(
        0.5, 1,
        if (position == 3) 0 else 0.5,
        1,
        unit = "mm"
      )
    )
}


plots_c <- Map(
  make_density_subplot,
  posterior_parameters,
  seq_along(posterior_parameters)
)

panel_c <- (
  wrap_plots(
    plots_c,
    ncol = 1,
    guides = "collect"
  ) &
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
    plot.margin = margin(0, 0.5, 0, 0.5, unit = "mm")  )
    )

# Prevent individual density subplots from receiving tags
panel_c_titled <- panel_c +
  plot_annotation(
    title = "ABC posteriors",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 9,
        hjust = 0.5
      )
    )
  )

panel_c <- panel_c_titled &
  theme(
    legend.position = "bottom",
    legend.box.spacing = grid::unit(0, "mm"),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm")
  )


panel_c <- wrap_elements(
  full = panel_c,
  clip = FALSE
)


#####################################################
## intermediate panel with info
#####################################################

box_colours <- c(
  red    = "#F6D6D6",
  blue   = "#D7E8F6",
  violet = "#E7DDF2",
  orange = "#F9E2CC"
)

segment_titles <- data.frame(
  x = c(0.5, 1.5, 2.5),
  y = 0.76,
  label = c(
    "All populations",
    "WSS populations",
    "WSS populations"
  )
)

population_boxes <- data.frame(
  x = c(
    0.125, 0.375, 0.625, 0.875,
    1.27, 1.73,
    2.27, 2.73
  ),
  y = 0.48,
  label = c(
    "25 Khoe-San", "50 RHGn", "11 eRHG", "17 wRHG",
    "25 Khoe-San/Nama", "50 RHGn/MSL",
    "5 Nama", "3 MSL"
  ),
  fill = unname(box_colours[c(
    "red", "blue", "violet", "orange",
    "red", "blue",
    "red", "blue"
  )])
)

summary_boxes <- data.frame(
  x = c(0.5, 1.5, 2.5),
  y = 0.21,
  label = c(
    "56 summary statistics",
    "20 summary statistics",
    "20 summary statistics"
  )
)

segment_dividers <- data.frame(
  x = c(1, 2),
  xend = c(1, 2),
  y = 0.04,
  yend = 0.96
)

panel_d_information <- ggplot() +
  geom_segment(
    data = segment_dividers,
    aes(x = x, xend = xend, y = y, yend = yend),
    colour = "grey80",
    linewidth = 0.3
  ) +
  geom_text(
    data = segment_titles,
    aes(x = x, y = y, label = label),
    fontface = "bold",
    size = 2.5
  ) +
  geom_label(
    data = population_boxes,
    aes(x = x, y = y, label = label, fill = fill),
    colour = "grey20",
    linewidth = 0.15,
    label.padding = unit(0.3, "mm"),
    label.r = unit(0.7, "mm"),
    size = 2.5,
    show.legend = FALSE
  ) +
  geom_label(
    data = summary_boxes,
    aes(x = x, y = y, label = label),
    fill = "grey92",
    colour = "grey20",
    linewidth = 0.15,
    label.padding = unit(0.3, "mm"),
    label.r = unit(0.7, "mm"),
    size = 2.5
  ) +
  scale_fill_identity() +
  scale_x_continuous(
    limits = c(0, 3),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(0, 0.5, 0, 0.5, unit = "mm")  )


panel_d_information <- wrap_elements(
  full = panel_d_information,
  clip = FALSE,
  ignore_tag = TRUE
)



##############################
## final figure
##############################


empty_panel <- function() {
  ggplot() +
    theme_void() +
    theme(plot.margin = margin(1, 1, 1, 1, unit = "mm"))
}

figure_layout <- c(
  area(t = 1, l = 1, b = 2, r = 3),  # A: 3/4 width
  area(t = 1, l = 4, b = 1, r = 4),  # B
  area(t = 2, l = 4, b = 2, r = 4),  # C
  area(t = 3, l = 1, b = 3, r = 4),  # information strip
  area(t = 4, l = 1, b = 4, r = 4),  # D
  area(t = 5, l = 1, b = 5, r = 4)   # E
)

d_height <- 0.95
e_height <- 1.25  # increased from 1.05; adjust if necessary


final_figure <- (
  panel_a +
    panel_b +
    panel_c +
    panel_d_information +
    wrap_elements(full = panel_d, clip = FALSE) +
    wrap_elements(full = panel_e, clip = FALSE) +
    plot_layout(
      design = figure_layout,
      widths = rep(1, 4),
      heights = c(
        2 / 3,
        4 / 3,
        d_height / 3,
        d_height,
        e_height
      )
    )
) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 11)
    )
  )


ggsave(
  "~/demog-afr/plots/ABC_composite_figure.pdf",
  final_figure,
  width = 180,
  height = 259,
  units = "mm",
  device = grDevices::cairo_pdf
)




#################################################################################
######### SI Figures on prediction error & posteriors
#################################################################################

######### prediction error
parameters_c <- colnames(cv_abc$true)

pe_unfiltered <- do.call(
  rbind,
  lapply(parameters_c, function(p) {
    extract_prediction_error(cv_abc, p)
  })
)

pe_all <- pe_unfiltered

# Preserve the parameter names/order stored in the CV object
tol_max <- 0.21

parameter_order <- colnames(cv_abc$true)

# Filter the original prediction-error table
pe_plot <- pe_all |>
  dplyr::transmute(
    parameter = as.character(parameter),
    tolerance = as.numeric(as.character(tolerance)),
    prediction_error = as.numeric(as.character(prediction_error))
  ) |>
  dplyr::filter(
    parameter %in% parameter_order,
    is.finite(tolerance),
    tolerance < tol_max,
    is.finite(prediction_error)
  ) |>
  dplyr::arrange(parameter, tolerance)

# Explicitly repeat all curves for every focal-parameter facet
all_curves <- tidyr::expand_grid(
  focus_parameter = parameter_order,
  pe_plot
) |>
  dplyr::mutate(
    focus_parameter = factor(
      focus_parameter,
      levels = parameter_order
    ),
    parameter = factor(
      parameter,
      levels = parameter_order
    ),
    curve_id = interaction(focus_parameter, parameter, drop = TRUE)
  ) |>
  dplyr::arrange(focus_parameter, parameter, tolerance)

# Red curve in each facet
focus_curves <- all_curves |>
  dplyr::filter(
    as.character(parameter) == as.character(focus_parameter)
  )

supp_prediction_errors <- ggplot(
  all_curves,
  aes(
    x = tolerance,
    y = prediction_error,
    group = curve_id
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.3,
    colour = "grey45"
  ) +
  geom_line(
    colour = "grey78",
    linewidth = 0.35
  ) +
  geom_line(
    data = focus_curves,
    colour = "#8B0000",
    linewidth = 0.75
  ) +
  facet_wrap(
    ~focus_parameter,
    ncol = 4,
    axes = "all",
    axis.labels = "all"
  ) +
  scale_x_continuous(
    breaks = sort(unique(pe_plot$tolerance)),
    labels = function(x) {
      format(
        x,
        scientific = FALSE,
        trim = TRUE,
        drop0trailing = TRUE
      )
    },
    guide = guide_axis(n.dodge = 2)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 1.14)) +
  labs(
    title = "Prediction errors for ABC parameters",
    x = "Tolerance",
    y = "Normalized prediction error"
  ) +
  theme_pub +
  theme(
    legend.position = "none",
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.ticks = element_line(colour = "black", linewidth = 0.25),
    axis.ticks.length = grid::unit(1, "mm"),
    axis.text = element_text(size = 5.5),
    axis.title = element_text(size = 7),
    strip.text = element_text(face = "bold", size = 7),
    panel.spacing = grid::unit(1.5, "mm"),
    plot.margin = margin(2, 2, 2, 2, unit = "mm"),
    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    )
  )


ggsave(
  "~/demog-afr/plots/SI_errors.pdf",
  supp_prediction_errors,
  width = 180,
  height = 275,
  units = "mm",
  device = grDevices::cairo_pdf
)


######### posteriors
# Set this to your original full ABC parameter table:
# prior_draws <- param_table_used_for_abc

# Parameter names and order
parameter_order <- colnames(myabc$adj.values)

# Adjusted ABC posterior draws
posterior_wide <- as.data.frame(
  myabc$adj.values[, parameter_order, drop = FALSE]
)

# Full set of prior draws used for the ABC simulations
prior_wide <- as.data.frame(
  attr(myabc, "components")$parameters[, parameter_order, drop = FALSE]
)

to_long <- function(x) {
  tidyr::pivot_longer(
    x,
    cols = dplyr::all_of(parameter_order),
    names_to = "parameter",
    values_to = "value"
  ) |>
    dplyr::transmute(
      parameter = factor(parameter, levels = parameter_order),
      value = as.numeric(value)
    ) |>
    dplyr::filter(is.finite(value))
}

prior_long <- to_long(prior_wide)
posterior_long <- to_long(posterior_wide)

make_density <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(unique(x)) < 2) {
    return(tibble(value = rep(x[1], 2), density = c(0, 0)))
  }
  
  d <- stats::density(
    x,
    n = 512,
    from = min(x),
    to = max(x)
  )
  
  tibble(value = d$x, density = d$y)
}

density_table <- function(dat) {
  dat |>
    group_by(parameter) |>
    group_modify(~ make_density(.x$value)) |>
    ungroup()
}

prior_density <- density_table(prior_long)
posterior_density <- density_table(posterior_long)

posterior_summary <- posterior_long |>
  group_by(parameter) |>
  summarise(
    lower_95 = quantile(value, 0.025, names = FALSE),
    median = median(value),
    upper_95 = quantile(value, 0.975, names = FALSE),
    .groups = "drop"
  )

supp_abc_posteriors <- ggplot() +
  # 95% posterior interval
  geom_rect(
    data = posterior_summary,
    aes(
      xmin = lower_95,
      xmax = upper_95,
      ymin = -Inf,
      ymax = Inf
    ),
    inherit.aes = FALSE,
    fill = "#F6D6D6",
    alpha = 0.7,
    colour = NA
  ) +
  
  # Prior
  geom_line(
    data = prior_density,
    aes(
      x = value,
      y = density,
      group = parameter,
      colour = "Prior"
    ),
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  
  # Adjusted ABC posterior
  geom_line(
    data = posterior_density,
    aes(
      x = value,
      y = density,
      group = parameter,
      colour = "Posterior"
    ),
    linewidth = 0.6
  ) +
  
  # Posterior median
  geom_vline(
    data = posterior_summary,
    aes(xintercept = median),
    inherit.aes = FALSE,
    colour = "#8B0000",
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  
  facet_wrap(
    ~ parameter,
    ncol = 4,
    scales = "free",
    axes = "all",
    axis.labels = "all"
  ) +
  
  expand_limits(y = 0) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_colour_manual(
    values = c(
      Prior = "grey55",
      Posterior = "#C43C39"
    ),
    breaks = c("Prior", "Posterior"),
    name = NULL
  ) +
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linetype = c("dashed", "solid"),
        linewidth = c(0.5, 0.7)
      )
    )
  ) +
  labs(
    title = "ABC prior and posterior ranges for all parameters",
    x = "Parameter value",
    y = "Density"
  ) +
  theme_pub +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.ticks = element_line(colour = "black", linewidth = 0.25),
    axis.ticks.length = grid::unit(1, "mm"),
    axis.text = element_text(size = 5.5),
    axis.title = element_text(size = 7),
    strip.text = element_text(face = "bold", size = 7),
    panel.spacing = grid::unit(1.5, "mm"),
    plot.margin = margin(2, 2, 2, 2, unit = "mm"),
    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    )
  )




ggsave(
  "~/demog-afr/plots/SI_posteriors.pdf",
  supp_abc_posteriors,
  width = 180,
  height = 275,
  units = "mm",
  device = grDevices::cairo_pdf
)
