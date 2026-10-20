# =============================================================================
# launch_ui_qc(): interactive quality control of ASCAT.sc single-cell results
# -----------------------------------------------------------------------------
#   qc <- getQCs(res)
#   launch_ui_qc(qc)       # tune the thresholds, then click "Done"; qc is updated:
#   qc$pass                # ids of the cells that pass
#   qc$fail                # ids of the cells that fail
#   qc$cells               # one row per cell: status, failed metrics, values
#   qc$QC_session$thresholds  # thresholds used, in original units
#   launch_ui_qc(qc)       # re-opens with the same thresholds and decisions
#
# The results are added to the object you pass in (nothing is removed), and
# it is also returned, so qc <- launch_ui_qc(qc) works the same way.
# =============================================================================

launch_ui_qc <- function(res, theme = c("light", "dark"), previous = res$QC_session, ...)
{
  # Name of the object passed in, so it can be updated where it lives
  caller   <- parent.frame()
  res_name <- if (is.name(substitute(res))) as.character(substitute(res)) else NULL
  theme <- match.arg(theme)
  needed <- c("bslib", "plotly", "DT", "htmltools")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("launch_ui_qc() needs ", paste(missing, collapse = ", "), ". Install with ",
         "install.packages(c(", paste0("\"", missing, "\"", collapse = ", "), "))",
         call. = FALSE)
  if (is.null(res) || is.null(res$QC_metrics))
    stop("res$QC_metrics not found. Run res <- getQCs(res) first.", call. = FALSE)
  # previous can also be an object from an earlier launch_ui_qc() call
  if (!is.null(previous$QC_session)) previous <- previous$QC_session

  qd <- .qc_prepare(res, previous)

  # Shared across browser sessions: a page refresh keeps the gates, and the
  # final state is what launch_ui_qc() returns.
  st <- new.env(parent = emptyenv())
  st$defaults <- .qc_initial_state(qd)
  st$state    <- .qc_initial_state(qd, previous = previous)

  label <- if (is.null(res_name)) "res" else res_name
  app <- shiny::shinyApp(
    ui     = function(req) .qc_ui(qd, st$state, theme),
    server = .qc_server(qd, res, st, theme, label)
  )
  tryCatch(shiny::runApp(app, ...), interrupt = function(e) invisible(NULL))

  res <- .qc_attach(res, .qc_result(qd, st$state))
  # Update the caller's object as well, so launch_ui_qc(qc) needs no assignment
  saved <- !is.null(res_name) && exists(res_name, envir = caller, inherits = FALSE) &&
    !inherits(try(assign(res_name, res, envir = caller), silent = TRUE), "try-error")
  message(sprintf("ASCAT.sc QC: %s of %s cells pass (%s).", format(length(res$pass), big.mark = ","),
                  format(qd$n, big.mark = ","), .qc_pct(length(res$pass), qd$n)))
  message(if (saved)
    sprintf("Saved in %1$s: %1$s$pass, %1$s$fail and %1$s$cells. Run launch_ui_qc(%1$s) again to continue where you left off.",
            res_name)
  else "Keep the result with res <- launch_ui_qc(res), then use res$pass, res$fail and res$cells.")
  invisible(res)
}

# ── UI ───────────────────────────────────────────────────────────────────────
# A segmented control that is a regular Shiny radio input.
.qc_segmented <- function(id, choices, selected = choices[[1]], label = NULL)
{
  shiny::div(
    id = id, class = "shiny-input-radiogroup shiny-input-container qc-seg",
    role = "radiogroup", `aria-label` = label,
    lapply(seq_along(choices), function(i) {
      v <- choices[[i]]
      shiny::tags$label(
        shiny::tags$input(type = "radio", name = id, value = v,
                          checked = if (identical(v, selected)) NA else NULL),
        shiny::span(names(choices)[i]))
    }))
}

.qc_ui <- function(qd, init, theme)
{
  tags <- shiny::tags; div <- shiny::div; span <- shiny::span; icon <- shiny::icon
  M   <- qd$metrics
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

  grouped <- lapply(split(M, factor(M$group, levels = unique(M$group))),
                    function(g) stats::setNames(g$metric, g$label))
  pick <- function(pref, avoid = character()) {
    hit <- setdiff(intersect(pref, M$metric), avoid)
    if (length(hit)) hit[1] else setdiff(M$metric, avoid)[1]
  }
  x0 <- pick(c("n_reads_total", "MAPD_gc_corrected"))
  y0 <- pick(c("MAPD_gc_corrected", "spikiness"), x0)
  z0 <- pick(c("MSRSI_non_integerness", "breakpoints", "spikiness"), c(x0, y0))
  if (is.na(z0)) z0 <- x0
  has_baf <- "distance_integer_BAF" %in% M$metric

  field <- function(label, ...) div(class = "qc-field", span(class = "qc-label", label), ...)
  card_head <- function(title, sub = NULL, right = NULL)
    div(class = "qc-card-head",
        div(class = "qc-card-head-left", span(class = "qc-card-title", title),
            if (!is.null(sub)) span(class = "qc-card-sub", sub)),
        right)
  legend <- function(...) div(class = "qc-legend", ...)
  tab_title <- function(ic, txt) shiny::tagList(icon(ic), span(class = "qc-tab-text", txt))

  # ---- Header ----
  dl_item <- function(id, ic, txt)
    tags$li(shiny::downloadLink(id, shiny::tagList(icon(ic), txt), class = "dropdown-item"))
  header <- tags$header(
    class = "qc-header",
    div(class = "qc-brand", .qc_logo(),
        div(class = "qc-title", "ASCAT.sc ", span("Quality control"))),
    div(class = "qc-pills",
        span(class = "qc-pill", fmt(qd$n), " cells"),
        span(class = "qc-pill", nrow(M), " metrics"),
        if (is.finite(qd$n_bins)) span(class = "qc-pill", fmt(qd$n_bins), " bins"),
        if (has_baf) span(class = "qc-pill", "Allele-specific")),
    div(class = "qc-header-actions",
        tags$button(type = "button", class = "qc-btn qc-btn-icon", `data-qc-toggle-sidebar` = NA,
                    title = "Show or hide the threshold panel", `aria-pressed` = "false",
                    icon("table-columns")),
        div(class = "dropdown",
            tags$button(type = "button", class = "qc-btn dropdown-toggle",
                        `data-bs-toggle` = "dropdown", `aria-expanded` = "false",
                        icon("arrow-up-from-bracket"), "Export"),
            tags$ul(class = "dropdown-menu dropdown-menu-end",
                    tags$li(tags$h6(class = "dropdown-header", "CSV")),
                    dl_item("dl_pass", "circle-check", "Passing cells"),
                    dl_item("dl_fail", "circle-xmark", "Failing cells"),
                    dl_item("dl_table", "table", "All cells with status"),
                    dl_item("dl_gates", "sliders", "Thresholds"))),
        bslib::input_dark_mode(id = "dark_mode", mode = theme),
        shiny::actionButton("done", "Done", icon = icon("check"), class = "qc-btn qc-btn-primary",
                            title = "Close the app and return the result to R")))

  # ---- Sidebar: gates ----
  gate_ui <- function(m) {
    sp <- M[m, ]
    div(class = "qc-gate", `data-metric` = m,
        div(class = "qc-gate-head",
            bslib::input_switch(paste0("on_", m), NULL, value = isTRUE(init$on[[m]])),
            span(class = "qc-gate-name", sp$label),
            if (!is.na(sp$help))
              bslib::tooltip(span(class = "qc-help", tabindex = "0", icon("circle-info")),
                             paste0(sp$help,
                                    if (sp$n_na) sprintf(" %d cells have no value and always pass.", sp$n_na)),
                             placement = "right")),
        div(class = "qc-gate-body",
            shiny::uiOutput(paste0("gs_", m)),
            shiny::sliderInput(paste0("gate_", m), NULL, min = sp$min, max = sp$max,
                               value = init$bounds[[m]], step = sp$step, ticks = FALSE,
                               width = "100%")))
  }
  groups <- unique(M$group)
  gate_sections <- lapply(groups, function(g) shiny::tagList(
    div(class = "qc-label qc-group-label", g),
    lapply(M$metric[M$group == g], gate_ui)))

  flag_ui <- if (nrow(qd$flags)) div(
    class = "qc-side-section",
    div(class = "qc-label qc-group-label", paste0("Flags ", .qc_s[["mid"]], " exclude cells that are")),
    lapply(qd$flags$flag, function(f) {
      fr <- qd$flags[f, ]
      div(class = "qc-flag",
          bslib::input_switch(paste0("flag_", f), NULL, value = isTRUE(init$flags_on[[f]])),
          span(fr$label),
          if (!is.na(fr$help))
            bslib::tooltip(span(class = "qc-help", tabindex = "0", icon("circle-info")), fr$help,
                           placement = "right"),
          span(class = "qc-flag-n", fmt(fr$n), " cells"))
    }))

  sidebar <- tags$aside(
    class = "qc-sidebar", `aria-label` = "Thresholds",
    div(class = "qc-side-section",
        div(class = "qc-label", "Thresholds"),
        div(class = "qc-auto",
            span(paste0("median ", .qc_s[["pm"]])),
            shiny::numericInput("k", NULL, value = init$k, min = 0.5, max = 10, step = 0.5,
                                width = "56px"),
            span(paste0(.qc_s[["times"]], " MAD")),
            div(class = "qc-auto-actions",
                shiny::actionButton("apply_k", "Apply", class = "qc-btn qc-btn-sm",
                                    title = "Recompute every threshold in use from the median and MAD")))),
    gate_sections,
    flag_ui,
    div(class = "qc-side-foot",
        shiny::actionButton("reset", "Reset thresholds", icon = icon("rotate-left"),
                            class = "qc-btn qc-btn-sm")),
    div(class = "qc-kbd-hints",
        div(tags$kbd(.qc_s[["larr"]]), " ", tags$kbd(.qc_s[["rarr"]]), " previous / next cell"),
        div(tags$kbd("lasso"), " select cells in Explore")))

  # ---- Tabs ----
  overview <- shiny::tagList(
    shiny::uiOutput("kpis"),
    div(class = "qc-grid-2",
        div(class = "qc-card",
            card_head("Pass and fail per metric", "cells within / outside each threshold"),
            shiny::uiOutput("impact")),
        div(class = "qc-card",
            card_head("Metric correlations", paste0("Spearman ", .qc_s[["rho"]], ", strength ", .qc_s[["mid"]], " click a square to plot the pair")),
            plotly::plotlyOutput("corr", height = "340px"))),
    div(class = "qc-card", style = "margin-top:16px",
        card_head("Distributions", "all cells, stacked by overall status",
                  legend(span(tags$i(class = "pass"), "Pass"), span(tags$i(class = "fail"), "Fail"),
                         span(tags$i(class = "line"), "threshold"))),
        div(class = "qc-dists", lapply(M$metric, function(m)
          div(class = "qc-dist",
              div(class = "qc-dist-head", span(class = "qc-dist-name", M[m, "label"]),
                  shiny::textOutput(paste0("distgate_", m), inline = TRUE,
                                    container = function(...) span(class = "qc-dist-gate", ...))),
              plotly::plotlyOutput(paste0("hist_", m), height = "128px"))))))

  axis_select <- function(id, sel)
    div(class = "qc-select-sm", shiny::selectInput(id, NULL, grouped, selected = sel, width = "190px"))
  explore <- shiny::tagList(
    div(class = "qc-toolbar",
        field("X", axis_select("xvar", x0)),
        field("Y", axis_select("yvar", y0)),
        shiny::conditionalPanel("input.dim === '3d'", field("Z", axis_select("zvar", z0))),
        field("Colour", div(class = "qc-select-sm",
                            shiny::selectInput("color_by", NULL,
                                               c(list(Status = c("Pass / fail" = "status")), grouped),
                                               selected = "status", width = "170px"))),
        div(class = "qc-spacer"),
        .qc_segmented("dim", c("3D" = "3d", "2D" = "2d"), "3d", "Plot dimensions")),
    shiny::uiOutput("selbar"),
    div(class = "qc-card",
        plotly::plotlyOutput("explore", height = "600px"),
        shiny::conditionalPanel(
          "input.dim === '3d'",
          div(class = "qc-hint", icon("cube"),
              paste0("Drag to rotate ", .qc_s[["mid"]], " scroll to zoom ", .qc_s[["mid"]], " right-drag to pan ", .qc_s[["mid"]], " click a cell to see its profile."),
              " The box marks the thresholds on the three axes; switch to 2D to select cells with the lasso.")),
        shiny::conditionalPanel(
          "input.dim === '2d'",
          div(class = "qc-hint", icon("hand-pointer"),
              "Click a cell to see its profile. Use the box or lasso tool to select several cells."))),
    shiny::conditionalPanel(
      "output.has_sel",
      div(class = "qc-card",
          div(class = "qc-card-head", shiny::uiOutput("explore_sel_head", class = "qc-explore-sel-head"),
              shiny::actionButton("to_inspector", "Open in Inspector", icon = icon("arrow-right"),
                                  class = "qc-btn qc-btn-sm")),
          div(class = "qc-cn-row", .qc_cn_legend()),
          plotly::plotlyOutput("explore_profile", height = "250px"))))

  cells <- shiny::tagList(
    div(class = "qc-toolbar",
        .qc_segmented("tbl_filter", c("All" = "all", "Passing" = "pass", "Failing" = "fail",
                                      "Manual" = "manual"), "all", "Filter cells"),
        div(class = "qc-spacer"),
        span(class = "qc-hint", icon("hand-pointer"), "Click a row to open the cell in the Inspector")),
    div(class = "qc-card", DT::DTOutput("tbl")))

  inspector <- shiny::tagList(
    div(class = "qc-toolbar",
        div(class = "qc-select-cell",
            shiny::selectizeInput("pick", NULL, choices = NULL, width = "260px",
                                  options = list(placeholder = paste0("Search cell", .qc_s[["ell"]])))),
        shiny::actionButton("prev", NULL, icon = icon("chevron-left"), class = "qc-btn qc-btn-icon",
                            title = paste0("Previous cell (", .qc_s[["larr"]], ")")),
        shiny::actionButton("nxt", NULL, icon = icon("chevron-right"), class = "qc-btn qc-btn-icon",
                            title = paste0("Next cell (", .qc_s[["rarr"]], ")")),
        shiny::uiOutput("insp_pos", inline = TRUE, class = "qc-pos"),
        div(class = "qc-spacer"),
        field("Browse", .qc_segmented("browse", c("All" = "all", "Failing" = "fail", "Passing" = "pass"),
                                      "all", "Cells to browse")),
        field("Decision", .qc_segmented("ovr", c("Auto" = "auto", "Keep" = "keep", "Exclude" = "exclude"),
                                        "auto", "Manual decision for this cell"))),
    div(class = "qc-card",
        shiny::uiOutput("insp_head", class = "qc-insp-head"),
        div(class = "qc-cn-row", .qc_cn_legend()),
        plotly::plotlyOutput("insp_profile", height = "340px")),
    div(class = "qc-card",
        card_head("Metrics", paste0("where this cell sits among all cells ", .qc_s[["mid"]], " the line marks this cell")),
        shiny::uiOutput("insp_tiles")))

  main <- tags$main(
    class = "qc-main",
    shiny::uiOutput("tab_status", class = "qc-tab-status"),
    bslib::navset_pill(
      id = "tab",
      bslib::nav_panel(tab_title("gauge-high", "Overview"), value = "overview", overview),
      bslib::nav_panel(tab_title("braille", "Explore"), value = "explore", explore),
      bslib::nav_panel(tab_title("table-list", "Cells"), value = "cells", cells),
      bslib::nav_panel(tab_title("dna", "Inspector"), value = "inspector", inspector)))

  bslib::page(
    title = paste0("ASCAT.sc ", .qc_s[["mid"]], " Quality control"), lang = "en", theme = .qc_bs_theme(),
    # Theme before first paint, so dark mode never flashes white.
    tags$head(tags$script(shiny::HTML(sprintf(
      "document.documentElement.setAttribute('data-bs-theme', '%s');", theme))),
      tags$link(rel = "icon", type = "image/webp", href = .qc_logo_uri("dark"))),
    # In <body>, after every stylesheet Shiny adds to <head>, so it wins ties.
    tags$style(shiny::HTML(.qc_css())),
    tags$script(shiny::HTML(.qc_js)),
    div(class = "qc-app", header, div(class = "qc-body", sidebar, main)))
}

# ── Server ───────────────────────────────────────────────────────────────────
.qc_server <- function(qd, res, st, theme, label = "res")
{
  M  <- qd$metrics
  ms <- M$metric
  fl <- qd$flags$flag
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  tags <- shiny::tags; div <- shiny::div; span <- shiny::span; icon <- shiny::icon

  # Profiles are extracted on demand and cached per cell. If extraction fails
  # (e.g. an unusual res layout) the reason is shown instead of the plot.
  prof_cache <- new.env(parent = emptyenv())
  prof_err <- new.env(parent = emptyenv())
  get_profile <- function(id) {
    if (!qd$has_profiles) return(NULL)
    if (exists(id, envir = prof_cache, inherits = FALSE)) return(get(id, envir = prof_cache))
    j <- qd$prof_idx[match(id, qd$ids)]
    sex <- if (!is.na(j) && length(res$sex) >= j) res$sex[j] else NULL
    pd <- tryCatch(.qc_profile_data(res, j, sex), error = function(e) {
      assign(id, conditionMessage(e), envir = prof_err); NULL })
    if (length(ls(prof_cache)) > 80) rm(list = ls(prof_cache), envir = prof_cache)
    assign(id, pd, envir = prof_cache)
    pd
  }
  no_profile_msg <- function(id) {
    err <- get0(id, envir = prof_err, inherits = FALSE)
    if (is.null(err)) "No copy-number profile is available for this cell."
    else paste("Could not build this cell's profile:", err)
  }
  status_pill <- function(s) {
    lab <- c(pass = "Pass", fail = "Fail", kept = paste0("Kept ", .qc_s[["mid"]], " manual"),
             excluded = paste0("Excluded ", .qc_s[["mid"]], " manual"))[[s]]
    ic <- if (s %in% c("pass", "kept")) "circle-check" else "circle-xmark"
    span(class = paste("qc-status", s), icon(ic), lab)
  }
  unescape <- function(x) {
    for (p in list(c("&lt;", "<"), c("&gt;", ">"), c("&quot;", "\""), c("&#39;", "'"), c("&amp;", "&")))
      x <- gsub(p[1], p[2], x, fixed = TRUE)
    x
  }

  ev <- function(event, source) plotly::event_data(event, source = source)

  function(input, output, session)
  {
    # The Explore plot is not drawn until its tab is opened, and plotly warns
    # (after each flush) about events of plots it has not seen yet. Register
    # them up front; the IDs follow plotly's "<event>-<source>" convention.
    session$userData$plotlyShinyEventIDs <- unique(c(
      session$userData$plotlyShinyEventIDs,
      "plotly_click-explore", "plotly_selected-explore", "plotly_deselect-explore",
      "plotly_click-corr"))

    pal <- shiny::reactive(.qc_plot_palette(
      if (input$dark_mode %in% c("dark", "light")) input$dark_mode else theme))

    # ---- Gate state --------------------------------------------------------
    gate_b <- function(m) {
      v <- input[[paste0("gate_", m)]]
      if (length(v) == 2) v else st$state$bounds[[m]]
    }
    gate_on <- function(m) {
      v <- input[[paste0("on_", m)]]
      if (is.null(v)) isTRUE(st$state$on[[m]]) else isTRUE(v)
    }
    bounds   <- shiny::reactive(stats::setNames(lapply(ms, gate_b), ms))
    on       <- shiny::reactive(stats::setNames(vapply(ms, gate_on, logical(1)), ms))
    flags_on <- shiny::reactive(stats::setNames(vapply(fl, function(f) {
      v <- input[[paste0("flag_", f)]]
      if (is.null(v)) isTRUE(st$state$flags_on[[f]]) else isTRUE(v)
    }, logical(1)), fl))
    overrides <- shiny::reactiveVal(st$state$overrides)
    k_val <- shiny::reactive({
      k <- suppressWarnings(as.numeric(input$k))
      if (length(k) != 1 || !is.finite(k) || k <= 0) st$state$k else k
    })

    cls <- shiny::reactive(.qc_classify(qd, bounds(), on(), flags_on(), overrides()))
    # Plots redraw once dragging pauses
    gates_d <- shiny::debounce(shiny::reactive(
      list(cls = cls(), bounds = bounds(), on = on())), 150)

    state_now <- shiny::reactive(list(bounds = bounds(), on = on(), flags_on = flags_on(),
                                      overrides = overrides(), k = k_val()))
    shiny::observe(st$state <- state_now())

    set_override <- function(ids, value) {
      ov <- overrides()
      if (value == "auto") ov <- ov[!names(ov) %in% ids] else ov[ids] <- value
      overrides(ov)
    }

    # ---- Sidebar gates -----------------------------------------------------
    lapply(ms, function(m) {
      sp <- M[m, ]
      med <- .qc_fmt(if (sp$log) 10^sp$median else sp$median)
      output[[paste0("gs_", m)]] <- shiny::renderUI({
        b <- gate_b(m); is_on <- gate_on(m)
        nfail <- if (is_on) sum(.qc_gate_fail(qd$S[, m], sp, b)) else 0L
        shiny::tagList(
          span(class = paste("qc-gate-badge", if (nfail == 0) "is-zero"),
               title = "Cells outside this threshold", paste(fmt(nfail), "fail")),
          div(class = "qc-gate-read",
              span(class = "qc-gate-limit", .qc_gate_text(sp, b)),
              span(class = "qc-gate-med", paste("median", med))),
          .qc_spark_svg(qd$hist[[m]], b[1], b[2], on = is_on))
      })
    })

    shiny::observeEvent(input$apply_k, {
      k <- k_val()
      for (m in ms) if (gate_on(m))
        shiny::updateSliderInput(session, paste0("gate_", m), value = .qc_default_bounds(M[m, ], k))
      shiny::showNotification(sprintf(paste0("Thresholds in use set to median ", .qc_s[["pm"]], " %s ", .qc_s[["times"]], " MAD."), k), duration = 3)
    })
    shiny::observeEvent(input$reset, {
      d <- st$defaults
      for (m in ms) {
        shiny::updateSliderInput(session, paste0("gate_", m), value = d$bounds[[m]])
        bslib::update_switch(paste0("on_", m), value = isTRUE(d$on[[m]]))
      }
      for (f in fl) bslib::update_switch(paste0("flag_", f), value = FALSE)
      shiny::updateNumericInput(session, "k", value = d$k)
      shiny::showNotification("Thresholds reset to their defaults. Manual decisions were kept.", duration = 3)
    })

    # ---- Summary ------------------------------------------------------------
    output$tab_status <- shiny::renderUI({
      np <- sum(cls()$pass)
      shiny::tagList(span(class = "qc-dot pass"),
                     span(tags$b(fmt(np)), " of ", fmt(qd$n), " cells pass"),
                     span(class = "qc-faint", .qc_pct(np, qd$n)))
    })

    output$kpis <- shiny::renderUI({
      cl <- cls(); n <- qd$n; np <- sum(cl$pass); nf <- n - np
      ov <- overrides(); nk <- sum(ov == "keep"); ne <- sum(ov == "exclude")
      kpi <- function(ic, label, value, small = NULL, sub = NULL, meter = NULL)
        div(class = "qc-kpi",
            div(class = "qc-kpi-label", ic, label),
            div(class = "qc-kpi-value", value, if (!is.null(small)) tags$small(small)),
            if (!is.null(sub)) div(class = "qc-kpi-sub", sub),
            if (!is.null(meter)) div(class = "qc-meter", div(style = sprintf("width:%.1f%%", 100 * meter))))
      manual <- if (nk + ne == 0) "No manual decisions" else shiny::tagList(
        sprintf(paste0("Manual: %d kept ", .qc_s[["mid"]], " %d excluded ", .qc_s[["mid"]], " "), nk, ne),
        shiny::actionLink("clear_ovr", "clear"))
      div(class = "qc-kpis",
          kpi(icon("layer-group"), "Cells", fmt(n),
              sub = if (is.finite(qd$n_bins)) paste(fmt(qd$n_bins), "bins per cell") else "in res$QC_metrics"),
          kpi(icon("circle-check", class = "icon-pass"), "Passing", fmt(np), .qc_pct(np, n),
              meter = if (n) np / n else 0),
          kpi(icon("circle-xmark", class = "icon-fail"), "Failing", fmt(nf), .qc_pct(nf, n),
              sub = sprintf("%s fail a metric or flag", fmt(sum(!cl$gate_pass)))),
          kpi(icon("sliders"), "Metrics in use", sum(on()), sprintf("of %d", length(ms)),
              sub = manual))
    })
    shiny::observeEvent(input$clear_ovr, {
      ov <- overrides()
      shiny::showModal(shiny::modalDialog(
        title = "Clear manual decisions?",
        sprintf("This removes %d keep and %d exclude decisions. Thresholds are not changed.",
                sum(ov == "keep"), sum(ov == "exclude")),
        easyClose = TRUE,
        footer = shiny::tagList(shiny::modalButton("Cancel"),
                                shiny::actionButton("clear_ovr_ok", "Clear", class = "qc-btn qc-btn-primary"))))
    })
    shiny::observeEvent(input$clear_ovr_ok, { overrides(stats::setNames(character(), character())); shiny::removeModal() })

    # Pass / fail per metric (and flag) in use; the last row combines them all.
    output$impact <- shiny::renderUI({
      cl <- cls(); o <- on(); fo <- flags_on()
      F <- cbind(cl$fails[, o, drop = FALSE], cl$flag_fails[, fo, drop = FALSE])
      if (ncol(F) == 0)
        return(div(class = "qc-empty", icon("sliders"), "No metric is in use. Switch one on in the panel."))
      n <- qd$n
      labs <- c(M$label[o], qd$flags$label[fo])
      n_fail <- colSums(F)
      pct <- function(x) .qc_pct(x, n)
      row <- function(label, nf, total = FALSE) {
        np <- n - nf
        tip <- sprintf("%s: %s pass (%s), %s fail (%s)", label, fmt(np), pct(np), fmt(nf), pct(nf))
        extra <- if (total) "is-total"
        shiny::tagList(
          div(class = paste("qc-pf-label", extra), title = tip, label),
          div(class = paste("qc-pf-bar", extra), title = tip,
              if (np > 0) div(class = "qc-pf-pass", style = sprintf("flex: %d 1 0", np)),
              if (nf > 0) div(class = "qc-pf-fail", style = sprintf("flex: %d 1 0", nf))),
          div(class = paste("qc-pf-n", extra), title = tip, fmt(np)),
          div(class = paste("qc-pf-n", extra, if (nf == 0) "is-zero"), title = tip, fmt(nf)))
      }
      shiny::tagList(
        div(class = "qc-pf",
            div(class = "qc-pf-h", "Metric"),
            div(class = "qc-pf-h"),
            div(class = "qc-pf-h qc-pf-hn", span(class = "qc-dot pass"), "Passing"),
            div(class = "qc-pf-h qc-pf-hn", span(class = "qc-dot fail"), "Failing"),
            lapply(order(-n_fail), function(j) row(labs[j], n_fail[j])),
            div(class = "qc-pf-sep"),
            row("All metrics", sum(rowSums(F) > 0), total = TRUE)),
        div(class = "qc-hint", style = "margin-top:12px",
            "All metrics: a cell passes only if it is within every threshold.",
            if (length(overrides())) " Manual keep / exclude decisions are not counted here."))
    })

    output$corr <- plotly::renderPlotly(.qc_plot_corr(qd, pal()))
    shiny::observeEvent(ev("plotly_click", "corr"), {
      d <- ev("plotly_click", "corr")
      xm <- M$metric[match(d$x, M$label)]; ym <- M$metric[match(d$y, M$label)]
      if (length(xm) != 1 || is.na(xm) || is.na(ym) || xm == ym) return()
      shiny::updateSelectInput(session, "xvar", selected = xm)
      shiny::updateSelectInput(session, "yvar", selected = ym)
      bslib::nav_select("tab", "explore")
    })

    lapply(ms, function(m) {
      output[[paste0("hist_", m)]] <- plotly::renderPlotly({
        g <- gates_d()
        .qc_plot_hist(qd, m, g$bounds[[m]], g$on[[m]], g$cls$pass, pal())
      })
      output[[paste0("distgate_", m)]] <- shiny::renderText(
        if (gate_on(m)) .qc_gate_text(M[m, ], gate_b(m)) else "not used")
    })

    # ---- Explore ------------------------------------------------------------
    sel <- shiny::reactiveVal(NULL)
    lasso <- shiny::reactiveVal(character())

    output$explore <- plotly::renderPlotly({
      g <- gates_d()
      shiny::req(input$xvar, input$yvar)
      dim <- if (identical(input$dim, "2d")) "2d" else "3d"
      if (dim == "3d") shiny::req(input$zvar)
      .qc_plot_scatter(qd, input$xvar, input$yvar, input$zvar, g$cls$status, pal(),
                       g$bounds, g$on, selected = sel(),
                       color_by = if (isTRUE(input$color_by %in% c("status", ms))) input$color_by else "status",
                       dim = dim)
    })
    shiny::observeEvent(ev("plotly_click", "explore"), {
      id <- ev("plotly_click", "explore")$customdata
      id <- if (length(id)) as.character(unlist(id))[1] else NULL
      if (!is.null(id) && id %in% qd$ids) sel(id)
    })
    shiny::observeEvent(ev("plotly_selected", "explore"), {
      d <- ev("plotly_selected", "explore")
      ids <- if (is.null(d) || !"customdata" %in% names(d)) character()
             else intersect(unique(as.character(unlist(d$customdata))), qd$ids)
      lasso(ids)
    }, ignoreNULL = FALSE)
    shiny::observeEvent(ev("plotly_deselect", "explore"), lasso(character()))

    output$selbar <- shiny::renderUI({
      ids <- lasso()
      if (!length(ids)) return(NULL)
      s <- cls()$status[ids]
      div(class = "qc-selbar",
          icon("object-group"), span(tags$b(fmt(length(ids))), " cells selected"),
          span(class = "qc-muted", sprintf(paste0("%s pass ", .qc_s[["mid"]], " %s fail"), fmt(sum(s %in% c("pass", "kept"))),
                                           fmt(sum(s %in% c("fail", "excluded"))))),
          div(class = "qc-spacer"),
          shiny::actionButton("sel_keep", "Keep", icon = icon("circle-check"), class = "qc-btn qc-btn-sm"),
          shiny::actionButton("sel_excl", "Exclude", icon = icon("circle-xmark"), class = "qc-btn qc-btn-sm qc-btn-danger"),
          shiny::actionButton("sel_auto", "Back to auto", icon = icon("rotate-left"), class = "qc-btn qc-btn-sm"),
          shiny::actionButton("sel_clear", NULL, icon = icon("xmark"), class = "qc-btn qc-btn-sm qc-btn-ghost",
                              title = "Clear the selection"))
    })
    shiny::observeEvent(input$sel_keep, set_override(lasso(), "keep"))
    shiny::observeEvent(input$sel_excl, set_override(lasso(), "exclude"))
    shiny::observeEvent(input$sel_auto, set_override(lasso(), "auto"))
    shiny::observeEvent(input$sel_clear, lasso(character()))

    output$has_sel <- shiny::reactive(!is.null(sel()))
    shiny::outputOptions(output, "has_sel", suspendWhenHidden = FALSE)

    cell_summary <- function(id) {
      i <- match(id, qd$ids); cl <- cls()
      fails <- .qc_failed_labels(list(fails = cl$fails[i, , drop = FALSE],
                                      flag_fails = cl$flag_fails[i, , drop = FALSE]), list(
                                        metrics = M, flags = qd$flags, n = 1))
      meta <- list()
      if (!is.null(qd$df$ploidy)) meta <- c(meta, list("ploidy ", tags$b(.qc_fmt(qd$df$ploidy[i]))))
      if (!is.null(qd$df$purity)) meta <- c(meta, list(paste0(" ", .qc_s[["mid"]], " purity "), tags$b(.qc_fmt(qd$df$purity[i]))))
      if (!is.null(qd$df$n_reads_total))
        meta <- c(meta, list(paste0(" ", .qc_s[["mid"]], " "), tags$b(.qc_fmt(qd$df$n_reads_total[i])), " reads"))
      list(id = id, status = cl$status[[i]], fails = fails, meta = meta)
    }

    output$explore_sel_head <- shiny::renderUI({
      shiny::req(sel())
      s <- cell_summary(sel())
      shiny::tagList(span(class = "qc-cell-id", s$id), status_pill(s$status),
                     span(class = "qc-cell-meta", s$meta),
                     if (nzchar(s$fails)) span(class = "qc-cell-fails", "Fails: ", tags$b(s$fails)))
    })
    output$explore_profile <- plotly::renderPlotly({
      shiny::req(sel())
      shiny::validate(shiny::need(qd$has_profiles,
        "Copy-number profiles need res$allTracks.processed, res$allProfiles and res$allSolutions."))
      pd <- get_profile(sel())
      shiny::validate(shiny::need(!is.null(pd), no_profile_msg(sel())))
      .qc_plot_profile(pd, pal(), compact = TRUE)
    })
    shiny::observeEvent(input$to_inspector, bslib::nav_select("tab", "inspector"))

    # ---- Cells table -------------------------------------------------------
    tbl_data <- shiny::reactive({
      cl <- cls(); ov <- overrides()
      keep <- switch(if (is.null(input$tbl_filter)) "all" else input$tbl_filter,
                     pass = cl$pass, fail = !cl$pass, manual = qd$ids %in% names(ov),
                     rep(TRUE, qd$n))
      lab <- c(pass = "Pass", fail = "Fail", kept = "Kept", excluded = "Excluded")
      # First two failed gates, "+N" for the rest, full list on hover
      fails <- strsplit(.qc_failed_labels(cl, qd), ", ", fixed = TRUE)
      fails_html <- vapply(fails, function(f) {
        if (!length(f)) return("")
        txt <- htmltools::htmlEscape(paste(utils::head(f, 2), collapse = ", "))
        if (length(f) > 2) txt <- sprintf('%s <span class="qc-more">+%d</span>', txt, length(f) - 2)
        sprintf('<span class="qc-fails" title="%s">%s</span>',
                htmltools::htmlEscape(paste(f, collapse = ", "), attribute = TRUE), txt)
      }, character(1))
      d <- data.frame(Cell = qd$ids,
                      Status = sprintf('<span class="qc-status %s">%s</span>', cl$status, lab[cl$status]),
                      `Failed metrics` = fails_html,
                      check.names = FALSE, stringsAsFactors = FALSE)
      for (m in ms) d[[M[m, "label"]]] <- qd$df[[m]]
      d[keep, , drop = FALSE]
    })
    output$tbl <- DT::renderDT({
      d <- shiny::isolate(tbl_data())
      dt <- DT::datatable(
        d, rownames = FALSE, escape = c(-2, -3), selection = "none", style = "bootstrap5",
        class = "table table-sm", width = "100%",
        options = list(pageLength = 25, lengthMenu = c(25, 50, 100, 500), scrollX = TRUE,
                       dom = "<'d-flex justify-content-between align-items-center mb-2'fl>t<'d-flex justify-content-between align-items-center mt-2'ip>",
                       order = list(), autoWidth = FALSE,
                       language = list(search = "", searchPlaceholder = paste0("Search cells", .qc_s[["ell"]]),
                                       lengthMenu = "_MENU_ per page",
                                       info = paste0("_START_", .qc_s[["ndash"]], "_END_ of _TOTAL_"),
                                       infoFiltered = "(of _MAX_)", zeroRecords = "No cells match")),
        callback = DT::JS(
          "table.on('click', 'tbody tr', function () {",
          "  var d = table.row(this).data();",
          "  if (d) Shiny.setInputValue('tbl_click', d[0], {priority: 'event'});",
          "});"))
      lab_int <- M$label[M$integer]
      lab_big <- M$label[M$log]
      lab_dec <- setdiff(M$label, c(lab_int, lab_big))
      if (length(lab_dec)) dt <- DT::formatSignif(dt, lab_dec, digits = 3)
      if (length(c(lab_int, lab_big))) dt <- DT::formatRound(dt, c(lab_int, lab_big), digits = 0)
      dt
    }, server = TRUE)
    tbl_proxy <- DT::dataTableProxy("tbl")
    shiny::observeEvent(tbl_data(), DT::replaceData(tbl_proxy, tbl_data(), resetPaging = FALSE,
                                                    rownames = FALSE), ignoreInit = TRUE)
    shiny::observeEvent(input$tbl_click, {
      id <- input$tbl_click
      if (!id %in% qd$ids) id <- unescape(id)
      if (id %in% qd$ids) { sel(id); bslib::nav_select("tab", "inspector") }
    })

    # ---- Inspector ---------------------------------------------------------
    browse_ids <- shiny::reactive({
      p <- cls()$pass
      switch(if (is.null(input$browse)) "all" else input$browse,
             fail = qd$ids[!p], pass = qd$ids[p], qd$ids)
    })
    insp_id <- shiny::reactive({
      id <- sel()
      if (!is.null(id)) return(id)
      b <- browse_ids()
      if (length(b)) b[1] else qd$ids[1]
    })
    go <- function(dir) {
      b <- browse_ids()
      if (!length(b)) return()
      i <- match(insp_id(), b)
      j <- if (is.na(i)) (if (dir > 0) 1 else length(b)) else max(1, min(length(b), i + dir))
      sel(b[j])
    }
    shiny::observeEvent(input$prev, go(-1))
    shiny::observeEvent(input$nxt, go(1))
    shiny::observeEvent(input$qc_nav, if (input$tab %in% c("inspector", "explore")) go(input$qc_nav))
    shiny::observeEvent(input$browse, {
      b <- browse_ids()
      if (length(b) && (is.null(sel()) || !sel() %in% b)) sel(b[1])
    }, ignoreInit = TRUE)
    # The cell shown in the Inspector is the selected cell everywhere else too
    shiny::observeEvent(input$tab, {
      if (identical(input$tab, "inspector") && is.null(sel())) sel(insp_id())
    })

    shiny::updateSelectizeInput(session, "pick", choices = qd$ids, server = TRUE,
                                selected = shiny::isolate(insp_id()))
    shiny::observeEvent(input$pick, {
      if (length(input$pick) && input$pick %in% qd$ids && !identical(input$pick, insp_id()))
        sel(input$pick)
    }, ignoreInit = TRUE)
    shiny::observeEvent(insp_id(), {
      if (!identical(input$pick, insp_id()))
        shiny::updateSelectizeInput(session, "pick", choices = qd$ids, server = TRUE,
                                    selected = insp_id())
    }, ignoreInit = TRUE)

    # Decision control mirrors the inspected cell's override
    shiny::observeEvent(list(insp_id(), overrides()), {
      cur <- overrides()[insp_id()]
      cur <- if (is.na(cur)) "auto" else unname(cur)
      if (!identical(input$ovr, cur)) shiny::updateRadioButtons(session, "ovr", selected = cur)
    })
    shiny::observeEvent(input$ovr, {
      cur <- overrides()[insp_id()]
      cur <- if (is.na(cur)) "auto" else unname(cur)
      if (!identical(cur, input$ovr)) set_override(insp_id(), input$ovr)
    }, ignoreInit = TRUE)

    output$insp_pos <- shiny::renderUI({
      b <- browse_ids(); i <- match(insp_id(), b)
      what <- switch(if (is.null(input$browse)) "all" else input$browse,
                     fail = "failing", pass = "passing", "cells")
      if (!length(b)) return(span(sprintf("no %s cells", what)))
      span(if (is.na(i)) sprintf("not in the %s list (%s)", what, fmt(length(b)))
           else sprintf("%s of %s %s", fmt(i), fmt(length(b)), what))
    })

    output$insp_head <- shiny::renderUI({
      s <- cell_summary(insp_id())
      pd <- get_profile(insp_id())
      shiny::tagList(
        div(class = "qc-cell-head",
            span(class = "qc-cell-id", s$id), status_pill(s$status),
            span(class = "qc-cell-meta", s$meta),
            if (!is.null(pd)) span(class = "qc-cell-meta qc-faint",
                                   sprintf(paste0(.qc_s[["mid"]], " %s segments"), fmt(nrow(pd$segs))))),
        div(class = "qc-cell-fails",
            if (nzchar(s$fails)) shiny::tagList("Fails: ", tags$b(s$fails))
            else if (s$status == "excluded") "Passes every metric, excluded manually."
            else "Passes every metric in use."))
    })
    output$insp_profile <- plotly::renderPlotly({
      shiny::validate(shiny::need(qd$has_profiles,
        "Copy-number profiles need res$allTracks.processed, res$allProfiles and res$allSolutions."))
      pd <- get_profile(insp_id())
      shiny::validate(shiny::need(!is.null(pd), no_profile_msg(insp_id())))
      .qc_plot_profile(pd, pal())
    })

    output$insp_tiles <- shiny::renderUI({
      i <- match(insp_id(), qd$ids); cl <- cls(); b <- bounds(); o <- on()
      tiles <- lapply(ms, function(m) {
        sp <- M[m, ]; s <- qd$S[i, m]
        pct <- if (is.finite(s)) round(100 * mean(qd$S[, m] <= s, na.rm = TRUE)) else NA
        failing <- isTRUE(cl$fails[i, m])
        st_ <- if (!o[[m]]) c("", "off", "not used")
               else if (failing) c("fail", "fail", paste("outside", .qc_gate_text(sp, b[[m]])))
               else c("pass", "pass", paste("within", .qc_gate_text(sp, b[[m]])))
        div(class = paste("qc-tile", if (o[[m]] && failing) "is-fail"),
            div(class = "qc-tile-label", title = sp$label, sp$label),
            div(class = "qc-tile-row",
                span(class = "qc-tile-value", .qc_fmt(qd$df[[m]][i])),
                if (!is.na(pct)) span(class = "qc-tile-pct", title = sprintf("Higher than or equal to %d%% of cells", pct),
                                      sprintf("P%d", pct))),
            .qc_spark_svg(qd$hist[[m]], b[[m]][1], b[[m]][2], on = o[[m]], marker = s, height = 22),
            div(class = paste("qc-tile-status", st_[1]), span(class = paste("qc-dot", st_[2])), st_[3]))
      })
      flag_tiles <- lapply(fl, function(f) {
        fr <- qd$flags[f, ]; v <- qd$df[[f]][i]
        hit <- isTRUE(v %in% fr$exclude_if); act <- isTRUE(flags_on()[[f]])
        st_ <- if (!act) c("", "off", "flag off") else if (hit) c("fail", "fail", "excluded by flag")
               else c("pass", "pass", "not flagged")
        div(class = paste("qc-tile", if (act && hit) "is-fail"),
            div(class = "qc-tile-label", title = fr$label, fr$label),
            div(class = "qc-tile-row", span(class = "qc-tile-value",
                                            if (is.na(v)) .qc_s[["dash"]] else if (isTRUE(v)) "Yes" else "No")),
            div(class = paste("qc-tile-status", st_[1]), style = "margin-top:6px",
                span(class = paste("qc-dot", st_[2])), st_[3]))
      })
      div(class = "qc-tiles", tiles, flag_tiles)
    })

    # ---- Export & finish ---------------------------------------------------
    result_now <- function() .qc_result(qd, state_now())
    stamp <- function() format(Sys.Date())
    output$dl_pass <- shiny::downloadHandler(
      filename = function() sprintf("ascatsc_qc_passing_%s.csv", stamp()),
      content = function(file) {
        t <- result_now()$cells
        utils::write.csv(t[t$qc_pass, , drop = FALSE], file, row.names = FALSE)
      })
    output$dl_fail <- shiny::downloadHandler(
      filename = function() sprintf("ascatsc_qc_failing_%s.csv", stamp()),
      content = function(file) {
        t <- result_now()$cells
        utils::write.csv(t[!t$qc_pass, , drop = FALSE], file, row.names = FALSE)
      })
    output$dl_table <- shiny::downloadHandler(
      filename = function() sprintf("ascatsc_qc_all_cells_%s.csv", stamp()),
      content = function(file) utils::write.csv(result_now()$cells, file, row.names = FALSE))
    output$dl_gates <- shiny::downloadHandler(
      filename = function() sprintf("ascatsc_qc_thresholds_%s.csv", stamp()),
      content = function(file) {
        r <- result_now()
        g <- data.frame(type = "metric", r$thresholds, exclude_if = NA, mad_k = r$k)
        f <- r$flags
        if (nrow(f)) {
          f <- data.frame(type = "flag", metric = f$flag, label = f$label, enabled = f$enabled,
                          direction = NA, min = NA, max = NA, data_min = NA, data_max = NA,
                          n_fail = f$n_fail, exclude_if = f$exclude_if, mad_k = r$k)
          g <- rbind(g, f)
        }
        utils::write.csv(g, file, row.names = FALSE)
      })

    shiny::observeEvent(input$done, {
      np <- sum(cls()$pass)
      shiny::insertUI("body", "beforeEnd", immediate = TRUE, div(
        class = "qc-done", role = "status",
        div(class = "qc-done-card", .qc_logo(),
            div(class = "qc-done-title", "Quality control finished"),
            div(class = "qc-done-text",
                sprintf("%s of %s cells pass (%s). The result is saved in %s in your R session; you can close this tab.",
                        fmt(np), fmt(qd$n), .qc_pct(np, qd$n), label)),
            tags$code(paste(paste0(label, c("$pass", "$fail", "$cells")),
                            collapse = paste0(" ", .qc_s[["mid"]], " "))))))
      shiny::stopApp()
    })
  }
}

# =============================================================================
# Internal helpers. Names start with "." so that exportPattern("^[[:alpha:]]+")
# in NAMESPACE keeps them internal (and hidden from ls() when sourced).
# =============================================================================

# ── Metric & flag metadata ───────────────────────────────────────────────────
# dir:  "upper" = high values are bad (gate has a max), "lower" = low values
#       are bad (gate has a min), "both" = two-sided.
# log:    gate and plot on a log10 scale.
# on:     gate enabled by default.
# robust: slider covers median +/- 20 MAD rather than the full data range, so
#         a few extreme cells do not squeeze the useful thresholds together.
.qc_metric_meta <- function()
{
  m <- function(label, group, dir, help, on = TRUE, log = FALSE, robust = group %in% c("Noise", "Fit"))
    list(label = label, group = group, dir = dir, help = help, on = on, log = log, robust = robust)
  list(
    n_reads_total = m("Total reads", "Coverage", "lower", log = TRUE,
      "Total number of reads across all bins. Low counts give noisy copy-number profiles."),
    MAPD_gc_corrected = m("MAPD", "Noise", "upper",
      "Median absolute pairwise difference of GC-corrected logR between neighbouring bins. Measures bin-level noise."),
    spikiness = m("Spikiness", "Noise", "upper",
      "Median second-order difference of logR. Detects isolated outlier bins. Lower is smoother."),
    MBRSM_dispersion = m("Segment dispersion", "Noise", "upper",
      "Scaled MAD of within-segment residuals. Measures intra-segment noise."),
    autocorr = m("Autocorrelation", "Noise", "upper",
      "Lag-1 autocorrelation of segment residuals. High values suggest wave-like bias (GC, replication timing)."),
    sd_noise = m("Noise (getFilters)", "Noise", "upper",
      "Median absolute difference of consecutive logR bins, as computed by getFilters()."),
    distance_integer_logR = m("LogR integer distance", "Fit", "upper",
      "Weighted distance of segment copy numbers from the nearest integer. Lower is a better fit."),
    distance_integer_BAF = m("BAF integer distance", "Fit", "upper",
      "Weighted distance of segment BAFs from the closest integer allelic state. Lower is cleaner."),
    MSRSI_non_integerness = m("Non-integerness", "Fit", "upper",
      "Scaled MAD of segment logR deviations from the integer-expected logR (MSRSI). Lower is more integer-like."),
    breakpoints = m("Breakpoints", "Fit", "upper",
      "Number of copy-number state changes across the genome. Very high counts usually mean noise is being segmented."),
    state_mode = m("Modal CN state", "Solution", "both", on = FALSE,
      "Most frequent copy-number state. Typically 2 for near-diploid cells."),
    ploidy = m("Ploidy", "Solution", "both", on = FALSE,
      "Fitted ploidy from the purity/ploidy grid search (res$allSolutions)."),
    purity = m("Purity", "Solution", "both", on = FALSE,
      "Fitted purity from the grid search (res$allSolutions). Single cells are expected close to 1.")
  )
}

.qc_group_order <- c("Coverage", "Noise", "Fit", "Solution", "Other")

# Logical columns become exclusion switches. exclude_if is the value that
# fails a cell when the switch is on.
.qc_flag_meta <- function()
{
  f <- function(label, exclude_if, help)
    list(label = label, exclude_if = exclude_if, help = help)
  list(
    ambiguous_flag = f("Ambiguous fit", TRUE,
      "searchGrid() could not settle on a purity/ploidy solution."),
    ambiguous = f("Ambiguous fit", TRUE,
      "searchGrid() could not settle on a purity/ploidy solution."),
    noise_filter_out = f("Noise outlier (getFilters)", TRUE,
      "Flagged by getFilters() as being among the noisiest cells."),
    n_reads_filter_out = f("Low reads (getFilters)", TRUE,
      "Flagged by getFilters() as being among the cells with fewest reads."),
    keep_cell_filter = f("Not kept by getFilters()", FALSE,
      "getFilters() did not keep this cell (keep_cell_filter is FALSE).")
  )
}

# ── Small utilities ──────────────────────────────────────────────────────────
.qc_or <- function(a, b) if (is.null(a)) b else a

# Typographic symbols used in labels. They are built from code points, not
# written into strings: R code must be ASCII to be portable, and in a C locale
# (common on clusters) R turns such characters into the text "<U+2264>".
.qc_s <- vapply(c(le = 0x2264, ge = 0x2265, dash = 0x2014, ndash = 0x2013,
                  pm = 0x00b1, times = 0x00d7, mid = 0x00b7, ell = 0x2026,
                  larr = 0x2190, rarr = 0x2192, uarr = 0x2191, darr = 0x2193,
                  lsaquo = 0x2039, rsaquo = 0x203a, rho = 0x03c1),
                intToUtf8, character(1))

# A "nice" slider step: 1, 2 or 5 times a power of ten.
.qc_nice_step <- function(x)
{
  if (!is.finite(x) || x <= 0) return(1)
  p <- 10^floor(log10(x))
  m <- x / p
  (if (m < 1.5) 1 else if (m < 3.5) 2 else if (m < 7.5) 5 else 10) * p
}

# Snap to the step grid; the rounding removes floating-point dust.
.qc_snap <- function(v, step, fn = round)
{
  digits <- max(0, -floor(log10(step))) + 2
  round(fn(v / step + 1e-9) * step, digits)
}

.qc_sig <- function(v, digits = 3)
{
  s <- formatC(signif(v, digits), digits = digits, format = "fg", flag = "#")
  sub("\\.$", "", trimws(s))
}

# Compact, human formatting: 1.29M, 65.3k, 0.243, 12. Vectorised, since hover
# labels format every cell.
.qc_fmt <- function(x, digits = 3)
{
  x <- as.numeric(x)
  out <- rep(.qc_s[["dash"]], length(x))
  ok <- is.finite(x)
  v <- x[ok]; a <- abs(v)
  mil <- a >= 1e6
  kil <- a >= 1e4 & !mil
  int <- !mil & !kil & v == round(v)
  oth <- !(mil | kil | int)
  r <- character(length(v))
  if (any(mil)) r[mil] <- paste0(.qc_sig(v[mil] / 1e6, digits), "M")
  if (any(kil)) r[kil] <- paste0(.qc_sig(v[kil] / 1e3, digits), "k")
  if (any(int)) r[int] <- format(v[int], big.mark = ",", scientific = FALSE, trim = TRUE)
  if (any(oth)) r[oth] <- .qc_sig(v[oth], digits)
  out[ok] <- r
  out
}

.qc_pct <- function(n, total)
  if (total == 0) "0%" else paste0(formatC(100 * n / total, format = "f", digits = 1), "%")

# Histogram counts over fixed breaks; values outside are clamped into the
# first / last bin.
.qc_hist_counts <- function(s, breaks)
{
  s <- s[is.finite(s)]
  s <- pmin(pmax(s, breaks[1]), breaks[length(breaks)])
  graphics::hist(s, breaks = breaks, plot = FALSE)$counts
}

# ── Data preparation ─────────────────────────────────────────────────────────
# Returns everything the app needs, computed once:
#   df      the QC table with a cell_id column (+ ploidy/purity if available)
#   ids     cell ids
#   S       numeric matrix of gate metrics on the slider scale (log10 if log)
#   metrics one row per gate metric: label, group, dir, log, slider range ...
#   flags   one row per logical flag column
#   hist    per-metric histogram (on the slider scale) for the sparklines
#   prof_idx index of each cell in res$allProfiles (NA if not found)
.qc_prepare <- function(res, previous = NULL)
{
  qc <- as.data.frame(res$QC_metrics, stringsAsFactors = FALSE)
  n  <- nrow(qc)
  if (n == 0) stop("res$QC_metrics has no rows.")

  rn <- rownames(qc)
  ids <- if (!is.null(rn) && !identical(rn, as.character(seq_len(n)))) rn
         else if (length(res$allProfiles) == n && !is.null(names(res$allProfiles))) names(res$allProfiles)
         else paste0("cell_", seq_len(n))
  ids <- make.unique(as.character(ids))

  # Map each QC row to its entry in res$allProfiles / allSolutions
  match_idx <- function(lst) {
    if (is.null(lst)) return(rep(NA_integer_, n))
    nm <- names(lst)
    if (!is.null(nm) && all(ids %in% nm)) return(match(ids, nm))
    if (length(lst) == n) return(seq_len(n))
    rep(NA_integer_, n)
  }
  prof_idx <- match_idx(res$allProfiles)
  sol_idx  <- match_idx(res$allSolutions)

  # Fitted ploidy / purity / ambiguity from the grid search, when available
  if (!all(is.na(sol_idx))) {
    sol_field <- function(f) vapply(sol_idx, function(j) {
      s <- if (is.na(j)) NULL else res$allSolutions[[j]]
      v <- if (is.list(s)) s[[f]] else NULL
      if (length(v) == 1 && !is.na(suppressWarnings(as.numeric(v)))) as.numeric(v) else NA_real_
    }, numeric(1))
    if (!"ploidy" %in% names(qc)) qc$ploidy <- sol_field("ploidy")
    if (!"purity" %in% names(qc)) qc$purity <- sol_field("purity")
    if (!any(c("ambiguous_flag", "ambiguous") %in% names(qc))) {
      amb <- sol_field("ambiguous")
      if (!all(is.na(amb))) qc$ambiguous <- as.logical(amb)
    }
  }

  meta_m <- .qc_metric_meta()
  meta_f <- .qc_flag_meta()

  # ---- gate metrics: numeric, not all NA, at least two distinct values ----
  is_num <- vapply(qc, function(x) is.numeric(x) && !is.logical(x), logical(1))
  cand <- names(qc)[is_num]
  cand <- cand[vapply(cand, function(col) {
    x <- qc[[col]]
    length(unique(x[is.finite(x)])) >= 2
  }, logical(1))]

  prev_g <- if (is.list(previous)) .qc_or(previous[["thresholds"]], previous[["gates"]]) else NULL
  if (!is.data.frame(prev_g)) prev_g <- NULL
  specs <- lapply(cand, function(col) {
    x  <- qc[[col]]
    mm <- .qc_or(meta_m[[col]],
                 list(label = gsub("_", " ", col), group = "Other", dir = "both",
                      help = NA_character_, on = TRUE, log = FALSE, robust = TRUE))
    ok  <- is.finite(x)
    log <- isTRUE(mm$log) && all(x[ok] > 0)
    s   <- if (log) log10(x) else x
    rng <- range(s[is.finite(s)])
    med <- stats::median(s, na.rm = TRUE)
    md  <- stats::mad(s, na.rm = TRUE)
    if (isTRUE(mm$robust) && is.finite(md) && md > 0) {
      if (mm$dir %in% c("lower", "both")) rng[1] <- max(rng[1], med - 20 * md)
      if (mm$dir %in% c("upper", "both")) rng[2] <- min(rng[2], med + 20 * md)
    }
    # Keep the finite limits of a previous run inside the slider range
    if (!is.null(prev_g) && col %in% prev_g$metric) {
      lim <- unlist(prev_g[prev_g$metric == col, c("min", "max")])
      lim <- lim[is.finite(lim) & (!log | lim > 0)]
      if (length(lim)) rng <- range(c(rng, if (log) log10(lim) else lim))
    }
    is_int <- !log && all(x[ok] == round(x[ok]))
    step <- .qc_nice_step(diff(rng) / 200)
    if (is_int) step <- max(1, round(step))
    lo <- .qc_snap(rng[1], step, floor)
    hi <- .qc_snap(rng[2], step, ceiling)
    if (hi <= lo) hi <- lo + step
    data.frame(metric = col, label = mm$label, group = mm$group, dir = mm$dir,
               help = .qc_or(mm$help, NA_character_), on = isTRUE(mm$on),
               log = log, integer = is_int, min = lo, max = hi, step = step,
               median = med, mad = md,
               n_below = sum(s < lo, na.rm = TRUE), n_above = sum(s > hi, na.rm = TRUE),
               n_na = sum(!ok), stringsAsFactors = FALSE)
  })
  metrics <- do.call(rbind, specs)
  if (is.null(metrics)) stop("No numeric QC metrics with variation were found in res$QC_metrics.")
  ord <- order(match(metrics$group, .qc_group_order), match(metrics$metric, names(meta_m)),
               na.last = TRUE)
  metrics <- metrics[ord, , drop = FALSE]
  rownames(metrics) <- metrics$metric

  S <- vapply(metrics$metric, function(col) {
    x <- qc[[col]]
    if (metrics[col, "log"]) suppressWarnings(log10(x)) else as.numeric(x)
  }, numeric(n))
  S <- matrix(S, nrow = n, dimnames = list(ids, metrics$metric))

  # ---- flags: logical columns ----
  lg <- names(qc)[vapply(qc, is.logical, logical(1))]
  lg <- lg[vapply(lg, function(col) any(!is.na(qc[[col]])), logical(1))]
  flags <- do.call(rbind, lapply(lg, function(col) {
    fm <- .qc_or(meta_f[[col]], list(label = gsub("_", " ", col), exclude_if = TRUE,
                                     help = NA_character_))
    data.frame(flag = col, label = fm$label, exclude_if = fm$exclude_if,
               help = .qc_or(fm$help, NA_character_),
               n = sum(qc[[col]] %in% fm$exclude_if), stringsAsFactors = FALSE)
  }))
  if (is.null(flags))
    flags <- data.frame(flag = character(), label = character(), exclude_if = logical(),
                        help = character(), n = integer(), stringsAsFactors = FALSE)
  rownames(flags) <- flags$flag

  # ---- sparkline histograms (fixed 36 bins over the slider range; cells
  #      beyond the range pile up in the end bins) ----
  hist <- lapply(stats::setNames(metrics$metric, metrics$metric), function(col) {
    sp <- metrics[col, ]
    br <- seq(sp$min, sp$max, length.out = 37)
    cnt <- .qc_hist_counts(S[, col], br)
    list(breaks = br, counts = cnt, mids = (br[-1] + br[-length(br)]) / 2)
  })

  df <- cbind(data.frame(cell_id = ids, stringsAsFactors = FALSE), qc)
  rownames(df) <- NULL

  has_profiles <- !is.null(res$allTracks.processed) && !is.null(res$allProfiles) &&
    !is.null(res$allSolutions) && !all(is.na(prof_idx))

  list(df = df, ids = ids, n = n, S = S, metrics = metrics, flags = flags,
       hist = hist, prof_idx = prof_idx, has_profiles = has_profiles,
       n_bins = if (has_profiles) tryCatch(sum(vapply(
         res$allTracks.processed[[prof_idx[which(!is.na(prof_idx))[1]]]]$lCTS,
         nrow, integer(1))), error = function(e) NA_integer_) else NA_integer_)
}

# Default gate on the slider scale: median +/- k * MAD, applied only on the
# side(s) where the metric's direction says values are bad.
.qc_default_bounds <- function(sp, k = 3)
{
  full <- c(sp$min, sp$max)
  if (!is.finite(sp$mad) || sp$mad == 0) return(full)
  lo <- if (sp$dir %in% c("lower", "both"))
    max(sp$min, .qc_snap(sp$median - k * sp$mad, sp$step, floor)) else sp$min
  hi <- if (sp$dir %in% c("upper", "both"))
    min(sp$max, .qc_snap(sp$median + k * sp$mad, sp$step, ceiling)) else sp$max
  c(lo, hi)
}

# Initial gate state: defaults, optionally overridden by a previous result.
# Gates that are off by default (ploidy, purity, ...) start fully open:
# those metrics are often multimodal, so median +/- MAD means little.
.qc_initial_state <- function(qd, k = 3, previous = NULL)
{
  M <- qd$metrics
  bounds <- lapply(stats::setNames(M$metric, M$metric), function(m)
    if (M[m, "on"]) .qc_default_bounds(M[m, ], k) else c(M[m, "min"], M[m, "max"]))
  on <- stats::setNames(M$on, M$metric)
  flags_on <- stats::setNames(rep(FALSE, nrow(qd$flags)), qd$flags$flag)
  overrides <- stats::setNames(character(), character())

  if (!is.null(previous)) {
    g <- .qc_or(previous[["thresholds"]], previous[["gates"]])   # $gates: older sessions
    if (is.data.frame(g) && all(c("metric", "enabled", "min", "max") %in% names(g))) {
      for (r in seq_len(nrow(g))) {
        m <- g$metric[r]
        if (!m %in% M$metric) next
        to_s <- function(v) if (M[m, "log"]) log10(v) else v
        lo <- if (is.finite(g$min[r])) to_s(g$min[r]) else M[m, "min"]
        hi <- if (is.finite(g$max[r])) to_s(g$max[r]) else M[m, "max"]
        bounds[[m]] <- c(max(M[m, "min"], min(lo, M[m, "max"])),
                         min(M[m, "max"], max(hi, M[m, "min"])))
        on[m] <- isTRUE(as.logical(g$enabled[r]))
      }
    }
    f <- previous$flags
    if (is.data.frame(f) && all(c("flag", "enabled") %in% names(f))) {
      keep <- f$flag %in% names(flags_on)
      flags_on[f$flag[keep]] <- as.logical(f$enabled[keep])
    }
    o <- previous$overrides
    if (length(o) && !is.null(names(o))) {
      o <- o[names(o) %in% qd$ids & o %in% c("keep", "exclude")]
      overrides <- stats::setNames(as.character(o), names(o))
    }
    if (is.numeric(previous$k) && length(previous$k) == 1) k <- previous$k
  }
  list(bounds = bounds, on = on, flags_on = flags_on, overrides = overrides, k = k)
}

# ── Classification ───────────────────────────────────────────────────────────
# A gate pinned to the end of its slider is open on that side, so cells
# beyond the rounded slider range are never failed by accident.
.qc_gate_fail <- function(s, sp, b)
{
  eps <- sp$step / 1000
  lo <- if (b[1] <= sp$min + eps) -Inf else b[1]
  hi <- if (b[2] >= sp$max - eps)  Inf else b[2]
  !is.na(s) & (s < lo - eps | s > hi + eps)
}

.qc_classify <- function(qd, bounds, on, flags_on, overrides)
{
  M <- qd$metrics
  fails <- vapply(M$metric, function(m) {
    if (!isTRUE(on[[m]])) return(logical(qd$n))
    .qc_gate_fail(qd$S[, m], M[m, ], bounds[[m]])
  }, logical(qd$n))
  fails <- matrix(fails, nrow = qd$n, dimnames = list(qd$ids, M$metric))
  ffl <- vapply(qd$flags$flag, function(f) {
    if (!isTRUE(flags_on[[f]])) return(logical(qd$n))
    qd$df[[f]] %in% qd$flags[f, "exclude_if"]
  }, logical(qd$n))
  ffl <- matrix(ffl, nrow = qd$n, dimnames = list(qd$ids, qd$flags$flag))

  gate_pass <- rowSums(fails) == 0 & rowSums(ffl) == 0
  ov <- overrides[qd$ids]
  status <- ifelse(!is.na(ov) & ov == "exclude", "excluded",
            ifelse(!is.na(ov) & ov == "keep" & !gate_pass, "kept",
            ifelse(gate_pass, "pass", "fail")))
  names(status) <- qd$ids
  list(fails = fails, flag_fails = ffl, gate_pass = gate_pass,
       pass = stats::setNames(status %in% c("pass", "kept"), qd$ids), status = status)
}

# Labels of the gates/flags each cell fails, e.g. "MAPD, Spikiness".
.qc_failed_labels <- function(cls, qd)
{
  lab <- c(qd$metrics$label, qd$flags$label)
  all_f <- cbind(cls$fails, cls$flag_fails)
  if (ncol(all_f) == 0) return(rep("", qd$n))
  apply(all_f, 1, function(r) paste(lab[r], collapse = ", "))
}

# Gates as a data frame on the original scale (for export and return value).
.qc_gates_table <- function(qd, bounds, on)
{
  M <- qd$metrics
  to_x <- function(m, v) if (M[m, "log"]) 10^v else v
  data.frame(
    metric = M$metric, label = M$label, enabled = as.logical(on[M$metric]),
    direction = M$dir,
    min = vapply(M$metric, function(m) {
      b <- bounds[[m]][1]; if (b <= M[m, "min"]) -Inf else to_x(m, b)
    }, numeric(1)),
    max = vapply(M$metric, function(m) {
      b <- bounds[[m]][2]; if (b >= M[m, "max"]) Inf else to_x(m, b)
    }, numeric(1)),
    data_min = vapply(M$metric, function(m) min(qd$df[[m]], na.rm = TRUE), numeric(1)),
    data_max = vapply(M$metric, function(m) max(qd$df[[m]], na.rm = TRUE), numeric(1)),
    n_fail = NA_integer_, stringsAsFactors = FALSE, row.names = NULL)
}

.qc_result <- function(qd, state)
{
  cls <- .qc_classify(qd, state$bounds, state$on, state$flags_on, state$overrides)
  gates <- .qc_gates_table(qd, state$bounds, state$on)
  gates$n_fail <- as.integer(colSums(cls$fails))
  flags <- data.frame(flag = qd$flags$flag, label = qd$flags$label,
                      enabled = as.logical(state$flags_on[qd$flags$flag]),
                      exclude_if = qd$flags$exclude_if,
                      n_fail = as.integer(colSums(cls$flag_fails)),
                      stringsAsFactors = FALSE, row.names = NULL)
  ov <- state$overrides[qd$ids]
  cells <- cbind(qd$df["cell_id"],
                 data.frame(qc_pass = unname(cls$pass), status = unname(cls$status),
                            override = ifelse(is.na(ov), "", unname(ov)),
                            failed_metrics = .qc_failed_labels(cls, qd),
                            stringsAsFactors = FALSE),
                 qd$df[setdiff(names(qd$df), "cell_id")])
  rownames(cells) <- qd$ids
  list(pass = qd$ids[cls$pass], fail = qd$ids[!cls$pass], cells = cells,
       thresholds = gates, flags = flags, overrides = state$overrides, k = state$k)
}

# Store a result in res, the way getQCs() stores QC_metrics: res$pass and
# res$fail (cell ids), res$cells (per-cell table) and res$QC_session
# (thresholds, flags and manual decisions, which launch_ui_qc() resumes from).
.qc_attach <- function(res, out)
{
  res$pass  <- out$pass
  res$fail  <- out$fail
  res$cells <- out$cells
  res$QC_session <- out[c("thresholds", "flags", "overrides", "k")]
  res
}

# ── Sparkline (inline SVG, coloured by CSS so theme switches are free) ───────
.qc_spark_svg <- function(h, lo, hi, on = TRUE, marker = NULL, height = 28,
                          class = "qc-spark")
{
  W <- 240
  n  <- length(h$counts)
  bw <- W / n
  mx <- max(h$counts, 1)
  ht <- ifelse(h$counts > 0, pmax(1.5, (height - 2) * h$counts / mx), 0)
  inside <- h$mids >= lo & h$mids <= hi
  cls <- if (!on) rep("off", n) else ifelse(inside, "in", "out")
  bars <- sprintf('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" rx="1" class="%s"/>',
                  (seq_len(n) - 1) * bw + 0.6, height - ht, max(bw - 1.2, 0.8), ht, cls)
  mk <- ""
  if (!is.null(marker) && is.finite(marker)) {
    rng <- range(h$breaks)
    mx_x <- W * (min(max(marker, rng[1]), rng[2]) - rng[1]) / diff(rng)
    mk <- sprintf('<line x1="%.2f" x2="%.2f" y1="0" y2="%d" class="mk"/>', mx_x, mx_x, height)
  }
  shiny::HTML(sprintf(
    '<svg class="%s" viewBox="0 0 %d %d" preserveAspectRatio="none" aria-hidden="true">%s<line x1="0" x2="%d" y1="%d" y2="%d" class="base"/>%s</svg>',
    class, W, height, paste(bars, collapse = ""), W, height, height, mk))
}

# Readout of a gate in words: "≤ 0.243", "≥ 65.3k", "0.1 – 0.4", "no limit".
.qc_gate_text <- function(sp, b)
{
  to_x <- function(v) if (sp$log) 10^v else v
  eps <- sp$step / 1000
  lo_open <- b[1] <= sp$min + eps
  hi_open <- b[2] >= sp$max - eps
  if (lo_open && hi_open) return("no limit")
  if (lo_open) return(paste(.qc_s[["le"]], .qc_fmt(to_x(b[2]))))
  if (hi_open) return(paste(.qc_s[["ge"]], .qc_fmt(to_x(b[1]))))
  paste(.qc_fmt(to_x(b[1])), .qc_s[["ndash"]], .qc_fmt(to_x(b[2])))
}

# ── Theme ────────────────────────────────────────────────────────────────────
.qc_font <- "-apple-system, BlinkMacSystemFont, 'SF Pro Text', 'Inter', 'Segoe UI', system-ui, sans-serif"
.qc_mono <- "ui-monospace, 'SF Mono', SFMono-Regular, Menlo, Consolas, monospace"

.qc_bs_theme <- function()
{
  bslib::bs_theme(
    version = 5, preset = "bootstrap",
    primary = "#0A84FF", success = "#34C759", danger = "#FF3B30",
    warning = "#FF9500", info = "#0A84FF",
    base_font = bslib::font_collection("-apple-system", "BlinkMacSystemFont",
                                       "SF Pro Text", "Inter", "Segoe UI",
                                       "system-ui", "sans-serif"),
    code_font = bslib::font_collection("ui-monospace", "SF Mono", "SFMono-Regular",
                                       "Menlo", "Consolas", "monospace"),
    "font-size-base" = "0.875rem",
    "border-radius" = "0.6rem",
    "border-radius-lg" = "1rem",
    "border-radius-sm" = "0.4rem"
  )
}

# Copy-number colours from ASCATsc_plot.R. CN 2 and 3 are darkened on light
# backgrounds, as in the ASCAT.sc Explorer web app.
.qc_cn_colors <- function(mode)
{
  cols <- get_total_cn_colors()
  if (mode == "light") { cols["2"] <- "#888888"; cols["3"] <- "#c49631" }
  cols
}

# Plot colours per theme. Green = pass and red = fail are the only strong
# colours; gate shading, the selected cell and magnitude scales are neutral.
# Pass/fail pairs clear the colour-blind separation target on each card
# surface (dataviz validator: dE 9.4 light, 11.2 dark; target >= 8).
.qc_plot_palette <- function(mode = c("dark", "light"))
{
  mode <- match.arg(mode)
  common <- list(mode = mode, font = .qc_font)
  if (mode == "dark")
    c(common, list(
      pass = "#12A37A", fail = "#FF453A",
      surface = "#242426", text = "#f5f5f7", text2 = "rgba(235,235,245,0.78)",
      muted = "rgba(235,235,245,0.55)", grid = "#323235", axis = "#48484a",
      band = "rgba(255,255,255,0.035)", raw = "rgba(235,235,245,0.30)",
      seg = "rgba(235,235,245,0.62)", ring = "#f5f5f7",
      gate_wash = "rgba(235,235,245,0.05)",
      hover_bg = "#2c2c2e", hover_border = "rgba(255,255,255,0.12)",
      seq = c("#636366", "#f2f2f7"), mag = c("#2c2c2e", "#f2f2f7"),
      cn = .qc_cn_colors("dark")))
  else
    c(common, list(
      pass = "#058A5E", fail = "#D92D20",
      surface = "#ffffff", text = "#1d1d1f", text2 = "rgba(60,60,67,0.78)",
      muted = "rgba(60,60,67,0.6)", grid = "#ececf0", axis = "#d1d1d6",
      band = "rgba(0,0,0,0.028)", raw = "rgba(60,60,67,0.28)",
      seg = "rgba(29,29,31,0.55)", ring = "#1d1d1f",
      gate_wash = "rgba(60,60,67,0.045)",
      hover_bg = "#ffffff", hover_border = "rgba(0,0,0,0.10)",
      seq = c("#c7c7cc", "#1d1d1f"), mag = c("#f2f2f7", "#1d1d1f"),
      cn = .qc_cn_colors("light")))
}

.qc_axis <- function(pal, title = NULL, ...)
  utils::modifyList(list(
    title = if (is.null(title)) NULL else list(text = title, font = list(size = 11, color = pal$muted)),
    color = pal$muted, gridcolor = pal$grid, linecolor = pal$axis, zeroline = FALSE,
    showline = TRUE, ticks = "", tickfont = list(size = 10, color = pal$muted),
    automargin = TRUE), list(...))

.qc_layout <- function(p, pal, ...)
  plotly::layout(p,
    font = list(family = pal$font, color = pal$text, size = 11),
    paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
    hoverlabel = list(bgcolor = pal$hover_bg, bordercolor = pal$hover_border,
                      font = list(family = pal$font, color = pal$text, size = 12)),
    ...)

# Log-scale tick labels for metrics gated on log10.
.qc_log_ticks <- function(rng)
{
  e <- seq(floor(rng[1] * 2) / 2, ceiling(rng[2] * 2) / 2, by = 0.5)
  if (length(e) > 7) e <- e[e == round(e)]
  list(tickvals = e, ticktext = .qc_fmt(10^e, 2))
}

# ── Per-cell copy-number profile ─────────────────────────────────────────────
# Bins are put on the copy-number scale exactly as plotSolution() does:
# centred by normaliseByPloidy()'s mean segment value, then converted with
# transform_bulk2tumour(), so they line up with the fitted segments.
.qc_profile_data <- function(res, j, sex = NULL)
{
  if (is.na(j)) return(NULL)
  tr   <- res$allTracks.processed[[j]]
  prof <- res$allProfiles[[j]]
  sol  <- res$allSolutions[[j]]
  if (!is.list(tr) || is.null(tr$lCTS) || !is.data.frame(prof) || !is.list(sol))
    return(NULL)
  rho <- as.numeric(.qc_or(sol$purity, 1))
  psi <- as.numeric(.qc_or(sol$ploidy, 2))
  lcts <- tr$lCTS
  chrs <- .qc_or(names(lcts), as.character(res$chr))
  if (length(chrs) != length(lcts)) return(NULL)

  bins <- do.call(rbind, lapply(seq_along(lcts), function(k) {
    x <- lcts[[k]]
    data.frame(chr = chrs[k], start = as.numeric(x$start), end = as.numeric(x$end),
               smoothed = as.numeric(x$smoothed), stringsAsFactors = FALSE)
  }))
  meanS <- if (!is.null(tr$lSegs)) mean(unlist(lapply(tr$lSegs, function(s)
    rep(s$output$seg.mean, s$output$num.mark))), na.rm = TRUE) else 0
  if (!is.finite(meanS)) meanS <- 0
  ismale <- identical(sex, "male")
  isX <- gsub("chr", "", bins$chr) %in% c("X", "23")
  bins$cn <- transform_bulk2tumour(bins$smoothed - meanS, rho, psi,
                                   ismale = ismale, isX = isX,
                                   isPON = isTRUE(res$isPON))

  segs <- data.frame(chr = as.character(prof$chromosome),
                     start = as.numeric(prof$start), end = as.numeric(prof$end),
                     cn = as.numeric(prof$total_copy_number),
                     cn_raw = as.numeric(prof$total_copy_number_logr),
                     stringsAsFactors = FALSE)

  # Genome coordinates (Mb), chromosomes in track order
  chr_len <- vapply(chrs, function(cc) max(c(bins$end[bins$chr == cc],
                                             segs$end[segs$chr == cc]), na.rm = TRUE), numeric(1))
  offs <- stats::setNames(c(0, cumsum(chr_len)[-length(chr_len)]), chrs)
  bins$x <- (offs[bins$chr] + (bins$start + bins$end) / 2) / 1e6
  segs <- segs[segs$chr %in% chrs, , drop = FALSE]
  segs$x0 <- (offs[segs$chr] + segs$start) / 1e6
  segs$x1 <- (offs[segs$chr] + segs$end) / 1e6

  # Segment of each bin (for one hover layer carrying both values)
  bins$seg_cn <- NA_real_
  for (cc in chrs) {
    bi <- which(bins$chr == cc); si <- which(segs$chr == cc)
    if (!length(bi) || !length(si)) next
    k <- findInterval((bins$start[bi] + bins$end[bi]) / 2, segs$start[si])
    k[k == 0] <- 1
    bins$seg_cn[bi] <- segs$cn[si][k]
  }

  chr_tab <- data.frame(chr = gsub("chr", "", chrs), x0 = offs / 1e6,
                        x1 = (offs + chr_len) / 1e6, stringsAsFactors = FALSE)
  chr_tab$mid <- (chr_tab$x0 + chr_tab$x1) / 2
  ymax <- min(10, max(5, ceiling(max(segs$cn, na.rm = TRUE)) + 1))
  list(bins = bins, segs = segs, chr = chr_tab, purity = rho, ploidy = psi,
       ymax = ymax, sex = sex)
}

.qc_plot_profile <- function(pd, pal, compact = FALSE)
{
  b <- pd$bins
  ymax <- pd$ymax
  yclip <- pmin(pmax(b$cn, 0), ymax + 0.3)
  hover <- sprintf("<b>CN %s</b>  segment<br>%s  bin<br>chr%s:%s Mb",
                   ifelse(is.na(b$seg_cn), .qc_s[["dash"]], b$seg_cn),
                   ifelse(is.finite(b$cn), formatC(b$cn, format = "f", digits = 2), .qc_s[["dash"]]),
                   gsub("chr", "", b$chr),
                   formatC(((b$start + b$end) / 2) / 1e6, format = "f", digits = 1))

  shapes <- lapply(seq_len(nrow(pd$chr)), function(k) {
    if (k %% 2 == 1) return(NULL)
    list(type = "rect", xref = "x", yref = "paper", x0 = pd$chr$x0[k], x1 = pd$chr$x1[k],
         y0 = 0, y1 = 1, fillcolor = pal$band, line = list(width = 0), layer = "below")
  })
  shapes <- Filter(Negate(is.null), shapes)

  p <- plotly::plot_ly()
  p <- plotly::add_trace(p, type = "scattergl", mode = "markers", x = b$x, y = yclip,
                         marker = list(size = if (compact) 2.5 else 3, color = pal$raw),
                         text = hover, hovertemplate = "%{text}<extra></extra>",
                         name = "Bins", showlegend = FALSE)
  s <- pd$segs
  # Non-integer segment means: the gap to the integer bar is the non-integerness
  ok <- is.finite(s$cn_raw)
  if (any(ok)) {
    yr <- pmin(pmax(s$cn_raw[ok], 0), ymax + 0.3)
    p <- plotly::add_trace(p, type = "scatter", mode = "lines",
                           x = c(rbind(s$x0[ok], s$x1[ok], NA)), y = c(rbind(yr, yr, NA)),
                           line = list(color = pal$seg, width = 1.5), hoverinfo = "skip",
                           name = "Segment mean", showlegend = FALSE)
  }
  for (st in sort(unique(s$cn[is.finite(s$cn)]))) {
    k <- which(s$cn == st)
    key <- as.character(min(6, max(0, st)))                  # 6+ share a colour
    col <- .qc_or(unname(pal$cn[key]), pal$seg)
    yv <- min(st, ymax + 0.3)
    p <- plotly::add_trace(p, type = "scatter", mode = "lines",
                           x = c(rbind(s$x0[k], s$x1[k], NA)), y = rep(c(yv, yv, NA), length(k)),
                           line = list(color = col, width = if (compact) 4 else 5),
                           hoverinfo = "skip", name = paste("CN", st), showlegend = FALSE)
  }
  p <- .qc_layout(p, pal,
    showlegend = FALSE, hovermode = "closest", shapes = shapes,
    margin = list(l = 44, r = 10, t = 6, b = 30),
    xaxis = .qc_axis(pal, tickvals = pd$chr$mid, ticktext = pd$chr$chr, showgrid = FALSE,
                     range = c(0, max(pd$chr$x1)), tickfont = list(size = if (compact) 9 else 10,
                     color = pal$muted), fixedrange = FALSE),
    yaxis = .qc_axis(pal, "Copy number", range = c(-0.25, ymax + 0.45), dtick = 1,
                     tick0 = 0, showline = FALSE))
  plotly::config(p, displaylogo = FALSE, responsive = TRUE,
                 modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d",
                                            "hoverClosestCartesian", "hoverCompareCartesian",
                                            "toggleSpikelines"),
                 toImageButtonOptions = list(format = "png", scale = 2))
}

# ── Explore scatter (2D / 3D) ────────────────────────────────────────────────
# Metrics are plotted on the slider scale (log10 for read counts), so gate
# boxes and axes share one coordinate system; hovers show real values.
.qc_plot_scatter <- function(qd, xm, ym, zm = NULL, status, pal, bounds, on,
                             selected = NULL, color_by = "status", dim = "2d")
{
  M  <- qd$metrics
  ms <- c(xm, ym, if (dim == "3d") zm)
  raw <- function(m) qd$df[[m]]
  hover <- sprintf("<b>%s</b><br>%s", qd$ids, do.call(paste, c(lapply(ms, function(m)
    sprintf("%s  %s", .qc_fmt(raw(m)), M[m, "label"])), sep = "<br>")))
  n  <- qd$n
  sz <- if (n <= 1000) 8 else if (n <= 5000) 6 else 4.5
  # 2D opens on the slider window when extreme cells lie beyond it (they are
  # counted in an annotation; Autoscale shows them).
  clipped <- function(m) M[m, "n_above"] + M[m, "n_below"] > 0
  ax <- function(m) {
    sp <- M[m, ]
    a <- .qc_axis(pal, sp$label)
    win <- if (dim == "2d" && clipped(m)) c(sp$min, sp$max) else range(qd$S[, m], na.rm = TRUE)
    if (dim == "2d" && clipped(m)) a$range <- win + c(-1, 1) * 0.04 * diff(win)
    if (sp$log) a <- c(a, .qc_log_ticks(win))
    a
  }
  data_rng <- function(m) range(qd$S[, m], na.rm = TRUE)
  # Gate window on the slider scale; open sides extend to the data range
  gate_win <- function(m) {
    sp <- M[m, ]; r <- data_rng(m)
    if (!isTRUE(on[[m]])) return(r)
    b <- bounds[[m]]; eps <- sp$step / 1000
    c(if (b[1] <= sp$min + eps) min(r[1], b[1]) else b[1],
      if (b[2] >= sp$max - eps) max(r[2], b[2]) else b[2])
  }
  groups <- list(
    list(key = "pass",     name = "Pass",              col = pal$pass, sym = "circle"),
    list(key = "fail",     name = "Fail",              col = pal$fail, sym = "circle"),
    list(key = "kept",     name = "Kept (manual)",     col = pal$pass, sym = "diamond"),
    list(key = "excluded", name = "Excluded (manual)", col = pal$fail, sym = "diamond"))

  p <- plotly::plot_ly(source = "explore")
  type <- if (dim == "3d") "scatter3d" else "scattergl"
  add_pts <- function(p, k, ...) {
    args <- list(p, type = type, mode = "markers",
                 x = qd$S[k, xm], y = qd$S[k, ym], customdata = qd$ids[k],
                 text = hover[k], hovertemplate = "%{text}<extra></extra>", ...)
    if (dim == "3d") args$z <- qd$S[k, zm]
    do.call(plotly::add_trace, args)
  }
  if (color_by == "status") {
    for (g in groups) {
      k <- which(status == g$key)
      if (!length(k)) next
      p <- add_pts(p, k, name = sprintf("%s  %s", g$name, format(length(k), big.mark = ",")),
                   marker = list(size = if (dim == "3d") sz * 0.55 else sz, color = g$col,
                                 symbol = g$sym, opacity = if (g$key %in% c("fail", "excluded")) 0.95 else 0.7,
                                 line = list(color = pal$surface, width = if (dim == "3d") 0 else 1)))
    }
  } else {
    cv <- qd$S[, color_by]
    p <- add_pts(p, seq_len(n), name = M[color_by, "label"], showlegend = FALSE,
                 marker = list(size = if (dim == "3d") sz * 0.55 else sz, color = cv,
                               colorscale = list(c(0, pal$seq[1]), c(1, pal$seq[2])),
                               opacity = 0.85,
                               line = list(color = pal$surface, width = if (dim == "3d") 0 else 1),
                               colorbar = c(list(title = list(text = M[color_by, "label"],
                                                              font = list(size = 11, color = pal$muted)),
                                                 thickness = 10, len = 0.6, outlinewidth = 0,
                                                 tickfont = list(size = 10, color = pal$muted)),
                                            if (M[color_by, "log"])
                                              .qc_log_ticks(range(cv, na.rm = TRUE)))))
  }
  if (!is.null(selected) && selected %in% qd$ids) {
    k <- match(selected, qd$ids)
    # a neutral ring around the cell, so its own pass/fail colour still shows
    p <- add_pts(p, k, name = "Selected", showlegend = FALSE,
                 marker = if (dim == "3d")
                   list(size = 11, symbol = "circle-open", color = pal$ring, opacity = 1)
                 else
                   list(size = 17, color = "rgba(0,0,0,0)", opacity = 1,
                        line = list(color = pal$ring, width = 2.5)))
  }

  # 3D: the gates on the three axes as a see-through (neutral) box
  if (dim == "3d" && any(vapply(c(xm, ym, zm), function(m) isTRUE(on[[m]]), logical(1)))) {
    gx <- gate_win(xm); gy <- gate_win(ym); gz <- gate_win(zm)
    vx <- gx[c(1, 2, 2, 1, 1, 2, 2, 1)]
    vy <- gy[c(1, 1, 2, 2, 1, 1, 2, 2)]
    vz <- gz[c(1, 1, 1, 1, 2, 2, 2, 2)]
    p <- plotly::add_trace(p, type = "mesh3d", x = vx, y = vy, z = vz,
                           i = c(0, 0, 4, 4, 0, 0, 3, 3, 0, 0, 1, 1),
                           j = c(1, 2, 5, 6, 1, 5, 2, 6, 3, 7, 2, 6),
                           k = c(2, 3, 6, 7, 5, 4, 6, 7, 7, 4, 6, 5),
                           # facecolor, not color: plotly's R wrapper remaps `color`
                           facecolor = rep(pal$ring, 12), opacity = 0.05, flatshading = TRUE,
                           hoverinfo = "skip", showlegend = FALSE, showscale = FALSE,
                           name = "Thresholds")
    e <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1), c(5, 6), c(6, 7), c(7, 8), c(8, 5),
               c(1, 5), c(2, 6), c(3, 7), c(4, 8))
    seg <- function(v) c(t(cbind(v[e[, 1]], v[e[, 2]], NA)))
    p <- plotly::add_trace(p, type = "scatter3d", mode = "lines", x = seg(vx), y = seg(vy),
                           z = seg(vz), line = list(color = pal$muted, width = 2, dash = "dash"),
                           hoverinfo = "skip", showlegend = TRUE, name = "Thresholds")
  }

  uirev <- paste(dim, xm, ym, zm, color_by)
  leg <- list(orientation = "h", x = 0, xanchor = "left", y = 1.02, yanchor = "bottom",
              itemsizing = "constant", font = list(size = 12, color = pal$text2),
              bgcolor = "rgba(0,0,0,0)")

  if (dim == "3d") {
    sax <- function(m) {
      a <- list(title = list(text = M[m, "label"], font = list(size = 11, color = pal$muted)),
                color = pal$muted, gridcolor = pal$grid, zerolinecolor = pal$grid,
                backgroundcolor = "rgba(0,0,0,0)", showbackground = FALSE,
                tickfont = list(size = 10, color = pal$muted))
      if (M[m, "log"]) a <- c(a, .qc_log_ticks(range(qd$S[, m], na.rm = TRUE)))
      a
    }
    p <- .qc_layout(p, pal, legend = leg, uirevision = uirev,
                    margin = list(l = 0, r = 0, t = 30, b = 0),
                    scene = list(xaxis = sax(xm), yaxis = sax(ym), zaxis = sax(zm),
                                 camera = list(eye = list(x = 1.3, y = 1.3, z = 0.85))))
  } else {
    shapes <- list()
    box <- c(gate_win(xm), gate_win(ym))
    if (isTRUE(on[[xm]]) || isTRUE(on[[ym]]))
      shapes[[1]] <- list(type = "rect", xref = "x", yref = "y", x0 = box[1], x1 = box[2],
                          y0 = box[3], y1 = box[4], fillcolor = pal$gate_wash,
                          line = list(width = 0), layer = "below")
    thr <- function(m, axis) {
      if (!isTRUE(on[[m]])) return(list())
      b <- bounds[[m]]; sp <- M[m, ]; out <- list()
      for (v in b[c(b[1] > sp$min + sp$step / 1000, b[2] < sp$max - sp$step / 1000)])
        out[[length(out) + 1]] <- if (axis == "x")
          list(type = "line", xref = "x", yref = "paper", x0 = v, x1 = v, y0 = 0, y1 = 1,
               line = list(color = pal$muted, width = 1, dash = "dash"), layer = "below")
        else
          list(type = "line", xref = "paper", yref = "y", x0 = 0, x1 = 1, y0 = v, y1 = v,
               line = list(color = pal$muted, width = 1, dash = "dash"), layer = "below")
      out
    }
    shapes <- c(shapes, thr(xm, "x"), thr(ym, "y"))
    ann <- list()
    note <- function(x, y, xa, ya, txt)
      list(xref = "paper", yref = "paper", x = x, y = y, xanchor = xa, yanchor = ya,
           showarrow = FALSE, text = txt, font = list(size = 11, color = pal$muted))
    fc <- function(v) format(v, big.mark = ",")
    if (M[ym, "n_above"]) ann <- c(ann, list(note(0.01, 0.99, "left", "top",
                                                   sprintf(paste0(.qc_s[["uarr"]], " %s more cells above"), fc(M[ym, "n_above"])))))
    if (M[ym, "n_below"]) ann <- c(ann, list(note(0.01, 0.01, "left", "bottom",
                                                   sprintf(paste0(.qc_s[["darr"]], " %s more cells below"), fc(M[ym, "n_below"])))))
    if (M[xm, "n_above"]) ann <- c(ann, list(note(0.99, 0.01, "right", "bottom",
                                                   sprintf(paste0("%s more cells to the right ", .qc_s[["rarr"]]), fc(M[xm, "n_above"])))))
    if (M[xm, "n_below"]) ann <- c(ann, list(note(0.01, 0.06, "left", "bottom",
                                                   sprintf(paste0(.qc_s[["larr"]], " %s more cells to the left"), fc(M[xm, "n_below"])))))
    if (!is.null(selected) && selected %in% qd$ids) {
      k <- match(selected, qd$ids)
      if (is.finite(qd$S[k, xm]) && is.finite(qd$S[k, ym]))
        ann[[length(ann) + 1]] <- list(x = qd$S[k, xm], y = qd$S[k, ym], text = selected, showarrow = TRUE,
                         arrowhead = 0, arrowwidth = 1, arrowcolor = pal$muted, ax = 18, ay = -22,
                         font = list(size = 11, color = pal$text), bgcolor = pal$hover_bg,
                         bordercolor = pal$hover_border, borderpad = 3)
    }
    p <- .qc_layout(p, pal, legend = leg, uirevision = uirev, shapes = shapes,
                    annotations = ann, hovermode = "closest", hoverdistance = 24,
                    dragmode = "zoom", margin = list(l = 10, r = 10, t = 36, b = 10),
                    xaxis = ax(xm), yaxis = ax(ym))
  }
  p <- plotly::config(p, displaylogo = FALSE, responsive = TRUE,
                      modeBarButtonsToRemove = c("hoverClosestCartesian", "hoverCompareCartesian",
                                                 "toggleSpikelines"),
                      toImageButtonOptions = list(format = "png", scale = 2))
  p <- plotly::event_register(p, "plotly_click")
  p <- plotly::event_register(p, "plotly_selected")
  plotly::event_register(p, "plotly_deselect")
}

# ── Distribution of one metric, stacked by overall pass / fail ───────────────
.qc_plot_hist <- function(qd, m, b, on, pass, pal)
{
  sp <- qd$metrics[m, ]
  s  <- qd$S[, m]
  br <- seq(sp$min, sp$max, length.out = 41)
  cp <- .qc_hist_counts(s[pass],  br)
  cf <- .qc_hist_counts(s[!pass], br)
  mids <- (br[-1] + br[-length(br)]) / 2
  w  <- br[2] - br[1]
  to_x <- function(v) if (sp$log) 10^v else v
  rng_txt <- sprintf(paste0("%s ", .qc_s[["ndash"]], " %s"), .qc_fmt(to_x(br[-length(br)])), .qc_fmt(to_x(br[-1])))

  p <- plotly::plot_ly()
  p <- plotly::add_bars(p, x = mids, y = cp, width = w, name = "Pass",
                        marker = list(color = pal$pass, line = list(color = pal$surface, width = 1)),
                        text = rng_txt, textposition = "none", hovertemplate = "<b>%{y}</b> pass<br>%{text}<extra></extra>")
  p <- plotly::add_bars(p, x = mids, y = cf, width = w, name = "Fail",
                        marker = list(color = pal$fail, line = list(color = pal$surface, width = 1)),
                        text = rng_txt, textposition = "none", hovertemplate = "<b>%{y}</b> fail<br>%{text}<extra></extra>")
  shapes <- list()
  eps <- sp$step / 1000
  if (isTRUE(on)) {
    if (b[1] > sp$min + eps)
      shapes <- c(shapes, list(
        list(type = "rect", xref = "x", yref = "paper", x0 = sp$min, x1 = b[1], y0 = 0, y1 = 1,
             fillcolor = pal$gate_wash, line = list(width = 0), layer = "below"),
        list(type = "line", xref = "x", yref = "paper", x0 = b[1], x1 = b[1], y0 = 0, y1 = 1,
             line = list(color = pal$text2, width = 1, dash = "dash"))))
    if (b[2] < sp$max - eps)
      shapes <- c(shapes, list(
        list(type = "rect", xref = "x", yref = "paper", x0 = b[2], x1 = sp$max, y0 = 0, y1 = 1,
             fillcolor = pal$gate_wash, line = list(width = 0), layer = "below"),
        list(type = "line", xref = "x", yref = "paper", x0 = b[2], x1 = b[2], y0 = 0, y1 = 1,
             line = list(color = pal$text2, width = 1, dash = "dash"))))
  }
  xa <- .qc_axis(pal, NULL, range = c(sp$min, sp$max), showgrid = FALSE, fixedrange = TRUE,
                 nticks = 5)
  if (sp$log) xa <- c(xa, .qc_log_ticks(c(sp$min, sp$max)))
  # Cells beyond the window are drawn in the end bins; say how many
  ann <- list()
  edge <- function(x, xa_, txt)
    list(xref = "paper", yref = "paper", x = x, y = 1, xanchor = xa_, yanchor = "top",
         showarrow = FALSE, text = txt, font = list(size = 10, color = pal$muted),
         bgcolor = pal$surface)
  if (sp$n_above) ann <- c(ann, list(edge(1, "right", sprintf(paste0("%d beyond %s ", .qc_s[["rsaquo"]]), sp$n_above,
                                                              .qc_fmt(to_x(sp$max))))))
  if (sp$n_below) ann <- c(ann, list(edge(0, "left", sprintf(paste0(.qc_s[["lsaquo"]], " %d below %s"), sp$n_below,
                                                             .qc_fmt(to_x(sp$min))))))
  p <- .qc_layout(p, pal, barmode = "stack", bargap = 0, showlegend = FALSE, shapes = shapes,
                  annotations = ann,
                  margin = list(l = 30, r = 6, t = 4, b = 22), hovermode = "closest",
                  xaxis = xa, yaxis = .qc_axis(pal, NULL, showline = FALSE, fixedrange = TRUE,
                                               nticks = 4))
  plotly::config(p, displayModeBar = FALSE, responsive = TRUE)
}

# ── Spearman correlation between gate metrics ────────────────────────────────
.qc_plot_corr <- function(qd, pal)
{
  M <- qd$metrics
  R <- suppressWarnings(stats::cor(qd$S, method = "spearman", use = "pairwise.complete.obs"))
  lab <- M$label
  p <- plotly::plot_ly(source = "corr")
  # Neutral strength scale (|rho|): red and green stay reserved for fail / pass.
  # The sign is in the tooltip.
  rho_txt <- matrix(formatC(R, format = "f", digits = 2), nrow(R))
  p <- plotly::add_heatmap(p, z = abs(R), x = lab, y = lab, text = rho_txt,
                           zmin = 0, zmax = 1, xgap = 2, ygap = 2,
                           colorscale = list(c(0, pal$mag[1]), c(1, pal$mag[2])),
                           colorbar = list(thickness = 8, len = 0.7, outlinewidth = 0,
                                           tickvals = c(0, 0.5, 1),
                                           tickfont = list(size = 10, color = pal$muted),
                                           title = list(text = paste0("|", .qc_s[["rho"]], "|"),
                                                        font = list(color = pal$muted))),
                           hovertemplate = paste0("<b>", .qc_s[["rho"]], " = %{text}</b><br>%{y}<br>%{x}<extra></extra>"))
  p <- .qc_layout(p, pal, margin = list(l = 10, r = 10, t = 6, b = 10),
                  xaxis = .qc_axis(pal, NULL, showgrid = FALSE, showline = FALSE, tickangle = -40,
                                   fixedrange = TRUE),
                  yaxis = .qc_axis(pal, NULL, showgrid = FALSE, showline = FALSE,
                                   autorange = "reversed", fixedrange = TRUE))
  p <- plotly::config(p, displayModeBar = FALSE, responsive = TRUE)
  plotly::event_register(p, "plotly_click")
}

# Copy-number colour key for the profile plots (colours come from the CSS, so
# they follow the theme)
.qc_cn_legend <- function()
  shiny::div(class = "qc-cn-legend", shiny::span(class = "qc-cn-title", "Copy number"),
             lapply(c(0:5, "6+"), function(k)
               shiny::span(shiny::tags$i(style = sprintf("background: var(--qc-cn-%s)", sub("+", "", k, fixed = TRUE))),
                           k)))

# ── ASCAT.sc hex logo ────────────────────────────────────────────────────────
# The image is set in the stylesheet: the dark hex in dark mode, the light hex in light
# mode, so a theme switch swaps it without a re-render.
.qc_logo <- function(class = NULL)
  shiny::span(class = paste(c("qc-logo-hex", class), collapse = " "),
              role = "img", `aria-label` = "ASCAT.sc")

# =============================================================================
# Embedded assets: stylesheet, client script and logos. They live here, not in
# separate files, so that this one file is the whole app.
# =============================================================================

.qc_logo_uri <- function(which = c("dark", "light"))
{
  b64 <- switch(match.arg(which), dark = .qc_logo_dark_b64, light = .qc_logo_light_b64)
  paste0("data:image/webp;base64,", gsub("[[:space:]]", "", b64))
}

.qc_css <- function()
{
  css <- gsub("@LOGO_DARK@", .qc_logo_uri("dark"), .qc_css_src, fixed = TRUE)
  gsub("@LOGO_LIGHT@", .qc_logo_uri("light"), css, fixed = TRUE)
}

# ── Stylesheet ───────────────────────────────────────────────────────────────
.qc_css_src <- r"---(
/* =============================================================================
   ASCAT.sc QC app (launch_ui_qc) - design system
   Same tokens as the ASCAT.sc Explorer web app: Apple system colours,
   glass chrome, flat opaque content cards. Every colour is a variable, so
   switching theme (data-bs-theme on <html>) never needs a re-render.
   ============================================================================= */

/* -- Tokens --------------------------------------------------------------- */
html[data-bs-theme="light"], html:not([data-bs-theme]) {
  color-scheme: light;
  --qc-bg: #f5f5f7;
  --qc-card: #ffffff;
  --qc-inset: #f5f5f7;
  --qc-item: #e8e8ed;
  --qc-seg-on: #ffffff;
  --qc-border: rgba(0, 0, 0, 0.08);
  --qc-border-2: rgba(0, 0, 0, 0.14);
  --qc-text: #1d1d1f;
  --qc-text-2: rgba(60, 60, 67, 0.78);
  --qc-text-3: rgba(60, 60, 67, 0.6);
  --qc-text-4: rgba(60, 60, 67, 0.4);
  --qc-track: #dcdce1;
  --qc-switch-off: #d1d1d6;
  /* Colour carries meaning only: green = pass, red = fail. Controls are
     neutral graphite. Pass/fail were checked for colour-blind separation. */
  --qc-pass: #058a5e;
  --qc-fail: #d92d20;
  --qc-good-text: #046c4a;
  --qc-bad-text: #b42318;
  --qc-pass-soft: rgba(5, 138, 94, 0.12);
  --qc-fail-soft: rgba(217, 45, 32, 0.10);
  --qc-fail-mid: rgba(217, 45, 32, 0.32);
  --qc-pass-line: rgba(5, 138, 94, 0.5);
  --qc-fail-line: rgba(217, 45, 32, 0.45);
  --qc-spark-in: rgba(60, 60, 67, 0.42);
  --qc-accent: #48484a;
  --qc-accent-soft: rgba(60, 60, 67, 0.08);
  --qc-accent-text: #1d1d1f;
  --qc-accent-rgb: 60, 60, 67;
  --qc-switch-on: #1d1d1f;
  --qc-primary-bg: #1d1d1f;
  --qc-primary-text: #ffffff;
  --qc-hover: rgba(0, 0, 0, 0.035);
  --qc-focus: rgba(60, 60, 67, 0.45);
  --qc-cn-2: #888888;
  --qc-cn-3: #c49631;
  --qc-glass: rgba(255, 255, 255, 0.72);
  --qc-glass-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.8), 0 4px 24px rgba(0, 0, 0, 0.05);
  --qc-menu: rgba(255, 255, 255, 0.94);
  --qc-menu-shadow: 0 10px 30px rgba(0, 0, 0, 0.14), 0 0 0 0.5px rgba(0, 0, 0, 0.08);
  --qc-scroll: #c7c7cc;
}
html[data-bs-theme="dark"] {
  color-scheme: dark;
  --qc-bg: #1c1c1e;
  --qc-card: #242426;
  --qc-inset: #1c1c1e;
  --qc-item: #3a3a3c;
  --qc-seg-on: #5a5a5e;
  --qc-border: rgba(255, 255, 255, 0.08);
  --qc-border-2: rgba(255, 255, 255, 0.14);
  --qc-text: #f5f5f7;
  --qc-text-2: rgba(235, 235, 245, 0.78);
  --qc-text-3: rgba(235, 235, 245, 0.55);
  --qc-text-4: rgba(235, 235, 245, 0.35);
  --qc-track: #48484a;
  --qc-switch-off: #39393d;
  /* Colour carries meaning only: green = pass, red = fail. Controls are
     space grey. Pass/fail were checked for colour-blind separation. */
  --qc-pass: #12a37a;
  --qc-fail: #ff453a;
  --qc-good-text: #3fd19f;
  --qc-bad-text: #ff7b72;
  --qc-pass-soft: rgba(18, 163, 122, 0.18);
  --qc-fail-soft: rgba(255, 69, 58, 0.16);
  --qc-fail-mid: rgba(255, 69, 58, 0.36);
  --qc-pass-line: rgba(18, 163, 122, 0.6);
  --qc-fail-line: rgba(255, 69, 58, 0.5);
  --qc-spark-in: rgba(235, 235, 245, 0.40);
  --qc-accent: #aeaeb2;
  --qc-accent-soft: rgba(174, 174, 178, 0.16);
  --qc-accent-text: #f5f5f7;
  --qc-accent-rgb: 174, 174, 178;
  --qc-switch-on: #e5e5ea;
  --qc-primary-bg: #e5e5ea;
  --qc-primary-text: #1c1c1e;
  --qc-hover: rgba(255, 255, 255, 0.05);
  --qc-focus: rgba(235, 235, 245, 0.45);
  --qc-cn-2: #c1c1c1;
  --qc-cn-3: #e3b55f;
  --qc-glass: rgba(28, 28, 30, 0.66);
  --qc-glass-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.06), 0 4px 24px rgba(0, 0, 0, 0.35);
  --qc-menu: rgba(44, 44, 46, 0.96);
  --qc-menu-shadow: 0 10px 30px rgba(0, 0, 0, 0.5), 0 0 0 0.5px rgba(255, 255, 255, 0.1);
  --qc-scroll: #48484a;
}
html[data-bs-theme] {
  /* copy-number colours of the profile plot (get_total_cn_colors, capped at 6+) */
  --qc-cn-0: #153570;
  --qc-cn-1: #577aba;
  --qc-cn-4: #d6804f;
  --qc-cn-5: #b3402e;
  --qc-cn-6: #821010;
  --qc-header-h: 52px;
  --qc-sidebar-w: 324px;
  --qc-mono: ui-monospace, "SF Mono", SFMono-Regular, Menlo, Consolas, monospace;
  /* Bootstrap colour-mode variables, so Bootstrap/DT/selectize parts follow */
  --bs-body-bg: var(--qc-bg);
  --bs-body-color: var(--qc-text);
  --bs-emphasis-color: var(--qc-text);
  --bs-secondary-color: var(--qc-text-3);
  --bs-tertiary-color: var(--qc-text-4);
  --bs-secondary-bg: var(--qc-item);
  --bs-tertiary-bg: var(--qc-inset);
  --bs-border-color: var(--qc-border);
  --bs-border-color-translucent: var(--qc-border);
  --bs-link-color: var(--qc-accent);
  --bs-link-hover-color: var(--qc-accent);
  --bs-link-color-rgb: var(--qc-accent-rgb);
  --bs-primary-rgb: var(--qc-accent-rgb);
  --bs-primary: var(--qc-accent);
}

/* -- Base ----------------------------------------------------------------- */
html, body { background: var(--qc-bg); }
body {
  color: var(--qc-text);
  -webkit-font-smoothing: antialiased;
  letter-spacing: -0.003em;
}
.qc-app { min-height: 100vh; }
.qc-app :focus-visible { outline: 2px solid var(--qc-focus); outline-offset: 2px; }
.qc-mono { font-family: var(--qc-mono); font-variant-numeric: tabular-nums; }
.qc-muted { color: var(--qc-text-3); }
.qc-faint { color: var(--qc-text-4); }

::-webkit-scrollbar { width: 9px; height: 9px; }
::-webkit-scrollbar-track { background: transparent; }
::-webkit-scrollbar-thumb {
  background: var(--qc-scroll); border-radius: 9px;
  border: 2px solid transparent; background-clip: content-box;
}
* { scrollbar-color: var(--qc-scroll) transparent; scrollbar-width: thin; }

/* Quick updates should not flash; slow widgets dim after a short delay. */
.qc-app .recalculating { opacity: 1; }
.qc-app .html-widget-output.recalculating {
  opacity: 0.6; transition: opacity 0.25s ease 0.25s;
}

/* -- Header (glass chrome, sticky over scrolling content) ----------------- */
.qc-header {
  position: sticky; top: 0; z-index: 1030;
  height: var(--qc-header-h);
  display: flex; align-items: center; gap: 12px;
  padding: 0 16px;
  background: var(--qc-glass);
  -webkit-backdrop-filter: blur(24px) saturate(180%);
  backdrop-filter: blur(24px) saturate(180%);
  border-bottom: 1px solid var(--qc-border);
  box-shadow: var(--qc-glass-shadow);
}
.qc-brand { display: flex; align-items: center; gap: 10px; min-width: 0; }
/* Hex logo: dark hex in dark mode, light hex in light mode (embedded WebP) */
.qc-logo-hex {
  display: inline-block; flex: none; width: 35px; height: 40px;
  background: url("@LOGO_LIGHT@") center / contain no-repeat;
}
/* The dark hex's fill is close to the dark surfaces; a hairline glow that
   follows its outline keeps the hexagon readable. */
html[data-bs-theme="dark"] .qc-logo-hex {
  background-image: url("@LOGO_DARK@");
  filter: drop-shadow(0 0 0.6px rgba(255, 255, 255, 0.55));
}
.qc-title { font-size: 14px; font-weight: 600; letter-spacing: -0.01em; white-space: nowrap; }
.qc-title span { color: var(--qc-text-3); font-weight: 400; }
.qc-pills { display: flex; align-items: center; gap: 6px; flex-wrap: nowrap; overflow: hidden; }
.qc-pill {
  display: inline-flex; align-items: center; gap: 6px;
  font-size: 11.5px; line-height: 1; padding: 5px 9px; border-radius: 999px;
  background: var(--qc-item); color: var(--qc-text-3); white-space: nowrap;
}
.qc-pill .qc-dot { width: 6px; height: 6px; }
.qc-header-actions { margin-left: auto; display: flex; align-items: center; gap: 8px; }

.qc-btn, .qc-app .btn.qc-btn {
  display: inline-flex; align-items: center; gap: 6px;
  font-size: 12px; font-weight: 500; line-height: 1;
  padding: 7px 11px; border-radius: 8px; border: 0;
  background: var(--qc-item); color: var(--qc-text-2);
  transition: background-color 0.15s ease, color 0.15s ease, filter 0.15s ease;
  white-space: nowrap;
}
.qc-btn:hover, .qc-app .btn.qc-btn:hover { color: var(--qc-text); filter: brightness(1.07); background: var(--qc-item); }
.qc-btn:active { filter: brightness(0.95); }
.qc-app .btn.qc-btn-primary, .qc-btn-primary { background: var(--qc-primary-bg); color: var(--qc-primary-text); }
.qc-app .btn.qc-btn-primary:hover, .qc-btn-primary:hover { background: var(--qc-primary-bg); color: var(--qc-primary-text); filter: brightness(1.06); }
.qc-btn-ghost { background: transparent; }
.qc-btn-sm { padding: 5px 9px; font-size: 11.5px; }
.qc-btn-icon { padding: 7px 8px; }
.qc-btn-danger { color: var(--qc-bad-text); }
.qc-app .dropdown-toggle::after { display: none; }

bslib-input-dark-mode {
  --text-1: var(--qc-text-2) !important;
  --text-2: var(--qc-text-3) !important;
  display: inline-flex; align-items: center; justify-content: center;
  width: 31px; height: 27px; border-radius: 8px; background: var(--qc-item);
}

/* Glass dropdown menus (export) */
.qc-app .dropdown-menu {
  --bs-dropdown-bg: var(--qc-menu);
  background: var(--qc-menu);
  border: 1px solid var(--qc-border); border-radius: 12px;
  box-shadow: var(--qc-menu-shadow);
  padding: 5px; min-width: 230px; font-size: 12.5px;
}
.qc-app .dropdown-header {
  font-size: 10.5px; text-transform: uppercase; letter-spacing: 0.07em;
  color: var(--qc-text-4); padding: 6px 10px 4px;
}
.qc-app .dropdown-item {
  display: flex; align-items: center; gap: 9px;
  border-radius: 8px; padding: 7px 10px; color: var(--qc-text-2);
}
.qc-app .dropdown-item .fa, .qc-app .dropdown-item svg { width: 14px; color: var(--qc-text-3); }
.qc-app .dropdown-item:hover, .qc-app .dropdown-item:focus {
  background: var(--qc-accent-soft); color: var(--qc-text);
}
.qc-app .dropdown-item small { margin-left: auto; color: var(--qc-text-4); }

/* -- Layout --------------------------------------------------------------- */
.qc-body {
  display: grid;
  grid-template-columns: var(--qc-sidebar-w) minmax(0, 1fr);
  transition: grid-template-columns 0.22s ease;
}
body.qc-sidebar-hidden .qc-body { grid-template-columns: 0 minmax(0, 1fr); }
.qc-sidebar {
  position: sticky; top: var(--qc-header-h); align-self: start;
  height: calc(100vh - var(--qc-header-h));
  overflow-y: auto; overflow-x: hidden;
  padding: 14px 12px 28px 14px;
  border-right: 1px solid var(--qc-border);
  background: var(--qc-glass);
}
body.qc-sidebar-hidden .qc-sidebar { visibility: hidden; padding: 0; border: 0; }
.qc-main { position: relative; min-width: 0; padding: 16px 22px 48px; }

/* -- Section labels ------------------------------------------------------- */
.qc-label {
  font-size: 10.5px; font-weight: 600; text-transform: uppercase;
  letter-spacing: 0.08em; color: var(--qc-text-3);
}
.qc-side-section { margin-bottom: 14px; }
.qc-side-section > .qc-label { display: flex; align-items: center; justify-content: space-between; margin: 0 2px 8px; }
.qc-group-label { margin: 14px 2px 7px; color: var(--qc-text-4); }

/* Auto-threshold row */
.qc-auto {
  display: flex; align-items: center; gap: 6px; flex-wrap: wrap;
  padding: 9px 10px; border-radius: 12px;
  background: var(--qc-card); border: 1px solid var(--qc-border);
  font-size: 12px; color: var(--qc-text-2);
}
.qc-auto .form-group { margin: 0; width: auto; }
.qc-auto input.form-control {
  width: 52px; height: 26px; padding: 2px 6px; text-align: center;
  font-family: var(--qc-mono); font-size: 12px;
  background: var(--qc-item); border: 0; border-radius: 7px; color: var(--qc-text);
}
.qc-auto .qc-auto-actions { margin-left: auto; display: flex; gap: 6px; }

/* -- Gate cards ----------------------------------------------------------- */
.qc-gate {
  position: relative;
  padding: 9px 10px 4px; margin-bottom: 7px;
  border-radius: 12px;
  background: var(--qc-card); border: 1px solid var(--qc-border);
  transition: border-color 0.15s ease;
}
.qc-gate:hover { border-color: var(--qc-border-2); }
.qc-gate-head {
  display: flex; align-items: center; gap: 8px;
  min-height: 20px; padding-right: 64px;
}
.qc-gate-head .shiny-input-container { width: auto !important; margin: 0; }
.qc-gate-name { font-size: 12.5px; font-weight: 500; color: var(--qc-text); line-height: 1.2; }
.qc-help { color: var(--qc-text-4); font-size: 11px; cursor: help; }
.qc-help:hover { color: var(--qc-text-2); }
.qc-gate-badge {
  position: absolute; top: 9px; right: 10px;
  font-size: 11px; line-height: 1; padding: 4px 7px; border-radius: 999px;
  font-variant-numeric: tabular-nums;
  display: inline-flex; align-items: center; gap: 5px;
  background: var(--qc-item); color: var(--qc-text-2);
}
.qc-gate-badge::before { content: ""; width: 6px; height: 6px; border-radius: 50%; background: var(--qc-fail); }
.qc-gate-badge.is-zero { color: var(--qc-text-4); }
.qc-gate-badge.is-zero::before { background: var(--qc-text-4); }
.qc-gate-read {
  display: flex; justify-content: space-between; align-items: baseline;
  margin: 6px 1px 3px; font-size: 11.5px;
}
.qc-gate-limit { font-family: var(--qc-mono); font-variant-numeric: tabular-nums; color: var(--qc-text-2); }
.qc-gate-med { font-size: 10.5px; color: var(--qc-text-4); }
.qc-gate .shiny-input-container:has(.irs) { margin: 0 0 2px; }
.qc-gate:has(.form-check-input:not(:checked)) .qc-gate-body { opacity: 0.42; }
.qc-gate:has(.form-check-input:not(:checked)) .irs--shiny .irs-bar { background: var(--qc-text-4); }
.qc-gate:has(.form-check-input:not(:checked)) .qc-gate-badge { display: none; }

/* Sparklines: neutral inside the gate, red for the cells it removes. */
.qc-spark { display: block; width: 100%; height: 28px; overflow: visible; }
.qc-spark rect.in { fill: var(--qc-spark-in); }
.qc-spark rect.out { fill: var(--qc-fail); opacity: 0.85; }
.qc-spark rect.off { fill: var(--qc-text-4); }
.qc-spark line.base { stroke: var(--qc-border-2); stroke-width: 1; vector-effect: non-scaling-stroke; }
.qc-spark line.mk { stroke: var(--qc-text); stroke-width: 2.5; vector-effect: non-scaling-stroke; }

/* Flags */
.qc-flag {
  display: flex; align-items: center; gap: 9px;
  padding: 8px 10px; margin-bottom: 6px; border-radius: 10px;
  background: var(--qc-card); border: 1px solid var(--qc-border);
  font-size: 12.5px; color: var(--qc-text-2);
}
.qc-flag .shiny-input-container { width: auto !important; margin: 0; }
.qc-flag-n { margin-left: auto; font-size: 11px; color: var(--qc-text-4); font-variant-numeric: tabular-nums; }

.qc-side-foot { display: flex; gap: 6px; margin-top: 12px; }
.qc-side-foot .btn { flex: 1; justify-content: center; }
.qc-kbd-hints { margin-top: 16px; font-size: 11px; color: var(--qc-text-4); line-height: 1.9; }
kbd, .qc-kbd {
  font-family: var(--qc-mono); font-size: 10.5px; padding: 1px 5px; border-radius: 5px;
  background: var(--qc-item); color: var(--qc-text-3); box-shadow: none;
}

/* -- iOS switch (bslib input_switch) -------------------------------------- */
.qc-app .form-switch { padding-left: 0; margin: 0; min-height: 0; display: flex; align-items: center; }
.qc-app .form-switch .form-check-input,
.qc-app .form-switch .form-check-input:focus,
.qc-app .form-switch .form-check-input:checked {
  --bs-form-switch-bg: url("data:image/svg+xml,%3csvg xmlns='http://www.w3.org/2000/svg' viewBox='-4 -4 8 8'%3e%3ccircle r='3.25' fill='%23fff'/%3e%3c/svg%3e");
}
.qc-app .form-switch .form-check-input {
  width: 30px; height: 18px; margin: 0; flex: none; float: none;
  border: 0; border-radius: 9px; cursor: pointer;
  background-color: var(--qc-switch-off);
  box-shadow: none;
  transition: background-position 0.18s ease, background-color 0.18s ease;
}
.qc-app .form-switch .form-check-input:checked { background-color: var(--qc-switch-on); }
/* Dark mode: the "on" track is bright, so its knob is dark, and on/off can never be confused */
html[data-bs-theme="dark"] .qc-app .form-switch .form-check-input:checked {
  --bs-form-switch-bg: url("data:image/svg+xml,%3csvg xmlns='http://www.w3.org/2000/svg' viewBox='-4 -4 8 8'%3e%3ccircle r='3.25' fill='%231c1c1e'/%3e%3c/svg%3e");
}
.qc-app .form-switch .form-check-input:focus { box-shadow: 0 0 0 3px var(--qc-accent-soft); }
.qc-app .form-switch .form-check-label:empty,
.qc-app .form-switch .form-check-label > span:empty { display: none; }

/* -- Sliders (ion.rangeSlider, labels replaced by our own readout) -------- */
.qc-app .irs--shiny { height: 22px; }
.qc-app .irs--shiny .irs-line {
  top: 9px; height: 4px; border: 0; border-radius: 4px; background: var(--qc-track);
}
.qc-app .irs--shiny .irs-bar {
  top: 9px; height: 4px; border: 0; background: var(--qc-accent);
}
.qc-app .irs--shiny .irs-handle {
  top: 2px; width: 18px; height: 18px; border: 0; border-radius: 50%;
  background: #fff; cursor: grab;
  box-shadow: 0 0.5px 4px rgba(0, 0, 0, 0.35), 0 0 0 0.5px rgba(0, 0, 0, 0.06);
}
.qc-app .irs--shiny .irs-handle:hover, .qc-app .irs--shiny .irs-handle.state_hover { background: #fff; }
.qc-app .irs--shiny .irs-handle:active { cursor: grabbing; }
.qc-app .irs--shiny .irs-min, .qc-app .irs--shiny .irs-max,
.qc-app .irs--shiny .irs-from, .qc-app .irs--shiny .irs-to,
.qc-app .irs--shiny .irs-single, .qc-app .irs--shiny .irs-grid { display: none; }

/* -- Segmented controls --------------------------------------------------- */
/* Tabs (bslib navset_pill) */
.qc-main > .nav-pills, .qc-main .qc-tabs > .nav-pills, .qc-main > .tabbable > .nav-pills {
  display: flex; width: fit-content; flex-wrap: nowrap; gap: 2px; padding: 2px;
  border-radius: 10px; background: var(--qc-item); margin: 0 auto 16px;
}
.qc-main .nav-pills .nav-link {
  display: inline-flex; align-items: center; gap: 7px;
  padding: 6px 14px; border-radius: 8px;
  font-size: 12.5px; font-weight: 500; color: var(--qc-text-3);
  transition: background-color 0.15s ease, color 0.15s ease;
}
.qc-main .nav-pills .nav-link:hover { color: var(--qc-text); }
.qc-main .nav-pills .nav-link.active, .qc-main .nav-pills .show > .nav-link {
  background: var(--qc-seg-on); color: var(--qc-text);
  box-shadow: 0 1px 2px rgba(0, 0, 0, 0.14), 0 0 0 0.5px rgba(0, 0, 0, 0.04);
}
.qc-main .nav-pills .nav-link .fa, .qc-main .nav-pills .nav-link svg { opacity: 0.8; font-size: 11.5px; }

/* Radio-based segmented control (.qc_segmented) */
.qc-seg {
  display: inline-flex; gap: 2px; padding: 2px; margin: 0;
  border-radius: 9px; background: var(--qc-item); vertical-align: middle;
}
.qc-seg.shiny-input-container { width: auto !important; margin-bottom: 0; }
.qc-seg label { margin: 0; position: relative; display: inline-flex; }
.qc-seg input { position: absolute; opacity: 0; width: 1px; height: 1px; pointer-events: none; }
.qc-seg span {
  display: inline-flex; align-items: center; gap: 6px;
  padding: 4px 11px; border-radius: 7px; cursor: pointer;
  font-size: 12px; font-weight: 500; color: var(--qc-text-3); white-space: nowrap;
  transition: background-color 0.15s ease, color 0.15s ease;
}
.qc-seg span:hover { color: var(--qc-text); }
.qc-seg input:checked + span {
  background: var(--qc-seg-on); color: var(--qc-text);
  box-shadow: 0 1px 2px rgba(0, 0, 0, 0.14), 0 0 0 0.5px rgba(0, 0, 0, 0.04);
}
.qc-seg input:focus-visible + span { outline: 2px solid var(--qc-focus); outline-offset: 1px; }

/* Status summary in the tab row */
.qc-tab-status {
  position: absolute; top: 18px; right: 22px;
  display: flex; align-items: center; gap: 10px;
  font-size: 12px; color: var(--qc-text-3);
}
.qc-tab-status b { color: var(--qc-text); font-weight: 600; }

/* -- Cards ---------------------------------------------------------------- */
.qc-card {
  background: var(--qc-card); border: 1px solid var(--qc-border);
  border-radius: 16px; padding: 14px 16px; margin-bottom: 16px; min-width: 0;
}
.qc-card-head {
  display: flex; align-items: center; justify-content: space-between;
  gap: 12px; margin-bottom: 10px; min-height: 22px;
}
.qc-card-head-left { display: flex; align-items: baseline; gap: 10px; min-width: 0; flex-wrap: wrap; }
.qc-card-title {
  font-size: 11px; font-weight: 600; text-transform: uppercase;
  letter-spacing: 0.07em; color: var(--qc-text-3); white-space: nowrap;
}
.qc-card-sub { font-size: 12px; color: var(--qc-text-4); }
.qc-grid-2 { display: grid; grid-template-columns: minmax(0, 5fr) minmax(0, 7fr); gap: 16px; }
.qc-grid-2 > .qc-card { margin-bottom: 0; }
.qc-stack > * + * { margin-top: 16px; }

/* Toolbar rows (filters sit in one row above what they scope) */
.qc-toolbar {
  display: flex; align-items: center; gap: 10px; flex-wrap: wrap; margin-bottom: 12px;
}
.qc-toolbar .form-group { margin: 0; }
.qc-toolbar .qc-field { display: flex; align-items: center; gap: 7px; }
.qc-toolbar .qc-field > .qc-label { white-space: nowrap; }
.qc-toolbar .qc-spacer { flex: 1; }
.qc-hint { font-size: 11.5px; color: var(--qc-text-4); }
.qc-hint .fa { margin-right: 4px; }

/* Legend chips (HTML legends above charts) */
.qc-legend { display: flex; align-items: center; gap: 14px; font-size: 11.5px; color: var(--qc-text-3); flex-wrap: wrap; }
.qc-legend i { display: inline-block; width: 10px; height: 10px; border-radius: 3px; margin-right: 6px; vertical-align: -1px; }
.qc-legend i.line { height: 0; width: 14px; border-radius: 0; border-top: 1.5px dashed var(--qc-text-2); vertical-align: 3px; }
.qc-legend i.pass { background: var(--qc-pass); }
.qc-legend i.fail { background: var(--qc-fail); }
.qc-legend i.fail-soft { background: var(--qc-fail-mid); }

/* -- Selects (selectize) -------------------------------------------------- */
.qc-app .selectize-control { margin: 0; }
.qc-app .selectize-input, .qc-app .selectize-input.full, .qc-app .selectize-input.focus {
  min-height: 30px; padding: 6px 28px 6px 10px !important;
  background: var(--qc-item); border: 0; border-radius: 8px; box-shadow: none;
  font-size: 12.5px; color: var(--qc-text); line-height: 18px;
}
.qc-app .selectize-input input { color: var(--qc-text); font-size: 12.5px; }
.qc-app .selectize-input input::placeholder { color: var(--qc-text-4); }
.qc-app .selectize-control.single .selectize-input::after { border-color: var(--qc-text-3) transparent transparent transparent; right: 11px; }
.qc-app .selectize-dropdown {
  background: var(--qc-menu); color: var(--qc-text-2);
  border: 1px solid var(--qc-border); border-radius: 10px;
  box-shadow: var(--qc-menu-shadow); margin-top: 4px; z-index: 1050; font-size: 12.5px;
}
.qc-app .selectize-dropdown-content { padding: 4px; max-height: 320px; }
.qc-app .selectize-dropdown .option, .qc-app .selectize-dropdown .optgroup-header {
  border-radius: 6px; padding: 6px 9px; color: var(--qc-text-2); background: transparent;
}
.qc-app .selectize-dropdown .optgroup-header {
  font-size: 10.5px; text-transform: uppercase; letter-spacing: 0.06em; color: var(--qc-text-4);
}
.qc-app .selectize-dropdown .active { background: var(--qc-accent-soft); color: var(--qc-text); }
.qc-app .selectize-dropdown .selected { color: var(--qc-text); font-weight: 500; }
.qc-select-sm .selectize-control { min-width: 170px; }
.qc-select-cell .selectize-control { min-width: 240px; }
.qc-select-cell .selectize-input { font-family: var(--qc-mono); }

/* -- KPI tiles ------------------------------------------------------------ */
.qc-kpis { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px; margin-bottom: 16px; }
.qc-kpi {
  background: var(--qc-card); border: 1px solid var(--qc-border);
  border-radius: 16px; padding: 13px 16px 14px; min-width: 0;
}
.qc-kpi-label { display: flex; align-items: center; gap: 6px; font-size: 12px; color: var(--qc-text-3); }
.qc-kpi-label .fa { font-size: 11px; }
.qc-kpi-value { margin-top: 4px; font-size: 26px; font-weight: 600; letter-spacing: -0.02em; line-height: 1.15; color: var(--qc-text); }
.qc-kpi-value small { font-size: 13px; font-weight: 500; color: var(--qc-text-3); letter-spacing: 0; margin-left: 4px; }
.qc-kpi-sub { margin-top: 3px; font-size: 12px; color: var(--qc-text-3); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.qc-kpi-sub a { color: var(--qc-text-2); text-decoration: underline; text-underline-offset: 2px; cursor: pointer; }
.qc-kpi-sub a:hover { text-decoration: underline; }
.qc-kpi .icon-pass { color: var(--qc-good-text); }
.qc-kpi .icon-fail { color: var(--qc-bad-text); }
.qc-meter { height: 6px; border-radius: 3px; background: var(--qc-pass-soft); overflow: hidden; margin-top: 10px; }
.qc-meter > div { height: 100%; background: var(--qc-pass); border-radius: 3px; transition: width 0.3s ease; }

/* -- Pass / fail per metric -------------------------------------------- */
.qc-pf {
  display: grid; grid-template-columns: minmax(110px, max-content) minmax(60px, 1fr) auto auto;
  gap: 10px 16px; align-items: center;
}
.qc-pf-h {
  display: flex; align-items: center; gap: 6px;
  font-size: 10.5px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.06em;
  color: var(--qc-text-4); padding-bottom: 7px; border-bottom: 1px solid var(--qc-border);
}
.qc-pf-hn { justify-content: flex-end; }
.qc-pf-label { font-size: 12.5px; color: var(--qc-text-2); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
/* 100% bar: green share passes, red share fails, a 2 px gap between them */
.qc-pf-bar { display: flex; gap: 2px; height: 8px; }
.qc-pf-pass { background: var(--qc-pass); opacity: 0.55; border-radius: 4px 0 0 4px; }
.qc-pf-fail { background: var(--qc-fail); min-width: 3px; border-radius: 0 4px 4px 0; }
.qc-pf-bar > :only-child { border-radius: 4px; }
.qc-pf-n { font-size: 12.5px; color: var(--qc-text); text-align: right; min-width: 4.5em; font-variant-numeric: tabular-nums; }
.qc-pf-n.is-zero { color: var(--qc-text-4); }
/* the combined row, under a full-width rule */
.qc-pf-sep { grid-column: 1 / -1; border-top: 1px solid var(--qc-border); }
.qc-pf-label.is-total, .qc-pf-n.is-total { font-weight: 600; color: var(--qc-text); }
.qc-empty { padding: 26px 8px; text-align: center; font-size: 12.5px; color: var(--qc-text-4); }
.qc-empty .fa { display: block; font-size: 18px; margin-bottom: 8px; opacity: 0.7; }

/* -- Distributions grid --------------------------------------------------- */
.qc-dists { display: grid; grid-template-columns: repeat(auto-fill, minmax(250px, 1fr)); gap: 10px; }
.qc-dist { background: var(--qc-inset); border-radius: 12px; padding: 9px 10px 4px; min-width: 0; }
.qc-dist-head { display: flex; justify-content: space-between; align-items: baseline; gap: 8px; }
.qc-dist-name { font-size: 12px; font-weight: 500; color: var(--qc-text-2); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.qc-dist-gate { font-size: 11px; color: var(--qc-text-4); font-family: var(--qc-mono); white-space: nowrap; }
.qc-dist.is-off .qc-dist-name::after { content: "not used"; margin-left: 6px; font-size: 10.5px; color: var(--qc-text-4); font-weight: 400; }

/* -- Selection bar (lasso) ------------------------------------------------ */
.qc-selbar {
  display: flex; align-items: center; gap: 10px; flex-wrap: wrap;
  padding: 8px 10px 8px 14px; margin-bottom: 16px; border-radius: 12px;
  background: var(--qc-accent-soft); border: 1px solid var(--qc-border-2);
  font-size: 12.5px; color: var(--qc-text);
}
.qc-selbar .qc-spacer { flex: 1; }

/* -- Status pills --------------------------------------------------------- */
.qc-status {
  display: inline-flex; align-items: center; gap: 5px;
  padding: 3px 9px; border-radius: 999px; font-size: 11.5px; font-weight: 600; line-height: 1.2;
  white-space: nowrap;
}
.qc-status.pass { background: var(--qc-pass-soft); color: var(--qc-good-text); }
.qc-status.fail { background: var(--qc-fail-soft); color: var(--qc-bad-text); }
.qc-status.kept { background: transparent; color: var(--qc-good-text); box-shadow: inset 0 0 0 1px var(--qc-pass-line); }
.qc-status.excluded { background: transparent; color: var(--qc-bad-text); box-shadow: inset 0 0 0 1px var(--qc-fail-line); }
.qc-dot { display: inline-block; width: 7px; height: 7px; border-radius: 50%; flex: none; }
.qc-dot.pass { background: var(--qc-pass); }
.qc-dot.fail { background: var(--qc-fail); }
.qc-dot.off { background: var(--qc-text-4); }

/* -- Inspector ------------------------------------------------------------ */
.qc-cell-head { display: flex; align-items: center; gap: 10px; flex-wrap: wrap; }
.qc-cell-id { font-family: var(--qc-mono); font-size: 14px; font-weight: 600; color: var(--qc-text); }
.qc-cell-meta { font-size: 12px; color: var(--qc-text-3); }
.qc-cell-meta b { color: var(--qc-text-2); font-weight: 600; font-variant-numeric: tabular-nums; }
.qc-cell-fails { font-size: 12px; color: var(--qc-text-3); margin-top: 6px; }
.qc-cell-fails b { color: var(--qc-bad-text); font-weight: 500; }
.qc-pos { font-size: 12px; color: var(--qc-text-3); font-variant-numeric: tabular-nums; white-space: nowrap; }

.qc-tiles { display: grid; grid-template-columns: repeat(auto-fill, minmax(178px, 1fr)); gap: 10px; }
.qc-tile { background: var(--qc-inset); border-radius: 14px; padding: 10px 12px 9px; min-width: 0; }
.qc-tile-label { font-size: 11.5px; color: var(--qc-text-3); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.qc-tile-row { display: flex; align-items: baseline; justify-content: space-between; gap: 6px; }
.qc-tile-value { font-size: 18px; font-weight: 600; color: var(--qc-text); letter-spacing: -0.01em; }
.qc-tile-pct { font-size: 11px; color: var(--qc-text-4); font-variant-numeric: tabular-nums; white-space: nowrap; }
.qc-tile .qc-spark { height: 22px; margin: 6px 0 6px; }
.qc-tile-status { display: flex; align-items: center; gap: 6px; font-size: 11.5px; color: var(--qc-text-3); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.qc-tile-status.pass { color: var(--qc-good-text); }
.qc-tile-status.fail { color: var(--qc-bad-text); }
.qc-tile.is-fail { box-shadow: inset 0 0 0 1px var(--qc-fail-line); }

/* Copy-number legend above the profile plots */
.qc-cn-row { display: flex; justify-content: flex-end; margin: -2px 0 2px; }
.qc-cn-legend { display: flex; align-items: center; gap: 10px; flex-wrap: wrap; font-size: 11px; color: var(--qc-text-3); }
.qc-cn-legend .qc-cn-title { color: var(--qc-text-4); }
.qc-cn-legend i { display: inline-block; width: 14px; height: 5px; border-radius: 2px; margin-right: 4px; vertical-align: 1px; }

.qc-explore-sel-head { display: flex; align-items: center; gap: 10px; flex-wrap: wrap; }

/* -- Tables (DT) ---------------------------------------------------------- */
.qc-app .dataTables_wrapper, .qc-app .dt-container { font-size: 12.5px; color: var(--qc-text-2); }
.qc-app table.dataTable { --bs-table-bg: transparent; --bs-table-color: var(--qc-text-2); --bs-table-striped-bg: transparent; border-collapse: collapse !important; }
.qc-app table.dataTable thead th {
  font-size: 10.5px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em;
  color: var(--qc-text-3); white-space: nowrap; background: var(--qc-card);
  border-bottom: 1px solid var(--qc-border-2) !important; border-top: 0 !important;
}
.qc-app table.dataTable tbody td {
  border-top: 1px solid var(--qc-border) !important; border-bottom: 0 !important;
  font-variant-numeric: tabular-nums; white-space: nowrap; padding: 7px 10px;
  color: var(--qc-text-2); background: transparent; box-shadow: none !important;
}
.qc-app table.dataTable tbody td:first-child { font-family: var(--qc-mono); color: var(--qc-text); }
.qc-app table.dataTable tbody tr { cursor: pointer; }
.qc-app table.dataTable tbody tr:hover > * { background: var(--qc-hover) !important; }
.qc-app table.dataTable.no-footer { border-bottom: 0; }
.qc-app .dataTables_filter input, .qc-app .dt-search input,
.qc-app .dataTables_length select, .qc-app .dt-length select {
  background-color: var(--qc-item); border: 0; border-radius: 8px; color: var(--qc-text);
  font-size: 12.5px; padding: 5px 10px; box-shadow: none;
}
.qc-app .dataTables_info, .qc-app .dt-info { color: var(--qc-text-4); font-size: 12px; }
.qc-app .pagination { --bs-pagination-bg: transparent; --bs-pagination-border-color: transparent; --bs-pagination-color: var(--qc-text-3); --bs-pagination-hover-bg: var(--qc-item); --bs-pagination-hover-color: var(--qc-text); --bs-pagination-hover-border-color: transparent; --bs-pagination-active-bg: var(--qc-primary-bg); --bs-pagination-active-color: var(--qc-primary-text); --bs-pagination-active-border-color: var(--qc-primary-bg); --bs-pagination-disabled-bg: transparent; --bs-pagination-disabled-color: var(--qc-text-4); --bs-pagination-disabled-border-color: transparent; --bs-pagination-focus-box-shadow: none; gap: 2px; }
.qc-app .page-link { border-radius: 7px !important; font-size: 12px; padding: 4px 10px; }
.qc-app .dataTables_processing, .qc-app .dt-processing { background: var(--qc-menu); color: var(--qc-text-2); border-radius: 10px; border: 1px solid var(--qc-border); }

/* -- Tooltips, popovers, notifications, modals ---------------------------- */
.tooltip {
  --bs-tooltip-bg: var(--qc-menu); --bs-tooltip-color: var(--qc-text-2);
  --bs-tooltip-max-width: 270px; --bs-tooltip-font-size: 12px; --bs-tooltip-padding-x: 10px; --bs-tooltip-padding-y: 7px;
  --bs-tooltip-border-radius: 9px; --bs-tooltip-opacity: 1;
}
.tooltip .tooltip-inner { box-shadow: var(--qc-menu-shadow); text-align: left; line-height: 1.45; }
#shiny-notification-panel { width: 320px; }
.shiny-notification {
  background: var(--qc-menu); color: var(--qc-text-2); border: 1px solid var(--qc-border);
  border-radius: 12px; box-shadow: var(--qc-menu-shadow); font-size: 12.5px; opacity: 1;
}
.modal-content {
  --bs-modal-bg: var(--qc-card); --bs-modal-color: var(--qc-text); --bs-modal-border-color: var(--qc-border);
  --bs-modal-border-radius: 16px; --bs-modal-header-border-color: var(--qc-border); --bs-modal-footer-border-color: var(--qc-border);
}
.modal-title { font-size: 15px; font-weight: 600; }

/* -- Responsive ----------------------------------------------------------- */
@media (max-width: 1280px) {
  .qc-grid-2 { grid-template-columns: minmax(0, 1fr); }
}
@media (max-width: 1180px) {
  .qc-tab-status { position: static; justify-content: center; margin: 0 0 10px; }
}
@media (max-width: 1100px) {
  .qc-kpis { grid-template-columns: repeat(2, minmax(0, 1fr)); }
  .qc-pills .qc-pill:not(:first-child) { display: none; }
}
@media (max-width: 900px) {
  .qc-body, body.qc-sidebar-hidden .qc-body { grid-template-columns: minmax(0, 1fr); }
  .qc-sidebar { position: static; height: auto; border-right: 0; border-bottom: 1px solid var(--qc-border); }
  body.qc-sidebar-hidden .qc-sidebar { display: none; }
  .qc-main { padding: 14px 14px 40px; }
  .qc-title span, .qc-pills { display: none; }
}
@media (max-width: 560px) {
  .qc-kpis { grid-template-columns: minmax(0, 1fr); }
  .qc-main .nav-pills .nav-link { padding: 6px 9px; }
  .qc-main .nav-pills .nav-link .qc-tab-text { display: none; }
}

/* -- "Done" overlay (shown while the app hands the result back to R) ------ */
.qc-done {
  position: fixed; inset: 0; z-index: 3000;
  display: flex; align-items: center; justify-content: center; padding: 16px;
  background: var(--qc-glass);
  -webkit-backdrop-filter: blur(24px) saturate(180%);
  backdrop-filter: blur(24px) saturate(180%);
}
body:has(.qc-done) #shiny-disconnected-overlay { display: none; }
.qc-done-card {
  max-width: 420px; text-align: center; padding: 28px 30px 24px; border-radius: 20px;
  background: var(--qc-card); border: 1px solid var(--qc-border); box-shadow: var(--qc-menu-shadow);
}
.qc-done-card .qc-logo-hex { display: block; width: 104px; height: 120px; margin: 0 auto 14px; }
.qc-done-title { font-size: 17px; font-weight: 600; margin-bottom: 6px; }
.qc-done-text { font-size: 13px; color: var(--qc-text-3); line-height: 1.5; margin-bottom: 14px; }
.qc-done-card code {
  display: inline-block; font-family: var(--qc-mono); font-size: 12px; padding: 5px 10px;
  border-radius: 8px; background: var(--qc-item); color: var(--qc-text-2);
}

/* Shiny gives non-empty uiOutput() divs display: contents; these need a box. */
.qc-insp-head { display: block; width: 100%; margin-bottom: 10px; }
.qc-fails { color: var(--qc-text-2); }
.qc-more {
  display: inline-block; margin-left: 3px; padding: 1px 6px; border-radius: 999px;
  font-size: 10.5px; background: var(--qc-item); color: var(--qc-text-3);
}
)---"

# ── Client script (keyboard navigation, sidebar toggle) ──────────────────────
.qc_js <- r"---(
// ASCAT.sc QC app (launch_ui_qc): the only client-side behaviour.
// Everything else is driven from R.
(function () {
  // Left / right arrows (or k / j) step through cells in the Inspector.
  document.addEventListener('keydown', function (e) {
    if (e.metaKey || e.ctrlKey || e.altKey || e.defaultPrevented) return;
    var t = e.target;
    if (t && t.closest && t.closest('input, textarea, select, [contenteditable="true"], .irs, .selectize-control, .dropdown-menu')) return;
    var dir = 0;
    if (e.key === 'ArrowRight' || e.key === 'j') dir = 1;
    else if (e.key === 'ArrowLeft' || e.key === 'k') dir = -1;
    if (!dir || !window.Shiny || !Shiny.setInputValue) return;
    e.preventDefault();
    Shiny.setInputValue('qc_nav', dir, { priority: 'event' });
  });

  // Show / hide the gate sidebar; plotly charts re-measure afterwards.
  document.addEventListener('click', function (e) {
    var b = e.target.closest && e.target.closest('[data-qc-toggle-sidebar]');
    if (!b) return;
    var hidden = document.body.classList.toggle('qc-sidebar-hidden');
    b.setAttribute('aria-pressed', hidden ? 'true' : 'false');
    setTimeout(function () { window.dispatchEvent(new Event('resize')); }, 240);
  });
})();
)---"

# ── Logos: the ASCAT.sc hex, dark and light (lossless WebP, 240 px, base64) ───
.qc_logo_dark_b64 <- r"---(
UklGRg4mAABXRUJQVlA4TAEmAAAvz8A7EDWH47aRHKlM/lG33927d0RMAP5IdUf7nFcs1XXCFUkAUEgSOpIYtQiAtEnOIm0ntiVx
QCoD4AbToguADik0briNVPt8pG0rRdpwE2u4lvZyNWf3Pniq3fz/HElTNofTd7v9+/177w6tA7j8w7rufWd6qqsKDYW0ta48xkOP
j4qjMV8LtWZZGwFay37IYy4EMtD0UxZa2YV6rQuAEOS+QZBBoSGEfsAmg/Zw10NMBDz47aI2glMeFp6sSJJtW+kzXdXncYEv0IMK
BCMAIbhwFEmSrVRVP/4cwb8WBCCEOwIkQJJk2pZt27Zt27Zt27Zt2za/bds2z947sCDZbdvgpdAFc+kEzIzxvi1R23ZcshM7aVR1
cffuwunqrrZt2z62rbaN2LZt27Zt2621qv5dtU9V986Yf+zkYLTj1IqdnXVGtXr1KHZqZXbc052MTuzkXRllx7bNoSu29cfpVbHf
UdaObVZO3L3i0X8ytJM/tvnFeYex1cNYoz+2R3+8M/Uf23ZSSWZOGLZtGya25GL/PzwCSVt49/47cNvGkeTt7YxcS2bLA6RIAOC4
zc1Iu6cZnVra3arM7cndZumR8ln7qVe9m92D6Zr8QvpCCjMMCq6igUZhFeYfvVdqmPoBfSDjBwRWHJZPmLaboLCKg0zF0h8QPyBN
n4lpXmCeRuu9QdAoCI4AgIWkjNHumd1Vtm2+4L5wlVHZtm376rVtvmM8BJI2294/8/vvQJEkRYqBgl7t07Uf0LW2/9jlNGN75vX/
0+aYm2Pbtu2PrzZOYGK7nBOIyllT52xSZuVT+E49VWx7r+lyAjbLQZW19mnENt84vW3btvtJZXPHZjfHENt2a6sbM3ZyCE7qQReb
O+dgtYNKVUprLbnatq1tY/3f9/2/pHAy9QRw69aR9k49o+7czSCGMDMzMzOZwxxDa4m2batNdfY+9aIRHCo4VNBariV97/+/B8vt
vwO3kRTJtcy8e/gA9xg62gvEbi8Q7gWihycxdkduGh9urPYCiXOFJy3oEWwDQcG5Ox5eZDkGftQeGF4MMug/H9ZLYK8+CPvHMLjr
AhVcBulTTjDkYVA2IL4qXPQd1cdhf4gUF+otWYAFxYZCCKlYGA0y8MJhsKiSRrtDyqBCXcsZhqVwRMkLwdPZposhVSVZHFIETlIL
Dcsw3BQWaJpKA9ovWvl+YHJjZ2nIfQMZIR76tiTABj4skY+z1KfFFuxcN9g4S4oN7XIfC8PCcqV1fDrxBGBoXUa1kiAnZf7p3ACE
UmiktD2OLMOHQaRiWBVyc2Ah1dD+zhWksab0esFv5ZF/r0BrzdJr6NwAF8H6Yo4+rKTCUJyx7kvxrksvO0IkmEo9Wy6I62VZzRi9
tN3DoFUz7JfFIE4cj9mnoECTzyQ8obsMfFZjHPmxRWDgUUL6WUMU6TRlOT/xStsPlG4YLAZx2AJIj3PwQaGKhUiQiTeL4lwzbAVh
0JN6ONdFlVKlYsOiILCsFp6Ua4aNIMXrKihtIGY/BeWapsyKycOmP1iHMQzuqrxpCsqB+3o7De9DK/JhlvqGR9zD/n1JbcI4UtZi
84C7m8+3oS2ldczESwWHdZ1kDwhFnW1LNuBzAMRSaFEJirggEufKIFMlWQJyw/vGA+sg4RpAWmFI/QOquv4HzG7sLABBoKBLIX2Y
gw8aaMsKPfU5u7DBdYPpEDwPB9r1eXk30iqD9UzRmdQMkyHwpEhbvqCN12XhjyPrmASximEsRAI7VgPxX93Ga6UpZeIPuYrKQF8x
jITcIUcM8fUct0E7+XBYQc6HLH6uYpiXVQZLaRfm5d2IRgzifhl0aoZZWcWQYpt5ej5taaWog4msfQgCBjG60b0SUAH9XUM0lk+Z
7TNW5mU3DIZ0I6NYCO1ZDj5orh0cfIJ4m8Gt7mBMgJTBiHRqzpdS0F7pUny6PNyQoFozdIfsPy8itJGnN2jjtWVsn86K2Nr233Br
zRhGd2Ogk6OA9n1OWkCUPRAQZAH0nYTEi7oyjuytUnyID8RB9MEorYlQIsh5lZ1dzdATUnBoE/fmzToFw5S+spkFFxjkqyQdD/do
Pb+LbXUjDFzBgP6Sam8MrJothsF9szwyvj3/3gdDBUTh9aIv4SI71w1aHe4J3ppXgnBxzMhg5xjM6j1Dm8M94nrORyAmL/zm6VIn
ECSqJE0O94hYXg4xlPEbUmkRGN3YtX+4hy4B9k6eg+dsgxW0C5nn4zIE1h1MuxAptqxXcvxjsIXSOmaWIwb9FhkVRERO6oIQbE5T
sIgm30sFuMHiiVDFaAny/P9w/a3PKCwjnzK2r6TcqwKNG4cWupEOLoL4so8+WEgnz5nJ7qHzrEsv1XuDFFOpZ/PyNt4qDLZteNKo
GSqzECSelht0I3Zh5I+BOm8MPBVD3eFeeZT8JlcQq8mHnvihsJTb3DCoOdxDF8D6+PWbbbCe1jGTrqM41ozVIQR90tERiAWVthtp
wLYMUawZq0H+Bxs4U1hUIbb0rwHCZdh/UMMq/5w4DO6B/3/UvHLe3iev7Jroy3W5vrSrpZ/Yhj/J+fqeW/Ihcd3QtKhC8ZB6V2E3
0tp1OcsLiZHAmOO7kbX0/xlnZQXiRXSWdVITyO9oWyoh6GxJqueqjur1s6pdx+hImoUt/Wucxn/ymjkM0k2SXietQ+nhHjqvHGpQ
PSGdEiCRclbFqpL9U0IDUMF8RupvQMUrl6zHyywpUj/MFUTdq8tKbY0qD6ItXhIcXxr0iMrdtrULWn6dhoztmRwht/z9kIzjgteL
wKPXx9OgUtFRtjqD/WKIyBSlcpWA42KGIfWaYZNuuDDcZdg7GlJa8aI5OWpVPON+vJmYvyp22gGDoLO4EPS6GxfEN0FP46QyViPm
RoWKYbohTaHqnTckNAiJ1wZ6bnGgPlcbTMyL8qWVDq4rrFp5L0MB9MLigD5XFsj9MKRCzYgthTvjtG64QTSCoQbRqHkZuclVJRSt
CBUBWg8o99d1XttFQw0KsufyFjouBQGd3rOhYaW2CMloNDtMabX7uRKQ+MZwDehiiVbPqJY9s0p1io5sPQjIF0fDKFXx1mo8LQV6
jXpd6mrfGgIQDIlBvb0Y+NgloIsyVD30hAzYZw5UyCqV7++O2/CZqcKtJj5fNWpnvVYC8QG67HFGUTp5mQTTSaaTptI9dgWoWd6q
XZsPae/6Lgt91uE0GhAAzdEtycq8MikeeClAex+LRRcAmkJ/hv4RVaLqa5wqEYXfAgeY4ypRZWFGhqoQx+pEKpkRAEGQeQkRApbI
lkF3UOcF4ERuVGMlpMKp0FSo0dHHHATfe1YaMvM14Gr9VWFwMyv0EUAAE9kQ2BVViEqgBbkhQipSYpjwZdIh0hHSGclsJz329ZGV
xHY5CP2gWPmVKQY8wCHUWCoSky4e9AusApvAo1Ip5mMeEzuuviEzRCZ5Wa1IBx51CYwB1EoucnBgEaqEaqFRVCb7sE38Va/uYfcS
dtfWUd1R2bIAtNAD6ZhKaQya5QAgnwur/nASxbVUdrR2iSQgBOrANQCa01f/tcKf1eyuH8/M2IIPIofnBqBZ+m/S5V1OJiFA4A5c
BCBIgoIIYDAXq6nsEXxlIBfrLE8l87w9Pr3JfEOGyHCAU8VMGqPGgNn1YlQ34v88CqLzixc17wO081WoGb3Eh0r0dzsxUCpLLLWH
rPCM4fNVIGqXJEFrf/4Lk6CL6sv6RCn83MVwAKJKi+wZiBF0ujEQTW4Lc6lUZDGJuRUKhX18HfcUkcUvb9OY7HcpVgD1qw8Ut1JW
LltuephxgKHKi6wGrNRqq+UmC1+J4m5U7Iw2162TW83Es8/np0BXxaLOVQU1uhplzS3X0S7Unc6CVWeRKIqHLY7jJFwCkomsBxxA
jSVgR26GLHwsI50+GHB1JegLVVr1uIHClqvXybLVex+4+V12XwroSGsVu2Jm91LdFxsw5Crdm4xkxQNdUW7OEt0m7QOsQVXboBGl
QizkcVNEf/ecKtUNDrtJ4f4hrZZbJpM5O6DwoZoV3yQEYIfNrL+cW7vIRNYyIQ3pgGKptBqVU+UGE4tYZSjAEFWH9Zg6WRyTkNv6
98wjK8g5K1ANXqTQKEjjArTyZewTKpefAvm+uzfqVI/aVVE5F7dAluKuWllymBmA6vDKa4+MAaBAbpZCOLEcOZgttOuR95CR21zW
F4ovnTrf3Fv1UqG7Q3yIcc8lTGYUJ5NxKI5dx9zWX1fhfFZD6LeKghu9H0a38sq7nMBhNVIC1UJ3iG2apSSc/lC2mHYBwaFZD9Nk
a5B7bnLMezytrw2JARMDSwyZzgTL+WzagVim+gc/IAKpTq9kVgcJTWAMQLukHwpbLhtvWke5bn+D7kVC/Zi8XihgNhWe7ziZQmGx
L/woCD3odC/wHZW96jEaSuz1yuKwfwmKIgLk3umR3RxwkaMzLyGYInqv+PfflshCKlDJyqzo1/v7kgXM6GJGu1tWpVPDUN6iin3h
U/UbGNAgKvYYowCdV7p1y6cnVcy12czGMpECh8SJDpEAc0RCeuV8qXUGKeCmIbf1ShRlAVHPZceuBcNUKoHw5HazjF+pT7YKn6s3
AE0qpTELkE+yDteuKn5WnZIBCIY6B5TOClhr8pvKBEOaEccIxtAsQ2qE+1IZgJ5gwegCxVqtIjdR+EozDhphD5pAJbALUDAfJPVK
sbxqUgmAQsYhTnqTI4YqRAWoKOpgvR9wSy6AKkuoJyohFoguApyWSm+f6v7chxqoTuYGIEq8GcdIDbjRFOEsBAoIiBr9FEeoqIbV
P7go4hwAtE4dczYlMVnvmcjSTJyDp0LJggkrGTJcO4YAFVYjAIZIJ9E2BnGXxjFXB1bxVFxowmpdUUjquLqjsUUsU1BPoosCJyaO
k10FMEfLZ9GJdQ7ofQXLeJq1/Z2yPBsJjFExesXoPnfzKuLeKSYiINxKu7pFUtObuoXTeQlvxco8b40/XqojJF0hytYlbVak5eWe
0URlO5tImXjyFuMYEGVZ0HeTyFLn7bjSZACQEqo0cGIJRK1DRXMbzSkc+6DCtZrIPSvzjW1kG3N1BBXFEq+HQ54EEQ/oYoM4ndcX
WVrA/j0UkBIcTrkbF5RllPZmHe0oXKPwrMouAtz2CVlg4FdtMBgg+Ip6SCPKCsUXLtoJqSyLNlETWS6NSpMWW060EVJRFq2ipq4e
zIrF1v1cdjOK0aEPFGVQtGJnZ2lhQaL7bC2Pi2ZTZWT5ll6JHiFaLI+HjKZKfbdUg4ucqTcQMpgs6QW+JrFFRfGywrZrArVWmONj
2iw7bdaN3CtYLcRaHnsrqv2gOp9VcYMzZf+lECFBuJbgRTBSFzzEW0FcCpUiD25XrtcFN3oDSAzuWBYM/Hqil8+EMxi6TmdnN6gg
4OR364IJwcmUChFfj4IEzhBc1J8IDu4goR4FLnCb3oBgZ0opH/M1KcXYW9B22Yd8fIJZXMwDF9UJ5m4YvdSi3Yk0Q5a+uRmXgywQ
QRHFLJCl9FAONyfZMUj2y8Z7Gu9o+GRmtXm+HFEO1D5uwwdDjXf+mlX5EznSvFj7BL/ejC1qQpD3WNCCs5Ql0BJUnRlLzai7ABZE
y4C0TfprR8ODDW/Mzl9e6F4OqX36v3Y0PtDw7uwCFwrmMfOQbfPqDTS+MieHu6wsrXaRTnGkHJVNQnK9iz2PKDHFyANtlPooPZvX
DgXaXjsJCCeywUjLYzE6WbsclxGgrpB04l6v3uWmp7ve2mJoR5KdOpKkXUlal6RthlPro05I09Ab6AJsNxqw66RtSdoEaT+eFt46
L7R9ul/hrbOCz+02UNqW0CbEDqOHKxXaPS20c7pSweeO8z/b/TeNiYkJ0w464ToR+/v/dwHP6A34uJ0ATTQrplDGKDfFx6OeLxt9
bdBEnwcZIyQde3AvR69Aj1zzA+51g/BK/I4+iuEZw2gghjY4GEcdtOvEC8AL0r0i8Qfm+YL/chg+hxHyhiqTZOJOjV5BPQhJ2POa
N/w/uxVyGGXecGWSSop6QdVdpL7BPmu3g0affRwDPLaou0VfdZXb/CBxEQPe7ipdN7UeMt8yug99coSpUYw62pAixPdDmjpaJpPq
r579tIws3s9zP+3G4iIkrki6oWfW9Eq/UPqGsbxlOsVUEYOenX5gxpruXpVXMLa1TB9+LULVXaSfA7YdDlj7bWF3YyG9o2fbYOJf
kkw9sFol8wsdMt/Qhs/S+fIcULyWyciMZl6feqX7T/sBEpEl9eCENoLcM3LPUqY3S6ewGgWgqJC0/dMNHnznYLPU4QoN+Qb6MSuE
hRqA1biw9xrjElzAKVvbO2Ns5V7x4DtTT/1gccoM0trC/31f8+2o7qpfu+dddpbSMiKRCXqcACVlls70W1nXMmpARQlIvOn5NQ3y
7TSDU0Yrj1HZTum4KygzsFrI1v4Z4GYHKhJP7qxcq6GYOU996mdVqnb1XYllxoDFPHWoccyDCregLKfG2e/dGWOqaYazuhGiYjIw
LiGjR6mHnnI1L/LVrfLUr9L1zyPungokReRYwuywTN7jMCWsQa3uwqa/dU91K+kSkLq0/ZMbOA6dVjREFRxZGicFXOGiiCpRJLWM
8m1JofX4+IBQwjfnY0pA4FPVe9KWvYxGwf696ypXx/NvQG0ydyFDdhOOuCk4sswpJS4BHLgAzgGQ84wxhgXDijqTZO6mru8gFYEp
66xZacXfpsksnWgqx3AYAu84ThJLV0SQeSVpiYkVEgFgSMduAKmDCXEz5/8zLYeFQj9VBNx2TKk/FID4mibz1IEGbXjutqDyz4x0
jp66MhYUMcYBAkA8Yn3vxKSAV2MRfRDNAS0yPCMGIKQ8YvvDlTDRlwGHKF95THxNeTIq0IQ2HGcrqOyz96qNnupth6sf/5dQlJQM
qwfkhhxrpRq4FSK4m1zSBnEBMPRYclIMnmgJt6lP5nDLEH2ibHlLdlyrSmmDF1u8ish3cDZUeG5Uor8uXKv+L8q8wwqrVN0aH8C7
cXoQJmevV97XrDJm/xquTsU1CQEZAqCOp1enpn/0Xblxu2K1iT7NZAxYXuA03G6P3ekcluprQIuKkQf+fznR3KorlMfn+68CK8Yo
CLSkGpp8IPONzLek5bfnNr6vKl7c6wjb940SciTmvytJ1A7TSAaEACNd38wuVM1MWsE2XJoncw3OjoZid9UdVq9yUpNKyrlkptcF
EUZACDjhHIFv9iyS2UknWSRTw1zzknc4XXlwXbQSWyN2yNm18223DHQhioy/wgrjX8JzGQwo9IiqldjayhrKhoZdvqVzgpAxhoEH
QDZUh+ff4XALHIAT6t625V3ZAHvXoQJmfq/OV3X1q5eNCIGrM6Pr4AMcuglKTABMBRj+Jmf17VH7Lr8qPzqp0pStPIzHQMbofZzn
DxjjootWVzgHv73x73Wzc/VcLzpS9cheh8XAxfy6HnAVIBF8pdVVfpYWyjZ5DxVCo1t+4jGQZ1i+jlSpHayqRB1mnZtExiAiM96z
MTMafi4aEzWzrEv0HPDRtqO9bjqojIE/wz3kh6vWmbamYyFTVYXJaPBh4sBlS9Jtek6K6l0rBQ39DMrZ7ZGb1v2pTY1nei3lvYFJ
moHJBzK/kPmWtB33ub4nwCuTycx88VgOq2T1V0LgLs7YPbcLVR98YBupzJMJ5wdB/hV9dOc1MmTsrBf+Ry7iWWCCCcReGafvm6ce
mh6Rvn8W6GYPUYTzpLioFV0jBkRZDo2HDfbSfbY1BKrWfEOdS3ie4wDNiVGJLl0KW/9Y11avegBaAjvT0oTnPkCXd1ZQq6yA8vh6
xmrbYjzqo0gA0CjyMUzZJI/vheUQADAEQAQwYLSkxflARMOOGqaG+VUQM9UoCLvbgerHAnYgBsD8dAbMpJqYjoEwXZlUw1GnRERa
kVek13ZRN6ghUfadxF8JV+y1qmAj47pSQ28kHm+IkBuJfkQWxe8EfycEuQZupKvvBCsGAib2XeOrcWL4N4bohC3qXrgmXfO+Lfa2
+NvCbzPLyneeJLq/fH1AVMpZY6s288jmmk2/Rl/ZYCfgdszBOqvx1TWHv4ojGAvtgoYfd423dDQehI4D8U1CQwBK3mspt+Y3olnd
oon7wY2oX6SjuWUZ8XMGePvljTOa/iOOGt++aOYWzuhsc7sR0WJunjDuc3Trr/5cwplq8QYBNcqrB9Cy/i37ZlvWbTiZAahRV09E
TW5EC/u3oNsyb3R/WsjNGxKio00xo6Y3LhaadMeJNXNYOS3Idv9lBXVaHbyWDragIH9oqAlBqv9jtnc4JXwaWbNbENmoVvcdj3cA
FK5nsczplqvW+fLktiLr0m50mO94reDPU1t6RsUsd92jbHdOGpYGjbJUANN3xpyRjYJoXwfap7zxR5ajB/ZVU8vI2jgaUFl65CPL
whiBIbQ2HL8CVFSXLQsVxX4aam/TfmqVLTt+Pw059D6Zy1rxYHT/Bd4sNKfsWwjIp7Sl34jyYXRLBf47c4rmN5pD9iyIrUXdLJYl
3LLWqNd2HdeoXumxpBvR9VyFAHCr+QQ1u833ZqFntNIsJt4o2N5G4D5KTbNlRjwb5TrQjYIxIkBoSlQW/oWom5exlq+irQ+so5nB
vw6sRi0PG+U6NSNvX0bavbCM54a59gr8RlzhacH9o/mN5hdrIavl389iTB7ORtSJjWKZ4QpSma5hdpz7YXh35DOu3Y1av2HM6zVa
G3sMu9R4RSpzXoxFFAPRsF3WYERvOJDKwtR12r02QZZlRjpTTGEqr1X6e+hl7nqlTc9qpW2oeq5vac0m0s7Fmu/wPMYRbeNK5m+m
v+GYD6S1mxhrQ5xQgitY6MyKaF+HvSTRWG5P8njHE212MdH0Zll0tMtTXrCzKK/RviJZuB1Mkda4ytrkZsvQnZUBgHGdaZASqbBP
I5iCJxREc8yiYtZBNI2yo20DQNiNnulDatM5BkF6ibjQytnyyWhB5Ti5v8nRupDWradWOoY1vHU0h5UoAEWUGGdpbBsOCteoLE/r
7oMrlK61D21jLQ9hNQJUjnyuIWewZUWzG/VQI9Hki2h+N7cCUA804czn1YgWlE0922rB5aloa5+q3ZrOjWgh2ctvjwDwQ+3iH8Ib
AZVYWuSNqFGNQPWKLEu9JSOXhf7TsGyk7QtGr5SOy9gzN8p1nYwz7bT9M7vI1j5UZW99dfgl6OU+bXu5V5C2S9c/t4ol1rEsc9eb
ynEKz3RnmO9hQt7hwsN2Uicj6xJuiSzK6rEoiH5rvxgqXja9jHQ091tvCHV/3+6sM6K5ZJrihui9nOxLDV8/WaxDRbFcZgmi/itY
1YPxE0/fUoYnLZ+ubyIIKCqij/Jv42w3wdztkDFXR0LWhIrqhLN3rOYGyNn8IWBay3fJr+lzn3Rc3/ahMvTJUdmGfDzXEOZ6I5p0
NqoRY9pi/1lpgnWjhwHbnXn9U8wWdvMliOiV9mae9Y1ozqUR6iRXJtWHKPuUMqLZh73Mdu9zCKKZlsNiejBJ7OWhTSRIqwuV7fPX
rFIOT6JuDYJd7iMC3iq+ILUmuF6bueudNC2jhziLw7JC1mgIQIKf71mkUqtEZvrUrbsRD5jDygurhTaxjRyWUjPJiBYesUpBo6IE
h8fGHcpmk4mFZWSz3D9tjXZNJyNacBaTYmhAzhmqSOs7N1LcJSSiJd2a3ejVyZDLktdEb2RZ9r8MkJm/vIPKMH3q569dPad0rMJ1
lOmPIR2fp3ERfW34KSZDFGjmafiYmjbL8AfnmuMBkwVV4rCYABcG8f0PG0dGNkvbxmGztfV5lS4Q9cXCbnKTViZjaLCdtj2OaP7P
Rr1kI1zjgjVo2fEYFzAAJj9ljFmjqr/+C1HWzSnI0rZx5mPUYqNa9Z2H+QkhSHDWrcu/gyri76GXdPI+X92CMt0xpOnlnh3ANrxN
MnmHyFEcaOVpWKgcQxtVuINzjZ7f+KWfjTjB6IEF3shCXAXRFLKIaKv7WvTtcxdIBsVrxVy/cV2xa6lhtS53JKsbo/vTGG+jtCuo
D9vnBWecER26I8Noa4YZBRFXTYLmlLkW+kXd3377DRRYlyAcUtZX6a/4f9M8kcZYHaTvnSkd0+BB0nYBb3dxQs7mBWkXjF5GWx/U
EkADIAEixPmBVTSzimUm2X4/XGhR7Vbp8unncZt5zCpmlc06a7OIzGo9dGI6UAkQ0eAjXyGhoaOwTmAJolPbyWst5GYJWsats27U
KzyqjUq9WJOzKuvYsmDN+vGwVlrArc0tBtxCm/kHWSy5K5TpAgR/g97R4y/x+B3IS4ACRKyZdOIhzcBTht5p4Ot9b5V/sPIHqxL0
ak9WvZF7lrJCqGjuq/QnzfCUNhylu/T9L+NjQsoaZ9o9UWawjZRW8Wx8LgwSFjzxjCyviF8b4mVheb0UtaUXhbtj7kHDWiaafbgw
lUaD7NA2Oglzf3U4Ft2/2c0yiepVv5/MzWqjQyeGMfL7VslmzVuR4GKbbD86fcxqo7VGMpMODwgYsSq6bYsu1uBiE+69SCgkENnb
PlhtaxOqMvROjfJdodQgUIOdH8iqC4VrScsGud7zWfYbfv+E1UunvuMsjf3ed/w99JKPPakc2yBny8eJOCcHKihdeFE3S7Pb1CLc
84XosVzfk+nYTW7W6WQOy3nlw9JuZNX94JesCI/FBhvvjbgf7CC6XCc12D6dNLuwWueRRV8qa83QzeYTVut8bzEKPNeNnnnEoB1L
KOZHqDBP+C1dzUTXQ3A1RtfDgr8gSnXo04wxygdJy8iqV1k737K1vcn9Yh2bSJV2YIpx8AA0U/crq2/lbkEZXl69UXjFoyrbnmcd
PeoK7+OwpwpDjbKORURfP5Dq6yFiBhlR40hXmhWaq65MXOYaxauDCvheBzxy0RlxNM/sIhuKuRzXLohog98xAYApso0bRFSmrniF
79TAQw8dYYFUd/b4ly424WIDnr8JkDAaKiPKxsVEjBFs9GAX3s7OaoE2XPreieR/BxPyNKaIUa6bsXs6CGMGVg82kcDqbqLw+2ew
IgDRlf0P8nWBCviRhnuHwk0vHQqvQiqMfSpcDx8viv0XQ0JrYxLbF6xq8eo3TyJrEGWAEaqvH7f3AyVE79kEahox9RpPdJOFhajv
lZKtIqugFTHuc61RSqQXPOIF9xCYlHzFQ8RGLiIh4ESCEwSAl9LnV/r6IoQCginPr2Dy3Au/fx538Trx1EPCmdvo65chzhKcAxN4
R/0UiRMsknns5UHi6fuEs7cR9y6NMz1YqR/qUv+nmgsPvLUoxocFUpbZXulT1kdEs8sYSSfBjwiJkBBUtsUoCA4/UGSL+KEKHHm8
NBYVUWxDG4pMdgAVAAWOOa7Rf0aEREN7BHR5OTo6OvB8XQAUPFOQIQDwz/MAFZUBFde1VoDVGiqqETUeAEi1XrQgal1WvRCYLALg
AF3i/4gXDm16ltb6ecQXLLF5CN6dzwtDyl4YY3574ESIXUdVVT9bVZB7O5geV1AiYwdDP83XJWQSgjqvHLXCa6dhl/w6Mq1jicB6
EcWyfr/o5daRD7Qc0teSIHlLP+QXetJrOWIqhyjLR3IUuQhJG8i52NBk35im/H2shkbDOca935BcODLFonGB2RVsy/XwlOYjbYiP
ZUHY2TSAGDbIrDfQOTq39wficZWRQ4RS2s95smaSTz3xw9tWHT5UXYknEoKTQrDxC5PsNhwK1nn9+8brXKdgZBfBDO7cMobFhnch
aDp3x6hmpBq6SNaXsoHtGEAUerh7KJ7OdaO6kXekWta/2WaMKfW3Dr+76hc0bhzU+n4RXHhJPEnWgkTY0METYec65SMjWexgVz8R
15aCsB0DDJy7Y2xjpCcsEet9Dt4+kNeTnsCC2hruWV15ELhobUSU/WSZoor6i1TW4vDVsWIIqME2bLJnpNMoCFKnMCQbDMdVzmBx
Y7vTn59YAwJ3HmBRQ9r3bQXKIn7ukwUYPmTiKylRd2gxXLq6EkGA2hcrhvEQ2D/W+jKwtTb8u1tXMAK6ogd9fIfBjHQSM9gShoJz
d2jlQIedH9ej2+OJud1IT7oOcNBteH43ujLQg4pIP/o0GQqBfRBIvjGQa2huoLpSghhxLObTZGBRRcTBOm5433hWNp/Q4i4owYR6
pt4FTdsb4LYZ1J27OOrrgIpXOOuLOYCQObE7z3B3sXk41416O9SSWk37Y0wH41MP95VHTj0SW3sHYfs7y7LOxc0uGQE5gw2oHFne
culVbDa0y3NcGwCBO1I6Ur710uthgZQUT3ub9WachE8Qn8AFOtcNJjmkIzS/I/tJ46KK+EupJQNtSq8yqBDXdN0z0nbzT6mTtWUG
Ex0G0m69Pu7YvBphMVe50+lmToOULvhpTl67MwriK1hUbTnDXIeObD2xYujVjbDVXfvtyarSLpWm5k6oe/ow0m4jBFIXGeRVurHU
2HwL9cF8rgVjHXu4ayz2Kj2u6W2ORqiA+D2nqX3IS1jfg5KadCPGMsogyjYSNz61XFRhQ6zt++8bdwNza0YzCg4j0qma0ebegG1T
gprKvcEU80dyhLA968O2tW6EdIfF3T7m8LrRXR+YQRXE33maWjHKQfoiR/Y1gdpK5v065w6GNG0mHm+T+jMK1v4iCCo2V2hYMYnO
CnaxZqiFgA4pejYzvzgM7g5KLOlNDl5hN0J8jC3Adh4/O+cQOGCNsP8VQ5VRG1LJW+znwXS8w7lClNhWnh5BWFq9GzmOsAkGcec6
S3qVRXGC3pAFVNitDIGdfVITiMUYt1FSiR8/sUoHswU/4rBwEU26EcsxDgMvrJOEP/ZTYwjW35/emD+rc511fb6CQxO2k+Ol1Kgb
yXAL/2iwsiQDzuuDfdg16Ebgrh5sdrb20DyM7oqfL1pCHuibGCxmbENPfPfxJkbIrMo4eCL8g7Ps04KiyrO1HnzgapBlXQbA4O9H
S690KWa4jX1PqgZ3I1odKeIKgnv6ibA9wQaZePu1Dx/MZ3ej+/Pn35PzU0FZTboRizN+SZz//4sqyMMZn+N3jO7hiIZK7jFz5AAA
)---"

.qc_logo_light_b64 <- r"---(
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
)---"
