############################################################
# Interactive isotope tracing figure app
# (pipeline logic from "# metabolic figure pipeline V2.r")
#
# Open this file in RStudio and click "Run App" (or run
# shiny::runApp() from this folder). Then:
#   1. Upload an input workbook (.xlsx) and, if it holds several tables
#      (e.g. raw peak areas + fractional incorporation), pick one
#   2. Pick metabolites and conditions — the order you pick them in is
#      the left-to-right order in the figure (remove + re-add to move one)
#   3. Download the figure as PDF / PNG, or the summary table as CSV
############################################################

library(shiny)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)

############################################################
# Defaults — edit to match your experiments
############################################################

# Short code used in the sample names -> display label, and default order.
# Conditions found in a file but missing here are still offered, using the
# code from the sample name as their label.
conditions <- tribble(
    ~code     , ~label                , ~order ,
    "UL"      , "Unlabelled"          ,      1 ,
    "N"       , "Normoxia"            ,      2 ,
    "N MTX"   , "Normoxia + MTX"      ,      3 ,
    "H4"      , "Hypoxia 4h"          ,      4 ,
    "H4 MTX"  , "Hypoxia 4h + MTX"    ,      5 ,
    "H6"      , "Hypoxia 6h"          ,      6 ,
    "H6 MTX"  , "Hypoxia 6h + MTX"    ,      7 ,
    "H48"     , "Hypoxia 48h"         ,      9 ,
    "H48 MTX" , "Hypoxia 48h + MTX"   ,     11 ,
    "0.4 mM"  , "0.4mM serine"        ,     12 ,
    "0.8 mM"  , "0.8mM serine"        ,     13 ,
    "UL KO"   , "Unlabelled HIF1B KO" ,     14 ,
    "N KO"    , "Normoxia KO"         ,     15 ,
    "H48 KO"  , "Hypoxia 48h KO"      ,     16
)

# Preferred facet order; metabolites not listed here go after these
METABOLITE_ORDER <- c(
    "AMP", "ADP", "ATP",
    "Serine", "dTMP", "dTDP", "dTTP",
    "SAM", "SAH", "Methionine"
)

# Isotopologue fill colours; heavier isotopologues get darker shades
iso_colors <- c(
    "M+0" = "#FFFFFF",
    "M+1" = "#F4ABAB",
    "M+2" = "#E03131",
    "M+3" = "#7A0000"
)

CONDITION_MAP <- setNames(conditions$label, conditions$code)
condition_order <- conditions$label[order(conditions$order)]

############################################################
# Pipeline
############################################################

# Find every table in a workbook. A table is a header row containing
# isotope columns ("...M+0") followed by sample rows up to the first blank
# row, so instrument exports with several stacked tables (raw peak areas,
# "Fractional incorporation", ...) are split up and offered separately.
find_tables <- function(path) {
    sheets <- excel_sheets(path)
    tables <- list()

    for (sheet in sheets) {
        cells <- as.data.frame(read_excel(
            path, sheet = sheet, col_names = FALSE, col_types = "text",
            .name_repair = "unique_quiet"
        ))
        if (nrow(cells) == 0) next
        first <- cells[[1]]
        is_header <- apply(cells, 1, function(r) any(grepl("M\\+[0-9]+\\s*$", r)))

        for (h in which(is_header)) {
            end <- h
            while (end < nrow(cells) && !is.na(first[end + 1]) && !is_header[end + 1]) {
                end <- end + 1
            }
            if (end == h) next

            header <- unlist(cells[h, ], use.names = FALSE)
            keep <- !is.na(header)
            keep[1] <- TRUE
            block <- cells[(h + 1):end, keep, drop = FALSE]
            colnames(block) <- header[keep]
            colnames(block)[1] <- "Sample"

            # Title = nearest text-only row above the header (blank rows skipped)
            title <- NA
            for (r in rev(seq_len(h - 1))) {
                filled <- sum(!is.na(unlist(cells[r, ])))
                if (filled == 0) next
                if (filled == 1 && !is.na(first[r])) title <- first[r]
                break
            }

            block <- normalise_columns(block)
            if (is.na(title)) {
                title <- if (isTRUE(attr(block, "from_areas"))) "Peak areas" else "Table"
            }
            label <- paste0(title, " (rows ", h, "–", end, ")")
            if (length(sheets) > 1) label <- paste0(sheet, ": ", label)

            tables[[length(tables) + 1]] <- list(label = label, data = block)
        }
    }
    tables
}

# Bring a table's columns to the "%<metabolite>_M+n" / "Total <metabolite>" form
normalise_columns <- function(block) {
    # Repair a "%" typed in the wrong place (e.g. "d%TTP_M+2" -> "%dTTP_M+2")
    misplaced <- grepl("^[^%].*%.*M\\+[0-9]+", colnames(block))
    colnames(block)[misplaced] <- paste0("%", str_remove_all(colnames(block)[misplaced], "%"))

    # Isotope columns without "%" hold raw peak areas. Mark them as "%" columns
    # (renormalisation later turns areas into % of total) and add a
    # "Total <metabolite>" column as the summed area if there isn't one.
    area_cols <- grep("^[^%].*M\\+[0-9]+\\s*$", colnames(block), value = TRUE)
    area_cols <- area_cols[!grepl("^Total", area_cols)]
    if (length(area_cols) > 0) {
        areas <- sapply(block[area_cols], function(x) suppressWarnings(as.numeric(x)))
        areas <- matrix(areas, ncol = length(area_cols))
        met <- str_remove(area_cols, "[_ ]*M\\+[0-9]+\\s*$")
        for (m in unique(met)) {
            total_col <- paste("Total", m)
            if (!total_col %in% colnames(block)) {
                block[[total_col]] <- rowSums(areas[, met == m, drop = FALSE], na.rm = TRUE)
            }
        }
        colnames(block)[match(area_cols, colnames(block))] <- paste0("%", area_cols)
        attr(block, "from_areas") <- TRUE
    }
    block
}

# Turn one table into long isotope and total-abundance tables.
# Metabolites are read from the column names ("%AMP_M+0", "Total AMP").
load_tracing <- function(data) {
    data <- data %>% filter(!is.na(Sample))

    # Isotope columns start with "%"; "NF" (not found) is treated as NA
    iso_cols <- grep("^%.*M\\+[0-9]+", colnames(data), value = TRUE)
    data[iso_cols] <- lapply(data[iso_cols], function(x) {
        x <- as.character(x)
        x[x == "NF"] <- NA
        as.numeric(x)
    })

    # "1. Normoxia 2" -> "Normoxia" -> display label
    data$Condition <- data$Sample %>%
        sub("^[0-9]+\\. ", "", .) %>%
        sub(" [0-9]+$", "", .) %>%
        recode(!!!CONDITION_MAP)

    iso <- data %>%
        pivot_longer(cols = all_of(iso_cols), names_to = "variable", values_to = "fraction") %>%
        mutate(
            metabolite = variable %>%
                str_remove("^%\\s*") %>%
                str_remove("[_ ]*M\\+[0-9]+\\s*$"),
            isotope = str_extract(variable, "M\\+[0-9]+")
        ) %>%
        select(Sample, Condition, metabolite, isotope, fraction) %>%
        # Within each sample x metabolite, rescale detected isotopes so they
        # sum to 100. NF (NA) isotopes remain NA and contribute nothing.
        group_by(Sample, metabolite) %>%
        mutate(fraction = fraction / sum(fraction, na.rm = TRUE) * 100) %>%
        ungroup()

    total <- data %>%
        select(Sample, Condition, starts_with("Total")) %>%
        pivot_longer(cols = starts_with("Total"), names_to = "metabolite", values_to = "abundance") %>%
        mutate(
            metabolite = str_remove(metabolite, "^Total\\s*"),
            abundance = suppressWarnings(as.numeric(abundance))
        ) %>%
        filter(!is.na(abundance))

    metabolites <- unique(c(iso$metabolite, total$metabolite))
    metabolites <- c(intersect(METABOLITE_ORDER, metabolites), setdiff(metabolites, METABOLITE_ORDER))

    found <- unique(data$Condition)
    conds <- c(intersect(condition_order, found), sort(setdiff(found, condition_order)))

    list(iso = iso, total = total, metabolites = metabolites, conditions = conds)
}

# Mean ± SEM per condition x metabolite x isotope, plus stack positions
summarise_isotopes <- function(iso, metabolites, conds) {
    isotope_levels <- unique(iso$isotope)
    isotope_levels <- isotope_levels[order(as.integer(str_extract(isotope_levels, "[0-9]+$")))]

    iso %>%
        filter(metabolite %in% metabolites, Condition %in% conds) %>%
        group_by(Condition, metabolite, isotope) %>%
        summarise(
            mean_fraction = mean(fraction, na.rm = TRUE),
            sem = sd(fraction, na.rm = TRUE) / sqrt(sum(!is.na(fraction))),
            n = sum(!is.na(fraction)),
            .groups = "drop"
        ) %>%
        mutate(
            isotope = factor(isotope, levels = isotope_levels),
            metabolite = factor(metabolite, levels = metabolites),
            Condition = factor(Condition, levels = conds)
        ) %>%
        arrange(metabolite, Condition, isotope) %>%
        group_by(metabolite, Condition) %>%
        mutate(cumsum_top = cumsum(replace_na(mean_fraction, 0))) %>%
        ungroup()
}

isotope_palette <- function(isotopes) {
    extra <- setdiff(isotopes, names(iso_colors))
    if (length(extra) == 0) {
        return(iso_colors)
    }
    ramp <- colorRampPalette(c(iso_colors[length(iso_colors)], "#1A0000"))
    c(iso_colors, setNames(ramp(length(extra) + 1)[-1], extra))
}

build_figure <- function(iso_summary, total, metabolites, conds, opts) {
    p_iso <- ggplot(iso_summary, aes(x = Condition, y = mean_fraction, fill = isotope)) +
        geom_bar(
            stat = "identity",
            color = "black",
            width = 0.7,
            position = position_stack(reverse = TRUE),
            na.rm = TRUE
        ) +
        facet_wrap(~metabolite, nrow = 1) +
        scale_fill_manual(values = isotope_palette(levels(iso_summary$isotope))) +
        labs(y = "Fraction of total", x = "") +
        theme_classic(base_size = opts$base_size) +
        theme(
            strip.text = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_blank(),
            legend.position = opts$legend
        )

    if (opts$error_bars) {
        p_iso <- p_iso +
            geom_errorbar(
                aes(y = cumsum_top, ymin = cumsum_top - sem, ymax = cumsum_top + sem),
                width = 0.18,
                linewidth = 0.55,
                colour = "black",
                na.rm = TRUE
            )
    }
    if (opts$legend == "bottom") {
        p_iso <- p_iso + guides(fill = guide_legend(nrow = 1))
    }

    total_sub <- total %>%
        filter(metabolite %in% metabolites, Condition %in% conds) %>%
        mutate(
            metabolite = factor(metabolite, levels = metabolites),
            Condition = factor(Condition, levels = conds)
        )

    if (!opts$show_total || nrow(total_sub) == 0) {
        return(p_iso)
    }

    p_total <- ggplot(total_sub, aes(x = Condition, y = abundance)) +
        geom_boxplot(
            width = 0.55,
            fill = "grey88",
            outlier.shape = NA,
            linewidth = 0.45
        ) +
        geom_jitter(width = 0.12, size = 1.6, colour = "black") +
        facet_wrap(~metabolite, nrow = 1, scales = "free_y") +
        labs(y = "Total abundance", x = NULL) +
        theme_classic(base_size = opts$base_size - 1) +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = opts$base_size - 1),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.line = element_line(colour = "black")
        )

    p_iso / p_total + plot_layout(heights = c(1.2, 1))
}

############################################################
# UI
############################################################

ui <- fluidPage(
    titlePanel("Isotope tracing figure"),
    sidebarLayout(
        sidebarPanel(
            width = 3,
            fileInput("file", "Input workbook (.xlsx)", accept = c(".xlsx", ".xls")),
            selectInput("table", "Table in workbook", choices = NULL),
            helpText("Raw peak-area tables are converted to % of each sample's total."),
            selectizeInput(
                "metabolites", "Metabolites (in facet order)",
                choices = NULL, multiple = TRUE,
                options = list(plugins = list("remove_button"))
            ),
            selectizeInput(
                "conditions", "Conditions (in x-axis order)",
                choices = NULL, multiple = TRUE,
                options = list(plugins = list("remove_button"))
            ),
            actionButton("all_conditions", "Select all conditions", class = "btn-sm"),
            hr(),
            checkboxInput("show_total", "Show total-abundance panel", TRUE),
            checkboxInput("error_bars", "Show error bars (SEM)", TRUE),
            radioButtons("legend", "Legend", c("Bottom" = "bottom", "Right" = "right"), inline = TRUE),
            sliderInput("base_size", "Font size", min = 8, max = 20, value = 14),
            hr(),
            fluidRow(
                column(6, numericInput("width", "Width (in)", 10, min = 3, max = 30)),
                column(6, numericInput("height", "Height (in)", 8, min = 3, max = 30))
            ),
            downloadButton("download_pdf", "PDF"),
            downloadButton("download_png", "PNG (600 dpi)"),
            downloadButton("download_csv", "Summary CSV")
        ),
        mainPanel(
            width = 9,
            textOutput("status"),
            plotOutput("figure", height = "750px")
        )
    )
)

############################################################
# Server
############################################################

server <- function(input, output, session) {
    tables <- reactive({
        req(input$file)
        found <- tryCatch(
            find_tables(input$file$datapath),
            error = function(e) {
                validate(need(FALSE, paste("Could not read this file:", conditionMessage(e))))
            }
        )
        validate(need(length(found) > 0, "No table with isotope columns (M+0, M+1, ...) found in this file."))
        found
    })

    # New file -> list its tables, defaulting to the last (final) one
    observeEvent(tables(), {
        t <- tables()
        choices <- setNames(seq_along(t), vapply(t, `[[`, "", "label"))
        updateSelectInput(session, "table", choices = choices, selected = length(t))
    })

    dataset <- reactive({
        req(input$table)
        i <- as.integer(input$table)
        req(i <= length(tables()))
        load_tracing(tables()[[i]]$data)
    })

    # New file -> offer what it contains, preselecting everything
    observeEvent(dataset(), {
        d <- dataset()
        updateSelectizeInput(session, "metabolites", choices = d$metabolites, selected = d$metabolites)
        updateSelectizeInput(session, "conditions", choices = d$conditions, selected = d$conditions)
    })

    observeEvent(input$all_conditions, {
        req(dataset())
        updateSelectizeInput(session, "conditions", selected = dataset()$conditions)
    })

    iso_summary <- reactive({
        req(dataset(), input$metabolites, input$conditions)
        summarise_isotopes(dataset()$iso, input$metabolites, input$conditions)
    })

    figure <- reactive({
        validate(need(nrow(iso_summary()) > 0, "No isotope data for this selection."))
        build_figure(
            iso_summary(),
            dataset()$total,
            input$metabolites,
            input$conditions,
            list(
                show_total = input$show_total,
                error_bars = input$error_bars,
                legend = input$legend,
                base_size = input$base_size
            )
        )
    })

    output$status <- renderText({
        if (is.null(input$file)) {
            return("Upload a workbook to start.")
        }
        d <- dataset()
        paste0(
            input$file$name, " — ", tables()[[as.integer(input$table)]]$label, ": ",
            n_distinct(d$iso$Sample), " samples, ",
            length(d$metabolites), " metabolites, ", length(d$conditions), " conditions"
        )
    })

    output$figure <- renderPlot(figure(), res = 96)

    file_stem <- function() {
        paste0(
            tools::file_path_sans_ext(input$file$name), "_",
            paste(input$metabolites, collapse = "_")
        )
    }

    output$download_pdf <- downloadHandler(
        filename = function() paste0(file_stem(), "_figure.pdf"),
        content = function(file) {
            ggsave(file, figure(), width = input$width, height = input$height, device = "pdf")
        }
    )

    output$download_png <- downloadHandler(
        filename = function() paste0(file_stem(), "_figure.png"),
        content = function(file) {
            ggsave(file, figure(), width = input$width, height = input$height, dpi = 600, device = "png")
        }
    )

    output$download_csv <- downloadHandler(
        filename = function() paste0(file_stem(), "_summary.csv"),
        content = function(file) {
            write.csv(select(iso_summary(), -cumsum_top), file, row.names = FALSE)
        }
    )
}

shinyApp(ui, server)
