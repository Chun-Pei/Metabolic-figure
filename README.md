# Metabolic Figure

R scripts that turn stable-isotope tracing data (LC-MS, exported from an
Excel workbook) into Prism-style figures: a stacked bar chart of
isotopologue fractions (M+0 … M+3) on top of a total-abundance boxplot,
faceted by metabolite and condition.

## Interactive app

**`app.R`** is a Shiny app built on the V2 pipeline. Open it in RStudio and
click **Run App**, or run `shiny::runApp()` from this folder (it needs
`install.packages("shiny")` once). Upload a workbook. If it holds several
stacked tables, as raw instrument exports do (peak areas, then
"Fractional incorporation", ...), choose one from **Table in workbook**.
A raw peak-area table is converted to % of each sample's total. Then pick the
metabolites and conditions to show. The order you pick them in is the
left-to-right order in the figure; to move one, remove it and add it
again. Metabolites are detected from the column names, so you don't need
to edit any list. Download the figure as PDF or 600-dpi PNG, or the
mean ± SEM table as CSV.

## Quick start

Use **`# metabolic figure pipeline - unified.r`** — it supersedes all the
other scripts in this folder. Its data pipeline is copied directly from
`# metabolic figure pipeline V2.r` (the most complete earlier version),
so it behaves exactly the same way, just split into 3 outputs. Open it
in RStudio and run top to bottom:

1. A file-picker dialog opens — select your input workbook (it must
   contain all the metabolites you want plotted; see below if your data
   is split across several files).
2. The script produces **three figures**, each an isotopologue panel +
   total-abundance panel pair, saved next to your input file:
   - `<name>_Purines_figure.pdf/.png` — AMP, ADP, ATP
   - `<name>_Thymidine_figure.pdf/.png` — Serine, dTMP, dTDP, dTTP
   - `<name>_SAM_SAH_figure.pdf/.png` — SAM, SAH, Methionine

   A figure is skipped automatically if none of its metabolites are
   present in the input file.

If your metabolites are spread across more than one workbook (as in
`# metabolic figure pipeline -2 plot.r`), merge them into a single sheet
first (e.g. with a `full_join(..., by = "Sample")` in R, or by copy-paste
in Excel) before running this script.

### Editing for your experiment

Everything you're likely to need to change lives in the `USER INPUT`
section at the top:

- **`conditions`** — maps the short codes in your sample names (e.g. `N`,
  `H48 MTX`) to display labels and controls their left-to-right order.
- **`KEEP_CONDITIONS`** — which conditions actually appear in the figures.
- **`PANELS`** — which metabolite columns go into the Purines / Thymidine
  / SAM_SAH figure. Serine is grouped with Thymidine by default because it
  feeds the folate one-carbon pool used for thymidylate synthesis; move it
  into `SAM_SAH` if you'd rather track it with the methionine cycle. Add
  or remove entries to get more or fewer figures; the entry name becomes
  the file suffix.
- **`iso_colors`**, **`FIG_WIDTH` / `FIG_HEIGHT` / `FIG_DPI`** — fill
  colours and output size.

The script tells you about input problems instead of silently dropping
data: it warns when a sample's condition isn't listed in `conditions`, or
when a `KEEP_CONDITIONS` entry doesn't match any label. It also repairs
isotope column names with a misplaced `%` (e.g. `d%TTP_M+2`).

Nothing below the `Auto-built` divider needs to change.

## Expected input format

An Excel workbook where:

- The first column is the sample name, formatted like `"1. Normoxia 2"`
  (a leading `N. ` index and a trailing replicate number are stripped to
  recover the condition — here, `Normoxia`).
- Isotopologue fraction columns start with `%`, e.g. `% AMP M+0`, and may
  contain `NF` ("not found") for undetected isotopes, which is treated as
  missing data. Within each sample × metabolite, detected isotopes are
  renormalised to sum to 100%.
- Total-abundance columns start with `Total `, e.g. `Total AMP`.

See the files in [summary/](summary/) for worked examples of input
workbooks and their output figures.

## Requirements

R with the following packages: `readxl`, `dplyr`, `tidyr`, `ggplot2`,
`stringr`, `patchwork`.

```r
install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "stringr", "patchwork"))
```

## Other scripts in this folder

These are earlier, single-purpose versions kept for reference; the
unified pipeline above covers everything they did:

| Script | What it did differently |
|---|---|
| `# metabolic figure.r` | Original version — hardcoded filenames, AMP/ADP/ATP only. |
| `# metabolic figure pipeline.r` | Generalised with a file picker and more metabolites, but no `NF` handling or condition filtering. |
| `# metabolic figure pipeline -text.r` | Same as V2 below, minus the bottom legend. |
| `# metabolic figure pipeline -2 plot.r` | Joins two input files (e.g. one per instrument run) into a single figure. |
| `# metabolic figure pipeline V2.r` | Full-featured: `NF` handling, isotope renormalisation, KO conditions, bottom legend. |
| `# metabolic figure pipeline V2 -simple.r` | V2 restricted to the dTMP/dTDP/dTTP trio, with panels swapped (total abundance on top). |

`APTX_function.svg` is a standalone reference diagram of Aprataxin (APTX)
function, unrelated to the plotting pipeline.
