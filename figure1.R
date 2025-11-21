library(tidyverse)
library(janitor)
library(cowplot)
library(ggrepel)
library(patchwork)
library(readxl)

# Load gene2pubmed data (or download if not available)
if (file.exists(here::here("data", "gene2pubmed.RData"))) {
  load(here::here("data", "gene2pubmed.RData"))
} else {
  # Download and process data (only runs if saved data doesn't exist)
  gene2pubmedurl <- "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz"
  gene2pubmed_raw <- read_tsv(gene2pubmedurl, col_names = TRUE) |>
    janitor::clean_names()

  # Read data from create_gene_summary.R
  # File that has gene ids -> gene names
  load(here::here("data", "gene_summary.RData"))

  gene2pubmed <- gene2pubmed_raw |>
    dplyr::filter(number_tax_id == 9606) |> #only the rows corresponding to humans (#tax_id = 9606)
    group_by(gene_id) |>
    count(sort = TRUE) |>
    left_join(gene_summary, by = c("gene_id" = "ncbi_gene_id")) |>
    select(gene_id, approved_symbol, n, approved_name) |>
    filter(!is.na(approved_symbol))

  # Save for future use
  save(gene2pubmed, file = here::here("data", "gene2pubmed.RData"))
}

# Read mitochondrial genes from MitoCarta
mitocarta <- read_excel(
  here::here("data", "Human.MitoCarta3.0.xls"),
  sheet = "A Human MitoCarta3.0"
)

# Extract mitochondrial gene symbols (third column)
mito_genes <- mitocarta[[3]] # Third column with gene symbols

# Create mitochondrial gene subset
gene2pubmed_mito <- gene2pubmed |>
  filter(approved_symbol %in% mito_genes)

# Read IEMbase metabolic disorder data
iembase_disorders <- read_csv(
  here::here("data", "iembase_disorders.csv"),
  col_types = cols(.default = "c")
)

# Extract genes associated with metabolic disorders
# Use hgnc_gene_sym as primary gene symbol
metabolic_genes <- iembase_disorders |>
  filter(!is.na(hgnc_gene_sym) & hgnc_gene_sym != "") |>
  pull(hgnc_gene_sym) |>
  unique()

# Create metabolic disorder-associated gene subset
gene2pubmed_metabolic <- gene2pubmed |>
  filter(approved_symbol %in% metabolic_genes)

# Panel A - Full width
panel_a <-
  ggplot() +
  geom_point(
    data = gene2pubmed,
    mapping = aes(x = fct_reorder(approved_symbol, n), y = n),
    alpha = 0.2
  ) +
  geom_text_repel(
    data = subset(gene2pubmed, n > 3000),
    mapping = aes(
      x = fct_reorder(approved_symbol, n),
      y = n,
      label = approved_symbol
    )
  ) +
  labs(
    x = "Human Gene",
    y = "Number of Publications",
    title = "Number of Publications per human gene"
  ) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  coord_cartesian(xlim = c(0, 33000))

# Panel B - Mitochondrial genes subset (half width)
panel_b <-
  ggplot() +
  geom_point(
    data = gene2pubmed_mito,
    mapping = aes(x = fct_reorder(approved_symbol, n), y = n),
    alpha = 0.5,
    color = "darkred"
  ) +
  geom_text_repel(
    data = subset(gene2pubmed_mito, n > 500),
    mapping = aes(
      x = fct_reorder(approved_symbol, n),
      y = n,
      label = approved_symbol
    ),
    size = 3
  ) +
  labs(
    x = "Mitochondrial Genes",
    y = "Publications"
  ) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  coord_cartesian(xlim = c(0, 1200))

# Panel C - Metabolic disorder genes (half width)
panel_c <-
  ggplot() +
  geom_point(
    data = gene2pubmed_metabolic,
    mapping = aes(x = fct_reorder(approved_symbol, n), y = n),
    alpha = 0.5,
    color = "darkgreen"
  ) +
  geom_text_repel(
    data = subset(gene2pubmed_metabolic, n > 500),
    mapping = aes(
      x = fct_reorder(approved_symbol, n),
      y = n,
      label = approved_symbol
    ),
    size = 3
  ) +
  labs(
    x = "Inborn Errors of Metabolism Genes",
    y = "Publications"
  ) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  coord_cartesian(xlim = c(0, 1800))

# Publication binning analysis
# Create function to bin publication counts
bin_publications <- function(data, dataset_name) {
  data |>
    mutate(
      pub_bin = case_when(
        n == 0 ~ "0",
        n >= 1 & n <= 10 ~ "1-10",
        n >= 11 & n <= 100 ~ "11-100",
        n >= 101 & n <= 1000 ~ "101-1000",
        n >= 1001 ~ "1001+",
        TRUE ~ NA_character_
      ),
      pub_bin = factor(pub_bin, levels = c("0", "1-10", "11-100", "101-1000", "1001+"))
    ) |>
    count(pub_bin, .drop = FALSE) |>
    mutate(dataset = dataset_name)
}

# Apply binning to all three datasets
pub_bins <- bind_rows(
  bin_publications(gene2pubmed, "All Human Genes"),
  bin_publications(gene2pubmed_mito, "Mitochondrial Genes"),
  bin_publications(gene2pubmed_metabolic, "Metabolic Disorder Genes")
)

# Panel D - Publication distribution bins (full width)
panel_d <- ggplot(pub_bins, aes(x = pub_bin, y = n, fill = dataset)) +
  geom_col() +
  facet_wrap(~dataset, scales = "free_y") +
  labs(
    x = "Publication Count Bins",
    y = "Number of Genes"
  ) +
  theme_cowplot() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c(
    "All Human Genes" = "grey50",
    "Mitochondrial Genes" = "darkred",
    "Metabolic Disorder Genes" = "darkgreen"
  ))

# Compose figure: A on top (full width), B and C in middle (split width), D on bottom (full width)
figure1 <- panel_a / (panel_b + panel_c) / panel_d +
  plot_annotation(tag_levels = 'A')

# Display the composed figure
figure1

# Save figure as publication-quality PDF
ggsave(
  filename = here::here("figure1.pdf"),
  plot = figure1,
  width = 8.5,
  height = 10,
  units = "in",
  dpi = 300
)
