library(tidyverse)
library(janitor)
library(cowplot)
library(ggrepel)
library(patchwork)
library(readxl)

gene2pubmedurl <- "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz"
gene2pubmed_raw <- read_tsv(gene2pubmedurl, col_names = TRUE) |>
  janitor::clean_names()

#read data from create_gene_summary.R
#file that has gene ids -> gene names
load(here::here("data", "gene_summary.RData"))

gene2pubmed <- gene2pubmed_raw |>
  dplyr::filter(number_tax_id == 9606) |> #only the rows corresponding to humans (#tax_id = 9606)
  group_by(gene_id) |>
  count(sort = TRUE) |>
  left_join(gene_summary, by = c("gene_id" = "ncbi_gene_id")) |>
  select(gene_id, approved_symbol, n, approved_name) |>
  filter(!is.na(approved_symbol))

# Read mitochondrial genes from MitoCarta
mitocarta <- read_excel(
  here::here("data", "Human.MitoCarta3.0.xls"),
  sheet = "A Human MitoCarta3.0"
)

# Extract mitochondrial gene symbols (third column)
mito_genes <- mitocarta[[3]]  # Third column with gene symbols

# Create mitochondrial gene subset
gene2pubmed_mito <- gene2pubmed |>
  filter(approved_symbol %in% mito_genes)

# Read OMIM disease-gene associations
mim2gene <- read_tsv(
  here::here("data", "mim2gene.txt"),
  comment = "#",
  col_names = c("mim_number", "mim_entry_type", "entrez_gene_id", "approved_symbol", "ensembl_gene_id"),
  col_types = cols(.default = "c")
)

# Extract genes associated with diseases
# Filter for entries marked as "gene" or "gene/phenotype" (disease-relevant genes in OMIM)
disease_genes <- mim2gene |>
  filter(str_detect(mim_entry_type, "gene")) |>
  filter(!is.na(approved_symbol) & approved_symbol != "") |>
  pull(approved_symbol) |>
  unique()

# Create disease-associated gene subset
gene2pubmed_disease <- gene2pubmed |>
  filter(approved_symbol %in% disease_genes)

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
    x = "Mitochondrial Gene",
    y = "Number of Publications",
    title = "Publications per Mitochondrial Gene"
  ) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Panel C - Disease-associated genes (half width)
panel_c <-
  ggplot() +
  geom_point(
    data = gene2pubmed_disease,
    mapping = aes(x = fct_reorder(approved_symbol, n), y = n),
    alpha = 0.5,
    color = "darkblue"
  ) +
  geom_text_repel(
    data = subset(gene2pubmed_disease, n > 1000),
    mapping = aes(
      x = fct_reorder(approved_symbol, n),
      y = n,
      label = approved_symbol
    ),
    size = 3
  ) +
  labs(
    x = "Disease-Associated Gene",
    y = "Number of Publications",
    title = "Publications per Disease-Associated Gene (OMIM)"
  ) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Compose figure: A on top (full width), B and C below (split width)
figure1 <- panel_a / (panel_b + panel_c) +
  plot_annotation(tag_levels = 'A')

# Display the composed figure
figure1

# Save figure as publication-quality PDF
ggsave(
  filename = here::here("figure1.pdf"),
  plot = figure1,
  width = 8.5,
  height = 11,
  units = "in",
  dpi = 300
)
