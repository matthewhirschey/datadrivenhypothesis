library(tidyverse)
library(janitor)
library(cowplot)
library(ggrepel)
library(patchwork)

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

# Panel B - Placeholder (half width)
panel_b <- ggplot() +
  geom_point(data = gene2pubmed, aes(x = 1, y = 1)) +
  labs(title = "Panel B - Placeholder") +
  theme_cowplot()

# Panel C - Placeholder (half width)
panel_c <- ggplot() +
  geom_point(data = gene2pubmed, aes(x = 1, y = 1)) +
  labs(title = "Panel C - Placeholder") +
  theme_cowplot()

# Compose figure: A on top (full width), B and C below (split width)
figure1 <- panel_a / (panel_b + panel_c) +
  plot_annotation(tag_levels = 'A')

# Display the composed figure
figure1
