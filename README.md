# Data-Driven Hypothesis

Analysis of publication patterns across human genes, with focus on mitochondrial genes and inborn errors of metabolism.

## Figure 1

**Publication distribution across human genes and disease-associated subsets.**

**(A)** Scatter plot showing the number of publications per human gene (n=~33,000 genes). Genes are ordered by publication count, with highly studied genes (>3,000 publications) labeled. Most genes have relatively few publications, while a small subset dominates the literature.

**(B)** Publication counts for genes associated with inborn errors of metabolism (n=1,766 genes, IEMbase). Green points highlight metabolic disorder genes, with genes exceeding 500 publications labeled. These disease-relevant genes show moderate to high publication rates.

**(C)** Publication counts for mitochondrial genes (n=1,136 genes, MitoCarta3.0). Red points show mitochondrial gene publication patterns, with genes exceeding 500 publications labeled. Mitochondrial genes represent a well-studied cellular compartment.

**(D)** Distribution of genes across publication bins (0, 1-10, 11-100, 101-1000, 1001+ publications) for all three datasets. Bar plots show that the majority of genes in all categories have low publication counts (1-100 publications), with progressively fewer genes in higher publication bins. This pattern is consistent across all human genes, metabolic disorder genes, and mitochondrial genes, though the specialized subsets show enrichment in higher publication bins.

## Data Sources

- **Gene2PubMed**: NCBI database linking genes to PubMed citations
- **MitoCarta3.0**: Inventory of 1,136 human mitochondrial genes
- **IEMbase**: Database of 2,026 inborn errors of metabolism disorders
- **HGNC Gene Symbols**: Official gene nomenclature

## Files

- `figure1.R` - Main analysis script generating Figure 1
- `download_iembase.R` - Script to download IEMbase metabolic disorder data via API
- `data/` - Cached data files (gene2pubmed, IEMbase disorders, MitoCarta)
- `figure1.pdf` - Publication-quality output figure
