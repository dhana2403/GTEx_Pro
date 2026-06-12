# FEED-Curated GTEx Dataset (v1 - Pilot Release)

## Overview

The **FEED-Curated GTEx Dataset (v1)** is a **pilot, research-ready subset of GTEx gene expression data** prepared by FEED AI.


The dataset is intended to demonstrate FEED AI’s curation and dataset-adaptive normalization framework, and will be expanded in future versions to include broader gene and tissue coverage.

The dataset reduces preprocessing overhead by providing standardized expression matrices, harmonized annotations, and quality-controlled data suitable for machine learning, visualization, and downstream bioinformatics analysis.

The dataset is **normalized and batch-corrected using dataset-adaptive logic**, ensuring consistency across tissues and experimental conditions while preserving biological signal.

This resource is designed to support reproducible academic research and data-driven biological discovery.

---


## Dataset Access

- **Zenodo DOI:** https://doi.org/10.5281/zenodo.20668195  
- **Direct download (Zenodo download link):**  
  https://zenodo.org/api/records/20668195/draft/files-archive  
- **Interactive analysis platform:**  
  https://appfeedai.bio/
---

## Data Description

This v1 release includes a curated subset of GTEx gene expression data with:

- Gene-level expression values (pilot subset)
- Harmonized gene identifiers
- Multi-tissue structured format (limited coverage in v1)
- Sample-level organization
- Processed expression matrix for downstream analysis

---

## Data Processing Pipeline


## Features

- Curated multi-tissue gene expression data derived from GTEx  
- Standardized gene identifiers and tissue annotations  
- Quality-controlled, research-ready format  
- Normalized and batch-corrected using dataset-adaptive logic  
- Optimized for machine learning, visualization, and downstream analysis  
- Compatible with GTEx Pro analysis workflows  
- Designed to improve reproducibility across studies  

> ⚠️ This is an initial pipeline version and will evolve in future releases.

---

## Intended Use

This dataset is suitable for:

- Multi-tissue gene expression analysis
- Machine learning model development and benchmarking
- Gene expression visualization
- Bioinformatics method prototyping
- Exploratory data analysis

---

## Limitations

- Pilot-scale gene coverage (not full GTEx transcriptome)
- Limited tissue representation in v1
- Not intended for clinical or diagnostic use

Future versions will expand gene and tissue coverage.

---

## License

This dataset is released under:

**Creative Commons Attribution 4.0 International (CC BY 4.0)**

You are free to:
- Share and redistribute
- Adapt and build upon
- Use commercially or academically

**Condition:** Proper attribution is required.

---

## Citation

If you use this dataset, please cite:

```bibtex
@dataset{feed_gtex_v1_2026,
  author    = {Dhanalakshmi Jothi},
  title     = {FEED-Curated GTEx Dataset (v1 - Pilot Release)},
  publisher = {FEED AI},
  year      = {2026},
  version   = {v1},
  doi       = {10.5281/zenodo.20668195}
}


