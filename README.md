# Tenebrio_Metaproteome

Analysis code for the metaproteomic study of a bacterial consortium derived from the
*Tenebrio molitor* gut microbiome grown on polyethylene.

This repository contains the **specific analysis applied to the dataset reported in the
associated publication**. It is not a general-purpose tool: the scripts implement one
concrete workflow, with the parameters and processing steps described in the Materials
and Methods section of the paper, so that the reported results can be reproduced.

## Contents

| File | Description |
|---|---|
| `Rscript_Metaproteome.R` | Full quantitative analysis in R: filtering and quality control, normalisation, handling of missing values, differential expression analysis, co-expression analysis, biological-relevance filtering, over-representation and gene-set enrichment analysis, random-forest validation, and generation of the protein universe FASTA for external annotation. |
| `annotation_pipeline.sh` | Functional annotation of the protein universe with EggNOG-mapper and InterProScan, and synchronisation of the annotation tables back to the R results folder. |

## Requirements

**R** (≥ 4.2) with the packages listed at the top of `Rscript_Metaproteome.R`. CRAN and
Bioconductor dependencies are installed automatically on first run if missing. The main
ones are `data.table`, `dplyr`, `ggplot2`, `limma`, `DEqMS`, `MSstats`, `imputeLCMD`,
`missForest`, `WGCNA`, `clusterProfiler`, `rstanarm`, `vegan`, `randomForest` and
`Biostrings`.

**Shell environment** (Linux, or WSL on Windows) with `python3`, Biopython,
[EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) and its database, and
[InterProScan](https://www.ebi.ac.uk/interpro/about/interproscan/). Neither tool is
bundled here; both must be installed separately and placed in the folder structure the
script creates.

## Input files

`Rscript_Metaproteome.R` expects the following in `project_dir`, located automatically
by filename pattern:

- the protein quantification matrix exported from the search engine,
- `experimental_metadata*.csv` — sample annotation,
- `experimental_contrasts*.csv` — contrasts to be tested,
- the protein sequence database (FASTA) used for the search.

## Usage

1. Set `project_dir` and `exp_id` at the top of `Rscript_Metaproteome.R` and run the
   script. Phases 1–11 produce the statistical results and write the protein universe
   FASTA; phase 12 copies that FASTA into the annotation workspace.
2. Run `annotation_pipeline.sh`. It reads every FASTA placed in
   `~/FASTA_Annotations/1_FASTA_Files`, runs EggNOG-mapper and InterProScan, and leaves
   the annotation tables in `~/FASTA_Annotations/5_Annotation_results`. Setting
   `WIN_DEST` to the R output folder copies them back automatically; leaving it empty
   skips that step.
3. Continue with the enrichment phases of the R script, which read the annotation
   tables.

Each phase writes to its own numbered subfolder, together with a log and an automated
traceability report recording the number of proteins at each filtering step, the
normalisation method applied and the corresponding test statistics.

## Data availability

The mass spectrometry proteomics data have been deposited in the ProteomeXchange
Consortium via the PRIDE partner repository with the dataset identifier **PXD******
*(to be completed on publication)*.

## Citation

If you use this code, please cite the associated publication and this repository.
A machine-readable citation is provided in `CITATION.cff`.

> *(publication details to be completed)*

## License

MIT — see `LICENSE`.
