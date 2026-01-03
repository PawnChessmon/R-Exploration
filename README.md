# Rexploration Omega

```mermaid
flowchart TD
  A["Input tables"] --> B["DESeq2 model"]
  B --> C["DE tables: step2a"]
  B --> D["Scaled matrix (VST + z-score)"]
  D --> E["PCA"]
  D --> F["Heatmap of top DE genes"]
  B --> G["Volcano plots"]
  B --> H["MA plots"]
  D --> I["Boxplots: top 10 DE genes"]
  B --> J["Top 10 up/down lists"]
```

## TL;DR
- Run with p-values:
  - `nextflow run main.nf --pcol pvalue --pthresh 0.01 --lfc 1`
- Run with adjusted p-values:
  - `nextflow run main.nf --pcol padj --pthresh 0.05 --lfc 1`

Outputs are written to `output_step2a` and either `output_step2b` (pvalue) or `output_step2c` (padj).


## Sample Outputs

  ### PCA
  ![PCA plot](docs/images/pca_samples.png)

  ### Heatmap
  ![Heatmap of top DE genes](docs/images/heatmap_top_degenes.png)

  ### Volcano (gut vs duct)
  ![Volcano plot](docs/images/volcano_de_gut_duct.png)

  ### MA (gut vs duct)
  ![MA plot](docs/images/ma_de_gut_duct.png)

  ### Boxplots (Top 10 DE genes)
  ![Boxplots of top DE genes](docs/images/boxplots_top10.png)



## Overview
This pipeline runs DESeq2 on the provided experiment inputs, filters non-finite values, produces DE tables, PCA, heatmap, volcano/MA plots, boxplots, and top-10 up/down tables. It enforces the sample group order `gut -> duct -> node`.

### Inputs
Expected under `input/`:
- `em.csv` (tab-delimited expression matrix, first column `ID`)
- `sample_sheet.csv` (tab-delimited: `SAMPLE`, `SAMPLE_GROUP`)
- `annotations.csv` (tab-delimited: `Gene ID`, `Associated Gene Name`)

### Outputs
Always:
- `output_step2a/de_gut_duct.tsv`
- `output_step2a/de_duct_node.tsv`
- `output_step2a/de_node_gut.tsv`

When `--pcol pvalue`:
- `output_step2b/` with PCA, heatmap, volcano, MA, boxplots, and top-10 tables

When `--pcol padj`:
- `output_step2c/` with the same plots/tables based on adjusted p-values

### Parameters
- `--pcol` : `pvalue` or `padj`
- `--pthresh` : significance threshold (default `0.01`)
- `--lfc` : absolute log2 fold-change threshold (default `1`)
- `--input_dir` : input directory (default `input`)
- `--outdir` : output directory (default `.`)




## Run locally
```bash
nextflow run main.nf --pcol pvalue --pthresh 0.01 --lfc 1
```

## Docker
If you prefer Docker, create a minimal container with R and required packages. Then run:
```bash
nextflow run main.nf -profile docker --pcol pvalue --pthresh 0.01 --lfc 1
```

You will need a Docker-enabled `nextflow.config` profile that sets `process.container` to an image containing:
- R (>= 4.2 recommended)
- Bioconductor `DESeq2`
- CRAN: `ggplot2`, `ggrepel`, `pheatmap`

## Conda
Create a conda env with R + packages and use Nextflow's conda profile:
```bash
nextflow run main.nf -profile conda --pcol pvalue --pthresh 0.01 --lfc 1
```

Example `nextflow.config` snippet for conda:
```
profiles {
  conda {
    process.conda = 'conda/renv.yml'
  }
  docker {
    process.container = 'your-docker-image:tag'
  }
}
```

## Notes
- The pipeline filters out NaN/Inf values before writing outputs.
- Sample group order is enforced as `gut`, `duct`, `node`.
- Labels on volcano/MA plots use `geom_label_repel`.
