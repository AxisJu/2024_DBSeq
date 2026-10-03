# 2024_DBSeq

Bioinformatics analysis and statistical figure code for studying brain-wide transcriptomic remodeling following anterior thalamic nucleus deep brain stimulation in a chronic epilepsy model, with comparisons across human drug-resistant epilepsy cohorts.

## Repository Architecture

```text
2024_DBSeq/
├── analysis/
│   ├── 100_data-preprocessing/       # Quality control, donor annotation, integration, and cell typing
│   ├── 200_TRG-identification/       # Differential expression and treatment-reverted gene identification
│   ├── 300_TRG-profiling-inATN/      # ATN subclustering, gene profiling, and spatial reference mapping
│   ├── 400_deconv/                   # Connectome metrics and regional molecular response analyses
│   ├── 500_GRN/                      # Gene regulatory networks and promoter motif analyses
│   └── 600_Others/
│       ├── 202607-Squidiff/          # Generative modeling and virtual perturbation
│       ├── 260909-DOLPHIN/           # Exon graphs, splicing, and donor-level analyses
│       └── 260921-CUT&TAG/           # Peak annotation, motif scanning, and signal quantification
├── figures/
│   ├── figure 1/ ... figure 7/       # Main-figure analysis and plotting code
│   ├── extended data figure 1/       # Extended Data Figure 1b
│   ├── extended data figure 2/       # Cell annotation and reference visualization
│   ├── extended data figure 5/ ... extended data figure 9/
│   ├── extended data figure 11/ ... extended data figure 14/
│   └── shared/                      # Shared figure utilities
└── README.md
```

The source collection contains 78 Python modules and 82 R scripts. Figure directories follow the manuscript figure and panel numbering.

## Input Convention

This source release presents computation and plotting routines for code review. Python modules expose `compute(inputs)` and R scripts expose `compute(supplied_inputs)`. The calling environment supplies preloaded in-memory data mappings, helper functions, and analysis options. Data loading is handled outside this repository; this release does not provide a standalone data-to-figure execution workflow.

## Prerequisites

The routines use Python and R, with dependencies varying by analysis stage. Major packages include:

- Single-cell analysis: `scanpy`, `anndata`, `scvi-tools`, `cell2location`, and `Seurat`.
- Numerical analysis and modeling: `numpy`, `pandas`, `scipy`, `torch`, `torch-geometric`, and `pyro-ppl`.
- Differential expression and gene regulation: `edgeR`, `AUCell`, and pySCENIC-related tools.
- Visualization and data manipulation: `matplotlib`, `ggplot2`, `ComplexHeatmap`, `dplyr`, and `qs`.

Squidiff and DOLPHIN workflows also use their respective model implementations. The imports in each module identify its stage-specific dependencies.

## Publication

Associated manuscript: **Connectome topology shapes brain-wide transcriptomic remodeling induced by deep brain stimulation**.

## Third-Party Code

The DOLPHIN adaptation files retain the original MIT license and copyright notices in their source headers.
