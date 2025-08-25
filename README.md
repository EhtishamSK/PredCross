# Prediction and Optimization of Plant Breeding Cross Performance with SimpleMating

## Overview

This R script uses the **SimpleMating** package (Peixoto et al., 2024) to plan crosses, predict progeny performance, and optimize genetic gain for plant breeding. It processes genotypic, phenotypic, and marker effect data to perform cross planning, estimate mid-parental values (MPV), calculate total genetic variance (TGV), predict usefulness, and select optimal crosses for a single trait (e.g., AUDPC).

The workflow includes:

- Converting HapMap genotypic data to a numeric format.
- Generating a genetic map from marker names.
- Planning crosses using the "half" mating design.
- Thinning populations based on relatedness.
- Estimating MPV, TGV, and usefulness components.
- Optimizing crosses for maximum genetic gain.

## Features

- Converts HapMap files to SimpleMating-compatible genotype matrices.
- Creates genetic maps with chromosome and position data.
- Plans crosses with "half" mating design (2145 crosses).
- Thins populations using relatedness thresholds.
- Estimates MPV, TGV, and usefulness for single-trait analysis.
- Optimizes crosses with constraints (max 200 crosses).
- Exports results to CSV files for further analysis.

## Input Requirements

### Genotypic Data (`NMRIL_1973_geno_diploid.hmp.txt`)

- HapMap format with marker IDs, alleles, metadata, and taxa genotypes.

### Phenotypic Data (`pheno.txt`)

- Tab-delimited file with taxa IDs and trait values (e.g., AUDPC).

### Marker Effects (`marker_effects.csv`)

- CSV with marker IDs as row names, additive and dominance effects.

## Output

The script produces:

- `geno.csv`: Numeric genotype matrix.
- `map_data.csv`: Genetic map with chromosome, position, and marker IDs.
- `SingleTrait_MPV.csv`: Mid-parental values for crosses.
- `SingleTrait_TGV.csv`: Total genetic variance for crosses.
- `Cross_Summaries.csv`: Cross summaries (mean, variance, SD, usefulness).
- `Optimization_Input.csv`: Input for cross optimization.
- `selected_crosses.csv`: Optimized cross selections.
- `crosses_data.csv`: Metrics for optimized crosses.

## Author

**Ehtisham Khokhar**  
New Mexico State University  
Email: ehtishamshakeel@gmail.com

## Reference

Peixoto, M. A., Amadeu, R.R., Bhering, L. L., Ferrão, L. F. V., Munoz, P. R., & Resende, M. F. R. (2024). *SimpleMating: R-package for prediction and optimization of breeding crosses using genomic selection.* The Plant Genome, e20533. [https://doi.org/10.1002/tpg2.20533]
For issues or questions, contact the repository maintainer or refer to the SimpleMating documentation: GitHub - Resende-Lab/SimpleMating.
