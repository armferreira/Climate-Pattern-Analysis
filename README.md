# Climate Pattern Analysis Using R

## Overview

This project involves a comprehensive analysis of climate patterns using various statistical and clustering techniques in R. The primary focus is to explore and understand climate data through missing value analysis, correlation matrices, PCA, hierarchical clustering, K-means, PAM clustering, and Gaussian Mixture Models.

## Libraries Used
- `dplyr`: Data manipulation
- `mice`: Imputation for missing data
- `lsr`: Functions for statistical analysis
- `confintr`: Confidence intervals
- `ggcorrplot`: Visualization of correlation matrices
- `psych`: Multivariate analysis and scale construction
- `factoextra`: Extract and visualize results of multivariate data analyses
- `gpairs`: Matrix of scatterplots
- `cluster`: Cluster analysis
- `mclust`: Gaussian finite mixture models
- `ggplot2`: Data visualization
- `sf`: Simple features for spatial data
- `rnaturalearth`: World map data

## Analysis Workflow

### 1. Data Preparation
- Reading climate data from a CSV file
- Missing value analysis

### 2. Data Analysis
- Creating and analyzing correlation matrices
- Handling missing data with MICE (Multivariate Imputation by Chained Equations)
- Performing Principal Component Analysis (PCA)
- Implementing hierarchical clustering, K-means, PAM clustering, and Gaussian Mixture Models
- Profiling geographical clusters and analyzing sustainable development indices

## Key Steps in the Analysis

1. **Missing Values Analysis**: Assess and handle missing data in the climate dataset.
2. **Correlation Matrix Calculation**: Develop a custom function to calculate correlation matrices for various types of variables.
3. **Data Imputation**: Use MICE to impute missing values.
4. **Principal Component Analysis (PCA)**: Standardize the dataset and extract principal components.
5. **Clustering Techniques**:
   - Hierarchical clustering using single-linkage, complete-linkage, and Ward.D2 methods.
   - K-means clustering and validation using silhouette analysis.
   - Partitioning Around Medoids (PAM) clustering.
   - Gaussian Mixture Models (GMM) for advanced clustering.
6. **Geographical Profiling**: Merge climate data with geographical data and visualize clusters on a world map.

## Usage

- Ensure all required libraries are installed in R.
- The dataset 'clima.csv' should be placed in the specified directory.
- Run the script step by step to reproduce the analysis.

## Note

This project is part of an academic course in Pattern Recognition and aims to provide insights into climate patterns using advanced statistical methods.
