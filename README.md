# InSe-zscan-XAI-shiny

An R Shiny dashboard for **dual-timescale (ns vs ps)** analysis and benchmarking of machine-learning models on **InSe thin-film optical datasets**.  
The app provides an end-to-end workflow for uploading experimental files, training multiple regressors, comparing predictive metrics, visualizing feature importance, generating parity plots, and inspecting residual error distributions.

## Features

- **Two-channel workflow:** independent processing for **nanosecond (ns)** and **picosecond (ps)** datasets
- **Unified model benchmarking:** LR, DT, RF, SVM (RBF), and XGBoost
- **Performance metrics:** R², RMSE, MAE (exportable to Excel)
- **Explainability (XAI):** Random Forest feature importance (IncNodePurity) for ns vs ps
- **Diagnostics:** parity plots + residual density comparisons
- **High-resolution exports:** plots downloadable in **600 DPI**

## Requirements

- R (recommended: ≥ 4.1)
- R packages (installed automatically via `pacman` in the app):
  - `shiny`, `shinydashboard`, `ggplot2`, `randomForest`, `xgboost`, `e1071`,  
    `rpart`, `rpart.plot`, `writexl`, `dplyr`, `tidyr`, `pacman`

> The app also sets: `options(shiny.maxRequestSize = 100 * 1024^2)` to allow large uploads.


## Citation

If you use this software in academic work, please cite the Zenodo archived release (DOI will be added after Zenodo integration).

## License

Add your license here (e.g., MIT or GPL-3.0).
If you have not selected a license yet, MIT is recommended for broad reuse.

## Contact

For issues, feature requests, or reproducibility questions, please open an issue in this repository.

## Installation

Clone the repository:

```bash
git clone https://github.com/ckrmustafa/InSe-zscan-XAI-shiny.git
cd InSe-zscan-XAI-shiny
