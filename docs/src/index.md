# BioDemuX.jl

![package_logo](assets/BioDemuX_logo.png)

[![Build Status](https://github.com/I-Mihara/BioDemuX.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/I-Mihara/BioDemuX.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build and Release](https://github.com/I-Mihara/BioDemuX.jl/actions/workflows/release.yml/badge.svg)](https://github.com/I-Mihara/BioDemuX.jl/actions/workflows/release.yml)
[![DOI](https://img.shields.io/badge/DOI-10.24433%2FCO.5417078.v2-blue)](https://doi.org/10.24433/CO.5417078.v2)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14191427.svg)](https://doi.org/10.5281/zenodo.14191427)

## Overview
BioDemuX.jl is a Julia package designed for demultiplexing reads based on barcodes. Given a set of sequencing reads in FASTQ format and a reference barcode in CSV or TSV format, each read is assigned to a FASTQ file corresponding to its barcode. This barcode assignment process is designed to be robust to barcode mutations and allows you to adjust the permitted level of mutation through parameters.

![conceptual_diagram](assets/conceptual_diagram.png)

### Key Features

*   **Ultara-Fast Semi-Global Alignment**: Unrivaled speed with robust handling of indels and mismatches.
*   **Extensive Flexibility**: Custom search ranges, dual barcoding, adapter trimming, and multi-mode matching (`:semiglobal`, `:hamming`, `:exact`).
*   **High-Performance Architecture**: Multi-threaded streaming for efficient memory usage.
*   **Visualize Analysis**: Optional interactive **HTML reports**.
*   **Zero-Dependency CLI**: Standalone binary for easy deployment.

### References
To be published.

## Installation
1. Open the Julia REPL (by typing `julia` in the terminal).
2. Press `]` to enter the Pkg REPL mode.
3. Run the following command to install `BioDemuX.jl`:
    ```julia
    add BioDemuX
    ```
4. Press the `Backspace` key to exit the Pkg mode.

## Dataset
The datasets used in the paper are publicly available.
- **In vitro dataset and test code:**

Link to the dataset and test code: https://doi.org/10.24433/CO.5417078.v1
- **In silico dataset:**

Link to the dataset: https://doi.org/10.5281/zenodo.14178228

## Support or Contact
If you encounter any issues or have requests, please provide feedback by posting a new GitHub issue on our repository. We appreciate your input and will do our best to assist you!
