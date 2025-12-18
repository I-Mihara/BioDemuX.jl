![package_logo](BioDemuX_logo.png)

[![Build Status](https://github.com/I-Mihara/BioDemuX.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/I-Mihara/BioDemuX.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build and Release](https://github.com/I-Mihara/BioDemuX.jl/actions/workflows/release.yml/badge.svg)](https://github.com/I-Mihara/BioDemuX.jl/actions/workflows/release.yml)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://I-Mihara.github.io/BioDemuX.jl/)
[![DOI](https://img.shields.io/badge/DOI-10.24433%2FCO.5417078.v2-blue)](https://doi.org/10.24433/CO.5417078.v2)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14191427.svg)](https://doi.org/10.5281/zenodo.14191427)

# BioDemuX.jl

**BioDemuX.jl** is a high-performance Julia package for demultiplexing sequencing reads with robust barcode mutation handling. It is designed to be the definitive tool for accurate and efficient barcode assignment.

## 🚀 Key Features

*   **⚡ Ultra-Fast Semi-Global Alignment**:
    *   The core engine is built on highly optimized semi-global alignment. It delivers **unrivaled speed** while robustly handling **insertions, deletions (indels)**, and mismatches—recovering reads that standard Hamming-based tools discard.

*   **🔧 Extensive Flexibility**:
    *   **Precise Targeting**: Define custom search ranges to focus on specific regions.
    *   Supports **Dual Barcoding** and variable-length barcodes.
    *   **Adapter Trimming**: Support for simultaneous adapter matching and trimming.
    *   **Multi-Mode Matching**: Choose between `:semiglobal` for robustness, or `:hamming`/`:exact` for speed.

*   **🚀 High-Performance Architecture**:
    *   Multi-threaded & Channel-based streaming.
    *   Efficient memory usage via streaming I/O.

*   **📊 Visualize Analysis**:
    *   Optional interactive **HTML reports** for distribution and score analysis.

*   **💻 Zero-Dependency CLI**:
    *   Available as a standalone binary for Linux, macOS, and Windows. No Julia installation required to run in production pipelines.

## 📚 Documentation

**[Full Documentation](https://I-Mihara.github.io/BioDemuX.jl/)**

Visit our documentation for comprehensive guides, API references, and detailed usage examples.

## 📦 Installation
```julia
using Pkg; Pkg.add("BioDemuX")
```

## 🏃 Quick Start

```julia
using BioDemuX
execute_demultiplexing("reads.fastq", "barcodes.csv", "output_dir")
```

**[CLI Users]** Download the binary from [Releases](https://github.com/I-Mihara/BioDemuX.jl/releases) and run:
```bash
./biodemux reads.fastq barcodes.csv output_dir
```

## 👥 Support
If you encounter any issues, please open an issue on [GitHub](https://github.com/I-Mihara/BioDemuX.jl/issues).
