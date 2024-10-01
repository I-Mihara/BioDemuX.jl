# Demux

[![Build Status](https://github.com/I-Mihara/Demux.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/I-Mihara/Demux.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Overview
Demux.jl is a Julia package designed for demultiplexing reads based on barcodes. Given a set of sequencing reads in FASTQ format and a reference barcode in TSV format, each read is assigned to a FASTQ file corresponding to its barcode. This barcode assignment process is designed to be robust to barcode mutations and allows you to adjust the permitted level of mutation through parameters.

### Package features
* Fast and accurate semi-global alignment
* Robust to barcode mutations
* No restrictions on barcode size or position
* Usable as julia package
* Supports parallel computing

### References
To be published.

## Table of Contents
- [Installation](#installation)
- [Basic Usage](#basic-usage)
  - [Input](#input)
  - [Output](#output)
- [Tips to Speed Up Demultiplexing](#tips-to-speed-up-demultiplexing)
- [Options](#options)
  - [Example: How Barcode Length and Option Values Affect Classification](#example-how-barcode-length-and-option-values-affect-classification)
- [Support or Contact](#support-or-contact)


## Dependencies
* Julia >= 1.10.5 (includes the `Distributed` standard library)
* DataFrames >= 1.7.0
* CSV >= 0.10.14

## Installation
1. Open the Julia REPL (by typing `julia` in the terminal).
2. Press `]` to enter the Pkg REPL mode.
3. Run the following command to install `Demux.jl`:
    ```julia
    add Demux
    ```
4. Press the `Backspace` key to exit the Pkg mode.

## Basic Usage

The primary function of this package is `execute_demultiplexing()`. It classifies sequences in a FASTQ file by aligning them with reference barcodes in barcode file. Usage is as follows:
```Julia
using Demux
execute_demultiplexing(FASTQ_file, barcode_file, output_directory)
```

### Input

#### FASTQ File
* There is no restriction on the sequence length in the FASTQ file.
* The function can take one or two FASTQ files as input. In the case of using two FASTQ files, the command can be executed as follows:
```julia
execute_demultiplexing(FASTQ_01, FASTQ_02, barcode_file, output_directory)
```
When using two FASTQ files, sequences in the FASTQ_02 file are classified based on the alignment of the FASTQ_01 sequences with the reference barcodes in the barcode file. Hence, the corresponding reads in both FASTQ files must be in the same order and present in equal numbers.

#### Barcode Reference File
* The reference file is expected to be a TSV file containing the following columns: `ID`, `Full_seq`, `Full_annotation`, as shown below:
```
ID  Full_seq	Full_annotation
001-barcode ACAGACUACAAA LLLBBBBBBB33
```
* In the `Full_seq` column, the region specified as `B` in the `Full_annotation` column is considerd as the barcode.

### Output

* All output files will be saved in the specified `output_directory`.
* The names of the output files are based on the `ID` values in the barcode reference file. For example, if the reference file contains IDs such as `001` and `002`, the resulting output files will be named `001.fastq`, `002.fastq`, and so on.
* Sequences that do not match any barcode in the reference file are saved in `unknown.fastq`. Sequences that have ambiguous classification (i.e., they match multiple barcodes with similar scores) are saved in `ambiguous_classification.fastq`.
* The function will throw an error if the specified `output_directory` already exists to prevent overwriting. A new directory is created to store the output files.

## Tips to Speed Up Demultiplexing

### 1. Parallel Computing
`Demux.jl` supports parallel computing, allowing faster processing of large datasets. To utilize parallel processing, follow the steps below:

#### Starting Julia with Multiple Threads
To enable parallel processing, you need to start Julia with multiple threads. Use the `-p` flag followed by the number of desired threads:
```bash
./julia -p [number_of_threads]
```
#### Adding Worker Processes After Starting Julia
Even after starting the Julia REPL, you can add more worker processes using the `Distributed` module:
```julia
using Distributed
addprocs(n) # 'n' is the number of desired workers
```
#### Runnning `execute_demutltiplexing` with Parallel Computing
Once the worker processes are set up. you can perform parallel computing using the `execute_demultiplexing` function. `Demux.jl` automatically divides the files based on the available worker processes for faster computation:
```julia
@everywhere using Demux
execute_demultiplexing(FASTQ_R1, FASTQ_R2, barcode_file, output_directory)
```

### 2. Setting Options
Demux.jl skips calculation of unnecessary path in DP matrix based on the settings of `max_error_rate`,`mismatch`. and `indel`. By setting lower max error rate or higher penalty, you can further increase computation speed.

## Options

The `execute_demultiplexing` function provides several optional parameters to control the demultiplexing process:

```julia
execute_demultiplexing(FASTQ_file, barcode_file, output_directory, max_error_rate=0.2, min_delta=0.1, mismatch=1, indel=2, classify_both=true, bc_complement=true)
```

- **`max_error_rate::Float64`** (default: `0.2`): 
  - This is the maximum allowed error rate for matching sequences to barcodes. It is multiplied by the barcode's length to calculate the total penalty score that can be tolerated. If the sequence's alignment penalty exceeds this limit for all barcodes, it will be saved in `unknown.fastq`.

- **`min_delta::Float64`** (default: `0`): 
  - This defines the minimum difference in penalty scores needed to confidently assign a sequence to a barcode. It is multiplied by the barcode's length to determine the score difference required to avoid ambiguity. If the difference between the best match's penalty score and the second-best match's score is less than this threshold, the sequence is considered ambiguous and saved in `ambiguous_classification.fastq`.
  
- **`mismatch::Int`** (default: `1`): 
  - The penalty score for mismatches during sequence alignment. A higher value makes the alignment more strict for mismatches.

- **`indel::Int`** (default: `1`): 
  - The penalty score for insertions and deletions (indels) during sequence alignment. A higher value makes the alignment more strict for insertions or deletions.

- **`classify_both::Bool`** (default: `false`): 
  - If set to `true`, the function will classify both R1 and R2 sequences and output separate files for each. Otherwise, it classifies only R2 sequences by default.

- **`bc_complement::Bool`** (default: `false`): 
  - If set to true, the barcodes in the reference file are converted to their complementary sequences before comparison. This is useful when the barcode sequences are provided in the complementary orientation.

### Example: How Barcode Length and Option Values Affect Classification

We assume the case where Barcode Length is 10, `max_error_rate ` is 0.2, `min_delta` is 0.2, `mismatch` is 1, `indel` is 2.


- **Maximum Allowed Penalty Score**:
  - With a `max_error_rate` of 0.2 and a barcode length of 10, the maximum allowed penalty score for a sequence to still match a barcode is `0.2 * 10 = 2`.
- **Minimum Allowed Penalty Difference**:
  - With `min_delta = 0.2` and a barcode length of 10, the minimum required difference in scores between the best and second-best barcode matches is `0.2 * 10 = 2`.

- **Penalty Settings**:
  - **`mismatch = 1`**: Each mismatch in the sequence alignment contributes a penalty of 1.
  - **`indel = 2`**: Each insertion or deletion (indel) contributes a penalty of 2, making indels more costly than mismatches.

With this settings, the classification works as follows:
1. **Allowed Error**:
   - Since the maximum allowed penalty score is 2 (`0.2 * 10`):
     - The sequence can have **up to 2 mismatches** (since each mismatch has a penalty of 1).
     - The sequence can have **up to 1 indel** (since each indel has a penalty of 2).
     - For example, when the sequence have a combination of 1 mismatch and 1 indel, the penalty score is `1 (mismatch) + 2 (indel) = 3`. Since this score exceeds the maximum allowed penalty score of 2, it would **not** be allowed.

2. **Matching Process**:
   - During the alignment, the sequence is compared to each barcode in the reference file. The total penalty score (based on mismatches and indels) is calculated for each alignment.
   - If a sequence's penalty score with a barcode exceeds 2, it **cannot** be classified under that barcode.

3. **Ambiguous Classification**:
   - If a sequence matches multiple barcodes and the penalty scores of the best match and the second-best match differ by **less than 2** (the minimum allowed penalty difference), the sequence will be classified into `ambiguous_classification.fastq`.
   - For example, if the best matching barcode has a penalty score of 1 and the second-best has a score of 2, the difference is `2 - 1 = 1`, which is **less** than the threshold of 2. Therefore, the sequence is ambiguous.

4. **Unknown Classification**:
   - If the sequence fails to match **any** barcode within the maximum allowed penalty score of 2, it is classified as `unknown.fastq`.


## Support or Contanct
If you encounter any issues or have requests, please provide feedback by posting a new GitHub issue on our repository. We appreciate your input and will do our best to assist you!