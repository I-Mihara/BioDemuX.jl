# Usage

## Dependencies
- **Julia** >= 1.10 (includes the `Distributed` standard library)
- **Julia packages**
  - DataFrames >= 1.7.0
  - CSV >= 0.10.14
  - CodecZlib >= 0.7.0
  - BufferedStreams >= 1.2.2

## Basic Usage

The primary function of this package is `execute_demultiplexing()`. It classifies sequences in a FASTQ file by aligning them with reference barcodes in barcode file. Usage is as follows:
```Julia
using BioDemuX
execute_demultiplexing(FASTQ_file, barcode_file, output_directory)
```

### Input

#### FASTQ File
* There is no restriction on the sequence length in the FASTQ file.
* The input can be gzipped; the function automatically detects and processes the files accordingly.
* The function can take one or two FASTQ files as input. In the case of using two FASTQ files, the command can be executed as follows:
```julia
execute_demultiplexing(FASTQ_file1, FASTQ_file2, barcode_file, output_directory)
```
When using two FASTQ files, sequences in the `FASTQ_file2` are classified based on the alignment of the `FASTQ_file1` sequences with the barcodes in the barcode reference file. Hence, the corresponding reads in both FASTQ files must be in the same order and present in equal numbers.

#### Barcode Reference File
* The reference file is expected to be a CSV or TSV file containing the following columns: `ID`, `Full_seq`, `Full_annotation`, as shown below:
```
ID  Full_seq	Full_annotation
001-barcode ACAGACUACAAA XXXBBBBBBBXX
```
* In the `Full_seq` column, the region specified as `B` in the `Full_annotation` column is considered as the barcode.

* Alternatively, a FASTA file of barcode sequences can be used as the reference. In this case, each sequence in the FASTA file is treated as a full barcode (the entire sequence is considered the barcode region) and the header line of each entry (without the `>` prefix) is used as its `ID`.

### Output

* All output files will be saved in the specified `output_directory`.
* The output is gzipped depending on the input FASTQ format, or can be specified using the `gzip_output` option.
* The names of the output files are based on the filename of the FASTQ file as the prefix and the `ID` values in the barcode reference file. For example, if the FASTQ filename is `sample.fastq` and the reference file contains IDs such as `001` and `002`, the resulting output files will be named `sample.001.fastq`, `sample.002.fastq`, and so on. You can freely change the prefix by specifying the `output_prefix` argument.
* Sequences that do not match any barcode in the reference file are saved in `unknown.fastq`. Sequences that have ambiguous classification (i.e., they match multiple barcodes with similar scores) are saved in `ambiguous_classification.fastq`. These FASTQ files also have prefix like `sample.unknown.fastq` and `sample.ambiguous_classification.fastq`.
* If the `output_directory` does not exist, a new directory is created to store the output files.
