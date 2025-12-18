# Options & Tips

## Tips to Speed Up Demultiplexing

### 1. Parallel Computing
`BioDemuX.jl` supports parallel computing using multi-threading and channel-based streaming IO, allowing faster processing of large datasets with efficient memory usage. To utilize parallel processing, follow the steps below:

#### Starting Julia with Multiple Threads
To enable parallel processing, you need to start Julia with multiple threads. Use the `-t` (or `--threads`) flag followed by the number of desired threads:
```bash
./julia -t [number_of_threads]
```
For example, to use 8 threads:
```bash
./julia -t 8
```
You can use the `-t auto` option to automatically detect the number of available threads and start Julia with that number of threads.

You can also set the `JULIA_NUM_THREADS` environment variable.
```bash
export JULIA_NUM_THREADS=8
```

#### Running `execute_demultiplexing` with Parallel Computing
Once Julia is started with multiple threads, `BioDemuX.jl` automatically utilizes them for reading, writing, and processing data in parallel. No additional setup is required.

### 2. Setting Options
BioDemuX.jl skips calculation of unnecessary path in DP matrix based on the settings of `max_error_rate`,`mismatch`. and `indel`. By setting lower max error rate or higher penalty, you can further increase computation speed.

## Options Reference

The `execute_demultiplexing` function provides several optional parameters to control the demultiplexing process:

```julia
execute_demultiplexing(FASTQ_file, barcode_file, output_directory; barcode_file2=nothing, output_prefix="", gzip_output=nothing, max_error_rate=0.2, min_delta=0.0, mismatch=1, indel=1, nindel=nothing, bc_complement=false, bc_rev=false, ref_search_range="1:end", barcode_start_range="1:end", barcode_end_range="1:end", ref_search_range2="1:end", barcode_start_range2="1:end", barcode_end_range2="1:end", chunk_size=4000, channel_capacity=64, matching_algorithm=:semiglobal, log=false)
```

- **`max_error_rate::Float64`** (default: `0.2`): 
  - This is the maximum allowed error rate for matching sequences to barcodes. It is multiplied by the barcode's length to calculate the total penalty score that can be tolerated. If the sequence's alignment penalty exceeds this limit for all barcodes, it will be saved in `unknown.fastq`.

- **`matching_algorithm::Symbol`** (default: `:semiglobal`):
  - Specifies the algorithm used for barcode matching.
    - `:semiglobal`: Uses semi-global alignment allowing for mismatches, insertions, and deletions. This is the default and most robust mode.
    - `:hamming`: Uses Hamming distance, allowing only for mismatches (no insertions or deletions). Faster than semiglobal.
    - `:exact`: Uses exact string matching, allowing no mismatches or indels. Fastest mode. Use only when data quality is perfect or strict filtering is desired.

- **`min_delta::Float64`** (default: `0.0`): 
  - This defines the minimum difference in penalty scores needed to confidently assign a sequence to a barcode. It is multiplied by the barcode's length to determine the score difference required to avoid ambiguity. If the difference between the best match's penalty score and the second-best match's score is less than this threshold, the sequence is considered ambiguous and saved in `ambiguous_classification.fastq`.
  
- **`mismatch::Int`** (default: `1`): 
  - The penalty score for mismatches during sequence alignment. A higher value makes the alignment more strict for mismatches.

- **`indel::Int`** (default: `1`): 
  - The penalty score for insertions and deletions (indels) during sequence alignment. A higher value makes the alignment more strict for insertions or deletions.

- **`nindel::Union{Int, Nothing}`** (default: `nothing`):
  - The penalty score for insertions and deletions (indels) when the base is 'N' (wildcard). This should only be specified when using 'N' in barcodes or reads. Note that enabling this option may slow down the demultiplexing process.

- **`classify_both::Bool`** (default: `false`): 
  - If set to `true`, the function will classify sequences in both `FASTQ_file1` and `FASTQ_file2` based on the alignment with sequences in `FASTQ_file1` and output separate files for each. Otherwise, it classifies only R2 sequences by default.

- **`bc_complement::Bool`** (default: `false`): 
  - If set to true, the barcodes in the reference file are converted to their complementary sequences before alignment.

- **`bc_rev::Bool`** (default: `false`):
  - If set to true, the barcodes in the reference file are reversed before alignment.

- **`output_prefix1::String`** (default: `""`):
  - Specifies the prefix for the first set of output files when processing two FASTQ files. If not provided, the prefix defaults to the name of the `FASTQ_file1`. By setting this option, you can customize the file names for the first set of outputs.

- **`output_prefix2::String`** (default: `""`):
  - Specifies the prefix for the second set of output files when processing two FASTQ files. If not provided, the prefix defaults to the name of the `FASTQ_file2`. By setting this option, you can customize the file names for the second set of outputs.

- **`output_prefix::String`** (default: `""`):
  - Specifies the prefix for the output files when processing a single FASTQ file. If not provided, the prefix defaults to the name of the FASTQ file.

- **`gzip_output::Bool`** (default: `auto-detect`):
  - Controls whether the output FASTQ files are compressed (gzipped). By default (when this option is not explicitly set), the output will be gzipped if the input FASTQ files have a `.gz` extension, and uncompressed otherwise. You can set `gzip_output=true` to force gzipped output files or `gzip_output=false` to ensure output files are not compressed.

- **`ref_search_range::String`** (default: `"1:end"`):
  - Specifies the range within the read sequence where the barcode search should be performed. The format is `"start:end"`, where `start` and `end` can be integers (1-based index) or relative to the end of the sequence using `end` (e.g., `"1:20"`, `"end-19:end"`). This allows you to restrict the search to a specific region, improving performance and accuracy if the barcode position is known.

- **`barcode_start_range::String`** (default: `"1:end"`):
  - Specifies the allowed range for the start position of the barcode alignment within the read. The format is the same as `ref_search_range`. If the aligned barcode starts outside this range, it will not be considered a valid match.

- **`barcode_end_range::String`** (default: `"1:end"`):
  - Specifies the allowed range for the end position of the barcode alignment within the read. The format is the same as `ref_search_range`. If the aligned barcode ends outside this range, it will not be considered a valid match.

- **`chunk_size::Int`** (default: `4000`):
  - Specifies the number of reads to process in a single chunk. Larger chunk sizes can reduce overhead but increase memory usage.

- **`channel_capacity::Int`** (default: `64`):
  - Specifies the capacity of the input/output channels. The recycle channel capacity is automatically set to `2 * channel_capacity`. Increasing this value can improve throughput on systems with high memory availability.

- **`barcode_file2::Union{String, Nothing}`** (default: `nothing`):
  - Specifies the path to the second barcode reference file for dual-index demultiplexing. If provided, the function performs a two-step alignment: first against `barcode_file`, then against `barcode_file2`.

- **`ref_search_range2::String`** (default: `"1:end"`):
  - Specifies the search range for the second barcode. Format is the same as `ref_search_range`.

- **`barcode_start_range2::String`** (default: `"1:end"`):
  - Specifies the allowed start range for the second barcode. Format is the same as `barcode_start_range`.

- **`barcode_end_range2::String`** (default: `"1:end"`):
  - Specifies the allowed end range for the second barcode. Format is the same as `barcode_end_range`.

- **`summary::Bool`** (default: `false`):
  - If set to `true`, a summary report is generated after demultiplexing. The report includes statistics on total reads, matched reads, unmatched reads, ambiguous reads, and detailed distributions of alignment scores, positions, and lengths for each barcode.

- **`summary_format::Symbol`** (default: `:html`):
  - Specifies the format of the summary report. Options are:
    - `:text`: A human-readable text file (`summary.txt`).
    - `:html`: An interactive HTML report (`summary.html`) with tabs and dropdowns for detailed visualization.
    - `:json`: A JSON file (`summary.json`) for programmatic parsing.
    - `:stdout`: Prints the summary to the standard output.

- **`trim_side::Union{Int, Nothing}`** (default: `nothing`):
  - Specifies which side of the read to trim relative to the barcode.
    - `3`: Trims the 3' side (right side) of the barcode match. Keeps the sequence *before* the barcode.
    - `5`: Trims the 5' side (left side) of the barcode match. Keeps the sequence *after* the barcode.
  - If `nothing`, no trimming is performed.

- **`trim_side2::Union{Int, Nothing}`** (default: `nothing`):
  - Specifies trimming for the second barcode in dual-barcode mode. Options are the same as `trim_side`.

- **`log::Bool`** (default: `false`):
  - If set to `true`, basic logging information (start configuration and end duration) is printed to the standard error (`stderr`). 

## Example: How Barcode Length and Option Values Affect Classification

We assume the case where barcode length is 10, `max_error_rate ` is 0.2, `min_delta` is 0.2, `mismatch` is 1, `indel` is 2.


- **Maximum Allowed Penalty Score**:
  - With a `max_error_rate` of 0.2 and a barcode length of 10, the maximum allowed penalty score for a sequence to still match a barcode is `0.2 * 10 = 2`.
- **Minimum Allowed Penalty Difference**:
  - With `min_delta = 0.2` and a barcode length of 10, the minimum required difference in scores between the best and second-best barcode matches is `0.2 * 10 = 2`.

- **Penalty Settings**:
  - **`mismatch = 1`**: Each mismatch in the sequence alignment contributes a penalty of 1.
  - **`indel = 2`**: Each insertion or deletion (indel) contributes a penalty of 2, making indels more costly than mismatches.

With these settings, the classification works as follows:
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
