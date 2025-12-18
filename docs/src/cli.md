# CLI Usage

BioDemuX can also be used as a standalone command-line tool without installing Julia explicitly, or by building it from source.

### 1. Download Pre-built Binary
Go to the **[Releases](https://github.com/I-Mihara/BioDemuX.jl/releases)** page and download the executable for your OS (Windows, macOS, or Linux).

### 2. Build from Source
If you want to build the executable yourself:
```bash
git clone https://github.com/I-Mihara/BioDemuX.jl.git
cd BioDemuX.jl
julia --project=. create_app.jl
```
This creates the executable at `biodemux/bin/biodemux` (or `BioDemuX` on some systems).

### 3. Running the Tool
> [!NOTE]
> Julia apps may have a short startup delay (~0.5-1s) compared to C/Rust tools. To minimize this impact when processing many files, use **Directory Mode** (see below).

#### Directory Mode (Recommended for Speed)
Pass a **directory path** instead of a single file to process all FASTQ files inside it efficiently (single startup).

**Single-end:**
```bash
./biodemux fastq_dir/ barcodes.csv output_dir
```

**Paired-end:**
```bash
./biodemux R1_dir/ barcodes.csv output_dir --fastq2 R2_dir/
```
*Note: Files in `R1_dir` and `R2_dir` are sorted alphabetically and matched 1-to-1. Please ensure file counts and order match.*

#### Parallel Execution
To enable multi-threading, set the environment variable `JULIA_NUM_THREADS` before running.
```bash
export JULIA_NUM_THREADS=auto  # Uses all available cores
./biodemux ...
```

#### Full Options
You can pass any of the configuration options as command-line arguments. Use `--help` to see the full list.

**Example with options:**
```bash
./biodemux input.fastq barcodes.csv output_dir \
    --max-error-rate 0.2 \
    --mismatch 1 \
    --summary \
    --gzip-output
```
*Note: Boolean flags like `--summary` do not require a value. Keyword arguments like `max_error_rate` require a value after the flag.*

```bash
./biodemux --help
```
