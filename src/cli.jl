using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "fastq1"
        help = "Path to the first FASTQ file (Read 1) OR directory containing FASTQ files"
        required = true
        "barcode_file"
        help = "Path to the barcode file (CSV/TSV)"
        required = true
        "output_directory"
        help = "Directory to save demultiplexed files"
        required = true
        "--fastq2"
        help = "Path to the second FASTQ file (Read 2) OR directory (paired with fastq1)"
        default = nothing
        "--barcode-file2", "-B"
        help = "Path to the second barcode file for dual indexing"
        default = nothing
        "--output-prefix1", "-p"
        help = "Output prefix for Read 1"
        default = ""
        "--output-prefix2", "-P"
        help = "Output prefix for Read 2"
        default = ""
        "--gzip-output", "-z"
        help = "Compress output with GZIP. If not set, inferred from input filenames."
        action = :store_true
        "--no-gzip-output"
        help = "Force disable GZIP output"
        action = :store_true
        "--max-error-rate", "-e"
        help = "Maximum error rate allowed"
        arg_type = Float64
        default = 0.2
        "--min-delta", "-d"
        help = "Minimum delta between best and second best match"
        arg_type = Float64
        default = 0.0
        "--match", "-m"
        help = "Match score"
        arg_type = Int
        default = 0
        "--mismatch", "-M"
        help = "Mismatch penalty"
        arg_type = Int
        default = 1
        "--indel", "-i"
        help = "Indel penalty"
        arg_type = Int
        default = 1
        "--nindel", "-I"
        help = "N-indel penalty"
        arg_type = Int
        default = nothing
        "--classify-both", "-c"
        help = "Classify both reads in paired-end mode"
        action = :store_true
        "--bc-complement", "-C"
        help = "Use complement of barcode sequences"
        action = :store_true
        "--bc-rev", "-r"
        help = "Use reverse of barcode sequences"
        action = :store_true
        "--ref-search-range"
        help = "Range to search in reference (e.g., '1:20')"
        default = "1:end"
        "--barcode-start-range"
        help = "Range for barcode start (e.g., '1:5')"
        default = "1:end"
        "--barcode-end-range"
        help = "Range for barcode end"
        default = "1:end"
        "--ref-search-range2"
        help = "Range to search in reference for barcode 2"
        default = "1:end"
        "--barcode-start-range2"
        help = "Range for barcode 2 start"
        default = "1:end"
        "--barcode-end-range2"
        help = "Range for barcode 2 end"
        default = "1:end"
        "--chunk-size"
        help = "Chunk size for processing"
        arg_type = Int
        default = 4000
        "--channel-capacity"
        help = "Channel capacity"
        arg_type = Int
        default = 64
        "--trim-side"
        help = "Trim side for Read 1 (3 or 5)"
        arg_type = Int
        default = nothing
        "--trim-side2"
        help = "Trim side for Read 2 (3 or 5)"
        arg_type = Int
        default = nothing
        "--summary"
        help = "Generate summary report"
        action = :store_true
        "--summary-format"
        help = "Summary format (html, csv, etc.)"
        default = "html"
        "--matching-algorithm"
        help = "Matching algorithm (semiglobal, exact)"
        default = "semiglobal"
        "--log", "-l"
        help = "Enable logging to stderr"
        action = :store_true
    end

    return parse_args(s)
end

function julia_main()::Cint
    try
        args = parse_commandline()

        # Handle gzip flag logic (3-state: true, false, nothing)
        gzip_val = nothing
        if args["gzip-output"]
            gzip_val = true
        elseif args["no-gzip-output"]
            gzip_val = false
        end

        # Handle nindel logic
        nindel_val = args["nindel"] # is actually Union{Int, Nothing} effectively

        fastq1 = args["fastq1"]
        fastq2 = args["fastq2"]
        barcode_file = args["barcode_file"]
        output_dir = args["output_directory"]

        # Convert some string args to symbols
        summary_fmt_sym = Symbol(args["summary-format"])
        matching_algo_sym = Symbol(args["matching-algorithm"])

        if isdir(fastq1)
            # --- Directory Mode ---
            files1 = sort(filter(f -> endswith(f, ".fastq") || endswith(f, ".fq") || endswith(f, ".fastq.gz") || endswith(f, ".fq.gz"), readdir(fastq1, join=true)))

            if isempty(files1)
                println(stderr, "Error: No FASTQ files found in directory: $fastq1")
                return 1
            end

            if !isnothing(fastq2)
                # Paired Directory Mode
                if !isdir(fastq2)
                    println(stderr, "Error: fastq1 is a directory but fastq2 is a file. Both must be directories or both must be files.")
                    return 1
                end

                files2 = sort(filter(f -> endswith(f, ".fastq") || endswith(f, ".fq") || endswith(f, ".fastq.gz") || endswith(f, ".fq.gz"), readdir(fastq2, join=true)))

                if length(files1) != length(files2)
                    println(stderr, "Error: File count mismatch between input directories.")
                    println(stderr, "  $(fastq1): $(length(files1)) files")
                    println(stderr, "  $(fastq2): $(length(files2)) files")
                    return 1
                end

                # Process matched pairs
                for (f1, f2) in zip(files1, files2)
                    if args["log"]
                        println(stderr, "Processing pair: $(basename(f1)) and $(basename(f2))")
                    end

                    execute_demultiplexing(
                        f1,
                        f2,
                        barcode_file,
                        output_dir;
                        barcode_file2=args["barcode-file2"],
                        output_prefix1=args["output-prefix1"],
                        output_prefix2=args["output-prefix2"],
                        gzip_output=gzip_val,
                        max_error_rate=args["max-error-rate"],
                        min_delta=args["min-delta"],
                        match=args["match"],
                        mismatch=args["mismatch"],
                        indel=args["indel"],
                        nindel=nindel_val,
                        classify_both=args["classify-both"],
                        bc_complement=args["bc-complement"],
                        bc_rev=args["bc-rev"],
                        ref_search_range=args["ref-search-range"],
                        barcode_start_range=args["barcode-start-range"],
                        barcode_end_range=args["barcode-end-range"],
                        ref_search_range2=args["ref-search-range2"],
                        barcode_start_range2=args["barcode-start-range2"],
                        barcode_end_range2=args["barcode-end-range2"],
                        chunk_size=args["chunk-size"],
                        channel_capacity=args["channel-capacity"],
                        trim_side=args["trim-side"],
                        trim_side2=args["trim-side2"],
                        summary=args["summary"],
                        summary_format=summary_fmt_sym,
                        matching_algorithm=matching_algo_sym,
                        log=args["log"]
                    )

                end
            else
                # Single Directory Mode
                for f1 in files1
                    if args["log"]
                        println(stderr, "Processing file: $(basename(f1))")
                    end

                    execute_demultiplexing(
                        f1,
                        barcode_file,
                        output_dir;
                        barcode_file2=args["barcode-file2"],
                        output_prefix=args["output-prefix1"],
                        gzip_output=gzip_val,
                        max_error_rate=args["max-error-rate"],
                        min_delta=args["min-delta"],
                        match=args["match"],
                        mismatch=args["mismatch"],
                        indel=args["indel"],
                        nindel=nindel_val,
                        bc_complement=args["bc-complement"],
                        bc_rev=args["bc-rev"],
                        ref_search_range=args["ref-search-range"],
                        barcode_start_range=args["barcode-start-range"],
                        barcode_end_range=args["barcode-end-range"],
                        ref_search_range2=args["ref-search-range2"],
                        barcode_start_range2=args["barcode-start-range2"],
                        barcode_end_range2=args["barcode-end-range2"],
                        chunk_size=args["chunk-size"],
                        channel_capacity=args["channel-capacity"],
                        trim_side=args["trim-side"],
                        trim_side2=args["trim-side2"],
                        summary=args["summary"],
                        summary_format=summary_fmt_sym,
                        matching_algorithm=matching_algo_sym,
                        log=args["log"]
                    )

                end
            end
        else
            # --- File Mode (Original) ---
            if !isnothing(fastq2)
                if isdir(fastq2)
                    println(stderr, "Error: fastq1 is a file but fastq2 is a directory. Both must be files.")
                    return 1
                end

                # Paired End File
                execute_demultiplexing(
                    fastq1,
                    fastq2,
                    barcode_file,
                    output_dir;
                    barcode_file2=args["barcode-file2"],
                    output_prefix1=args["output-prefix1"],
                    output_prefix2=args["output-prefix2"],
                    gzip_output=gzip_val,
                    max_error_rate=args["max-error-rate"],
                    min_delta=args["min-delta"],
                    match=args["match"],
                    mismatch=args["mismatch"],
                    indel=args["indel"],
                    nindel=nindel_val,
                    classify_both=args["classify-both"],
                    bc_complement=args["bc-complement"],
                    bc_rev=args["bc-rev"],
                    ref_search_range=args["ref-search-range"],
                    barcode_start_range=args["barcode-start-range"],
                    barcode_end_range=args["barcode-end-range"],
                    ref_search_range2=args["ref-search-range2"],
                    barcode_start_range2=args["barcode-start-range2"],
                    barcode_end_range2=args["barcode-end-range2"],
                    chunk_size=args["chunk-size"],
                    channel_capacity=args["channel-capacity"],
                    trim_side=args["trim-side"],
                    trim_side2=args["trim-side2"],
                    summary=args["summary"],
                    summary_format=summary_fmt_sym,
                    matching_algorithm=matching_algo_sym,
                    log=args["log"]
                )
            else
                # Single End File
                execute_demultiplexing(
                    fastq1,
                    barcode_file,
                    output_dir;
                    barcode_file2=args["barcode-file2"],
                    output_prefix=args["output-prefix1"], # Use prefix1 as the main prefix for single end
                    gzip_output=gzip_val,
                    max_error_rate=args["max-error-rate"],
                    min_delta=args["min-delta"],
                    match=args["match"],
                    mismatch=args["mismatch"],
                    indel=args["indel"],
                    nindel=nindel_val,
                    bc_complement=args["bc-complement"],
                    bc_rev=args["bc-rev"],
                    ref_search_range=args["ref-search-range"],
                    barcode_start_range=args["barcode-start-range"],
                    barcode_end_range=args["barcode-end-range"],
                    ref_search_range2=args["ref-search-range2"],
                    barcode_start_range2=args["barcode-start-range2"],
                    barcode_end_range2=args["barcode-end-range2"],
                    chunk_size=args["chunk-size"],
                    channel_capacity=args["channel-capacity"],
                    trim_side=args["trim-side"],
                    trim_side2=args["trim-side2"],
                    summary=args["summary"],
                    summary_format=summary_fmt_sym,
                    matching_algorithm=matching_algo_sym,
                    log=args["log"]
                )
            end
        end

        return 0
    catch e
        showerror(stderr, e, catch_backtrace())
        return 1
    end
end
