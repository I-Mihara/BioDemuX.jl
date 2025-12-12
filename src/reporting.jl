function merge_stats(stats_list::Vector{DemuxStats})
    merged = DemuxStats()
    for s in stats_list
        merged.total_reads += s.total_reads
        merged.matched_reads += s.matched_reads
        merged.unmatched_reads += s.unmatched_reads
        merged.ambiguous_reads += s.ambiguous_reads

        merge!(+, merged.sample_counts, s.sample_counts)

        merge!(+, merged.bc1_pos_counts, s.bc1_pos_counts)
        merge!(+, merged.bc1_len_counts, s.bc1_len_counts)
        merge!(+, merged.bc1_score_counts, s.bc1_score_counts)

        for (bc_idx, counts) in s.bc1_per_bc_score_counts
            if !haskey(merged.bc1_per_bc_score_counts, bc_idx)
                merged.bc1_per_bc_score_counts[bc_idx] = Dict{Float64,Int}()
            end
            merge!(+, merged.bc1_per_bc_score_counts[bc_idx], counts)
        end
        for (bc_idx, counts) in s.bc1_per_bc_pos_counts
            if !haskey(merged.bc1_per_bc_pos_counts, bc_idx)
                merged.bc1_per_bc_pos_counts[bc_idx] = Dict{Int,Int}()
            end
            merge!(+, merged.bc1_per_bc_pos_counts[bc_idx], counts)
        end
        for (bc_idx, counts) in s.bc1_per_bc_len_counts
            if !haskey(merged.bc1_per_bc_len_counts, bc_idx)
                merged.bc1_per_bc_len_counts[bc_idx] = Dict{Int,Int}()
            end
            merge!(+, merged.bc1_per_bc_len_counts[bc_idx], counts)
        end

        merge!(+, merged.bc2_pos_counts, s.bc2_pos_counts)
        merge!(+, merged.bc2_len_counts, s.bc2_len_counts)
        merge!(+, merged.bc2_score_counts, s.bc2_score_counts)

        for (bc_idx, counts) in s.bc2_per_bc_score_counts
            if !haskey(merged.bc2_per_bc_score_counts, bc_idx)
                merged.bc2_per_bc_score_counts[bc_idx] = Dict{Float64,Int}()
            end
            merge!(+, merged.bc2_per_bc_score_counts[bc_idx], counts)
        end
        for (bc_idx, counts) in s.bc2_per_bc_pos_counts
            if !haskey(merged.bc2_per_bc_pos_counts, bc_idx)
                merged.bc2_per_bc_pos_counts[bc_idx] = Dict{Int,Int}()
            end
            merge!(+, merged.bc2_per_bc_pos_counts[bc_idx], counts)
        end
        for (bc_idx, counts) in s.bc2_per_bc_len_counts
            if !haskey(merged.bc2_per_bc_len_counts, bc_idx)
                merged.bc2_per_bc_len_counts[bc_idx] = Dict{Int,Int}()
            end
            merge!(+, merged.bc2_per_bc_len_counts[bc_idx], counts)
        end
    end
    return merged
end

function write_text_report(stats::DemuxStats, config::DemuxConfig, output_dir::String, fastq_path::String, bc_path::String, fastq_path2::Union{String,Nothing}=nothing, bc_path2::Union{String,Nothing}=nothing; bc_complement::Bool=false, bc_rev::Bool=false, trim_side::Union{Int,Nothing}=nothing, trim_side2::Union{Int,Nothing}=nothing, n_threads::Union{Int,Nothing}=nothing, duration::Union{Dates.Period,Nothing}=nothing)
    output_file = joinpath(output_dir, "summary.txt")
    mode = isfile(output_file) ? "a" : "w"

    open(output_file, mode) do io
        if mode == "a"
            println(io, "\n==================================================\n")
        end

        println(io, "BioDemuX Summary Report")
        println(io, "=======================")
        println(io, "Run Information:")
        println(io, "  Date: ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
        println(io, "  Input FASTQ(s): ", fastq_path, !isnothing(fastq_path2) ? ", " * fastq_path2 : "")
        println(io, "  Barcode File: ", bc_path)
        if !isnothing(bc_path2)
            println(io, "  Barcode File 2: ", bc_path2)
        end
        if !isnothing(n_threads)
            println(io, "  Threads: ", n_threads)
        end
        if !isnothing(duration)
            println(io, "  Duration: ", string(Dates.canonicalize(duration)))
        end
        println(io, "  Parameters:")
        println(io, "    Max Error Rate: ", config.max_error_rate)
        println(io, "    Min Delta: ", config.min_delta)
        println(io, "    Match: ", config.match, ", Mismatch: ", config.mismatch, ", Indel: ", config.indel)
        if !isnothing(config.nindel)
            println(io, "    N Indel: ", config.nindel)
        end
        if config.classify_both
            println(io, "    Classify Both: true")
        end
        if config.gzip_output
            println(io, "    Gzip Output: true")
        end
        if bc_complement
            println(io, "    BC Complement: true")
        end
        if bc_rev
            println(io, "    BC Reverse: true")
        end
        if !isnothing(trim_side)
            println(io, "    Trim Side: ", trim_side)
        end
        if !isnothing(trim_side2)
            println(io, "    Trim Side 2: ", trim_side2)
        end

        println(io, "")
        println(io, "Total Reads: ", stats.total_reads)
        println(io, "Matched Reads: ", stats.matched_reads, " (", round(stats.matched_reads / stats.total_reads * 100, digits=2), "%)")
        println(io, "Unmatched Reads: ", stats.unmatched_reads, " (", round(stats.unmatched_reads / stats.total_reads * 100, digits=2), "%)")
        println(io, "Ambiguous Reads: ", stats.ambiguous_reads, " (", round(stats.ambiguous_reads / stats.total_reads * 100, digits=2), "%)")
        println(io, "")
        println(io, "Barcode Counts:")

        # Sort by count descending
        sorted_counts = sort(collect(stats.sample_counts), by=x -> x[2], rev=true)

        for ((bc1_idx, bc2_idx), count) in sorted_counts
            name1 = config.ids[bc1_idx]
            name = name1
            if bc2_idx > 0
                name2 = config.ids2[bc2_idx]
                name = name1 * "." * name2
            end
            println(io, name, "\t", count, "\t", round(count / stats.total_reads * 100, digits=2), "%")
        end
    end
end

function write_json_report(stats::DemuxStats, config::DemuxConfig, output_dir::String, fastq_path::String, bc_path::String, fastq_path2::Union{String,Nothing}=nothing, bc_path2::Union{String,Nothing}=nothing; bc_complement::Bool=false, bc_rev::Bool=false, trim_side::Union{Int,Nothing}=nothing, trim_side2::Union{Int,Nothing}=nothing, n_threads::Union{Int,Nothing}=nothing, duration::Union{Dates.Period,Nothing}=nothing)
    # Simple JSON construction
    function dict_to_json(d)
        items = String[]
        for (k, v) in d
            key_str = isa(k, Tuple) ? string(k) : string(k)
            val_str = isa(v, Dict) ? dict_to_json(v) : string(v)
            push!(items, "\"$key_str\": $val_str")
        end
        return "{" * join(items, ", ") * "}"
    end

    # Build run info object
    run_info = """
    "run_info": {
        "date": "$(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))",
        "input_fastq": "$fastq_path",
        "input_fastq2": $(isnothing(fastq_path2) ? "null" : "\"$fastq_path2\""),
        "barcode_file": "$bc_path",
        "barcode_file2": $(isnothing(bc_path2) ? "null" : "\"$bc_path2\""),
        "threads": $(isnothing(n_threads) ? "null" : n_threads),
        "duration": $(isnothing(duration) ? "null" : "\"$(string(Dates.canonicalize(duration)))\""),
        "parameters": {
            "max_error_rate": $(config.max_error_rate),
            "min_delta": $(config.min_delta),
            "match": $(config.match),
            "mismatch": $(config.mismatch),
            "indel": $(config.indel),
            "nindel": $(isnothing(config.nindel) ? "null" : config.nindel),
            "classify_both": $(config.classify_both),
            "gzip_output": $(config.gzip_output),
            "bc_complement": $bc_complement,
            "bc_rev": $bc_rev,
            "trim_side": $(isnothing(trim_side) ? "null" : trim_side),
            "trim_side2": $(isnothing(trim_side2) ? "null" : trim_side2)
        }
    }
    """

    json_str = """
    {
        $run_info,
        "total_reads": $(stats.total_reads),
        "matched_reads": $(stats.matched_reads),
        "unmatched_reads": $(stats.unmatched_reads),
        "ambiguous_reads": $(stats.ambiguous_reads),
        "sample_counts": $(dict_to_json(stats.sample_counts)),
        "bc1_pos_counts": $(dict_to_json(stats.bc1_pos_counts)),
        "bc1_len_counts": $(dict_to_json(stats.bc1_len_counts)),
        "bc1_score_counts": $(dict_to_json(stats.bc1_score_counts)),
        "bc1_per_bc_score_counts": $(dict_to_json(stats.bc1_per_bc_score_counts)),
        "bc1_per_bc_pos_counts": $(dict_to_json(stats.bc1_per_bc_pos_counts)),
        "bc1_per_bc_len_counts": $(dict_to_json(stats.bc1_per_bc_len_counts)),
        "bc2_pos_counts": $(dict_to_json(stats.bc2_pos_counts)),
        "bc2_len_counts": $(dict_to_json(stats.bc2_len_counts)),
        "bc2_score_counts": $(dict_to_json(stats.bc2_score_counts)),
        "bc2_per_bc_score_counts": $(dict_to_json(stats.bc2_per_bc_score_counts)),
        "bc2_per_bc_pos_counts": $(dict_to_json(stats.bc2_per_bc_pos_counts)),
        "bc2_per_bc_len_counts": $(dict_to_json(stats.bc2_per_bc_len_counts))
    }
    """

    output_file = joinpath(output_dir, "summary.json")

    if isfile(output_file)
        # Append mode: Read existing, parse (simplistic), and append
        existing_content = read(output_file, String)
        existing_content = strip(existing_content)

        if startswith(existing_content, "[")
            # Already a list, remove closing ']' and append
            new_content = existing_content[1:end-1] * ", " * json_str * "]"
        else
            # Single object (legacy), convert to list
            new_content = "[" * existing_content * ", " * json_str * "]"
        end

        open(output_file, "w") do io
            write(io, new_content)
        end
    else
        # New file: Initialize as a list containing one object for consistency
        open(output_file, "w") do io
            write(io, "[" * json_str * "]")
        end
    end
end

function write_stdout_report(stats::DemuxStats, config::DemuxConfig, fastq_path::String, bc_path::String, fastq_path2::Union{String,Nothing}=nothing, bc_path2::Union{String,Nothing}=nothing; bc_complement::Bool=false, bc_rev::Bool=false, trim_side::Union{Int,Nothing}=nothing, trim_side2::Union{Int,Nothing}=nothing, n_threads::Union{Int,Nothing}=nothing, duration::Union{Dates.Period,Nothing}=nothing)
    println(stdout, "BioDemuX Summary Report")
    println(stdout, "=======================")
    println(stdout, "Run Information:")
    println(stdout, "  Date: ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
    println(stdout, "  Input FASTQ(s): ", fastq_path, !isnothing(fastq_path2) ? ", " * fastq_path2 : "")
    println(stdout, "  Barcode File: ", bc_path)
    if !isnothing(bc_path2)
        println(stdout, "  Barcode File 2: ", bc_path2)
    end
    if !isnothing(n_threads)
        println(stdout, "  Threads: ", n_threads)
    end
    if !isnothing(duration)
        println(stdout, "  Duration: ", string(Dates.canonicalize(duration)))
    end
    println(stdout, "  Parameters:")
    println(stdout, "    Max Error Rate: ", config.max_error_rate)
    println(stdout, "    Min Delta: ", config.min_delta)
    println(stdout, "    Match: ", config.match, ", Mismatch: ", config.mismatch, ", Indel: ", config.indel)
    if !isnothing(config.nindel)
        println(stdout, "    N Indel: ", config.nindel)
    end
    if config.classify_both
        println(stdout, "    Classify Both: true")
    end
    if config.gzip_output
        println(stdout, "    Gzip Output: true")
    end
    if bc_complement
        println(stdout, "    BC Complement: true")
    end
    if bc_rev
        println(stdout, "    BC Reverse: true")
    end
    if !isnothing(trim_side)
        println(stdout, "    Trim Side: ", trim_side)
    end
    if !isnothing(trim_side2)
        println(stdout, "    Trim Side 2: ", trim_side2)
    end

    println(stdout, "")
    println(stdout, "Total Reads: ", stats.total_reads)
    println(stdout, "Matched Reads: ", stats.matched_reads, " (", round(stats.matched_reads / stats.total_reads * 100, digits=2), "%)")
    println(stdout, "Unmatched Reads: ", stats.unmatched_reads, " (", round(stats.unmatched_reads / stats.total_reads * 100, digits=2), "%)")
    println(stdout, "Ambiguous Reads: ", stats.ambiguous_reads, " (", round(stats.ambiguous_reads / stats.total_reads * 100, digits=2), "%)")
    println(stdout, "")
    println(stdout, "Barcode Counts:")

    # Sort by count descending
    sorted_counts = sort(collect(stats.sample_counts), by=x -> x[2], rev=true)

    for ((bc1_idx, bc2_idx), count) in sorted_counts
        name1 = config.ids[bc1_idx]
        name = name1
        if bc2_idx > 0
            name2 = config.ids2[bc2_idx]
            name = name1 * "." * name2
        end
        println(stdout, name, "\t", count, "\t", round(count / stats.total_reads * 100, digits=2), "%")
    end
end

function write_html_report(stats::DemuxStats, config::DemuxConfig, output_dir::String, fastq_path::String, bc_path::String, fastq_path2::Union{String,Nothing}=nothing, bc_path2::Union{String,Nothing}=nothing; bc_complement::Bool=false, bc_rev::Bool=false, trim_side::Union{Int,Nothing}=nothing, trim_side2::Union{Int,Nothing}=nothing, n_threads::Union{Int,Nothing}=nothing, duration::Union{Dates.Period,Nothing}=nothing)
    # Rich HTML report with CSS/SVG

    # Helper to generate SVG bar chart (moved to top level)


    # Generate unique Run ID once
    run_id = string(Dates.now().instant.periods.value)

    # Build Parameter List
    param_list_html = """
        <ul style="margin-top: 8px; padding-left: 20px; list-style-type: disc;">
            <li>Max Error Rate: $(config.max_error_rate)</li>
            <li>Min Delta: $(config.min_delta)</li>
            <li>Match: $(config.match), Mismatch: $(config.mismatch), Indel: $(config.indel)</li>
    """

    # Conditional Parameters
    if !isnothing(config.nindel)
        param_list_html *= """<li>N Indel: $(config.nindel)</li>"""
    end
    if config.classify_both
        param_list_html *= """<li>Classify Both: true</li>"""
    end
    if config.gzip_output
        param_list_html *= """<li>Gzip Output: true</li>"""
    end
    if bc_complement
        param_list_html *= """<li>BC Complement: true</li>"""
    end
    if bc_rev
        param_list_html *= """<li>BC Reverse: true</li>"""
    end
    if !isnothing(trim_side)
        param_list_html *= """<li>Trim Side: $trim_side</li>"""
    end
    if !isnothing(trim_side2)
        param_list_html *= """<li>Trim Side 2: $trim_side2</li>"""
    end

    param_list_html *= """</ul>"""

    # Generate Report Content (The part that gets appended)
    report_body = """
            <div class="run-section">
            <div class="run-info">
                <h3>Run Information</h3>
                <ul>
                    <li><strong>Date:</strong> $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))</li>
                    <li><strong>Input FASTQ(s):</strong> $fastq_path $(!isnothing(fastq_path2) ? ", " * fastq_path2 : "")</li>
                    <li><strong>Barcode File:</strong> $bc_path</li>
                    $(!isnothing(bc_path2) ? "<li><strong>Barcode File 2:</strong> $bc_path2</li>" : "")
                    $(!isnothing(n_threads) ? "<li><strong>Threads:</strong> $n_threads</li>" : "")
                    $(!isnothing(duration) ? "<li><strong>Duration:</strong> $(string(Dates.canonicalize(duration)))</li>" : "")
                    <li><strong>Parameters:</strong>
                        $param_list_html
                    </li>
                </ul>
            </div>
            
            <div class="summary-grid">
                <div class="summary-card">
                    <h3>Total Reads</h3>
                    <div class="value">$(stats.total_reads)</div>
                </div>
                <div class="summary-card">
                    <h3>Matched</h3>
                    <div class="value" style="color: #2ea44f;">$(stats.matched_reads)</div>
                    <div class="percentage">$(round(stats.matched_reads / stats.total_reads * 100, digits=2))%</div>
                </div>
                <div class="summary-card">
                    <h3>Unmatched</h3>
                    <div class="value" style="color: #cf222e;">$(stats.unmatched_reads)</div>
                    <div class="percentage">$(round(stats.unmatched_reads / stats.total_reads * 100, digits=2))%</div>
                </div>
                <div class="summary-card">
                    <h3>Ambiguous</h3>
                    <div class="value" style="color: #d29922;">$(stats.ambiguous_reads)</div>
                    <div class="percentage">$(round(stats.ambiguous_reads / stats.total_reads * 100, digits=2))%</div>
                </div>
            </div>
            
            <h2>Barcode Statistics</h2>
            <table>
                <tr>
                    <th>Barcode</th>
                    <th>Count</th>
                    <th>Percentage</th>
                    <th>Distribution</th>
                </tr>
    """

    # Sort by count descending
    sorted_counts = sort(collect(stats.sample_counts), by=x -> x[2], rev=true)
    max_count = isempty(sorted_counts) ? 0 : sorted_counts[1][2]

    for ((bc1_idx, bc2_idx), count) in sorted_counts
        name1 = config.ids[bc1_idx]
        name = name1
        if bc2_idx > 0
            name2 = config.ids2[bc2_idx]
            name = name1 * "." * name2
        end

        pct = round(count / stats.total_reads * 100, digits=2)
        bar_width = round(count / max_count * 100, digits=1)

        report_body *= """
        <tr>
            <td>$name</td>
            <td>$count</td>
            <td>$pct%</td>
            <td style="width: 40%;">
                <div class="progress-bar-container">
                    <div class="progress-bar" style="width: $(bar_width)%"></div>
                    <div class="progress-label">$count</div>
                </div>
            </td>
        </tr>
        """
    end

    report_body *= """
            </table>
            
            <h2>Detailed Distributions</h2>
            
            <div class="tabs">
                <div class="tab active" onclick="openTab(event, 'Global_$run_id')">Global Stats</div>
                <div class="tab" onclick="openTab(event, 'PerBarcode_$run_id')">Per-Barcode Stats</div>
            </div>
            
            <div id="Global_$run_id" class="tab-content active">
                <h3>Global Stats</h3>
                <div class="summary-grid">
                    <div class="chart-container" style="grid-column: span 3;">
                        $(generate_histogram(stats.bc1_pos_counts, "Barcode 1 Start Position Distribution", "Start Position"))
                    </div>
                    <div class="chart-container" style="grid-column: span 3;">
                        $(generate_histogram(stats.bc1_len_counts, "Barcode 1 Length Distribution", "Length"))
                    </div>
                    <div class="chart-container" style="grid-column: span 3;">
                        $(generate_histogram(stats.bc1_score_counts, "Barcode 1 Score Distribution", "Score"))
                    </div>
                </div>
    """

    if config.is_dual
        report_body *= """
                <h3>Global Stats (Barcode 2)</h3>
                <div class="summary-grid">
                    <div class="chart-container" style="grid-column: span 3;">
                        $(generate_histogram(stats.bc2_pos_counts, "Barcode 2 Start Position Distribution", "Start Position"))
                    </div>
                    <div class="chart-container" style="grid-column: span 3;">
                        $(generate_histogram(stats.bc2_len_counts, "Barcode 2 Length Distribution", "Length"))
                    </div>
                    <div class="chart-container" style="grid-column: span 3;">
                        $(generate_histogram(stats.bc2_score_counts, "Barcode 2 Score Distribution", "Score"))
                    </div>
                </div>
        """
    end

    report_body *= """
            </div>
            
            <div id="PerBarcode_$run_id" class="tab-content">
                <div style="margin-bottom: 20px;">
                    <label for="barcode-select-$run_id" style="font-weight: 600; margin-right: 10px; color: #24292f;">Select Barcode:</label>
                    <select id="barcode-select-$run_id" onchange="showBarcodeStats(this, '$run_id')">
                        <option value="">-- Select a Barcode --</option>
    """

    unique_bc1s = sort(unique([x[1] for x in keys(stats.sample_counts)]))
    for bc1_idx in unique_bc1s
        bc_name = config.ids[bc1_idx]
        report_body *= """<option value="bc1_$bc1_idx">$bc_name</option>"""
    end

    report_body *= """
                    </select>
                </div>
    """

    for bc1_idx in unique_bc1s
        bc_name = config.ids[bc1_idx]

        # Extract stats for this barcode
        score_data = get(stats.bc1_per_bc_score_counts, bc1_idx, Dict{Float64,Int}())
        pos_data = get(stats.bc1_per_bc_pos_counts, bc1_idx, Dict{Int,Int}())
        len_data = get(stats.bc1_per_bc_len_counts, bc1_idx, Dict{Int,Int}())

        report_body *= """
            <div id="stats-bc1_$bc1_idx-$run_id" class="barcode-stats-$run_id" style="display: none;">
                <h3>Stats for $bc_name</h3>
                <div class="summary-grid">
                    <div class="chart-container" style="margin-top: 10px; min-width: 0;">
                        $(generate_histogram(score_data, "Score", "Score"))
                    </div>
                    <div class="chart-container" style="margin-top: 10px; min-width: 0;">
                        $(generate_histogram(pos_data, "Start Position", "Start Position"))
                    </div>
                    <div class="chart-container" style="margin-top: 10px; min-width: 0;">
                        $(generate_histogram(len_data, "Length", "Length"))
                    </div>
                </div>
                </div>
        """
    end

    if config.is_dual
        report_body *= """
        <div style="margin-top: 40px; margin-bottom: 20px;">
           <label for="barcode2-select-$run_id" style="font-weight: 600; margin-right: 10px; color: #24292f;">Select Barcode 2:</label>
           <select id="barcode2-select-$run_id" onchange="showBarcodeStats(this, '$run_id')">
               <option value="">-- Select a Barcode 2 --</option>
        """
        unique_bc2s = sort(unique([x[2] for x in keys(stats.sample_counts) if x[2] > 0]))
        for bc2_idx in unique_bc2s
            bc_name = config.ids2[bc2_idx]
            report_body *= """<option value="bc2_$bc2_idx">$bc_name</option>"""
        end
        report_body *= """</select></div>"""

        for bc2_idx in unique_bc2s
            bc_name = config.ids2[bc2_idx]
            score_data = get(stats.bc2_per_bc_score_counts, bc2_idx, Dict{Float64,Int}())
            pos_data = get(stats.bc2_per_bc_pos_counts, bc2_idx, Dict{Int,Int}())
            len_data = get(stats.bc2_per_bc_len_counts, bc2_idx, Dict{Int,Int}())

            report_body *= """
                <div id="stats-bc2_$bc2_idx-$run_id" class="barcode-stats-$run_id" style="display: none;">
                    <h3>Stats for $bc_name (Barcode 2)</h3>
                    <div class="summary-grid">
                        <div class="chart-container" style="margin-top: 10px; min-width: 0;">
                            $(generate_histogram(score_data, "Score", "Score"))
                        </div>
                        <div class="chart-container" style="margin-top: 10px; min-width: 0;">
                            $(generate_histogram(pos_data, "Start Position", "Start Position"))
                        </div>
                        <div class="chart-container" style="margin-top: 10px; min-width: 0;">
                            $(generate_histogram(len_data, "Length", "Length"))
                        </div>
                    </div>
                </div>
            """
        end
    end

    report_body *= """
            </div>
            </div>
            <hr style="margin: 60px 0; border: 0; border-top: 2px dashed #e1e4e8;">
    """

    output_file = joinpath(output_dir, "summary.html")

    if isfile(output_file)
        # Append mode
        existing_content = read(output_file, String)

        wrapped_body = """
        <div class="container" style="margin-top: 40px;">
        $report_body
        </div>
        """
        new_content = replace(existing_content, "</body>" => wrapped_body * "\n</body>")

        open(output_file, "w") do io
            write(io, new_content)
        end
    else
        # New file mode
        open(output_file, "w") do io
            write(io, HTML_REPORT_TEMPLATE_START)
            write(io, report_body)
            write(io, HTML_REPORT_TEMPLATE_END)
        end
    end
end

function generate_summary_report(stats::DemuxStats, config::DemuxConfig, output_dir::String, fastq_path::String, bc_path::String, fastq_path2::Union{String,Nothing}=nothing, bc_path2::Union{String,Nothing}=nothing; bc_complement::Bool=false, bc_rev::Bool=false, trim_side::Union{Int,Nothing}=nothing, trim_side2::Union{Int,Nothing}=nothing, n_threads::Union{Int,Nothing}=nothing, duration::Union{Dates.Period,Nothing}=nothing)
    if config.summary_format == :json
        write_json_report(stats, config, output_dir, fastq_path, bc_path, fastq_path2, bc_path2; bc_complement=bc_complement, bc_rev=bc_rev, trim_side=trim_side, trim_side2=trim_side2, n_threads=n_threads, duration=duration)
    elseif config.summary_format == :html
        write_html_report(stats, config, output_dir, fastq_path, bc_path, fastq_path2, bc_path2; bc_complement=bc_complement, bc_rev=bc_rev, trim_side=trim_side, trim_side2=trim_side2, n_threads=n_threads, duration=duration)
    elseif config.summary_format == :stdout
        write_stdout_report(stats, config, fastq_path, bc_path, fastq_path2, bc_path2; bc_complement=bc_complement, bc_rev=bc_rev, trim_side=trim_side, trim_side2=trim_side2, n_threads=n_threads, duration=duration)
    else
        write_text_report(stats, config, output_dir, fastq_path, bc_path, fastq_path2, bc_path2; bc_complement=bc_complement, bc_rev=bc_rev, trim_side=trim_side, trim_side2=trim_side2, n_threads=n_threads, duration=duration)
    end
end
