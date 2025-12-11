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
        # Append mode: Read existing, parse (manually if needed, but here we assume structure), and append
        # Since we don't want to depend on JSON.jl, we'll do simple string manipulation.
        # If it starts with '[', it's a list. If '{', it's a single object.
        existing_content = read(output_file, String)
        existing_content = strip(existing_content)

        if startswith(existing_content, "[")
            # Already a list, remove closing ']' and append
            new_content = existing_content[1:end-1] * ", " * json_str * "]"
        else
            # Single object, convert to list
            new_content = "[" * existing_content * ", " * json_str * "]"
        end

        open(output_file, "w") do io
            write(io, new_content)
        end
    else
        # New file: Write as a list containing one object for consistency? 
        # Or write as single object and convert later?
        # The user approved "list of objects". So let's write as a list even for the first one?
        # Or maybe keep single object for backward compatibility if only one run?
        # The prompt said "When appending... root structure will change".
        # So first run = single object. Second run = list of objects.
        open(output_file, "w") do io
            write(io, json_str)
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

    # Helper to generate SVG bar chart
    function generate_histogram(data::Dict{T,Int}, title::String, x_label::String) where T
        if isempty(data)
            return "<p>No data for $title</p>"
        end

        sorted_keys = sort(collect(keys(data)))
        values = [data[k] for k in sorted_keys]
        max_val = maximum(values)

        # SVG dimensions
        width = 600
        height = 300
        margin_left = 60
        margin_bottom = 60
        bar_width = (width - margin_left) / length(values) * 0.8

        # Color palette for bars (gradient-like)
        bar_color = "#3498db"

        svg = """<div class="chart-container"><h3>$title</h3><svg width="$width" height="$height" viewBox="0 0 $width $height">"""

        # Grid lines
        for i in 0:5
            y = height - margin_bottom - (i / 5) * (height - margin_bottom - 20)
            svg *= """<line x1="$margin_left" y1="$y" x2="$width" y2="$y" stroke="#eee" stroke-width="1" />"""
            # Y-axis labels
            label_val = round(Int, (i / 5) * max_val)
            svg *= """<text x="$(margin_left - 10)" y="$(y + 5)" font-size="10" text-anchor="end" fill="#7f8c8d">$label_val</text>"""
        end

        # Y-axis title
        svg *= """<text transform="rotate(-90)" x="$(-(height - margin_bottom)/2)" y="15" font-size="12" text-anchor="middle" fill="#2c3e50" font-weight="bold">Read Count</text>"""

        # Y-axis line
        svg *= """<line x1="$margin_left" y1="0" x2="$margin_left" y2="$(height - margin_bottom)" stroke="#ccc" />"""
        # X-axis line
        svg *= """<line x1="$margin_left" y1="$(height - margin_bottom)" x2="$width" y2="$(height - margin_bottom)" stroke="#ccc" />"""

        for (i, val) in enumerate(values)
            bar_height = (val / max_val) * (height - margin_bottom - 20)
            x = margin_left + (i - 1) * ((width - margin_left) / length(values)) + ((width - margin_left) / length(values) - bar_width) / 2
            y = height - margin_bottom - bar_height

            # Bar with hover effect (using simple title for now, CSS hover can be added)
            svg *= """<rect x="$x" y="$y" width="$bar_width" height="$bar_height" fill="$bar_color" class="bar"><title>$(sorted_keys[i]): $val</title></rect>"""

            # X-axis label
            if length(values) <= 20 || i % (div(length(values), 20) + 1) == 0
                label_x = x + bar_width / 2
                label_y = height - margin_bottom + 15
                svg *= """<text x="$label_x" y="$label_y" font-size="10" text-anchor="middle" fill="#7f8c8d">$(sorted_keys[i])</text>"""
            end
        end

        # X-axis title
        svg *= """<text x="$(width/2 + margin_left/2)" y="$(height - 10)" font-size="12" text-anchor="middle" fill="#2c3e50" font-weight="bold">$x_label</text>"""

        svg *= "</svg></div>"
        return svg
    end

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
                <div class="tab active" onclick="openTab(event, 'Global_$(Dates.now().instant.periods.value)')">Global Stats</div>
                <div class="tab" onclick="openTab(event, 'PerBarcode_$(Dates.now().instant.periods.value)')">Per-Barcode Stats</div>
            </div>
            
            <div id="Global_$(Dates.now().instant.periods.value)" class="tab-content active">
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
            
            <div id="PerBarcode_$(Dates.now().instant.periods.value)" class="tab-content">
                <div style="margin-bottom: 20px;">
                    <label for="barcode-select-$(Dates.now().instant.periods.value)" style="font-weight: 600; margin-right: 10px; color: #24292f;">Select Barcode:</label>
                    <select id="barcode-select-$(Dates.now().instant.periods.value)" onchange="showBarcodeStats(this, '$(Dates.now().instant.periods.value)')">
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
            <div id="stats-bc1_$bc1_idx-$(Dates.now().instant.periods.value)" class="barcode-stats-$(Dates.now().instant.periods.value)" style="display: none;">
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
           <label for="barcode2-select-$(Dates.now().instant.periods.value)" style="font-weight: 600; margin-right: 10px; color: #24292f;">Select Barcode 2:</label>
           <select id="barcode2-select-$(Dates.now().instant.periods.value)" onchange="showBarcodeStats(this, '$(Dates.now().instant.periods.value)')">
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
                <div id="stats-bc2_$bc2_idx-$(Dates.now().instant.periods.value)" class="barcode-stats-$(Dates.now().instant.periods.value)" style="display: none;">
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
        # Insert before closing body tag
        # We need to make sure we don't duplicate headers or scripts if they are static.
        # But our script is simple.
        # To support multiple runs, we need unique IDs for tabs and selects.
        # I added timestamp/ID to IDs above.

        # We also need to update the script to handle dynamic IDs or be generic.
        # The script functions openTab and showBarcodeStats need to be generic.
        # I'll update the script in the header to be generic.

        # Insert before </body>
        new_content = replace(existing_content, "</body>" => report_body * "\n</body>")

        open(output_file, "w") do io
            write(io, new_content)
        end
    else
        # New file mode
        html_header = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>BioDemuX Report</title>
            <style>
                body { font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background-color: #f0f2f5; color: #1a1a1a; margin: 0; padding: 40px; }
                .container { max_width: 1200px; margin: 0 auto; background: white; padding: 40px; border-radius: 12px; box-shadow: 0 4px 20px rgba(0,0,0,0.05); }
                h1 { color: #1a1a1a; border-bottom: 2px solid #f0f2f5; padding-bottom: 20px; margin-bottom: 30px; font-weight: 700; letter-spacing: -0.5px; }
                h2 { color: #2c3e50; margin-top: 40px; font-weight: 600; font-size: 1.5em; }
                h3 { color: #57606a; font-size: 1em; font-weight: 600; text-transform: uppercase; letter-spacing: 0.5px; }
                
                .summary-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 24px; margin-bottom: 40px; }
                .summary-card { background: white; padding: 24px; border-radius: 12px; border: 1px solid #e1e4e8; box-shadow: 0 2px 8px rgba(0,0,0,0.02); transition: transform 0.2s; }
                .summary-card:hover { transform: translateY(-2px); box-shadow: 0 4px 12px rgba(0,0,0,0.05); }
                .summary-card h3 { margin: 0 0 12px 0; font-size: 0.85em; color: #6e7781; }
                .summary-card .value { font-size: 2.5em; font-weight: 700; color: #1a1a1a; line-height: 1; margin-bottom: 8px; }
                .summary-card .percentage { font-size: 1em; color: #57606a; font-weight: 500; }
                
                table { border-collapse: separate; border-spacing: 0; width: 100%; margin-top: 20px; border-radius: 8px; overflow: hidden; border: 1px solid #e1e4e8; }
                th, td { padding: 16px 20px; text-align: left; border-bottom: 1px solid #e1e4e8; }
                th { background-color: #f6f8fa; font-weight: 600; color: #24292f; font-size: 0.9em; text-transform: uppercase; letter-spacing: 0.5px; }
                tr:last-child td { border-bottom: none; }
                tr:hover td { background-color: #f8f9fa; }
                
                .progress-bar-container { width: 100%; background-color: #e1e4e8; border-radius: 100px; height: 8px; position: relative; overflow: hidden; }
                .progress-bar { height: 100%; background-color: #3498db; border-radius: 100px; }
                .progress-label { display: none; } /* Hide label inside bar, show count in table cell */
                
                .chart-container { margin-top: 24px; padding: 24px; background: white; border: 1px solid #e1e4e8; border-radius: 12px; overflow-x: auto; }
                .bar:hover { opacity: 0.8; cursor: pointer; }
                
                .run-info { background: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 30px; border: 1px solid #e1e4e8; }
                .run-info h3 { margin-top: 0; color: #2c3e50; }
                .run-info ul { list-style: none; padding: 0; margin: 0; }
                .run-info li { margin-bottom: 8px; font-size: 0.95em; color: #57606a; }
                .run-info li strong { color: #24292f; font-weight: 600; margin-right: 8px; }

                .tabs { display: flex; border-bottom: 1px solid #e1e4e8; margin-top: 40px; gap: 8px; }
                .tab { padding: 12px 24px; cursor: pointer; background: transparent; border: none; border-bottom: 2px solid transparent; font-weight: 500; color: #57606a; transition: all 0.2s; }
                .tab:hover { color: #24292f; }
                .tab.active { border-bottom: 2px solid #3498db; color: #3498db; font-weight: 600; }
                .tab-content { display: none; padding: 30px 0; animation: fadeIn 0.3s ease; }
                .tab-content.active { display: block; }
                
                @keyframes fadeIn { from { opacity: 0; transform: translateY(10px); } to { opacity: 1; transform: translateY(0); } }

                select { padding: 8px 12px; border: 1px solid #e1e4e8; border-radius: 6px; font-size: 14px; color: #24292f; background-color: white; min-width: 200px; }
            </style>
            <script>
                function openTab(evt, tabName) {
                    // Find the parent container of the tabs to scope the search?
                    // Or just use unique IDs. I used unique IDs.
                    
                    // We need to handle multiple tab groups on one page.
                    // The 'tabName' is unique (e.g., Global_12345).
                    // But we need to close other tabs in the SAME group.
                    // How to identify the group?
                    // The tab buttons are siblings.
                    
                    var tabButton = evt.currentTarget;
                    var tabContainer = tabButton.parentElement;
                    var allTabs = tabContainer.getElementsByClassName("tab");
                    
                    // Remove active from all tabs in this container
                    for (var i = 0; i < allTabs.length; i++) {
                        allTabs[i].classList.remove("active");
                    }
                    
                    // Add active to clicked tab
                    tabButton.classList.add("active");
                    
                    // Hide all tab contents that correspond to this group
                    // This is tricky. We can infer the group ID from the tabName.
                    // tabName format: Name_ID
                    var parts = tabName.split("_");
                    var id = parts[parts.length-1];
                    
                    // Hide Global_ID and PerBarcode_ID
                    document.getElementById("Global_" + id).style.display = "none";
                    document.getElementById("PerBarcode_" + id).style.display = "none";
                    document.getElementById("Global_" + id).classList.remove("active");
                    document.getElementById("PerBarcode_" + id).classList.remove("active");
                    
                    // Show selected
                    document.getElementById(tabName).style.display = "block";
                    document.getElementById(tabName).classList.add("active");
                }

                function showBarcodeStats(select, runId) {
                    var bcId = select.value;
                    // Hide all stats in this run
                    var statsDivs = document.getElementsByClassName("barcode-stats-" + runId);
                    for (var i = 0; i < statsDivs.length; i++) {
                        statsDivs[i].style.display = "none";
                    }
                    if (bcId) {
                        document.getElementById("stats-" + bcId + "-" + runId).style.display = "block";
                    }
                }
            </script>
        </head>
        <body>
            <div class="container">
                <h1>BioDemuX Report</h1>
                $report_body
            </div>
        </body>
        </html>
        """

        open(output_file, "w") do io
            write(io, html_header)
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

function determine_filename_and_stats(seq::String, config::DemuxConfig, ws::SemiGlobalWorkspace, stats::DemuxStats)
    n = ncodeunits(seq)
    stats.total_reads += 1

    # --- Pass 1: Barcode 1 ---
    ref_search_range1 = resolve(config.ref_search_range, n)
    barcode_start_range1 = resolve(config.barcode_start_range, n)
    barcode_end_range1 = resolve(config.barcode_end_range, n)

    start_j1 = max(first(ref_search_range1), first(barcode_start_range1), 1)
    end_j1 = min(last(ref_search_range1), last(barcode_end_range1), n)
    max_start_pos1 = last(barcode_start_range1)
    min_end_pos1 = first(barcode_end_range1)

    if start_j1 > end_j1 || start_j1 > max_start_pos1 || end_j1 < min_end_pos1
        stats.unmatched_reads += 1
        return "unknown.fastq", -1, -1
    end

    ref_search_range1 = start_j1:end_j1

    min_score_bc1, score1, delta1, best_start1, best_end1 = find_best_matching_bc(seq, config.bc_seqs, config.bc_lengths_no_N, config, ws, ref_search_range1, max_start_pos1, min_end_pos1, config.trim_side, true)

    if min_score_bc1 == 0
        stats.unmatched_reads += 1
        return "unknown.fastq", -1, -1
    elseif delta1 < config.min_delta
        stats.ambiguous_reads += 1
        return "ambiguous_classification.fastq", -1, -1
    end

    # --- Pass 2: Barcode 2 (if dual) ---
    min_score_bc2 = 0
    score2 = 0.0
    best_start2 = -1
    best_end2 = -1

    if config.is_dual

        ref_search_range2 = resolve(config.ref_search_range2, n)
        barcode_start_range2 = resolve(config.barcode_start_range2, n)
        barcode_end_range2 = resolve(config.barcode_end_range2, n)

        start_j2 = max(first(ref_search_range2), first(barcode_start_range2), 1)
        end_j2 = min(last(ref_search_range2), last(barcode_end_range2), n)
        max_start_pos2 = last(barcode_start_range2)
        min_end_pos2 = first(barcode_end_range2)

        if start_j2 > end_j2 || start_j2 > max_start_pos2 || end_j2 < min_end_pos2
            stats.unmatched_reads += 1
            return "unknown.fastq", -1, -1
        end

        ref_search_range2 = start_j2:end_j2

        min_score_bc2, score2, delta2, best_start2, best_end2 = find_best_matching_bc(seq, config.bc_seqs2, config.bc_lengths_no_N2, config, ws, ref_search_range2, max_start_pos2, min_end_pos2, config.trim_side2, true)

        if min_score_bc2 == 0
            stats.unmatched_reads += 1
            return "unknown.fastq", -1, -1
        elseif delta2 < config.min_delta
            stats.ambiguous_reads += 1
            return "ambiguous_classification.fastq", -1, -1
        end

        # Combine IDs
        output_filename = string(config.ids[min_score_bc1]) * "." * string(config.ids2[min_score_bc2]) * ".fastq"
    else
        # Single mode
        output_filename = string(config.ids[min_score_bc1]) * ".fastq"
    end

    # If we reached here, it's a match
    stats.matched_reads += 1

    # Update stats
    key = (min_score_bc1, min_score_bc2)
    stats.sample_counts[key] = get(stats.sample_counts, key, 0) + 1

    stats.bc1_pos_counts[best_start1] = get(stats.bc1_pos_counts, best_start1, 0) + 1
    len1 = best_end1 - best_start1 + 1
    stats.bc1_len_counts[len1] = get(stats.bc1_len_counts, len1, 0) + 1
    stats.bc1_score_counts[score1] = get(stats.bc1_score_counts, score1, 0) + 1

    if !haskey(stats.bc1_per_bc_score_counts, min_score_bc1)
        stats.bc1_per_bc_score_counts[min_score_bc1] = Dict{Float64,Int}()
        stats.bc1_per_bc_pos_counts[min_score_bc1] = Dict{Int,Int}()
        stats.bc1_per_bc_len_counts[min_score_bc1] = Dict{Int,Int}()
    end
    stats.bc1_per_bc_score_counts[min_score_bc1][score1] = get(stats.bc1_per_bc_score_counts[min_score_bc1], score1, 0) + 1
    stats.bc1_per_bc_pos_counts[min_score_bc1][best_start1] = get(stats.bc1_per_bc_pos_counts[min_score_bc1], best_start1, 0) + 1
    stats.bc1_per_bc_len_counts[min_score_bc1][len1] = get(stats.bc1_per_bc_len_counts[min_score_bc1], len1, 0) + 1

    if config.is_dual
        stats.bc2_pos_counts[best_start2] = get(stats.bc2_pos_counts, best_start2, 0) + 1
        len2 = best_end2 - best_start2 + 1
        stats.bc2_len_counts[len2] = get(stats.bc2_len_counts, len2, 0) + 1
        stats.bc2_score_counts[score2] = get(stats.bc2_score_counts, score2, 0) + 1

        if !haskey(stats.bc2_per_bc_score_counts, min_score_bc2)
            stats.bc2_per_bc_score_counts[min_score_bc2] = Dict{Float64,Int}()
            stats.bc2_per_bc_pos_counts[min_score_bc2] = Dict{Int,Int}()
            stats.bc2_per_bc_len_counts[min_score_bc2] = Dict{Int,Int}()
        end
        stats.bc2_per_bc_score_counts[min_score_bc2][score2] = get(stats.bc2_per_bc_score_counts[min_score_bc2], score2, 0) + 1
        stats.bc2_per_bc_pos_counts[min_score_bc2][best_start2] = get(stats.bc2_per_bc_pos_counts[min_score_bc2], best_start2, 0) + 1
        stats.bc2_per_bc_len_counts[min_score_bc2][len2] = get(stats.bc2_per_bc_len_counts[min_score_bc2], len2, 0) + 1
    end

    trim_range = -1:-1

    # Calculate trim range for BC1
    trim_range1 = -1:-1
    if !isnothing(config.trim_side)
        if config.trim_side == 3
            safe_start = max(1, best_start1)
            trim_range1 = 1:(safe_start-1)
        elseif config.trim_side == 5
            trim_range1 = (best_end1+1):n
        end
    end

    # Calculate trim range for BC2 (if dual)
    trim_range2 = -1:-1
    if config.is_dual && !isnothing(config.trim_side2)
        if config.trim_side2 == 3
            safe_start2 = max(1, best_start2)
            trim_range2 = 1:(safe_start2-1)
        elseif config.trim_side2 == 5
            trim_range2 = (best_end2+1):n
        end
    end

    # Combine trim ranges
    if first(trim_range1) != -1 && first(trim_range2) != -1
        trim_range = intersect(trim_range1, trim_range2)
    elseif first(trim_range1) != -1
        trim_range = trim_range1
    elseif first(trim_range2) != -1
        trim_range = trim_range2
    end

    if config.gzip_output
        return output_filename * ".gz", first(trim_range), last(trim_range)
    else
        return output_filename, first(trim_range), last(trim_range)
    end
end
