struct SemiGlobalWorkspace
    DP::Vector{Int}
    origin::Union{Vector{Int},Nothing}
end
SemiGlobalWorkspace(max_m::Int) = SemiGlobalWorkspace(Vector{Int}(undef, max_m), nothing)
SemiGlobalWorkspace(max_m::Int, trim::Bool) = SemiGlobalWorkspace(Vector{Int}(undef, max_m), trim ? Vector{Int}(undef, max_m) : nothing)
const INF_INT = typemax(Int) ÷ 4

struct DynamicRange
    start_offset::Int
    start_from_end::Bool
    end_offset::Int
    end_from_end::Bool
end

Base.@kwdef struct DemuxConfig
    max_error_rate::Float64 = 0.2
    min_delta::Float64 = 0.0
    match::Int = 0
    mismatch::Int = 1
    indel::Int = 1
    nindel::Union{Int,Nothing} = nothing
    classify_both::Bool = false
    gzip_output::Bool = false

    # Range 1
    ref_search_range::DynamicRange = parse_dynamic_range("1:end")
    barcode_start_range::DynamicRange = parse_dynamic_range("1:end")
    barcode_end_range::DynamicRange = parse_dynamic_range("1:end")
    bc_seqs::Vector{String}
    bc_lengths_no_N::Vector{Int}
    ids::Vector{String}

    # Dual Index Support
    is_dual::Bool = false
    ref_search_range2::DynamicRange = parse_dynamic_range("1:end")
    barcode_start_range2::DynamicRange = parse_dynamic_range("1:end")
    barcode_end_range2::DynamicRange = parse_dynamic_range("1:end")
    bc_seqs2::Vector{String} = String[]
    bc_lengths_no_N2::Vector{Int} = Int[]
    ids2::Vector{String} = String[]

    # Trimming
    trim_side::Union{Int,Nothing} = nothing
    trim_side2::Union{Int,Nothing} = nothing

    # Summary
    summary::Bool = false
    summary_format::Symbol = :txt # :txt, :json, :html
end

function parse_part(s::String)
    s = strip(s)
    from_end = occursin("end", s)
    if from_end
        s = replace(s, "end" => "0")
    end

    # Handle simple arithmetic
    val = 0
    if occursin("-", s)
        p = split(s, '-')
        val = parse(Int, strip(p[1])) - parse(Int, strip(p[2]))
    elseif occursin("+", s)
        p = split(s, '+')
        val = parse(Int, strip(p[1])) + parse(Int, strip(p[2]))
    else
        val = parse(Int, s)
    end

    return val, from_end
end

function parse_dynamic_range(range_str::String)
    parts = split(range_str, ':')
    if length(parts) != 2
        error("Invalid range format: $range_str. Expected 'start:end'.")
    end


    start_offset, start_from_end = parse_part(String(parts[1]))
    end_offset, end_from_end = parse_part(String(parts[2]))

    return DynamicRange(start_offset, start_from_end, end_offset, end_from_end)
end

function resolve(dr::DynamicRange, len::Int)
    s = dr.start_from_end ? len + dr.start_offset : dr.start_offset
    e = dr.end_from_end ? len + dr.end_offset : dr.end_offset
    return max(1, s):min(len, e)
end

"""
This function aligns `query` to `ref`, using semiglobal alignment algorithm. 
# Returns
An alignment score as a float, where lower values indicate better alignment.
"""
abstract type AbstractScoring end

struct SimpleScoring <: AbstractScoring
    match::Int
    mismatch::Int
    indel::Int
end

struct NScoring <: AbstractScoring
    match::Int
    mismatch::Int
    indel::Int
    nindel::Int
end

abstract type AbstractOutputPolicy end

struct ScoreOnly <: AbstractOutputPolicy end

struct TracebackOutput <: AbstractOutputPolicy
    trim_side::Union{Int,Nothing}
end

@inline function init_result(::ScoreOnly, INF_INT)
    return INF_INT
end

@inline function init_result(::TracebackOutput, INF_INT)
    return (INF_INT, -1, -1)
end

@inline function update_result(::ScoreOnly, current_best, new_score, j)
    return min(current_best, new_score)
end

@inline function update_result(output::TracebackOutput, current_best, new_score, j, start_pos)
    (best_score, best_start, best_end) = current_best

    if new_score < best_score
        return (new_score, start_pos, j)
    elseif new_score == best_score
        if !isnothing(output.trim_side) && output.trim_side == 3 && start_pos > best_start
            return (new_score, start_pos, j)
        end
    end
    return current_best
end

@inline function finalize_result(::ScoreOnly, result, normalization)
    if result >= INF_INT
        return Inf
    end
    return result / normalization
end

@inline function finalize_result(::TracebackOutput, result, normalization)
    (score, start, end_pos) = result
    if score >= INF_INT
        return (Inf, start, end_pos)
    end
    return (score / normalization, start, end_pos)
end

@inline function max_indel_steps_for(scoring::SimpleScoring, allowed_error::Int)
    return div(allowed_error, scoring.indel)
end

@inline function max_indel_steps_for(scoring::NScoring, allowed_error::Int)
    return div(allowed_error, min(scoring.indel, scoring.nindel))
end

Base.@propagate_inbounds function step_scores_main(scoring::SimpleScoring, q::Base.CodeUnits{UInt8,String}, r::Base.CodeUnits{UInt8,String}, i::Int, j::Int, previous_score::Int, DP::Vector{Int}, m::Int, INF_INT::Int)
    indel = scoring.indel
    match = scoring.match
    mismatch = scoring.mismatch

    insertion_score = DP[i] + indel
    deletion_score = previous_score + indel
    substitution_score = DP[i-1] + (q[i] == r[j] ? match : mismatch)

    return insertion_score, deletion_score, substitution_score
end

Base.@propagate_inbounds function step_scores_main(scoring::NScoring, q::Base.CodeUnits{UInt8,String}, r::Base.CodeUnits{UInt8,String}, i::Int, j::Int, previous_score::Int, DP::Vector{Int}, m::Int, INF_INT::Int)
    indel = scoring.indel
    nindel = scoring.nindel
    match = scoring.match
    mismatch = scoring.mismatch

    is_N_q = q[i] == UInt8('N')
    cost = is_N_q ? nindel : indel

    insertion_score = DP[i] + cost
    deletion_score = previous_score + cost

    is_match = (q[i] == r[j] || is_N_q)
    substitution_score = DP[i-1] + (is_match ? match : mismatch)

    return insertion_score, deletion_score, substitution_score
end

Base.@propagate_inbounds function step_scores(scoring::SimpleScoring, q::Base.CodeUnits{UInt8,String}, r::Base.CodeUnits{UInt8,String}, i::Int, j::Int, previous_score::Int, DP::Vector{Int}, m::Int, INF_INT::Int)
    indel = scoring.indel
    match = scoring.match
    mismatch = scoring.mismatch

    insertion_score = (i == m ? INF_INT : DP[i] + indel)
    deletion_score = previous_score + indel
    substitution_score = (i == 1 ? 0 : DP[i-1]) + (q[i] == r[j] ? match : mismatch)

    return insertion_score, deletion_score, substitution_score
end

Base.@propagate_inbounds function step_scores(scoring::NScoring, q::Base.CodeUnits{UInt8,String}, r::Base.CodeUnits{UInt8,String}, i::Int, j::Int, previous_score::Int, DP::Vector{Int}, m::Int, INF_INT::Int)
    indel = scoring.indel
    nindel = scoring.nindel
    match = scoring.match
    mismatch = scoring.mismatch

    is_N_q = q[i] == UInt8('N')
    cost = is_N_q ? nindel : indel

    insertion_score = DP[i] + (i == m ? INF_INT : cost)
    deletion_score = previous_score + cost

    is_match = (q[i] == r[j] || is_N_q)
    substitution_score = (i == 1 ? 0 : DP[i-1]) + (is_match ? match : mismatch)

    return insertion_score, deletion_score, substitution_score
end

function semiglobal_alignment_core(
    ws::SemiGlobalWorkspace,
    q, r, m::Int, n::Int,
    max_error::Float64,
    scoring::S,
    output::O,
    ref_search_range::UnitRange{Int},
    max_start_pos::Int,
    min_end_pos::Int,
    normalization_length::Int
) where {S<:AbstractScoring,O<:AbstractOutputPolicy}

    if m == 0 || n == 0
        return finalize_result(output, init_result(output, INF_INT), normalization_length)
    end

    allowed_error = floor(Int, max_error * normalization_length)
    result = init_result(output, INF_INT)

    max_indel_steps = max_indel_steps_for(scoring, allowed_error)

    min_valid_start = min_end_pos - (m + max_indel_steps) + 1

    if min_valid_start > max_start_pos
        return finalize_result(output, result, normalization_length)
    end

    # Shrink ref_search_range start
    if min_valid_start > first(ref_search_range)
        ref_search_range = max(first(ref_search_range), min_valid_start):last(ref_search_range)
    end

    band_offset = max(m - n - max_indel_steps, -max_start_pos - max_indel_steps)

    DP = ws.DP

    # Origin tracking setup
    is_traceback = O <: TracebackOutput
    origin = ws.origin

    @inbounds for i in 1:m # Initialize the DP vector.
        DP[i] = scoring.indel * i
        if is_traceback
            origin[i] = 1 - i
        end
    end

    # Run DP column by column.
    lact = min(allowed_error + 1, m)
    @inbounds for j in ref_search_range
        previous_score_origin = j
        if j + band_offset >= 1
            fact = j + band_offset
            previous_score = allowed_error
        else
            fact = 1
            previous_score = 0
        end

        if fact > lact
            return finalize_result(output, result, normalization_length)
        end

        if fact <= lact
            # 1. First iteration (i = fact)
            insertion_score, deletion_score, substitution_score = step_scores(scoring, q, r, fact, j, previous_score, DP, m, INF_INT)

            # Determine origin for this step
            if is_traceback
                del_origin = previous_score_origin
                sub_origin = (fact == 1 ? j : origin[fact-1])
                ins_origin = origin[fact]
                best_score = deletion_score
                best_origin = del_origin

                if substitution_score < best_score
                    best_score = substitution_score
                    best_origin = sub_origin
                end

                if insertion_score < best_score
                    best_score = insertion_score
                    best_origin = ins_origin
                end

                current_origin = best_origin
            end

            if fact != 1
                DP[fact-1] = previous_score
                if is_traceback
                    origin[fact-1] = previous_score_origin
                end
            end
            previous_score = min(insertion_score, deletion_score, substitution_score)
            if is_traceback
                previous_score_origin = current_origin
            end

            # 2. Main Loop (fact+1 to min(lact, m-1))
            limit = (lact == m) ? m - 1 : lact

            for i in (fact+1):limit
                insertion_score, deletion_score, substitution_score = step_scores_main(scoring, q, r, i, j, previous_score, DP, m, INF_INT)

                if is_traceback
                    del_origin = previous_score_origin
                    sub_origin = origin[i-1]
                    ins_origin = origin[i]

                    best_score = deletion_score
                    best_origin = del_origin

                    if substitution_score < best_score
                        best_score = substitution_score
                        best_origin = sub_origin
                    end

                    if insertion_score < best_score
                        best_score = insertion_score
                        best_origin = ins_origin
                    end

                    current_origin = best_origin
                end

                DP[i-1] = previous_score
                if is_traceback
                    origin[i-1] = previous_score_origin
                end

                previous_score = min(insertion_score, deletion_score, substitution_score)
                if is_traceback
                    previous_score_origin = current_origin
                end
            end

            # 3. Last iteration (i = m) if needed
            if lact == m && lact > fact
                insertion_score, deletion_score, substitution_score = step_scores(scoring, q, r, m, j, previous_score, DP, m, INF_INT)

                if is_traceback
                    del_origin = previous_score_origin
                    ins_origin = origin[m]
                    sub_origin = origin[m-1]

                    best_score = deletion_score
                    best_origin = del_origin

                    if substitution_score < best_score
                        best_score = substitution_score
                        best_origin = sub_origin
                    end

                    if insertion_score < best_score
                        best_score = insertion_score
                        best_origin = ins_origin
                    end

                    current_origin = best_origin
                end

                DP[m-1] = previous_score
                if is_traceback
                    origin[m-1] = previous_score_origin
                end

                previous_score = min(insertion_score, deletion_score, substitution_score)
                if is_traceback
                    previous_score_origin = current_origin
                end
            end
        end

        DP[lact] = previous_score
        if is_traceback
            origin[lact] = previous_score_origin
        end

        if lact == m && previous_score <= allowed_error
            lact -= 1
            if j >= min_end_pos
                if previous_score == 0
                    do_early_exit = !is_traceback || (is_traceback && output.trim_side == 5)

                    if do_early_exit
                        if is_traceback
                            return finalize_result(output, (0.0, previous_score_origin, j), normalization_length)
                        else
                            return finalize_result(output, 0.0, normalization_length)
                        end
                    end
                end

                if is_traceback
                    result = update_result(output, result, previous_score, j, previous_score_origin)
                else
                    result = update_result(output, result, previous_score, j)
                end
            end
        end
        while lact > 0 && DP[lact] > allowed_error
            lact -= 1
        end
        lact += 1
    end
    return finalize_result(output, result, normalization_length)
end

function semiglobal_alignment(ws::SemiGlobalWorkspace, query::String, ref::String, max_error::Float64, match::Int, mismatch::Int, indel::Int, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int, trim_side::Union{Int,Nothing}=nothing, need_traceback::Bool=false)
    m = ncodeunits(query)
    n = ncodeunits(ref)
    q = codeunits(query)
    r = codeunits(ref)
    scoring = SimpleScoring(match, mismatch, indel)

    if isnothing(trim_side) && !need_traceback
        output = ScoreOnly()
    else
        output = TracebackOutput(trim_side)
    end

    return semiglobal_alignment_core(ws, q, r, m, n, max_error, scoring, output, ref_search_range, max_start_pos, min_end_pos, m)
end

function semiglobal_alignment_N(ws::SemiGlobalWorkspace, query::String, ref::String, max_error::Float64, match::Int, mismatch::Int, indel::Int, nindel::Int, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int, non_N_m::Int, trim_side::Union{Int,Nothing}=nothing, need_traceback::Bool=false)
    m = ncodeunits(query)
    n = ncodeunits(ref)
    q = codeunits(query)
    r = codeunits(ref)
    scoring = NScoring(match, mismatch, indel, nindel)

    if isnothing(trim_side) && !need_traceback
        output = ScoreOnly()
    else
        output = TracebackOutput(trim_side)
    end

    return semiglobal_alignment_core(ws, q, r, m, n, max_error, scoring, output, ref_search_range, max_start_pos, min_end_pos, non_N_m)
end

"""
Calculate and compare the similarity of a given sequence seq with the sequences in the given DataFrame bc_df.
# Returns
A tuple `(min_score_bc, min_score, delta, best_start, best_end)`.
"""
function find_best_matching_bc_no_delta(seq::String, bc_seqs::Vector{String}, bc_lengths_no_N::Vector{Int}, ws::SemiGlobalWorkspace, max_error_rate::Float64, match::Int, mismatch::Int, indel::Int, nindel::Union{Int,Nothing}, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int, trim_side::Union{Int,Nothing}, need_traceback::Bool)
    min_score = Inf
    min_score_bc = 0
    best_start = -1
    best_end = -1

    @inbounds for (i, bc) in pairs(bc_seqs)
        if isnothing(nindel)
            alignment_result = semiglobal_alignment(ws, bc, seq, max_error_rate, match, mismatch, indel, ref_search_range, max_start_pos, min_end_pos, trim_side, need_traceback)
        else
            alignment_result = semiglobal_alignment_N(ws, bc, seq, max_error_rate, match, mismatch, indel, nindel, ref_search_range, max_start_pos, min_end_pos, bc_lengths_no_N[i], trim_side, need_traceback)
        end

        if isnothing(trim_side) && !need_traceback
            score = alignment_result
            s, e = -1, -1
        else
            (score, s, e) = alignment_result
        end

        if score <= max_error_rate && score < min_score
            min_score = score
            min_score_bc = i
            max_error_rate = min(max_error_rate, min_score)
            best_start = s
            best_end = e
        end
    end
    return min_score_bc, min_score, Inf, best_start, best_end
end

function find_best_matching_bc_with_delta(seq::String, bc_seqs::Vector{String}, bc_lengths_no_N::Vector{Int}, ws::SemiGlobalWorkspace, max_error_rate::Float64, match::Int, mismatch::Int, indel::Int, nindel::Union{Int,Nothing}, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int, trim_side::Union{Int,Nothing}, need_traceback::Bool)
    min_score = Inf
    sub_min_score = Inf
    min_score_bc = 0
    best_start = -1
    best_end = -1

    @inbounds for (i, bc) in pairs(bc_seqs)
        if isnothing(nindel)
            alignment_result = semiglobal_alignment(ws, bc, seq, max_error_rate, match, mismatch, indel, ref_search_range, max_start_pos, min_end_pos, trim_side, need_traceback)
        else
            alignment_result = semiglobal_alignment_N(ws, bc, seq, max_error_rate, match, mismatch, indel, nindel, ref_search_range, max_start_pos, min_end_pos, bc_lengths_no_N[i], trim_side, need_traceback)
        end

        if isnothing(trim_side) && !need_traceback
            score = alignment_result
            s, e = -1, -1
        else
            (score, s, e) = alignment_result
        end

        if score <= max_error_rate
            if score < min_score
                sub_min_score = min_score
                min_score = score
                min_score_bc = i
                max_error_rate = min(max_error_rate, sub_min_score)
                best_start = s
                best_end = e
            elseif score < sub_min_score
                sub_min_score = score
                max_error_rate = min(max_error_rate, sub_min_score)
            end

        end
    end
    delta = sub_min_score - min_score
    return min_score_bc, min_score, delta, best_start, best_end
end

"""
Calculate and compare the similarity of a given sequence seq with the sequences in the given DataFrame bc_df.
# Returns
A tuple `(min_score_bc, min_score, delta, best_start, best_end)`.
"""


function find_best_matching_bc(seq::String, bc_seqs::Vector{String}, bc_lengths_no_N::Vector{Int}, config::DemuxConfig, ws::SemiGlobalWorkspace, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int, trim_side::Union{Int,Nothing}, need_traceback::Bool=false)
    if config.min_delta == 0.0
        return find_best_matching_bc_no_delta(seq, bc_seqs, bc_lengths_no_N, ws, config.max_error_rate, config.match, config.mismatch, config.indel, config.nindel, ref_search_range, max_start_pos, min_end_pos, trim_side, need_traceback)
    else
        return find_best_matching_bc_with_delta(seq, bc_seqs, bc_lengths_no_N, ws, config.max_error_rate, config.match, config.mismatch, config.indel, config.nindel, ref_search_range, max_start_pos, min_end_pos, trim_side, need_traceback)
    end
end



function determine_filename(seq::String, config::DemuxConfig, ws::SemiGlobalWorkspace)
    n = ncodeunits(seq)

    # --- Pass 1: Barcode 1 ---
    ref_search_range1 = resolve(config.ref_search_range, n)
    barcode_start_range1 = resolve(config.barcode_start_range, n)
    barcode_end_range1 = resolve(config.barcode_end_range, n)

    start_j1 = max(first(ref_search_range1), first(barcode_start_range1), 1)
    end_j1 = min(last(ref_search_range1), last(barcode_end_range1), n)
    max_start_pos1 = last(barcode_start_range1)
    min_end_pos1 = first(barcode_end_range1)

    if start_j1 > end_j1 || start_j1 > max_start_pos1 || end_j1 < min_end_pos1
        return "unknown.fastq", -1, -1
    end

    ref_search_range1 = start_j1:end_j1

    min_score_bc1, score1, delta1, best_start1, best_end1 = find_best_matching_bc(seq, config.bc_seqs, config.bc_lengths_no_N, config, ws, ref_search_range1, max_start_pos1, min_end_pos1, config.trim_side)

    if min_score_bc1 == 0
        return "unknown.fastq", -1, -1
    elseif delta1 < config.min_delta
        return "ambiguous_classification.fastq", -1, -1
    end

    # --- Pass 2: Barcode 2 (if dual) ---
    if config.is_dual

        ref_search_range2 = resolve(config.ref_search_range2, n)
        barcode_start_range2 = resolve(config.barcode_start_range2, n)
        barcode_end_range2 = resolve(config.barcode_end_range2, n)

        start_j2 = max(first(ref_search_range2), first(barcode_start_range2), 1)
        end_j2 = min(last(ref_search_range2), last(barcode_end_range2), n)
        max_start_pos2 = last(barcode_start_range2)
        min_end_pos2 = first(barcode_end_range2)

        if start_j2 > end_j2 || start_j2 > max_start_pos2 || end_j2 < min_end_pos2
            return "unknown.fastq", -1, -1
        end

        ref_search_range2 = start_j2:end_j2

        min_score_bc2, score2, delta2, best_start2, best_end2 = find_best_matching_bc(seq, config.bc_seqs2, config.bc_lengths_no_N2, config, ws, ref_search_range2, max_start_pos2, min_end_pos2, config.trim_side2)

        if min_score_bc2 == 0
            return "unknown.fastq", -1, -1
        elseif delta2 < config.min_delta
            return "ambiguous_classification.fastq", -1, -1
        end

        # Combine IDs
        output_filename = string(config.ids[min_score_bc1]) * "." * string(config.ids2[min_score_bc2]) * ".fastq"
    else
        # Single mode
        output_filename = string(config.ids[min_score_bc1]) * ".fastq"
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



mutable struct DemuxStats
    total_reads::Int
    matched_reads::Int
    unmatched_reads::Int
    ambiguous_reads::Int

    # Key: (bc1_idx, bc2_idx). bc2_idx=0 if single.
    sample_counts::Dict{Tuple{Int,Int},Int}

    bc1_pos_counts::Dict{Int,Int}
    bc1_len_counts::Dict{Int,Int}
    bc1_score_counts::Dict{Float64,Int}
    bc1_per_bc_score_counts::Dict{Int,Dict{Float64,Int}}
    bc1_per_bc_pos_counts::Dict{Int,Dict{Int,Int}}
    bc1_per_bc_len_counts::Dict{Int,Dict{Int,Int}}

    bc2_pos_counts::Dict{Int,Int}
    bc2_len_counts::Dict{Int,Int}
    bc2_score_counts::Dict{Float64,Int}
    bc2_per_bc_score_counts::Dict{Int,Dict{Float64,Int}}
    bc2_per_bc_pos_counts::Dict{Int,Dict{Int,Int}}
    bc2_per_bc_len_counts::Dict{Int,Dict{Int,Int}}
end

function DemuxStats()
    DemuxStats(
        0, 0, 0, 0,
        Dict{Tuple{Int,Int},Int}(),
        Dict{Int,Int}(), Dict{Int,Int}(), Dict{Float64,Int}(), Dict{Int,Dict{Float64,Int}}(), Dict{Int,Dict{Int,Int}}(), Dict{Int,Dict{Int,Int}}(),
        Dict{Int,Int}(), Dict{Int,Int}(), Dict{Float64,Int}(), Dict{Int,Dict{Float64,Int}}(), Dict{Int,Dict{Int,Int}}(), Dict{Int,Dict{Int,Int}}()
    )
end
