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



    ref_search_range::DynamicRange = parse_dynamic_range("1:end")
    barcode_start_range::DynamicRange = parse_dynamic_range("1:end")
    barcode_end_range::DynamicRange = parse_dynamic_range("1:end")
    bc_seqs::Vector{String}
    bc_lengths_no_N::Vector{Int}
    ids::Vector{String}



    is_dual::Bool = false
    ref_search_range2::DynamicRange = parse_dynamic_range("1:end")
    barcode_start_range2::DynamicRange = parse_dynamic_range("1:end")
    barcode_end_range2::DynamicRange = parse_dynamic_range("1:end")
    bc_seqs2::Vector{String} = String[]
    bc_lengths_no_N2::Vector{Int} = Int[]
    ids2::Vector{String} = String[]



    trim_side::Union{Int,Nothing} = nothing
    trim_side2::Union{Int,Nothing} = nothing



    summary::Bool = false
    summary_format::Symbol = :txt # :txt, :json, :html



    matching_algorithm::Symbol = :semiglobal # :semiglobal, :hamming
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
Aligns `query` to `ref` using exact string matching (byte-by-byte comparison).
No indels or mismatches allowed.
We scan `query` across `ref_search_range` in `ref`.
Returns `(0.0, start_pos, end_pos)` if found, otherwise `(Inf, -1, -1)`.
"""
function exact_align(query::String, ref::String, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int, trim_side::Union{Int,Nothing})
    m = ncodeunits(query)
    n = ncodeunits(ref)

    # Reference range constraint
    start_range_first = max(first(ref_search_range), 1)
    start_range_last = min(last(ref_search_range), max_start_pos, n - m + 1)

    if start_range_last < start_range_first
        return (Inf, -1, -1)
    end

    # If trim_side == 3, we prefer rightmost -> findprev
    # Search backward from the LAST possible character of the match (start_range_last + m - 1)
    if !isnothing(trim_side) && trim_side == 3
        # findprev(query, ref, last_idx) finds match starting at or before last_idx
        # We need the match to start at or before start_range_last.
        # However, findprev returns range of the match.
        # If we use `start_range_last + m - 1` as the limit, findprev guarantees the match ends at or before that limit (if searching strict match).

        range = findprev(query, ref, start_range_last + m - 1)
        if !isnothing(range)
            s = first(range)
            if s >= start_range_first
                # Check min_end_pos constraint
                if (s + m - 1) >= min_end_pos
                    return (0.0, s, s + m - 1)
                end
            end
        end
        return (Inf, -1, -1)

    else
        # Prefer leftmost (default) -> findnext
        range = findnext(query, ref, start_range_first)
        if !isnothing(range)
            s = first(range)
            if s <= start_range_last
                # Check min_end_pos constraint
                if (s + m - 1) >= min_end_pos
                    return (0.0, s, s + m - 1)
                else
                    # Match found but ends too early. Continue scan.
                    next_search = s + 1
                    while next_search <= start_range_last
                        range = findnext(query, ref, next_search)
                        if isnothing(range)
                            break
                        end
                        s = first(range)
                        if s > start_range_last
                            break
                        end
                        if (s + m - 1) >= min_end_pos
                            return (0.0, s, s + m - 1)
                        end
                        next_search = s + 1
                    end
                end
            end
        end
        return (Inf, -1, -1)
    end
end

"""
Aligns `query` to `ref` using Hamming distance.
Indels are NOT allowed (distance = infinity if lengths differ in an alignment context, but here we scan).
We scan `query` across `ref_search_range` in `ref`.
`N` in `query` matches anything. `N` in `ref` matches only `N` in `query`.
Returns `(score, start_pos, end_pos)`.
"""
function hamming_align(query::String, ref::String, max_error_rate::Float64, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int, trim_side::Union{Int,Nothing})
    m = ncodeunits(query)
    n = ncodeunits(ref)

    # Pre-calculate best score tracking
    best_score = Inf
    best_start = -1
    best_end = -1

    # allowed max errors count (floor)
    allowed_errors = floor(Int, max_error_rate * m)

    # Reference range constraint
    start_range_first = max(first(ref_search_range), 1)
    start_range_last = min(last(ref_search_range), max_start_pos, n - m + 1)

    if start_range_last < start_range_first
        # No valid start position
        return (Inf, -1, -1)
    end

    q_bytes = codeunits(query)
    r_bytes = codeunits(ref)

    @inbounds for j in start_range_first:start_range_last
        # Check end pos constraint
        end_pos = j + m - 1
        if end_pos < min_end_pos
            continue
        end

        current_mismatches = 0
        match_failed = false

        # Inner loop - scan the barcode
        for k in 0:(m-1)
            q_char = q_bytes[k+1]
            r_char = r_bytes[j+k]

            # Match if equal OR query is 'N' (wildcard)
            if q_char != r_char && q_char != 0x4E # 'N'
                current_mismatches += 1
                if current_mismatches > allowed_errors
                    match_failed = true
                    break
                end
            end
        end

        if !match_failed
            score = current_mismatches / m

            if score < best_score
                best_score = score
                best_start = j
                best_end = end_pos
            elseif score == best_score
                if !isnothing(trim_side) && trim_side == 3
                    if j > best_start
                        best_start = j
                        best_end = end_pos
                    end
                end
            end
        end
    end

    return (best_score, best_start, best_end)
end

"""
Calculate and compare the similarity of a given sequence seq with the sequences in the given DataFrame bc_df.
# Returns
A tuple `(min_score_bc, min_score, delta, best_start, best_end)`.
"""
function find_best_matching_bc_no_delta(seq::String, bc_seqs::Vector{String}, bc_lengths_no_N::Vector{Int}, ws::SemiGlobalWorkspace, matching_algorithm::Symbol, max_error_rate::Float64, match::Int, mismatch::Int, indel::Int, nindel::Union{Int,Nothing}, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int, trim_side::Union{Int,Nothing}, need_traceback::Bool)
    min_score = Inf
    min_score_bc = 0
    best_start = -1
    best_end = -1

    @inbounds for (i, bc) in pairs(bc_seqs)
        if matching_algorithm == :hamming
            alignment_result = hamming_align(bc, seq, max_error_rate, ref_search_range, max_start_pos, min_end_pos, trim_side)
        elseif matching_algorithm == :exact
            alignment_result = exact_align(bc, seq, ref_search_range, max_start_pos, min_end_pos, trim_side)
        else
            if isnothing(nindel)
                alignment_result = semiglobal_alignment(ws, bc, seq, max_error_rate, match, mismatch, indel, ref_search_range, max_start_pos, min_end_pos, trim_side, need_traceback)
            else
                alignment_result = semiglobal_alignment_N(ws, bc, seq, max_error_rate, match, mismatch, indel, nindel, ref_search_range, max_start_pos, min_end_pos, bc_lengths_no_N[i], trim_side, need_traceback)
            end
        end

        if isnothing(trim_side) && !need_traceback && matching_algorithm != :hamming && matching_algorithm != :exact
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

function find_best_matching_bc_with_delta(seq::String, bc_seqs::Vector{String}, bc_lengths_no_N::Vector{Int}, ws::SemiGlobalWorkspace, matching_algorithm::Symbol, max_error_rate::Float64, match::Int, mismatch::Int, indel::Int, nindel::Union{Int,Nothing}, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int, trim_side::Union{Int,Nothing}, need_traceback::Bool)
    min_score = Inf
    sub_min_score = Inf
    min_score_bc = 0
    best_start = -1
    best_end = -1

    @inbounds for (i, bc) in pairs(bc_seqs)
        if matching_algorithm == :hamming
            (score, s, e) = hamming_align(bc, seq, max_error_rate, ref_search_range, max_start_pos, min_end_pos, trim_side)
        elseif matching_algorithm == :exact
            (score, s, e) = exact_align(bc, seq, ref_search_range, max_start_pos, min_end_pos, trim_side)
        else
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
        return find_best_matching_bc_no_delta(seq, bc_seqs, bc_lengths_no_N, ws, config.matching_algorithm, config.max_error_rate, config.match, config.mismatch, config.indel, config.nindel, ref_search_range, max_start_pos, min_end_pos, trim_side, need_traceback)
    else
        return find_best_matching_bc_with_delta(seq, bc_seqs, bc_lengths_no_N, ws, config.matching_algorithm, config.max_error_rate, config.match, config.mismatch, config.indel, config.nindel, ref_search_range, max_start_pos, min_end_pos, trim_side, need_traceback)
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

"""
    match_barcode_pass(seq::String, n::Int, config::DemuxConfig, ws::SemiGlobalWorkspace, is_pass2::Bool, stats::Union{DemuxStats,Nothing}=nothing)

Helper function to run a single barcode matching pass.
Returns `(status, barcode_index, start_pos, end_pos)`.
Status can be: `:match`, `:unknown`, `:ambiguous`.
"""
function match_barcode_pass(seq::String, n::Int, config::DemuxConfig, ws::SemiGlobalWorkspace, is_pass2::Bool, stats::Union{DemuxStats,Nothing}=nothing)
    # Select ranges and params based on pass
    if !is_pass2
        ref_sr = config.ref_search_range
        bc_start_sr = config.barcode_start_range
        bc_end_sr = config.barcode_end_range
        bc_seqs = config.bc_seqs
        bc_lens = config.bc_lengths_no_N
        trim_side = config.trim_side
    else
        ref_sr = config.ref_search_range2
        bc_start_sr = config.barcode_start_range2
        bc_end_sr = config.barcode_end_range2
        bc_seqs = config.bc_seqs2
        bc_lens = config.bc_lengths_no_N2
        trim_side = config.trim_side2
    end

    # Resolve dynamic ranges
    r_search = resolve(ref_sr, n)
    b_start = resolve(bc_start_sr, n)
    b_end = resolve(bc_end_sr, n)

    start_j = max(first(r_search), first(b_start), 1)
    end_j = min(last(r_search), last(b_end), n)
    max_start_pos = last(b_start)
    min_end_pos = first(b_end)

    # Sanity check ranges
    if start_j > end_j || start_j > max_start_pos || end_j < min_end_pos
        return :unknown, 0, -1, -1
    end

    final_search_range = start_j:end_j

    # Find best match
    need_tb = !isnothing(trim_side) || !isnothing(stats)

    min_bc, score, delta, s, e = find_best_matching_bc(
        seq, bc_seqs, bc_lens, config, ws,
        final_search_range, max_start_pos, min_end_pos,
        trim_side, need_tb
    )

    if min_bc == 0
        return :unknown, 0, -1, -1
    elseif delta < config.min_delta
        return :ambiguous, 0, -1, -1
    end

    # Update stats if provided
    if !isnothing(stats)
        if !is_pass2
            # BC1 stats
            stats.bc1_pos_counts[s] = get(stats.bc1_pos_counts, s, 0) + 1

            len = e - s + 1
            stats.bc1_len_counts[len] = get(stats.bc1_len_counts, len, 0) + 1

            r_score = round(score, digits=2)
            stats.bc1_score_counts[r_score] = get(stats.bc1_score_counts, r_score, 0) + 1

            d_s = get!(() -> Dict{Float64,Int}(), stats.bc1_per_bc_score_counts, min_bc)
            d_s[r_score] = get(d_s, r_score, 0) + 1

            d_p = get!(() -> Dict{Int,Int}(), stats.bc1_per_bc_pos_counts, min_bc)
            d_p[s] = get(d_p, s, 0) + 1

            d_l = get!(() -> Dict{Int,Int}(), stats.bc1_per_bc_len_counts, min_bc)
            d_l[len] = get(d_l, len, 0) + 1
        else
            # BC2 stats
            stats.bc2_pos_counts[s] = get(stats.bc2_pos_counts, s, 0) + 1

            len = e - s + 1
            stats.bc2_len_counts[len] = get(stats.bc2_len_counts, len, 0) + 1

            r_score = round(score, digits=2)
            stats.bc2_score_counts[r_score] = get(stats.bc2_score_counts, r_score, 0) + 1

            d_s = get!(() -> Dict{Float64,Int}(), stats.bc2_per_bc_score_counts, min_bc)
            d_s[r_score] = get(d_s, r_score, 0) + 1

            d_p = get!(() -> Dict{Int,Int}(), stats.bc2_per_bc_pos_counts, min_bc)
            d_p[s] = get(d_p, s, 0) + 1

            d_l = get!(() -> Dict{Int,Int}(), stats.bc2_per_bc_len_counts, min_bc)
            d_l[len] = get(d_l, len, 0) + 1
        end
    end

    return :match, min_bc, s, e
end


function determine_filename(seq::String, config::DemuxConfig, ws::SemiGlobalWorkspace)
    n = ncodeunits(seq)

    # --- Pass 1: Barcode 1 ---
    status1, bc1_idx, start1, end1 = match_barcode_pass(seq, n, config, ws, false, nothing)

    suffix = config.gzip_output ? ".fastq.gz" : ".fastq"

    if status1 == :unknown
        return "unknown" * suffix, -1, -1
    elseif status1 == :ambiguous
        return "ambiguous_classification" * suffix, -1, -1
    end

    # --- Pass 2: Barcode 2 (if dual) ---
    bc2_idx = 0
    if config.is_dual
        status2, idx, start2, end2 = match_barcode_pass(seq, n, config, ws, true, nothing)

        if status2 == :unknown
            return "unknown" * suffix, -1, -1
        elseif status2 == :ambiguous
            return "ambiguous_classification" * suffix, -1, -1
        end
        bc2_idx = idx

        output_filename = string(config.ids[bc1_idx]) * "." * string(config.ids2[bc2_idx]) * suffix
    else
        output_filename = string(config.ids[bc1_idx]) * suffix
    end

    # Calculate trim range for BC1
    trim_start = -1
    trim_end = -1

    # We want the KEEP range.
    keep_start = 1
    keep_end = n

    if !isnothing(config.trim_side)
        if config.trim_side == 3
            # Keep 1 to start-1
            # But safe_start = max(1, start1)
            keep_end = max(1, start1) - 1
        elseif config.trim_side == 5
            # Keep end1+1 to n
            keep_start = end1 + 1
        end
    end

    if config.is_dual && !isnothing(config.trim_side2)
        if config.trim_side2 == 3
            # Trim everything after start2
            keep_end = min(keep_end, max(1, start2) - 1)
        elseif config.trim_side2 == 5
            # Trim everything before end2
            keep_start = max(keep_start, end2 + 1)
        end
    end

    # Check validity
    if keep_start > keep_end
        # Should result in empty string (represented as 1:0 for internal logic if needed, but here returns 1, 0)
        return output_filename, 1, 0
    end

    return output_filename, keep_start, keep_end
end

function determine_filename_and_stats(seq::String, config::DemuxConfig, ws::SemiGlobalWorkspace, stats::DemuxStats)
    n = ncodeunits(seq)
    stats.total_reads += 1

    # --- Pass 1: Barcode 1 ---
    status1, bc1_idx, start1, end1 = match_barcode_pass(seq, n, config, ws, false, stats)

    suffix = config.gzip_output ? ".fastq.gz" : ".fastq"

    if status1 == :unknown
        stats.unmatched_reads += 1
        return "unknown" * suffix, -1, -1
    elseif status1 == :ambiguous
        stats.ambiguous_reads += 1
        return "ambiguous_classification" * suffix, -1, -1
    end

    # --- Pass 2: Barcode 2 (if dual) ---
    bc2_idx = 0
    if config.is_dual
        status2, idx, start2, end2 = match_barcode_pass(seq, n, config, ws, true, stats)

        if status2 == :unknown
            stats.unmatched_reads += 1
            return "unknown" * suffix, -1, -1
        elseif status2 == :ambiguous
            stats.ambiguous_reads += 1
            return "ambiguous_classification" * suffix, -1, -1
        end
        bc2_idx = idx
        output_filename = string(config.ids[bc1_idx]) * "." * string(config.ids2[bc2_idx]) * suffix
    else
        output_filename = string(config.ids[bc1_idx]) * suffix
    end

    # Update Match Counts
    stats.matched_reads += 1
    key = (bc1_idx, bc2_idx)
    stats.sample_counts[key] = get(stats.sample_counts, key, 0) + 1

    # Calculate trim range (Same deduplicated logic)
    keep_start = 1
    keep_end = n

    if !isnothing(config.trim_side)
        if config.trim_side == 3
            keep_end = max(1, start1) - 1
        elseif config.trim_side == 5
            keep_start = end1 + 1
        end
    end

    if config.is_dual && !isnothing(config.trim_side2)
        if config.trim_side2 == 3
            keep_end = min(keep_end, max(1, start2) - 1)
        elseif config.trim_side2 == 5
            keep_start = max(keep_start, end2 + 1)
        end
    end

    if keep_start > keep_end
        return output_filename, 1, 0
    end

    return output_filename, keep_start, keep_end
end
