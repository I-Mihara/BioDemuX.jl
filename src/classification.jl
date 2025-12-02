struct SemiGlobalWorkspace
	DP::Vector{Int}
end
SemiGlobalWorkspace(max_m::Int) = SemiGlobalWorkspace(Vector{Int}(undef, max_m))
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
	nindel::Union{Int, Nothing} = nothing
	classify_both::Bool = false
	gzip_output::Bool = false
	ref_search_range::DynamicRange = parse_dynamic_range("1:end")
	barcode_start_range::DynamicRange = parse_dynamic_range("1:end")
	barcode_end_range::DynamicRange = parse_dynamic_range("1:end")
	bc_seqs::Vector{String}
	bc_lengths_no_N::Vector{Int}
	ids::Vector{String}
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

@inline function init_result(::ScoreOnly, INF_INT)
    return INF_INT
end

@inline function update_result(::ScoreOnly, current_best, new_score, j)
    return min(current_best, new_score)
end

@inline function finalize_result(::ScoreOnly, result, normalization)
    return result / normalization
end

@inline function max_indel_steps_for(scoring::SimpleScoring, allowed_error::Int)
    return div(allowed_error, scoring.indel)
end

@inline function max_indel_steps_for(scoring::NScoring, allowed_error::Int)
    return div(allowed_error, min(scoring.indel, scoring.nindel))
end

Base.@propagate_inbounds function step_scores_main(scoring::SimpleScoring, q::Base.CodeUnits{UInt8, String}, r::Base.CodeUnits{UInt8, String}, i::Int, j::Int, previous_score::Int, DP::Vector{Int}, m::Int, INF_INT::Int)
    indel = scoring.indel
    match = scoring.match
    mismatch = scoring.mismatch
    
    insertion_score = DP[i] + indel
    deletion_score = previous_score + indel
    substitution_score = DP[i-1] + (q[i] == r[j] ? match : mismatch)
    
    return insertion_score, deletion_score, substitution_score
end

Base.@propagate_inbounds function step_scores_main(scoring::NScoring, q::Base.CodeUnits{UInt8, String}, r::Base.CodeUnits{UInt8, String}, i::Int, j::Int, previous_score::Int, DP::Vector{Int}, m::Int, INF_INT::Int)
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

Base.@propagate_inbounds function step_scores(scoring::SimpleScoring, q::Base.CodeUnits{UInt8, String}, r::Base.CodeUnits{UInt8, String}, i::Int, j::Int, previous_score::Int, DP::Vector{Int}, m::Int, INF_INT::Int)
    indel = scoring.indel
    match = scoring.match
    mismatch = scoring.mismatch
    
    insertion_score = (i == m ? INF_INT : DP[i] + indel)
    deletion_score = previous_score + indel
    substitution_score = (i == 1 ? 0 : DP[i-1]) + (q[i] == r[j] ? match : mismatch)
    
    return insertion_score, deletion_score, substitution_score
end

Base.@propagate_inbounds function step_scores(scoring::NScoring, q::Base.CodeUnits{UInt8, String}, r::Base.CodeUnits{UInt8, String}, i::Int, j::Int, previous_score::Int, DP::Vector{Int}, m::Int, INF_INT::Int)
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
) where {S<:AbstractScoring, O<:AbstractOutputPolicy}
    
    if m == 0 || n == 0
        return Inf
    end
    
    allowed_error = floor(Int, max_error * normalization_length)
    result = init_result(output, INF_INT)
    
    max_indel_steps = max_indel_steps_for(scoring, allowed_error)
    band_offset = max(m - n - max_indel_steps, - max_start_pos - max_indel_steps)
    
    DP = ws.DP
    @inbounds for i in 1:m # Initialize the DP vector.
        DP[i] = scoring.indel * i
    end
    # Run DP column by column.
    lact = min(allowed_error + 1, m)
    @inbounds for j in ref_search_range
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
            if fact != 1
                DP[fact-1] = previous_score
            end
            previous_score = min(insertion_score, deletion_score, substitution_score)
            
            # 2. Main Loop (fact+1 to min(lact, m-1))
            limit = (lact == m) ? m - 1 : lact
            
            for i in (fact + 1):limit
                insertion_score, deletion_score, substitution_score = step_scores_main(scoring, q, r, i, j, previous_score, DP, m, INF_INT)
                DP[i-1] = previous_score
                previous_score = min(insertion_score, deletion_score, substitution_score)
            end
            
            # 3. Last iteration (i = m) if needed
            if lact == m && lact > fact
                insertion_score, deletion_score, substitution_score = step_scores(scoring, q, r, m, j, previous_score, DP, m, INF_INT)
                DP[m-1] = previous_score
                previous_score = min(insertion_score, deletion_score, substitution_score)
            end
        end
        
        DP[lact] = previous_score
        if lact == m && previous_score <= allowed_error
            lact -= 1
            if j >= min_end_pos
                if previous_score == 0
                    return finalize_result(output, 0.0, normalization_length)
                end
                result = update_result(output, result, previous_score, j)
            end
        end
        while lact > 0 && DP[lact] > allowed_error
            lact -= 1
        end
        lact += 1
    end
    return finalize_result(output, result, normalization_length)
end

"""
This function aligns `query` to `ref`, using semiglobal alignment algorithm. 
# Returns
An alignment score as a float, where lower values indicate better alignment.
"""
function semiglobal_alignment(ws::SemiGlobalWorkspace, query::String, ref::String, max_error::Float64, match::Int, mismatch::Int, indel::Int, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int)
    m = ncodeunits(query)
    n = ncodeunits(ref)
    q = codeunits(query)
    r = codeunits(ref)
    scoring = SimpleScoring(match, mismatch, indel)
    output = ScoreOnly()
    return semiglobal_alignment_core(ws, q, r, m, n, max_error, scoring, output, ref_search_range, max_start_pos, min_end_pos, m)
end


function semiglobal_alignment_N(ws::SemiGlobalWorkspace, query::String, ref::String, max_error::Float64, match::Int, mismatch::Int, indel::Int, nindel::Int, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int, non_N_m::Int)
    m = ncodeunits(query)
    n = ncodeunits(ref)
    q = codeunits(query)
    r = codeunits(ref)
    scoring = NScoring(match, mismatch, indel, nindel)
    output = ScoreOnly()
    return semiglobal_alignment_core(ws, q, r, m, n, max_error, scoring, output, ref_search_range, max_start_pos, min_end_pos, non_N_m)
end

"""
Calculate and compare the similarity of a given sequence seq with the sequences in the given DataFrame bc_df.
# Returns
A tuple `(min_score_bc, delta)`, where `min_score_bc` is the index of the best matching sequence in `bc_df`, and `delta` is the difference between the lowest and second-lowest scores.
"""
function find_best_matching_bc_no_delta(seq::String, bc_seqs::Vector{String}, bc_lengths_no_N::Vector{Int}, ws::SemiGlobalWorkspace, max_error_rate::Float64, match::Int, mismatch::Int, indel::Int, nindel::Union{Int, Nothing}, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int)
	min_score = Inf
	min_score_bc = 0
	
	@inbounds for (i, bc) in pairs(bc_seqs)
		if isnothing(nindel)
			alignment_score = semiglobal_alignment(ws, bc, seq, max_error_rate, match, mismatch, indel, ref_search_range, max_start_pos, min_end_pos)
		else
			alignment_score = semiglobal_alignment_N(ws, bc, seq, max_error_rate, match, mismatch, indel, nindel, ref_search_range, max_start_pos, min_end_pos, bc_lengths_no_N[i])
		end

		if alignment_score <= max_error_rate && alignment_score < min_score
			min_score = alignment_score
			min_score_bc = i
			max_error_rate = min(max_error_rate, min_score)
		end
	end
	return min_score_bc, Inf
end

function find_best_matching_bc_with_delta(seq::String, bc_seqs::Vector{String}, bc_lengths_no_N::Vector{Int}, ws::SemiGlobalWorkspace, max_error_rate::Float64, match::Int, mismatch::Int, indel::Int, nindel::Union{Int, Nothing}, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int)
	min_score = Inf
	sub_min_score = Inf
	min_score_bc = 0
	
	@inbounds for (i, bc) in pairs(bc_seqs)
		if isnothing(nindel)
			alignment_score = semiglobal_alignment(ws, bc, seq, max_error_rate, match, mismatch, indel, ref_search_range, max_start_pos, min_end_pos)
		else
			alignment_score = semiglobal_alignment_N(ws, bc, seq, max_error_rate, match, mismatch, indel, nindel, ref_search_range, max_start_pos, min_end_pos, bc_lengths_no_N[i])
		end
		if alignment_score <= max_error_rate
			if alignment_score < min_score
				sub_min_score = min_score
				min_score = alignment_score
				min_score_bc = i
				max_error_rate = min(max_error_rate, sub_min_score)
			elseif alignment_score < sub_min_score
				sub_min_score = alignment_score
				max_error_rate = min(max_error_rate, sub_min_score)
			end
			
		end
	end
	delta = sub_min_score - min_score
	return min_score_bc, delta
end

"""
Calculate and compare the similarity of a given sequence seq with the sequences in the given DataFrame bc_df.
# Returns
A tuple `(min_score_bc, delta)`, where `min_score_bc` is the index of the best matching sequence in `bc_df`, and `delta` is the difference between the lowest and second-lowest scores.
"""


function find_best_matching_bc(seq::String, config::DemuxConfig, ws::SemiGlobalWorkspace, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int)
	if config.min_delta == 0.0
		return find_best_matching_bc_no_delta(seq, config.bc_seqs, config.bc_lengths_no_N, ws, config.max_error_rate, config.match, config.mismatch, config.indel, config.nindel, ref_search_range, max_start_pos, min_end_pos)
	else
		return find_best_matching_bc_with_delta(seq, config.bc_seqs, config.bc_lengths_no_N, ws, config.max_error_rate, config.match, config.mismatch, config.indel, config.nindel, ref_search_range, max_start_pos, min_end_pos)
	end
end



function determine_filename(seq::String, config::DemuxConfig, ws::SemiGlobalWorkspace)
	n = ncodeunits(seq)
	ref_search_range = resolve(config.ref_search_range, n)
	barcode_start_range = resolve(config.barcode_start_range, n)
	barcode_end_range = resolve(config.barcode_end_range, n)
	
	start_j = max(first(ref_search_range), first(barcode_start_range), 1)
	end_j   = min(last(ref_search_range),  last(barcode_end_range), n)
	max_start_pos = last(barcode_start_range)
	min_end_pos = first(barcode_end_range)
	
	if start_j > end_j || start_j > max_start_pos || end_j < min_end_pos
		# Invalid range, treat as unknown.
		return "unknown.fastq"
	end
	
	ref_search_range = start_j:end_j
	
	min_score_bc, delta = find_best_matching_bc(seq, config, ws, ref_search_range, max_start_pos, min_end_pos)
	output_filename = ""
	if min_score_bc == 0
		output_filename = "unknown.fastq"
	elseif delta < config.min_delta
		output_filename = "ambiguous_classification.fastq"
	else
		output_filename = string(config.ids[min_score_bc]) * ".fastq"
	end
	if config.gzip_output
		return output_filename * ".gz"
	else
		return output_filename
	end
end
