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
	classify_both::Bool = false
	gzip_output::Bool = false
	ref_search_range::DynamicRange = parse_dynamic_range("1:end")
	barcode_start_range::DynamicRange = parse_dynamic_range("1:end")
	barcode_end_range::DynamicRange = parse_dynamic_range("1:end")
	bc_seqs::Vector{String}
	ids::Vector{String}
	ws::SemiGlobalWorkspace
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
function semiglobal_alignment(ws::SemiGlobalWorkspace, query::String, ref::String, max_error::Float64, match::Int, mismatch::Int, indel::Int, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int)
	m = ncodeunits(query)
	n = ncodeunits(ref)
	q = codeunits(query)
	r = codeunits(ref)
	
	if m == 0 || n == 0
		return Inf
	end
	
	allowed_error = floor(Int, max_error * m)
	score = INF_INT
	
	max_indel_steps = div(allowed_error, indel)
	band_offset = max(m - n - max_indel_steps, - max_start_pos - max_indel_steps)
	
	DP = ws.DP
    @inbounds for i in 1:m # Initialize the DP vector.
        DP[i] = indel * i
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
			return score / m
		end
		@inbounds for i in fact:lact
			insertion_score = (i == m ? INF_INT : DP[i]+indel)#→
			deletion_score = previous_score + indel#↓
			substitution_score = (i == 1 ? 0 : DP[i-1]) + (q[i] == r[j] ? match : mismatch)#↘︎
			if i != 1
				DP[i-1] = previous_score
			end
			previous_score = min(insertion_score, deletion_score, substitution_score)
		end
		DP[lact] = previous_score
		if lact == m && previous_score <= allowed_error
			lact -= 1
			if j >= min_end_pos
				if previous_score == 0
					return 0.0
				end
				score = min(score, previous_score)
			end
		end
		while lact > 0 && DP[lact] > allowed_error
			lact -= 1
		end
		lact += 1
	end
	return score / m
end

"""
Calculate and compare the similarity of a given sequence seq with the sequences in the given DataFrame bc_df.
# Returns
A tuple `(min_score_bc, delta)`, where `min_score_bc` is the index of the best matching sequence in `bc_df`, and `delta` is the difference between the lowest and second-lowest scores.
"""
function find_best_matching_bc_no_delta(seq::String, bc_seqs::Vector{String}, ws::SemiGlobalWorkspace, max_error_rate::Float64, match::Int, mismatch::Int, indel::Int, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int)
	min_score = Inf
	min_score_bc = 0
	
	@inbounds for (i, bc) in pairs(bc_seqs)
		alignment_score = semiglobal_alignment(ws, bc, seq, max_error_rate, match, mismatch, indel, ref_search_range, max_start_pos, min_end_pos)

		if alignment_score <= max_error_rate && alignment_score < min_score
			min_score = alignment_score
			min_score_bc = i
			max_error_rate = min(max_error_rate, min_score)
		end
	end
	return min_score_bc, Inf
end

function find_best_matching_bc_with_delta(seq::String, bc_seqs::Vector{String}, ws::SemiGlobalWorkspace, max_error_rate::Float64, match::Int, mismatch::Int, indel::Int, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int)
	min_score = Inf
	sub_min_score = Inf
	min_score_bc = 0
	
	@inbounds for (i, bc) in pairs(bc_seqs)
		alignment_score = semiglobal_alignment(ws, bc, seq, max_error_rate, match, mismatch, indel, ref_search_range, max_start_pos, min_end_pos)	
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


function find_best_matching_bc(seq::String, config::DemuxConfig, ref_search_range::UnitRange{Int}, max_start_pos::Int, min_end_pos::Int)
	if config.min_delta == 0.0
		return find_best_matching_bc_no_delta(seq, config.bc_seqs, config.ws, config.max_error_rate, config.match, config.mismatch, config.indel, ref_search_range, max_start_pos, min_end_pos)
	else
		return find_best_matching_bc_with_delta(seq, config.bc_seqs, config.ws, config.max_error_rate, config.match, config.mismatch, config.indel, ref_search_range, max_start_pos, min_end_pos)
	end
end



function determine_filename(seq::String, config::DemuxConfig)
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
	
	min_score_bc, delta = find_best_matching_bc(seq, config, ref_search_range, max_start_pos, min_end_pos)
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

const _fastq_ios = Dict{String, IO}()

function get_fastq_io(path::String)
    get!(_fastq_ios, path) do
        raw = open(path, "a")
        io  = endswith(lowercase(path), ".gz") ?
              GzipCompressorStream(raw) : raw
        return BufferedOutputStream(io)
    end
end

function write_fastq_entry(filepath, header, seq, plus, quality)
	io=get_fastq_io(filepath)
	write(io, header * "\n" * seq * "\n" * plus * "\n" * quality * "\n")
end

function close_all_fastq_ios()
	for io in values(_fastq_ios)
		flush(io)
		close(io)
	end
	empty!(_fastq_ios)
end


"""
Compare each sequence in the FASTQ_file1 file with the sequences in bc_df, and classify the sequences of the specified file based on that comparison.
"""
function classify_sequences(FASTQ_file1::String, FASTQ_file2::String, output_dir::String, output_prefix1::String, output_prefix2::String, config::DemuxConfig)
	if config.classify_both
		read_fastq(FASTQ_file1) do primary_file
			read_fastq(FASTQ_file2) do secondary_file
				header, seq, plus, quality_score = "", "", "", ""
				header2, seq2, plus2, quality_score2 = "", "", "", ""
				mode = "header"
				filename = ""
				for line1 in eachline(primary_file)
					line2 = readline(secondary_file)
					if line1[1] == '@' && mode == "header"
						header = line1
						header2 = line2
						mode = "seq"
					elseif mode == "seq"
						seq = line1
						seq2 = line2
						mode = "plus"
					elseif mode == "plus"
						plus = line1
						plus2 = line2
						mode = "quality_score"
					elseif mode == "quality_score"
						quality_score = line1
						quality_score2 = line2
						filename = determine_filename(seq, config)
						write_fastq_entry(output_dir * "/" * output_prefix1 * "." * filename, header, seq, plus, quality_score)
						write_fastq_entry(output_dir * "/" * output_prefix2 * "." * filename, header2, seq2, plus2, quality_score2)
						mode = "header"
					end
				end
			end
		end
	else
		read_fastq(FASTQ_file1) do primary_file
			read_fastq(FASTQ_file2) do secondary_file
				header2, seq2, plus2, quality_score2 = "", "", "", ""
				mode = "header"
				filename = ""
				for line1 in eachline(primary_file)
					line2 = readline(secondary_file)
					if line1[1] == '@' && mode == "header"
						header2 = line2
						mode = "seq"
					elseif mode == "seq"
						filename = determine_filename(line1, config)
						seq2 = line2
						mode = "plus"
					elseif mode == "plus"
						plus2 = line2
						mode = "quality_score"
					elseif mode == "quality_score"
						quality_score2 = line2
						write_fastq_entry(output_dir * "/" * output_prefix2 * "." * filename, header2, seq2, plus2, quality_score2)
						mode = "header"
					end
				end
			end
		end
	end
	close_all_fastq_ios()
end

function classify_sequences(FASTQ_file1::String, output_dir::String, output_prefix::String, config::DemuxConfig)
	read_fastq(FASTQ_file1) do file
		header, seq, plus, quality_score = "", "", "", ""
		mode = "header"
		filename = ""
		for line in eachline(file)
			if line[1] == '@' && mode == "header"
				header = line
				mode = "seq"
			elseif mode == "seq"
				seq = line
				mode = "plus"
			elseif mode == "plus"
				plus = line
				mode = "quality_score"
			elseif mode == "quality_score"
				quality_score = line
				filename = determine_filename(seq, config)
				write_fastq_entry(output_dir * "/" * output_prefix * "." * filename, header, seq, plus, quality_score)
				mode = "header"
			end
		end
	end
	close_all_fastq_ios()
end