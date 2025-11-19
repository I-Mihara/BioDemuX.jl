struct SemiGlobalWorkspace
	DP::Vector{Int}
end
SemiGlobalWorkspace(max_m::Int) = SemiGlobalWorkspace(Vector{Int}(undef, max_m))
const INF_INT = typemax(Int) ÷ 4
"""
This function aligns `query` to `ref`, using semiglobal alignment algorithm. 
# Returns
An alignment score as a float, where lower values indicate better alignment.
"""
function semiglobal_alignment(ws::SemiGlobalWorkspace, query::String, ref::String, max_error::Float64; match::Int = 0, mismatch::Int = 1, indel::Int = 1)
	m = ncodeunits(query)
	n = ncodeunits(ref)
	q = codeunits(query)
	r = codeunits(ref)
	if m == 0 || n == 0
		return Inf
	end
	score = INF_INT
	allowed_error = floor(Int, max_error * m)
	DP = ws.DP
    @inbounds for i in 1:m # Initialize the DP vector.
        DP[i] = indel * i
    end
	# Run DP column by column.
	lact = min(allowed_error + 1, m)
	@inbounds for j in 1:n
		if m - div(allowed_error, indel) - n + j >= 1#
			fact = m - div(allowed_error, indel) - n + j#
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
			if previous_score == 0
				return 0.0
			end
			lact -= 1
			score = min(score, previous_score)
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
function find_best_matching_bc(seq::String, bc_seqs::Vector{String}, ws::SemiGlobalWorkspace, max_error_rate::Float64, min_delta::Float64, mismatch::Int, indel::Int)
	min_score = Inf
	sub_min_score = Inf
	min_score_bc = 0
	if min_delta == 0.0 #to be modified
		for (i, bc) in enumerate(bc_seqs)
			alignment_score = semiglobal_alignment(ws, bc, seq, max_error_rate, mismatch = mismatch, indel = indel)

			if alignment_score <= max_error_rate && alignment_score < min_score
				min_score = alignment_score
				min_score_bc = i
				max_error_rate = min(max_error_rate, min_score)
			end
		end
	else	
		if alignment_score <= max_error_rate
			if alignment_score < min_score
				sub_min_score = min_score
				min_score = alignment_score
				min_score_bc = i
			elseif alignment_score < sub_min_score
				sub_min_score = alignment_score
			end
		end
	end
	delta = sub_min_score - min_score
	return min_score_bc, delta
end

function determine_filename(seq::String, bc_seqs::Vector{String}, ids::Vector{String}, ws::SemiGlobalWorkspace, max_error_rate::Float64, min_delta::Float64, mismatch::Int, indel::Int, gzip_output::Bool)
	min_score_bc, delta = find_best_matching_bc(seq, bc_seqs, ws, max_error_rate, min_delta, mismatch, indel)
	output_filename = ""
	if min_score_bc == 0
		output_filename = "unknown.fastq"
	elseif delta < min_delta
		output_filename = "ambiguous_classification.fastq"
	else
		output_filename = string(ids[min_score_bc]) * ".fastq"
	end
	if gzip_output
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
function classify_sequences(FASTQ_file1::String, FASTQ_file2::String, bc_df::DataFrame, output_dir::String, output_prefix1::String, output_prefix2::String, max_error_rate::Float64, min_delta::Float64, mismatch::Int, indel::Int, classify_both::Bool, gzip_output::Bool)
	bc_seqs = Vector{String}(bc_df.Full_seq)
    ids     = bc_df.ID
    max_m   = maximum(ncodeunits.(bc_seqs))
    ws      = SemiGlobalWorkspace(max_m)
	if classify_both
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
						filename = determine_filename(seq, bc_seqs, ids, ws, max_error_rate, min_delta, mismatch, indel, gzip_output)
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
						filename = determine_filename(line1, bc_seqs, ids, ws, max_error_rate, min_delta, mismatch, indel, gzip_output)
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

function classify_sequences(FASTQ_file1::String, bc_df::DataFrame, output_dir::String, output_prefix::String, max_error_rate::Float64, min_delta::Float64, mismatch::Int, indel::Int, gzip_output::Bool)
	bc_seqs = Vector{String}(bc_df.Full_seq)
    ids     = bc_df.ID
    max_m   = maximum(ncodeunits.(bc_seqs))
    ws      = SemiGlobalWorkspace(max_m)
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
				filename, best_start, best_end = determine_filename_trim(seq, bc_seqs, ids, ws, max_error_rate, min_delta, mismatch, indel, trim_end, gzip_output)
				write_fastq_entry(output_dir * "/" * output_prefix * "." * filename, header, seq, plus, quality_score)
				mode = "header"
			end
		end
	end
	close_all_fastq_ios()
end