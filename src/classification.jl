"""
This function aligns `query` to `ref`, using semiglobal alignment algorithm. 
# Returns
An alignment score as a float, where lower values indicate better alignment.
"""
function semiglobal_alignment(query::String, ref::String, max_error::Float64; match::Int = 0, mismatch::Int = 1, indel::Int = 1)
	m = Base.length(query)
	n = Base.length(ref)
	if m == 0 || n == 0
		return Inf
	end
	score = Inf
	allowed_error = floor(max_error * m) |> Int
	DP = [indel * i for i in 1:m] # Initialize the DP vector.
	# Run DP column by column.
	lact = min(allowed_error + 1, m)
	for j in 1:n
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
		for i in fact:lact
			insertion_score = (i == m ? Inf : DP[i]+indel)#→
			deletion_score = previous_score + indel#↓
			substitution_score = (i == 1 ? 0 : DP[i-1]) + (query[i] == ref[j] ? match : mismatch)#↘︎
			if i != 1
				DP[i-1] = previous_score
			end
			previous_score = min(insertion_score, deletion_score, substitution_score)
		end
		DP[lact] = previous_score
		while lact > 0 && DP[lact] > allowed_error
			lact -= 1
		end
		if lact == m
			score=min(score,previous_score)
		else
			lact += 1
		end
	end
	return score / m
end

"""
Calculate and compare the similarity of a given sequence seq with the sequences in the given DataFrame bc_df.
# Returns
A tuple `(min_score_bc, delta)`, where `min_score_bc` is the index of the best matching sequence in `bc_df`, and `delta` is the difference between the lowest and second-lowest scores.
"""
function find_best_matching_bc(seq::String, bc_df::DataFrame, max_error_rate::Float64, mismatch::Int, indel::Int)
	min_score = Inf
	sub_min_score = Inf
	min_score_bc = 0

	for (i, row) in enumerate(eachrow(bc_df))
		alignment_score = semiglobal_alignment(row.Full_seq, seq, max_error_rate, mismatch = mismatch, indel = indel)

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

function determine_filename(seq::String, bc_df::DataFrame, max_error_rate::Float64, min_delta::Float64, mismatch::Int, indel::Int)
	min_score_bc, delta = find_best_matching_bc(seq, bc_df, max_error_rate, mismatch, indel)

	if min_score_bc == 0
		return "unknown.fastq"
	elseif delta < min_delta
		return "ambiguous_classification.fastq"
	else
		return string(bc_df.ID[min_score_bc]) * ".fastq"
	end
end

function write_fastq_entry(filepath, header, seq, plus, quality)
	open(filepath, "a") do outputfile
		write(outputfile, header * "\n" * seq * "\n" * plus * "\n" * quality * "\n")
	end
end

"""
Compare each sequence in the FASTQ_file1 file with the sequences in bc_df, and classify the sequences of the specified file based on that comparison.
"""
function classify_sequences(FASTQ_file1::String, FASTQ_file2::String, bc_df::DataFrame, output_dir::String, output_prefix1::String, output_prefix2::String, max_error_rate::Float64, min_delta::Float64, mismatch::Int, indel::Int, classify_both::Bool)
	if classify_both
		open(FASTQ_file1, "r") do primary_file
			open(FASTQ_file2, "r") do secondary_file
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
						filename = determine_filename(seq, bc_df, max_error_rate, min_delta, mismatch, indel)
						write_fastq_entry(output_dir * "/" * output_prefix1 * "." * filename, header, seq, plus, quality_score)
						write_fastq_entry(output_dir * "/" * output_prefix2 * "." * filename, header2, seq2, plus2, quality_score2)
						mode = "header"
					end
				end
			end
		end
	else
		open(FASTQ_file1, "r") do primary_file
			open(FASTQ_file2, "r") do secondary_file
				header2, seq2, plus2, quality_score2 = "", "", "", ""
				mode = "header"
				filename = ""
				for line1 in eachline(primary_file)
					line2 = readline(secondary_file)
					if line1[1] == '@' && mode == "header"
						header2 = line2
						mode = "seq"
					elseif mode == "seq"
						filename = determine_filename(line1, bc_df, max_error_rate, min_delta, mismatch, indel)
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
end

function classify_sequences(FASTQ_file1::String, bc_df::DataFrame, output_dir::String, output_prefix::String, max_error_rate::Float64, min_delta::Float64, mismatch::Int, indel::Int)
	open(FASTQ_file1, "r") do file
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
				filename = determine_filename(seq, bc_df, max_error_rate, min_delta, mismatch, indel)
				write_fastq_entry(output_dir * "/" * output_prefix * "." * filename, header, seq, plus, quality_score)
				mode = "header"
			end
		end
	end
end