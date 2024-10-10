"""
Orchestrates the entire demultiplexing process for FASTQ files.
Handles the preprocessing, dividing, demultiplexing, and merging of files.
"""
function execute_demultiplexing(FASTQ_file1::String, FASTQ_file2::String, bc_file::String, output_dir::String; output_prefix1::String = "", output_prefix2::String = "", max_error_rate::Float64 = 0.2, min_delta::Float64 = 0.1, mismatch::Int = 1, indel::Int = 1, classify_both::Bool = false, bc_complement::Bool = false, bc_rev::Bool = false)
	if !isdir(output_dir)
		mkdir(output_dir)
	end
	workers = nworkers()
	FASTQ_file1, filename1 = preprocess_fastq(FASTQ_file1)
	FASTQ_file2, filename2 = preprocess_fastq(FASTQ_file2)
	if output_prefix1 == ""
		output_prefix1 = filename1
	end
	if output_prefix2 == ""
		output_prefix2 = filename2
	end
	bc_df = preprocess_bc_file(bc_file, bc_complement, bc_rev)
	if workers == 1
		classify_sequences(FASTQ_file1, FASTQ_file2, bc_df, output_dir, output_prefix1, output_prefix2, max_error_rate, min_delta, mismatch, indel, classify_both)
	else
		divided_dir = divide_fastq(FASTQ_file1, FASTQ_file2, output_dir, workers)
		pmap(x -> multi_demultiplex(x, bc_df, divided_dir, output_prefix1, output_prefix2, max_error_rate, min_delta, mismatch, indel, classify_both), 1:workers)

		paths = []
		for (root, dirs, files) in walkdir(divided_dir)
			for file in files
				push!(paths, joinpath(root, file))
			end
		end
		paths = filter(x -> occursin(r"thread", x), paths)
		if classify_both
			paths_file1 = filter(x -> occursin(output_prefix1, x), paths)
			paths_file2 = filter(x -> occursin(output_prefix2, x), paths)
			merge_fastq_files(paths_file1, bc_df, output_dir, output_prefix1)
			merge_fastq_files(paths_file2, bc_df, output_dir, output_prefix2)
		else
			merge_fastq_files(paths, bc_df, output_dir, output_prefix2)
		end
		rm(divided_dir,recursive=true)
	end
end

function execute_demultiplexing(FASTQ_file::String, bc_file::String, output_dir::String; output_prefix::String = "", max_error_rate::Float64 = 0.2, min_delta::Float64 = 0.1, mismatch::Int = 1, indel::Int = 1, bc_complement::Bool = false, bc_rev::Bool = false)
	if !isdir(output_dir)
		mkdir(output_dir)
	end
	workers = nworkers()
	FASTQ_file, filename = preprocess_fastq(FASTQ_file)
	if output_prefix == ""
		output_prefix = filename
	end
	bc_df = preprocess_bc_file(bc_file, bc_complement, bc_rev)
	if workers == 1
		classify_sequences(FASTQ_file, bc_df, output_dir, output_prefix, max_error_rate, min_delta, mismatch, indel)
	else
		divided_dir = divide_fastq(FASTQ_file, output_dir, workers)
		pmap(x -> multi_demultiplex(x, bc_df, divided_dir, output_prefix, max_error_rate, min_delta, mismatch, indel), 1:workers)
		paths = []
		for (root, dirs, files) in walkdir(divided_dir)
			for file in files
				push!(paths, joinpath(root, file))
			end
		end
		paths = filter(x -> occursin(r"thread", x), paths)
		merge_fastq_files(paths, bc_df, output_dir, output_prefix)
		rm(divided_dir,recursive=true)
	end
end