"""
Orchestrates the entire demultiplexing process for FASTQ files.
Handles the preprocessing, dividing, demultiplexing, and merging of files.
"""
function execute_demultiplexing(FASTQ_file1::String, FASTQ_file2::String, bc_file::String, output_dir::String; output_prefix1::String = "", output_prefix2::String = "", gzip_output::Union{Nothing, Bool} = nothing, max_error_rate::Float64 = 0.2, min_delta::Float64 = 0.1, mismatch::Int = 1, indel::Int = 1, classify_both::Bool = false, bc_complement::Bool = false, bc_rev::Bool = false)
	if !isdir(output_dir)
		mkdir(output_dir)
	end
	workers = nworkers()
	should_gzip = endswith(lowercase(FASTQ_file1), ".gz")
	gzip_output = isnothing(gzip_output) ? should_gzip : gzip_output
	if output_prefix1 == ""
		output_prefix1 = replace(basename(FASTQ_file1), r"\.fastq(\.gz)?$|\.fq(\.gz)?$" => "")
	end
	if output_prefix2 == ""
		output_prefix2 = replace(basename(FASTQ_file2), r"\.fastq(\.gz)?$|\.fq(\.gz)?$" => "")
	end
	bc_df = preprocess_bc_file(bc_file, bc_complement, bc_rev)
	if workers == 1
		classify_sequences(FASTQ_file1, FASTQ_file2, bc_df, output_dir, output_prefix1, output_prefix2, max_error_rate, min_delta, mismatch, indel, classify_both, gzip_output)
	else
		divided_dir = divide_fastq(FASTQ_file1, FASTQ_file2, output_dir, workers, gzip_output)
		pmap(x -> multi_demultiplex(x, bc_df, divided_dir, output_prefix1, output_prefix2, max_error_rate, min_delta, mismatch, indel, classify_both, gzip_output), 1:workers)

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
			merge_fastq_files(paths_file1, bc_df, output_dir, output_prefix1, gzip_output)
			merge_fastq_files(paths_file2, bc_df, output_dir, output_prefix2, gzip_output)
		else
			merge_fastq_files(paths, bc_df, output_dir, output_prefix2, gzip_output)
		end
		rm(divided_dir,recursive=true)
	end
end

function execute_demultiplexing(FASTQ_file::String, bc_file::String, output_dir::String; output_prefix::String = "", gzip_output::Union{Nothing, Bool} = nothing, max_error_rate::Float64 = 0.2, min_delta::Float64 = 0.1, mismatch::Int = 1, indel::Int = 1, bc_complement::Bool = false, bc_rev::Bool = false)
	if !isdir(output_dir)
		mkdir(output_dir)
	end
	workers = nworkers()
	should_gzip = endswith(lowercase(FASTQ_file), ".gz")
	gzip_output = isnothing(gzip_output) ? should_gzip : gzip_output
	if output_prefix == ""
		output_prefix = replace(basename(FASTQ_file), r"\.fastq(\.gz)?$|\.fq(\.gz)?$" => "")
	end
	bc_df = preprocess_bc_file(bc_file, bc_complement, bc_rev)
	if workers == 1
		classify_sequences(FASTQ_file, bc_df, output_dir, output_prefix, max_error_rate, min_delta, mismatch, indel, gzip_output)
	else
		divided_dir = divide_fastq(FASTQ_file, output_dir, workers, gzip_output)
		pmap(x -> multi_demultiplex(x, bc_df, divided_dir, output_prefix, max_error_rate, min_delta, mismatch, indel, gzip_output), 1:workers)
		paths = []
		for (root, dirs, files) in walkdir(divided_dir)
			for file in files
				push!(paths, joinpath(root, file))
			end
		end
		paths = filter(x -> occursin(r"thread", x), paths)
		merge_fastq_files(paths, bc_df, output_dir, output_prefix, gzip_output)
		rm(divided_dir,recursive=true)
	end
end