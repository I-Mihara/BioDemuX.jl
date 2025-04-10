"""
Preprocesses the barcode file by modifying sequences based on specific criteria.
"""
function preprocess_bc_file(bc_file::String, complement::Bool, rev::Bool)
	if endswith(lowercase(bc_file), ".fasta")
        sequences = String[]
        ids = String[]
        current_seq = ""
        current_id = ""
        open(bc_file, "r") do io
            for line in eachline(io)
                if startswith(line, '>')
                    if !isempty(current_seq)
                        push!(sequences, current_seq)
                        current_seq = ""
                    end
                    current_id = replace(strip(line[2:end]), r"\s.*$" => "")
                    push!(ids, current_id)
                else
                    current_seq *= strip(line)
                end
            end
            if !isempty(current_seq)
                push!(sequences, current_seq)
            end
        end
        bc_df = DataFrame(
            ID = ids,
            Full_seq = sequences,
            Full_annotation = ["B" ^ length(seq) for seq in sequences]
        )
    else
        delim = endswith(lowercase(bc_file), ".csv") ? ',' : '\t'
        bc_df = CSV.read(bc_file, DataFrame, delim = delim)
    end

	for i in 1:nrow(bc_df)
		prefix_region = 1:findfirst('B', bc_df.Full_annotation[i])-1
		suffix_region = findlast('B', bc_df.Full_annotation[i])+1:length(bc_df.Full_annotation[i])
		prefix = SubString(bc_df.Full_seq[i], prefix_region)
		suffix = SubString(bc_df.Full_seq[i], suffix_region)
		bc_df.Full_seq[i] = SubString(bc_df.Full_seq[i], length(prefix)+1)
		bc_df.Full_seq[i] = SubString(bc_df.Full_seq[i], 1, length(bc_df.Full_seq[i])-length(suffix))
	end
	bc_df.Full_seq = uppercase.(bc_df.Full_seq)
	bc_df.Full_seq = replace.(bc_df.Full_seq, "U" => "T")
	if complement == true
		bc_df.Full_seq = replace.(bc_df.Full_seq, "A" => "T", "T" => "A", "G" => "C", "C" => "G")
	end
	if rev == true
		bc_df.Full_seq = reverse.(bc_df.Full_seq)
	end
	return bc_df
end

"""
Open a fastq file and apply a function to it with the correct decompressor.
"""
function read_fastq(f::Function, filepath::String)
    if endswith(lowercase(filepath), ".gz")
        open(GzipDecompressorStream, filepath, "r") do io
            f(io)
        end
    else
        open(filepath, "r") do io
            f(io)
        end
    end
end

"""
Open a fastq file and apply a function to it with the correct compressor.
"""
function write_fastq(f::Function, filepath::String)
    if endswith(lowercase(filepath), ".gz")
        open(GzipCompressorStream, filepath, "a") do io
            f(io)
        end
    else
        open(filepath, "a") do io
            f(io)
        end
    end
end

function count_fastq_lines(file::String)
    total_lines = 0
    if endswith(lowercase(file), ".gz")
        open(GzipDecompressorStream, file, "r") do io
            for _ in eachline(io)
                total_lines += 1
            end
        end
    else
        total_lines = countlines(file)
    end
    return total_lines
end

"""
Divides a pair of FASTQ files into smaller parts for parallel processing.
It calculates the number of reads per worker and uses the split command to divide the files.
"""
function divide_fastq(FASTQ_file1::String, FASTQ_file2::String, output_dir::String, workers::Int, gzip_output::Bool)
	divided_dir = mktempdir(output_dir)
	total_lines = count_fastq_lines(FASTQ_file1)
	total_reads = total_lines ÷ 4
	reads_per_worker = cld(total_reads, workers)
	lines_per_worker = reads_per_worker * 4
	read_fastq(FASTQ_file1) do in_io1
		for i in 0:workers-1
			out_file1 = "$divided_dir/file1_$(lpad(i,5,"0")).fastq" * (gzip_output ? ".gz" : "")
			write_fastq(out_file1) do out_io1
				for _ in 1:lines_per_worker
					if eof(in_io1) break end
					write(out_io1, readline(in_io1, keep=true))
				end
			end
		end
	end
	read_fastq(FASTQ_file2) do in_io2
		for i in 0:workers-1
			out_file2 = "$divided_dir/file2_$(lpad(i,5,"0")).fastq" * (gzip_output ? ".gz" : "")
			write_fastq(out_file2) do out_io2
				for _ in 1:lines_per_worker
					if eof(in_io2) break end
					write(out_io2, readline(in_io2, keep=true))
				end
			end
		end
	end
	return divided_dir
end

"""
Divides a single FASTQ file for parallel processing.
"""
function divide_fastq(FASTQ_file::String, output_dir::String, workers::Int, gzip_output::Bool)
	divided_dir = mktempdir(output_dir)
	total_lines = count_fastq_lines(FASTQ_file)
	total_reads = total_lines ÷ 4
	reads_per_worker = cld(total_reads, workers)
	lines_per_worker = reads_per_worker * 4
	read_fastq(FASTQ_file) do in_io
		for i in 0:workers-1
			out_file = "$divided_dir/$(lpad(i,5,"0")).fastq" * (gzip_output ? ".gz" : "")
			write_fastq(out_file) do out_io
				for _ in 1:lines_per_worker
					if eof(in_io) break end
					write(out_io, readline(in_io, keep=true))
				end
			end
		end
	end
	return divided_dir
end

function merge_fastq_files(paths::Vector, bc_df::DataFrame, output_dir::String, prefix::String, gzip_output::Bool)
	paths_unknown = filter(x -> occursin(r".*unknown.fastq", x), paths)
	if paths_unknown != []
		output_file = "$output_dir/$prefix.unknown.fastq" * (gzip_output ? ".gz" : "")
		write_fastq(output_file) do out_io
			for path in paths_unknown
				read_fastq(path) do in_io
					write(out_io, read(in_io))
				end
			end
		end
	end

	paths_ambiguous_classification = filter(x -> occursin(r".*ambiguous_classification.fastq", x), paths)
	if paths_ambiguous_classification != []
		output_file = "$output_dir/$prefix.ambiguous_classification.fastq" * (gzip_output ? ".gz" : "")
		write_fastq(output_file) do out_io
			for path in paths_ambiguous_classification
				read_fastq(path) do in_io
					write(out_io, read(in_io))
				end
			end
		end
	end

	for i in 1:nrow(bc_df)
		regex = string(bc_df.ID[i]) * ".fastq"
		paths_matched = filter(x -> occursin(regex, x), paths)
		if paths_matched != []
			output_file = "$output_dir/$prefix.$(bc_df.ID[i]).fastq" * (gzip_output ? ".gz" : "")
			write_fastq(output_file) do out_io
				for path in paths_matched
					read_fastq(path) do in_io
						write(out_io, read(in_io))
					end
				end
			end
		end
	end
end

function multi_demultiplex(thread_num::Int, bc_df::DataFrame, divided_dir::String, output_prefix1::String, output_prefix2::String, max_error_rate::Float64, min_delta::Float64, mismatch::Int, indel::Int, classify_both::Bool, gzip_output::Bool)
	FASTQ_file1 = divided_dir * "/file1_" * lpad((thread_num - 1), 5, "0") * ".fastq" * (gzip_output ? ".gz" : "")
	FASTQ_file2 = divided_dir * "/file2_" * lpad((thread_num - 1), 5, "0") * ".fastq" * (gzip_output ? ".gz" : "")
	mkdir(divided_dir * "/thread" * string(thread_num))
	output_dir = divided_dir * "/thread" * string(thread_num)
	classify_sequences(FASTQ_file1, FASTQ_file2, bc_df, output_dir, output_prefix1, output_prefix2, max_error_rate, min_delta, mismatch, indel, classify_both, gzip_output)
end

function multi_demultiplex(thread_num::Int, bc_df::DataFrame, divided_dir::String, output_prefix::String, max_error_rate::Float64, min_delta::Float64, mismatch::Int, indel::Int, gzip_output::Bool)
	FASTQ_file = divided_dir * "/" * lpad((thread_num - 1), 5, "0") * ".fastq" * (gzip_output ? ".gz" : "")
	mkdir(divided_dir * "/thread" * string(thread_num))
	output_dir = divided_dir * "/thread" * string(thread_num)
	classify_sequences(FASTQ_file, bc_df, output_dir, output_prefix, max_error_rate, min_delta, mismatch, indel, gzip_output)
end