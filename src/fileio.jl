"""
Preprocesses the barcode file by modifying sequences based on specific criteria.
"""
function preprocess_bc_file(bc_file::String, complement::Bool, rev::Bool)
	sequences = String[]
	ids = String[]
	annotations = String[]

	if endswith(lowercase(bc_file), ".fasta")
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
		annotations = ["B" ^ length(seq) for seq in sequences]
	else
		delim = endswith(lowercase(bc_file), ".csv") ? ',' : '\t'
		df = CSV.read(bc_file, DataFrame, delim = delim)
		sequences = Vector{String}(df.Full_seq)
		ids = Vector{String}(df.ID)
		annotations = Vector{String}(df.Full_annotation)
	end

	for i in 1:length(sequences)
		if length(sequences[i]) != length(annotations[i])
			error("Length mismatch between sequence and annotation for ID: $(ids[i])")
		end
		sequences[i] = String([c for (c, a) in zip(sequences[i], annotations[i]) if a == 'B'])
	end
	
	sequences = uppercase.(sequences)
	sequences = replace.(sequences, "U" => "T")
	if complement == true
		sequences = replace.(sequences, "A" => "T", "T" => "A", "G" => "C", "C" => "G")
	end
	if rev == true
		sequences = reverse.(sequences)
	end
	
	bc_lengths_no_N = [count(c -> c != 'N', bc) for bc in sequences]
	
	return sequences, bc_lengths_no_N, ids
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
