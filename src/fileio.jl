"""
    preprocess_bc_file(bc_file::String, complement::Bool, rev::Bool)

Preprocesses the barcode file (FASTA or CSV/TSV).
Returns a tuple: `(sequences, lengths_no_N, ids)`.
"""
function preprocess_bc_file(bc_file::String, complement::Bool, rev::Bool)
    sequences = String[]
    ids = String[]

    if endswith(lowercase(bc_file), ".fasta") || endswith(lowercase(bc_file), ".fa")
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
        # For FASTA, assume annotation is all 'B's (all bases used)
        annotations = ["B"^length(seq) for seq in sequences]
    else
        # CSV or TSV
        delim = endswith(lowercase(bc_file), ".csv") ? ',' : '\t'
        df = CSV.read(bc_file, DataFrame, delim=delim)
        sequences = Vector{String}(df.Full_seq)
        ids = Vector{String}(df.ID)
        annotations = Vector{String}(df.Full_annotation)
    end



    for i in 1:length(sequences)
        if length(sequences[i]) != length(annotations[i])
            error("Length mismatch between sequence and annotation for ID: $(ids[i])")
        end
        # Filter for 'B' annotated bases only
        sequences[i] = String([c for (c, a) in zip(sequences[i], annotations[i]) if a == 'B'])
    end



    sequences = uppercase.(sequences)
    sequences = replace.(sequences, "U" => "T")

    if complement
        complement_dict = Dict(
            'A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G',
            'a' => 't', 't' => 'a', 'g' => 'c', 'c' => 'g',
            'N' => 'N', 'n' => 'n'
        )
        sequences = map(seq -> map(c -> get(complement_dict, c, c), seq), sequences)
    end
    if rev
        sequences = reverse.(sequences)
    end

    bc_lengths_no_N = [count(c -> c != 'N', bc) for bc in sequences]

    return sequences, bc_lengths_no_N, ids
end

"""
Helper to open files with gzip support.
"""
function smart_open(f::Function, filepath::String, mode::String)
    is_gzip = endswith(lowercase(filepath), ".gz")

    if is_gzip
        if mode == "r"
            open(GzipDecompressorStream, filepath, mode) do io
                f(io)
            end
        else
            open(GzipCompressorStream, filepath, mode) do io
                f(io)
            end
        end
    else
        open(filepath, mode) do io
            f(io)
        end
    end
end

"""
    read_fastq(f::Function, filepath::String)

Open a FASTQ file (potentially gzipped) for reading.
"""
function read_fastq(f::Function, filepath::String)
    smart_open(f, filepath, "r")
end

"""
    write_fastq(f::Function, filepath::String)

Open a FASTQ file (potentially gzipped) for appending.
"""
function write_fastq(f::Function, filepath::String)
    smart_open(f, filepath, "a")
end
