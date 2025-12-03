"""
Orchestrates the entire demultiplexing process for FASTQ files using a channel-based streaming approach.
"""

struct FastqChunk
    headers::Vector{String}
    seqs::Vector{String}
    pluses::Vector{String}
    quals::Vector{String}
end

function FastqChunk(capacity::Int)
    FastqChunk(
        sizehint!(String[], capacity),
        sizehint!(String[], capacity),
        sizehint!(String[], capacity),
        sizehint!(String[], capacity)
    )
end

function Base.empty!(chunk::FastqChunk)
    empty!(chunk.headers)
    empty!(chunk.seqs)
    empty!(chunk.pluses)
    empty!(chunk.quals)
end

function Base.push!(chunk::FastqChunk, h, s, p, q)
    push!(chunk.headers, h)
    push!(chunk.seqs, s)
    push!(chunk.pluses, p)
    push!(chunk.quals, q)
end

struct Chunk
    id::Int
    data::FastqChunk
    data2::Union{FastqChunk, Nothing} # For paired-end
end

function reader_task(input_channel::Channel{Chunk}, recycle_channel::Channel{Chunk}, file1::String, file2::Union{String, Nothing}, chunk_size::Int)
    read_fastq(file1) do io1
        if !isnothing(file2)
            read_fastq(file2) do io2
                chunk_id = 1
                while !eof(io1) && !eof(io2)
                    # Try to recycle a chunk, or create a new one
                    local c1::FastqChunk
                    local c2::FastqChunk
                    if isready(recycle_channel)
                        chunk = take!(recycle_channel)
                        # We reuse the FastqChunk objects but create a new Chunk wrapper with updated ID
                        c1 = chunk.data
                        c2 = chunk.data2
                        empty!(c1)
                        empty!(c2)
                    else
                        c1 = FastqChunk(chunk_size)
                        c2 = FastqChunk(chunk_size)
                    end
                    
                    for _ in 1:chunk_size
                        if eof(io1) || eof(io2) break end
                        
                        # Read 4 lines for file 1
                        push!(c1, readline(io1), readline(io1), readline(io1), readline(io1))
                        
                        # Read 4 lines for file 2
                        push!(c2, readline(io2), readline(io2), readline(io2), readline(io2))
                    end
                    
                    if !isempty(c1.headers)
                        put!(input_channel, Chunk(chunk_id, c1, c2))
                        chunk_id += 1
                    end
                end
            end
        else
            chunk_id = 1
            while !eof(io1)
                local c1::FastqChunk
                if isready(recycle_channel)
                    chunk = take!(recycle_channel)
                    c1 = chunk.data
                    empty!(c1)
                else
                    c1 = FastqChunk(chunk_size)
                end
                
                for _ in 1:chunk_size
                    if eof(io1) break end
                    push!(c1, readline(io1), readline(io1), readline(io1), readline(io1))
                end
                
                if !isempty(c1.headers)
                    put!(input_channel, Chunk(chunk_id, c1, nothing))
                    chunk_id += 1
                end
            end
        end
    end
end

struct ResultChunk
    chunk::Chunk
    filenames::Vector{String}
end

function worker_task(input_channel::Channel{Chunk}, output_channel::Channel{ResultChunk}, config::DemuxConfig)
    # Create thread-local workspace
    max_m = maximum(ncodeunits.(config.bc_seqs))
    ws = SemiGlobalWorkspace(max_m)
    
    for chunk in input_channel
        n = length(chunk.data.headers)
        filenames = Vector{String}(undef, n)
        
        for i in 1:n
            seq = chunk.data.seqs[i]
            filenames[i] = determine_filename(seq, config, ws)
        end
        
        put!(output_channel, ResultChunk(chunk, filenames))
    end
end

function writer_task(output_channel::Channel{ResultChunk}, recycle_channel::Channel{Chunk}, output_dir::String, output_prefix1::String, output_prefix2::String, config::DemuxConfig)
    # Dictionary to keep track of open file handles
    file_handles = Dict{String, BufferedOutputStream}()
    
    function get_handle(filename::String)
        if !haskey(file_handles, filename)
            path = joinpath(output_dir, filename)
            raw = open(path, "a")
            io = endswith(lowercase(path), ".gz") ? GzipCompressorStream(raw) : raw
            file_handles[filename] = BufferedOutputStream(io)
        end
        return file_handles[filename]
    end
    
    function write_entry(io, h, s, p, q)
        write(io, h, "\n", s, "\n", p, "\n", q, "\n")
    end

    # Buffer for reordering chunks
    pending_chunks = Dict{Int, ResultChunk}()
    next_chunk_id = 1
    
    for result_chunk in output_channel
        pending_chunks[result_chunk.chunk.id] = result_chunk
        
        while haskey(pending_chunks, next_chunk_id)
            rc = pending_chunks[next_chunk_id]
            delete!(pending_chunks, next_chunk_id)
            
            chunk = rc.chunk
            filenames = rc.filenames
            n = length(chunk.data.headers)
            
            for i in 1:n
                filename = filenames[i]
                
                if config.classify_both && !isnothing(chunk.data2)
                    # Paired end, classify both
                    # Output 1
                    fname1 = output_prefix1 * "." * filename
                    io1 = get_handle(fname1)
                    write_entry(io1, chunk.data.headers[i], chunk.data.seqs[i], chunk.data.pluses[i], chunk.data.quals[i])
                    
                    # Output 2
                    fname2 = output_prefix2 * "." * filename
                    io2 = get_handle(fname2)
                    write_entry(io2, chunk.data2.headers[i], chunk.data2.seqs[i], chunk.data2.pluses[i], chunk.data2.quals[i])
                elseif !isnothing(chunk.data2)
                    # Paired end, classify only read 2 based on read 1
                    fname = output_prefix2 * "." * filename
                    io = get_handle(fname)
                    write_entry(io, chunk.data2.headers[i], chunk.data2.seqs[i], chunk.data2.pluses[i], chunk.data2.quals[i])
                else
                    # Single end
                    fname = output_prefix1 * "." * filename
                    io = get_handle(fname)
                    write_entry(io, chunk.data.headers[i], chunk.data.seqs[i], chunk.data.pluses[i], chunk.data.quals[i])
                end
            end
            
            # Recycle the chunk
            # Empty the FastqChunk objects within the Chunk before recycling
            empty!(chunk.data)
            if !isnothing(chunk.data2)
                empty!(chunk.data2)
            end
            
            # Only put back if there is space, otherwise let GC handle it
            # This prevents deadlock if the system created more chunks than recycle capacity
            lock(recycle_channel)
            try
                if Base.n_avail(recycle_channel) < recycle_channel.sz_max
                    put!(recycle_channel, chunk)
                end
            finally
                unlock(recycle_channel)
            end
            
            next_chunk_id += 1
        end
    end
    
    # Close all handles
    for io in values(file_handles)
        flush(io)
        close(io)
    end
end

function execute_demultiplexing(FASTQ_file1::String, FASTQ_file2::String, bc_file::String, output_dir::String; output_prefix1::String = "", output_prefix2::String = "", gzip_output::Union{Nothing, Bool} = nothing, max_error_rate::Float64 = 0.2, min_delta::Float64 = 0.0, mismatch::Int = 1, indel::Int = 1, nindel::Union{Int, Nothing} = nothing, classify_both::Bool = false, bc_complement::Bool = false, bc_rev::Bool = false, ref_search_range::String = "1:end", barcode_start_range::String = "1:end", barcode_end_range::String = "1:end", chunk_size::Int = 4000, channel_capacity::Int = 64)
	if !isdir(output_dir)
		mkdir(output_dir)
	end
	
	should_gzip = endswith(lowercase(FASTQ_file1), ".gz")
	gzip_output = isnothing(gzip_output) ? should_gzip : gzip_output
	
	if output_prefix1 == ""
		output_prefix1 = replace(basename(FASTQ_file1), r"\.fastq(\.gz)?$|\.fq(\.gz)?$" => "")
	end
	if output_prefix2 == ""
		output_prefix2 = replace(basename(FASTQ_file2), r"\.fastq(\.gz)?$|\.fq(\.gz)?$" => "")
	end
	
	if indel <= 0
		error("indel must be > 0")
	end
	if !isnothing(nindel) && nindel <= 0
		error("nindel must be > 0")
	end
	
	bc_seqs, bc_lengths_no_N, ids = preprocess_bc_file(bc_file, bc_complement, bc_rev)
	
    # We do NOT create ws here, it is created in worker_task
	config = DemuxConfig(max_error_rate, min_delta, 0, mismatch, indel, nindel, classify_both, gzip_output, parse_dynamic_range(ref_search_range), parse_dynamic_range(barcode_start_range), parse_dynamic_range(barcode_end_range), bc_seqs, bc_lengths_no_N, ids)

    # Channel setup
    input_channel = Channel{Chunk}(channel_capacity)
    output_channel = Channel{ResultChunk}(channel_capacity)
    recycle_channel = Channel{Chunk}(2 * channel_capacity) # Buffer for recycling
    
    # Start tasks
    
    # Reader
    reader = Threads.@spawn begin
        try
            reader_task(input_channel, recycle_channel, FASTQ_file1, FASTQ_file2, chunk_size)
        finally
            close(input_channel)
        end
    end
    
    # Writer
    writer = Threads.@spawn begin
        writer_task(output_channel, recycle_channel, output_dir, output_prefix1, output_prefix2, config)
    end
    
    # Workers
    workers = [Threads.@spawn worker_task(input_channel, output_channel, config) for _ in 1:Threads.nthreads()]
    
    # Monitor workers
    Threads.@spawn begin
        for w in workers
            wait(w)
        end
        close(output_channel)
    end
    
    wait(writer)
end

function execute_demultiplexing(FASTQ_file::String, bc_file::String, output_dir::String; output_prefix::String = "", gzip_output::Union{Nothing, Bool} = nothing, max_error_rate::Float64 = 0.2, min_delta::Float64 = 0.0, mismatch::Int = 1, indel::Int = 1, nindel::Union{Int, Nothing} = nothing, bc_complement::Bool = false, bc_rev::Bool = false, ref_search_range::String = "1:end", barcode_start_range::String = "1:end", barcode_end_range::String = "1:end", chunk_size::Int = 4000, channel_capacity::Int = 64)
	if !isdir(output_dir)
		mkdir(output_dir)
	end
	
	should_gzip = endswith(lowercase(FASTQ_file), ".gz")
	gzip_output = isnothing(gzip_output) ? should_gzip : gzip_output
	
	if output_prefix == ""
		output_prefix = replace(basename(FASTQ_file), r"\.fastq(\.gz)?$|\.fq(\.gz)?$" => "")
	end
	
	if indel <= 0
		error("indel must be > 0")
	end
	if !isnothing(nindel) && nindel <= 0
		error("nindel must be > 0")
	end
	
	bc_seqs, bc_lengths_no_N, ids = preprocess_bc_file(bc_file, bc_complement, bc_rev)
	
	config = DemuxConfig(max_error_rate, min_delta, 0, mismatch, indel, nindel, false, gzip_output, parse_dynamic_range(ref_search_range), parse_dynamic_range(barcode_start_range), parse_dynamic_range(barcode_end_range), bc_seqs, bc_lengths_no_N, ids)

    # Channel setup
    input_channel = Channel{Chunk}(channel_capacity)
    output_channel = Channel{ResultChunk}(channel_capacity)
    recycle_channel = Channel{Chunk}(2 * channel_capacity)
    
    reader = Threads.@spawn begin
        try
            reader_task(input_channel, recycle_channel, FASTQ_file, nothing, chunk_size)
        finally
            close(input_channel)
        end
    end
    
    # For single file, we pass output_prefix as output_prefix1 to writer_task
    writer = Threads.@spawn begin
        writer_task(output_channel, recycle_channel, output_dir, output_prefix, "", config)
    end
    
    workers = [Threads.@spawn worker_task(input_channel, output_channel, config) for _ in 1:Threads.nthreads()]
    
    Threads.@spawn begin
        for w in workers
            wait(w)
        end
        close(output_channel)
    end
    
    wait(writer)
end