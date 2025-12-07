__precompile__()

module BioDemuX

export
    semiglobal_alignment,
    find_best_matching_bc,
    determine_filename, preprocess_bc_file,
    read_fastq,
    write_fastq,
    execute_demultiplexing,
    DemuxConfig


using DataFrames, CSV
using CodecZlib
using Dates
using BufferedStreams
include("classification.jl")
include("fileio.jl")
include("reporting.jl")
include("core.jl")

end#module
