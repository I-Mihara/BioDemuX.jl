"""
BioDemuX: High-performance, flexible demultiplexing for FASTQ files.

Exports the main `execute_demultiplexing` function and configuration structures.
"""
module BioDemuX

export
    # Main entry point
    execute_demultiplexing,

    # Configuration
    DemuxConfig,

    # Low-level API (exposed for advanced usage or testing)
    semiglobal_alignment,
    find_best_matching_bc,
    determine_filename,
    preprocess_bc_file,
    read_fastq,
    write_fastq

using DataFrames, CSV
using CodecZlib
using Dates
using BufferedStreams

include("classification.jl")
include("fileio.jl")
include("templates.jl")
include("reporting.jl")
include("core.jl")

end#module
