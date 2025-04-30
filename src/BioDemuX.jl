__precompile__()

module BioDemuX

export
	semiglobal_alignment,
	find_best_matching_bc,
	determine_filename,
	get_fastq_io,
	close_all_fastq_ios,
	write_fastq_entry,
	classify_sequences,

	preprocess_bc_file,
	read_fastq,
	write_fastq,
	divide_fastq,
	multi_demultiplex,
	merge_fastq_files,
	execute_demultiplexing


using DataFrames, CSV, Distributed, CodecZlib, BufferedStreams
include("classification.jl")
include("fileio.jl")
include("core.jl")

end#module
