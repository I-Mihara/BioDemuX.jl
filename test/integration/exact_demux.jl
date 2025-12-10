
@testset "Exact Demultiplexing Integration" begin
    mktempdir() do output_dir
        # Create barcodes
        bc_file = joinpath(output_dir, "barcodes.csv")
        open(bc_file, "w") do io
            println(io, "ID,Full_seq,Full_annotation")
            # Strict barcodes
            println(io, "BC1,ACGTAC,BBBBBB")
            println(io, "BC2,CCCCCC,BBBBBB")
        end

        # Create Reads
        fastq_file = joinpath(output_dir, "reads.fastq")
        open(fastq_file, "w") do io
            write(io, "@read1\nACGTAC\n+\nIIIIII\n") # Perfect BC1
            write(io, "@read2\nCCCCCC\n+\nIIIIII\n") # Perfect BC2
            write(io, "@read3\nACATAC\n+\nIIIIII\n") # 1 mismatch BC1 (A vs G). Should be UNKNOWN in Exact mode.
            write(io, "@read4\nACGTAG\n+\nIIIIII\n") # 1 mismatch BC1 (G vs C). Should be UNKNOWN.
            write(io, "@read_indel\nACGGTAC\n+\nIIIIIII\n") # Insertion. Should be UNKNOWN.
        end

        BioDemuX.execute_demultiplexing(
            fastq_file, bc_file, output_dir,
            matching_algorithm=:exact,
            max_error_rate=0.0 # Strict
        )

        @test isfile(joinpath(output_dir, "reads.BC1.fastq"))
        @test isfile(joinpath(output_dir, "reads.BC2.fastq"))
        @test isfile(joinpath(output_dir, "reads.unknown.fastq"))

        # Count lines
        count_lines(f) = count(_ -> true, eachline(f))

        # BC1: read1 -> 4 lines
        @test count_lines(joinpath(output_dir, "reads.BC1.fastq")) == 4
        # BC2: read2 -> 4 lines
        @test count_lines(joinpath(output_dir, "reads.BC2.fastq")) == 4

        # Unknown: read3, read4, read_indel -> 12 lines
        @test count_lines(joinpath(output_dir, "reads.unknown.fastq")) == 12
    end
end
