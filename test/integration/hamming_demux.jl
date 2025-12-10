@testset "Hamming Distance Integration" begin
    # Test Single End, Hamming Mode
    @testset "Single End Hamming" begin
        mktempdir() do output_dir
            # Create barcodes
            bc_file = joinpath(output_dir, "barcodes.csv")
            open(bc_file, "w") do io
                println(io, "ID,Full_seq,Full_annotation")
                println(io, "BC1,ACGTAC,BBBBBB")
                println(io, "BC2,CCCCCC,BBBBBB")
            end

            # Create Reads
            fastq_file = joinpath(output_dir, "reads.fastq")
            open(fastq_file, "w") do io
                write(io, "@read1\nACGTAC\n+\nIIIIII\n") # Perfect BC1
                write(io, "@read2\nCCCCCC\n+\nIIIIII\n") # Perfect BC2
                write(io, "@read3\nACATAC\n+\nIIIIII\n") # 1 mismatch BC1 (A vs G) -> Score 1/6 = 0.166 <= 0.2. Match.
                write(io, "@read4\nACGTAG\n+\nIIIIII\n") # 1 mismatch BC1 (G vs C) -> Score 1/6 = 0.166 <= 0.2. Match.
                write(io, "@read_indel\nACGGTAC\n+\nIIIIIII\n") # Insertion G. "ACGTAC" vs "ACGGTAC".
                # Pos 1: A-A(m), C-C(m), G-G(m), T-G(mis), A-T(mis), C-A(mis). 3 mismatches. Score 0.5 > 0.2. Fail.
            end

            BioDemuX.execute_demultiplexing(
                fastq_file, bc_file, output_dir,
                matching_algorithm=:hamming,
                max_error_rate=0.2
            )

            @test isfile(joinpath(output_dir, "reads.BC1.fastq"))
            @test isfile(joinpath(output_dir, "reads.BC2.fastq"))
            @test isfile(joinpath(output_dir, "reads.unknown.fastq"))

            # Count lines
            count_lines(f) = count(_ -> true, eachline(f))

            # BC1: read1, read3, read4 -> 12 lines
            @test count_lines(joinpath(output_dir, "reads.BC1.fastq")) == 12
            # BC2: read2 -> 4 lines
            @test count_lines(joinpath(output_dir, "reads.BC2.fastq")) == 4
            # Unknown: read_indel -> 4 lines
            @test count_lines(joinpath(output_dir, "reads.unknown.fastq")) == 4
        end
    end

    # Test Paired End, Hamming Mode
    @testset "Paired End Hamming" begin
        mktempdir() do output_dir
            bc_file = joinpath(output_dir, "barcodes.csv")
            open(bc_file, "w") do io
                println(io, "ID,Full_seq,Full_annotation")
                println(io, "BC1,AAAA,BBBB")
            end

            # Create Dummy Paired Reads
            fq1 = joinpath(output_dir, "R1.fastq")
            fq2 = joinpath(output_dir, "R2.fastq")

            open(fq1, "w") do io
                write(io, "@seq1\nAAAA\n+\nIIII\n")
            end
            open(fq2, "w") do io
                write(io, "@seq1\nGGGG\n+\nIIII\n")
            end

            BioDemuX.execute_demultiplexing(
                fq1, fq2, bc_file, output_dir,
                matching_algorithm=:hamming,
                max_error_rate=0.0,
                classify_both=true
            )

            @test isfile(joinpath(output_dir, "R1.BC1.fastq"))
            @test isfile(joinpath(output_dir, "R2.BC1.fastq"))
        end
    end
end
