using BioDemuX
using Test
using CodecZlib

@testset "Dual Index Demultiplexing" begin
    # Create temporary files
    tmp_dir = mktempdir()

    # 1. Create Barcode Files
    # BC1: AAAA, CCCC
    # BC2: TTTT, GGGG
    bc_file1 = joinpath(tmp_dir, "bc1.tsv")
    open(bc_file1, "w") do io
        println(io, "Full_seq\tID\tFull_annotation")
        println(io, "AAAA\tID1_A\tBBBB")
        println(io, "CCCC\tID1_C\tBBBB")
    end

    bc_file2 = joinpath(tmp_dir, "bc2.tsv")
    open(bc_file2, "w") do io
        println(io, "Full_seq\tID\tFull_annotation")
        println(io, "TTTT\tID2_T\tBBBB")
        println(io, "GGGG\tID2_G\tBBBB")
    end

    # 2. Create FASTQ File
    # Read 1: [BC1] [Linker] [BC2] [Seq]
    # BC1 (4bp), Linker (4bp), BC2 (4bp)
    # R1: AAAA TATA TTTT ACGT... -> ID1_A.ID2_T
    # R2: CCCC TATA GGGG ACGT... -> ID1_C.ID2_G
    # R3: AAAA TATA GGGG ACGT... -> ID1_A.ID2_G
    # R4: AAAA TATA AAAA ACGT... -> Unknown (BC2 mismatch)

    fastq_file = joinpath(tmp_dir, "test.fastq")
    open(fastq_file, "w") do io
        # Read 1
        println(io, "@read1")
        println(io, "AAAATATATTTTACGT")
        println(io, "+")
        println(io, "IIIIIIIIIIIIIIII")

        # Read 2
        println(io, "@read2")
        println(io, "CCCCTATAGGGGACGT")
        println(io, "+")
        println(io, "IIIIIIIIIIIIIIII")

        # Read 3
        println(io, "@read3")
        println(io, "AAAATATAGGGGACGT")
        println(io, "+")
        println(io, "IIIIIIIIIIIIIIII")

        # Read 4 (Unknown)
        println(io, "@read4")
        println(io, "AAAATATAAAAAACGT")
        println(io, "+")
        println(io, "IIIIIIIIIIIIIIII")
    end

    output_dir = joinpath(tmp_dir, "output")
    mkdir(output_dir)

    # 3. Execute Demultiplexing
    # BC1: 1:4
    # BC2: 9:12 (4+4+1 = 9)

    BioDemuX.execute_demultiplexing(
        fastq_file,
        bc_file1,
        output_dir;
        barcode_file2=bc_file2,
        ref_search_range="1:4",
        ref_search_range2="9:12",
        max_error_rate=0.0,
        chunk_size=100
    )

    # 4. Verify Outputs
    # Expected files: test.ID1_A.ID2_T.fastq, test.ID1_C.ID2_G.fastq, test.ID1_A.ID2_G.fastq, test.unknown.fastq

    @test isfile(joinpath(output_dir, "test.ID1_A.ID2_T.fastq"))
    @test isfile(joinpath(output_dir, "test.ID1_C.ID2_G.fastq"))
    @test isfile(joinpath(output_dir, "test.ID1_A.ID2_G.fastq"))
    @test isfile(joinpath(output_dir, "test.unknown.fastq"))

    # Check content of ID1_A.ID2_T.fastq (should contain read1)
    content = read(joinpath(output_dir, "test.ID1_A.ID2_T.fastq"), String)
    @test occursin("@read1", content)

    # Check content of unknown.fastq (should contain read4)
    content_unknown = read(joinpath(output_dir, "test.unknown.fastq"), String)
    @test occursin("@read4", content_unknown)

    rm(tmp_dir, recursive=true)
end
