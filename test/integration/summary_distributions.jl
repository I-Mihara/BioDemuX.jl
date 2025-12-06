using Test
using BioDemuX

@testset "Summary Mode Distributions" begin
    mktempdir() do temp_dir
        # Setup
        bc_file = joinpath(temp_dir, "barcodes.fasta")
        open(bc_file, "w") do io
            write(io, ">BC1\nAAAA\n")
        end

        fastq_file = joinpath(temp_dir, "reads.fastq")
        open(fastq_file, "w") do io
            write(io, "@read1\nAAAA\n+\nIIII\n")       # Pos 1, Len 4
            write(io, "@read2\nNAAAA\n+\nIIIII\n")      # Pos 2, Len 4
            write(io, "@read3\nNNAAAA\n+\nIIIIII\n")    # Pos 3, Len 4
            write(io, "@read4\nAAAT\n+\nIIII\n")        # Pos 1, Len 4 (1 mismatch)
            write(io, "@read5\nAAA\n+\nIII\n")          # Pos 1, Len 3 (deletion)
        end

        output_dir = joinpath(temp_dir, "output_dist")

        BioDemuX.execute_demultiplexing(
            fastq_file, bc_file, output_dir;
            max_error_rate=0.3,
            summary=true,
            summary_format=:json
        )

        @test isfile(joinpath(output_dir, "summary.json"))
        content = read(joinpath(output_dir, "summary.json"), String)

        # Verify Global Stats Keys
        @test occursin("\"bc1_pos_counts\": {", content)
        @test occursin("\"bc1_len_counts\": {", content)
        @test occursin("\"bc1_score_counts\": {", content)

        # Verify Specific Counts
        # Positions: 1 (3 reads), 2 (1 read), 3 (1 read)
        @test occursin("\"1\":", content)
        @test occursin("\"2\":", content)
        @test occursin("\"3\":", content)

        # Lengths: 3 (1 read), 4 (4 reads)
        @test occursin("\"3\":", content)
        @test occursin("\"4\":", content)

        # Verify Per-Barcode Stats Keys
        @test occursin("\"bc1_per_bc_pos_counts\": {", content)
        @test occursin("\"bc1_per_bc_len_counts\": {", content)
        @test occursin("\"bc1_per_bc_score_counts\": {", content)
    end
end
