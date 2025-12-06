using Test
using BioDemuX
using BioDemuX: SemiGlobalWorkspace, semiglobal_alignment, determine_filename, DemuxConfig

@testset "Trimming Tests" begin
    # Setup
    ws = SemiGlobalWorkspace(100, true) # Enable trimming in workspace

    @testset "3' Trimming" begin
        # Read: PREFIX + BARCODE + SUFFIX
        # Keep PREFIX

        read = "AAAAATTTTTCCCCC"
        bc = "TTTTT"
        # Match at 6:10
        # 3' Trim: Keep 1:5 ("AAAAA")

        # Test semiglobal_alignment directly
        (score, start_pos, end_pos) = semiglobal_alignment(ws, bc, read, 0.0, 0, 1, 1, 1:length(read), 100, 1, 3)
        @test score == 0.0
        @test start_pos == 6
        @test end_pos == 10

        # Test determine_filename
        config = DemuxConfig(
            bc_seqs=[bc], bc_lengths_no_N=[5], ids=["id1"],
            trim_side=3,
            ref_search_range=BioDemuX.parse_dynamic_range("1:end"),
            barcode_start_range=BioDemuX.parse_dynamic_range("1:end"),
            barcode_end_range=BioDemuX.parse_dynamic_range("1:end")
        )

        fname, t_start, t_end = determine_filename(read, config, ws)
        trim_range = t_start:t_end
        @test fname == "id1.fastq"
        @test trim_range == 1:5
        @test read[trim_range] == "AAAAA"
    end

    @testset "5' Trimming" begin
        # Read: PREFIX + BARCODE + SUFFIX
        # Keep SUFFIX

        read = "AAAAATTTTTCCCCC"
        bc = "TTTTT"
        # Match at 6:10
        # 5' Trim: Keep 11:15 ("CCCCC")

        (score, start_pos, end_pos) = semiglobal_alignment(ws, bc, read, 0.0, 0, 1, 1, 1:length(read), 100, 1, 5)
        @test score == 0.0
        @test start_pos == 6
        @test end_pos == 10

        config = DemuxConfig(
            bc_seqs=[bc], bc_lengths_no_N=[5], ids=["id1"],
            trim_side=5,
            ref_search_range=BioDemuX.parse_dynamic_range("1:end"),
            barcode_start_range=BioDemuX.parse_dynamic_range("1:end"),
            barcode_end_range=BioDemuX.parse_dynamic_range("1:end")
        )

        fname, t_start, t_end = determine_filename(read, config, ws)
        trim_range = t_start:t_end
        @test fname == "id1.fastq"
        @test trim_range == 11:15
        @test read[trim_range] == "CCCCC"
    end

    @testset "Tie Breaking (3' Trim)" begin
        # Case where multiple alignments have same score.
        # We want the one starting more to the right.

        # Read: ACGTACGT
        # BC: ACGT
        # Match 1: 1:4
        # Match 2: 5:8
        # Both score 0.

        read = "ACGTACGT"
        bc = "ACGT"

        (score, start_pos, end_pos) = semiglobal_alignment(ws, bc, read, 0.0, 0, 1, 1, 1:length(read), 100, 1, 3)
        @test score == 0.0
        @test start_pos == 5
        @test end_pos == 8

        # Case: AAAA vs AA
        read = "AAAA"
        bc = "AA"
        (score, start_pos, end_pos) = semiglobal_alignment(ws, bc, read, 0.0, 0, 1, 1, 1:length(read), 100, 1, 3)
        @test score == 0.0
        @test start_pos == 3
        @test end_pos == 4
    end

    @testset "Dual Trimming" begin
        # Read: PREFIX + BC1 + MIDDLE + BC2 + SUFFIX
        # BC1: 5' trim (keep after)
        # BC2: 3' trim (keep before)
        # Result: MIDDLE

        # Read: AAAAA TTTTT CCCCC GGGGG TTTTT
        # BC1: TTTTT (at 6:10) -> Trim 5' -> Keep 11:end (CCCCC GGGGG TTTTT)
        # BC2: GGGGG (at 16:20) -> Trim 3' -> Keep 1:15 (AAAAA TTTTT CCCCC)
        # Intersect: 11:15 (CCCCC)

        read = "AAAAATTTTTCCCCCGGGGGTTTTT"
        bc1 = "TTTTT"
        bc2 = "GGGGG"

        config = DemuxConfig(
            bc_seqs=[bc1], bc_lengths_no_N=[5], ids=["id1"],
            is_dual=true,
            bc_seqs2=[bc2], bc_lengths_no_N2=[5], ids2=["id2"],
            trim_side=5, # Trim 5' of BC1
            trim_side2=3, # Trim 3' of BC2
            ref_search_range=BioDemuX.parse_dynamic_range("1:end"),
            barcode_start_range=BioDemuX.parse_dynamic_range("1:end"),
            barcode_end_range=BioDemuX.parse_dynamic_range("1:end"),
            ref_search_range2=BioDemuX.parse_dynamic_range("1:end"),
            barcode_start_range2=BioDemuX.parse_dynamic_range("1:end"),
            barcode_end_range2=BioDemuX.parse_dynamic_range("1:end")
        )

        fname, t_start, t_end = determine_filename(read, config, ws)
        trim_range = t_start:t_end
        @test fname == "id1.id2.fastq"
        @test trim_range == 11:15
        @test read[trim_range] == "CCCCC"
    end

    @testset "Performance / No Trim" begin
        # Ensure no regression (basic check)
        ws_no_trim = SemiGlobalWorkspace(100) # trim=false by default
        @test ws_no_trim.origin === nothing

        read = "AAAAATTTTTCCCCC"
        bc = "TTTTT"

        score = semiglobal_alignment(ws_no_trim, bc, read, 0.0, 0, 1, 1, 1:length(read), 100, 1)
        @test score == 0.0
        @test isa(score, Float64)
    end
end
