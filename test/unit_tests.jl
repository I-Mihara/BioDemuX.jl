using BioDemuX
using Test

@testset "Unit Tests" begin
    @testset "N Support in Alignment" begin
        # Mock workspace
        ws = BioDemuX.SemiGlobalWorkspace(100)
        
        # Barcode with Ns: ANNC (length 4, non-N length 2)
        bc = "ANNC"
        non_N_m = 2
        
        # Config
        max_error = 0.5
        match_score = 0
        mismatch_score = 1
        indel_score = 1
        nindel_score = 1 # Must be > 0
        ref_search_range = 1:4
        max_start_pos = 1
        min_end_pos = 4
        
        # 1. Perfect match (ignoring Ns)
        # Query: ANNC (barcode)
        # Read: ATTC (matches ANNC where N matches T)
        seq_match = "ATTC"
        score = BioDemuX.semiglobal_alignment_N(ws, bc, seq_match, max_error, match_score, mismatch_score, indel_score, nindel_score, ref_search_range, max_start_pos, min_end_pos, non_N_m)
        @test score == 0.0
        
        # 2. Mismatch
        # Query: ANNC (barcode)
        # Read: ATTG (matches ATT, mismatch G vs C at pos 4)
        # Cost: 1 mismatch. Score: 1 / 2 = 0.5
        seq_mismatch = "ATTG"
        score = BioDemuX.semiglobal_alignment_N(ws, bc, seq_mismatch, max_error, match_score, mismatch_score, indel_score, nindel_score, ref_search_range, max_start_pos, min_end_pos, non_N_m)
        @test score == 0.5
        
        # 3. Indel
        # Query: ANNC (barcode)
        # Read: ATT (deletion of C)
        # Cost: 1 indel. Score: 1 / 2 = 0.5
        seq_indel = "ATT"
        # min_end_pos must be <= read length (3) for this to match
        # ref_search_range must be within 1:3
        score = BioDemuX.semiglobal_alignment_N(ws, bc, seq_indel, max_error, match_score, mismatch_score, indel_score, nindel_score, 1:3, max_start_pos, 3, non_N_m)
        @test score == 0.5
    end

    @testset "Position Restriction Logic" begin
        # Test parse_dynamic_range
        dr1 = BioDemuX.parse_dynamic_range("1:10")
        @test dr1.start_offset == 1
        @test dr1.start_from_end == false
        @test dr1.end_offset == 10
        @test dr1.end_from_end == false
        
        dr2 = BioDemuX.parse_dynamic_range("1:end")
        @test dr2.start_offset == 1
        @test dr2.start_from_end == false
        @test dr2.end_offset == 0
        @test dr2.end_from_end == true
        
        dr3 = BioDemuX.parse_dynamic_range("end-5:end")
        @test dr3.start_offset == -5
        @test dr3.start_from_end == true
        @test dr3.end_offset == 0
        @test dr3.end_from_end == true
        
        # Test resolve
        len = 100
        r1 = BioDemuX.resolve(dr1, len)
        @test r1 == 1:10
        
        r2 = BioDemuX.resolve(dr2, len)
        @test r2 == 1:100
        
        r3 = BioDemuX.resolve(dr3, len)
        @test r3 == 95:100
    end
end
