using BioDemuX
using Test
using DataFrames

@testset "Feature Tests" begin

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
        score = BioDemuX.semiglobal_alignment_N(ws, bc, seq_indel, max_error, match_score, mismatch_score, indel_score, nindel_score, ref_search_range, max_start_pos, min_end_pos, non_N_m)
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

    @testset "Integration with N and Ranges" begin
        mktempdir() do temp_dir
            # Create dummy barcode file (FASTA)
            bc_file = joinpath(temp_dir, "barcodes.fasta")
            open(bc_file, "w") do io
                write(io, ">BC1\nANNC\n") # Matches AC
                write(io, ">BC2\nTTTT\n") # Matches TTTT
            end
            
            # Create dummy FASTQ file
            fastq_file = joinpath(temp_dir, "reads.fastq")
            open(fastq_file, "w") do io
                write(io, "@read1\nATTC\n+\nIIII\n")       # Should match BC1 (ANNC)
                write(io, "@read2\nTTTT\n+\nIIII\n")   # Should match BC2
                write(io, "@read3\nATTG\n+\nIIII\n")       # BC1 mismatch (score 0.5). If max_error=0.6, match.
                write(io, "@read4\nGGGG\n+\nIIII\n")   # No match (Score 1.0)
            end
            
            output_dir = joinpath(temp_dir, "output")
            
            # Run demultiplexing
            BioDemuX.execute_demultiplexing(
                fastq_file, bc_file, output_dir;
                max_error_rate=0.6,
                nindel=1,
                ref_search_range="1:end",
                barcode_start_range="1:end",
                barcode_end_range="1:end"
            )
            
            # Check output files
            # Prefix is "reads"
            @test isfile(joinpath(output_dir, "reads.BC1.fastq"))
            @test isfile(joinpath(output_dir, "reads.BC2.fastq"))
            @test isfile(joinpath(output_dir, "reads.unknown.fastq"))
            
            # Check content counts
            count_lines(f) = count(l -> true, eachline(f))
            
            # BC1: read1 and read3 -> 8 lines
            c1 = count_lines(joinpath(output_dir, "reads.BC1.fastq"))
            @test c1 == 8
            
            # BC2: read2 -> 4 lines
            c2 = count_lines(joinpath(output_dir, "reads.BC2.fastq"))
            @test c2 == 4
            
            # Unknown: read4 -> 4 lines
            @test count_lines(joinpath(output_dir, "reads.unknown.fastq")) == 4
        end
    end


    @testset "Range Restrictions Integration" begin
        mktempdir() do temp_dir
            # Create dummy barcode file
            bc_file = joinpath(temp_dir, "barcodes.fasta")
            open(bc_file, "w") do io
                write(io, ">BC1\nAAAA\n")
            end
            
            # Create dummy FASTQ file
            fastq_file = joinpath(temp_dir, "reads.fastq")
            open(fastq_file, "w") do io
                write(io, "@read1_start\nAAAATTTT\n+\nIIIIIIII\n") # Barcode at 1:4
                write(io, "@read2_end\nTTTTAAAA\n+\nIIIIIIII\n")   # Barcode at 5:8
                write(io, "@read3_mid\nTTAAAATT\n+\nIIIIIIII\n")   # Barcode at 3:6
            end
            
            output_dir_1 = joinpath(temp_dir, "output_1")
            
            # Test 1: Search only first 4 bases
            BioDemuX.execute_demultiplexing(
                fastq_file, bc_file, output_dir_1;
                max_error_rate=0.0,
                nindel=1,
                ref_search_range="1:4",
                barcode_start_range="1:end",
                barcode_end_range="1:end"
            )
            
            @test isfile(joinpath(output_dir_1, "reads.BC1.fastq"))
            @test isfile(joinpath(output_dir_1, "reads.unknown.fastq"))
            
            # Check counts
            count_lines(f) = count(l -> true, eachline(f))
            # read1_start should match (4 lines)
            @test count_lines(joinpath(output_dir_1, "reads.BC1.fastq")) == 4
            # read2_end and read3_mid should be unknown (8 lines)
            @test count_lines(joinpath(output_dir_1, "reads.unknown.fastq")) == 8
            
            output_dir_2 = joinpath(temp_dir, "output_2")
            
            # Test 2: Search only last 4 bases (5:8)
            # Note: ref_search_range is relative to read length.
            # "5:8" or "end-3:end"
            BioDemuX.execute_demultiplexing(
                fastq_file, bc_file, output_dir_2;
                max_error_rate=0.0,
                nindel=1,
                ref_search_range="5:8",
                barcode_start_range="1:end",
                barcode_end_range="1:end"
            )
            
            # read2_end should match
            @test count_lines(joinpath(output_dir_2, "reads.BC1.fastq")) == 4
            # read1_start and read3_mid should be unknown
            @test count_lines(joinpath(output_dir_2, "reads.unknown.fastq")) == 8

             output_dir_3 = joinpath(temp_dir, "output_3")
            
            # Test 3: Barcode must start at 1
            BioDemuX.execute_demultiplexing(
                fastq_file, bc_file, output_dir_3;
                max_error_rate=0.0,
                nindel=1,
                ref_search_range="1:end",
                barcode_start_range="1:1",
                barcode_end_range="1:end"
            )
            
            # read1_start should match (starts at 1)
            @test count_lines(joinpath(output_dir_3, "reads.BC1.fastq")) == 4
             # read2_end (starts at 5) and read3_mid (starts at 3) should be unknown
            @test count_lines(joinpath(output_dir_3, "reads.unknown.fastq")) == 8
        end
    end
end

