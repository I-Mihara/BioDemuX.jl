@testset "Integration Tests" begin
    @testset "1 FASTQ file (plain text)" begin
        mktempdir() do output_dir
            ideal_dir = "results/demo1_R1"
            files = readdir("FASTQ_files/demo1_R1")
            for file in files
                execute_demultiplexing("FASTQ_files/demo1_R1/$file", "reference_files/demo1.tsv", output_dir,)
            end
            check_output_files(output_dir, ideal_dir)
        end
    end

    @testset "2 FASTQ files (paired, plain text)" begin
        mktempdir() do output_dir
            ideal_dir = "results/demo1_R2"
            files_R1 = readdir("FASTQ_files/demo1_R1")
            files_R2 = readdir("FASTQ_files/demo1_R2")
            for (file_R1, file_R2) in zip(files_R1, files_R2)
                execute_demultiplexing("FASTQ_files/demo1_R1/$file_R1", "FASTQ_files/demo1_R2/$file_R2", "reference_files/demo1.tsv", output_dir)
            end
            check_output_files(output_dir, ideal_dir)
        end
    end

    @testset "2 FASTQ file (paired, gzipped, with options)" begin
        mktempdir() do output_dir
            ideal_dir = "results/demo2"
            files_R1 = readdir("FASTQ_files/demo2_R1")
            files_R2 = readdir("FASTQ_files/demo2_R2")
            for (file_R1, file_R2) in zip(files_R1, files_R2)
                filename_R1 = replace(basename(file_R1), r"\.fastq\.gz$" => "")
                filename_R2 = replace(basename(file_R2), r"\.fastq\.gz$" => "")
                execute_demultiplexing("FASTQ_files/demo2_R1/$file_R1",
                                       "FASTQ_files/demo2_R2/$file_R2",
                                       "reference_files/demo2.csv", output_dir;
                                       max_error_rate=0.25, min_delta=0.15,
                                       mismatch=1, indel=2,
                                       classify_both=true, bc_complement=true, bc_rev=true,
                                       output_prefix1="test_prefix1.$filename_R1",
                                       output_prefix2="test_prefix2.$filename_R2",
                                       gzip_output=false)
            end
            check_output_files(output_dir, ideal_dir)
        end
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
