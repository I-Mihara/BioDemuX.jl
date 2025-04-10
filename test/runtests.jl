using BioDemuX
using Test
using Distributed
using CodecZlib

function check_output_files(output_dir::String, ideal_dir::String)
    for ideal_file in readdir(ideal_dir)
        output_file_path = joinpath(output_dir, ideal_file)
        ideal_file_path = joinpath(ideal_dir, ideal_file)
        @test isfile(output_file_path)
        @test isfile(ideal_file_path)
        
        if endswith(lowercase(ideal_file), ".gz")
            output_content = read(GzipDecompressorStream(open(output_file_path)), String)
            ideal_content = read(GzipDecompressorStream(open(ideal_file_path)), String)
        else
            output_content = read(output_file_path, String)
            ideal_content = read(ideal_file_path, String)
        end
        
        @test output_content == ideal_content
    end
end

@testset "Single Thread Tests" begin
    @testset "Single thread with 1 FASTQ file (plain text)" begin
        mktempdir() do output_dir
            ideal_dir = "results/demo1_R1"
            files = readdir("FASTQ_files/demo1_R1")
            for file in files
                execute_demultiplexing("FASTQ_files/demo1_R1/$file", "reference_files/demo1.tsv", output_dir,)
            end
            check_output_files(output_dir, ideal_dir)
        end
    end

    @testset "Single thread with 2 FASTQ files (paired, plain text)" begin
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

    @testset "Single thread with 2 FASTQ file (paired, gzipped, with options)" begin
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
end

addprocs(3)
@everywhere using BioDemuX
@testset "Multi Thread Tests" begin
    @testset "Multi-thread with 1 FASTQ file (plain text)" begin
        mktempdir() do output_dir
            ideal_dir = "results/demo1_R1"
            files = readdir("FASTQ_files/demo1_R1")
            for file in files
                execute_demultiplexing("FASTQ_files/demo1_R1/$file", "reference_files/demo1.tsv", output_dir)
            end
            check_output_files(output_dir, ideal_dir)
        end
    end

    @testset "Multi-thread with 2 FASTQ files (paired, plain text)" begin
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

    @testset "Multi-thread with 2 FASTQ files (paired, gzip, with options)" begin
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
end