using BioDemuX
using Test
using Distributed

function check_output_files(output_dir::String, ideal_dir::String)
    ideal_files = readdir(ideal_dir)

    for ideal_file in ideal_files
        output_file_path = joinpath(output_dir, ideal_file)
        ideal_file_path = joinpath(ideal_dir, ideal_file)

        @test isfile(output_file_path)
        @test isfile(ideal_file_path)

        diff_result = readchomp(`diff $output_file_path $ideal_file_path`)
        @test diff_result == ""
    end
end

@testset "Single Thread Tests" begin
    @testset "Single thread with 1 FASTQ file" begin
        mktempdir() do output_dir
            ideal_dir = "results/demo1_R1"

            files = readdir("FASTQ_files/demo1_R1")
            for file in files
                execute_demultiplexing("FASTQ_files/demo1_R1/$file", "reference_files/demo1.tsv", output_dir)
            end
            check_output_files(output_dir, ideal_dir)
        end
    end

    @testset "Single thread with 2 FASTQ files" begin
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
    # @testset "Single thread with 2 FASTQ files with options" begin
    #     mktempdir() do output_dir
    #         ideal_dir = "results/demo2"

    #         files_R1 = readdir("FASTQ_files/demo2_R1")
    #         files_R2 = readdir("FASTQ_files/demo2_R2")
    #         for (file_R1, file_R2) in zip(files_R1, files_R2)
    #             filename = replace(basename(file_R1), r"\.fastq.gz$" => "")
    #             filename2 = replace(basename(file_R2), r"\.fastq.gz$" => "")
    #             execute_demultiplexing("FASTQ_files/demo2_R1/$file_R1", "FASTQ_files/demo2_R2/$file_R2", "reference_files/demo2.csv", output_dir,
    #                 max_error_rate=0.25, min_delta=0.15, mismatch=1, indel=2, classify_both=true, bc_complement=true, bc_rev=true,
    #                 output_prefix1="test_prefix1." * filename, output_prefix2="test_prefix2." * filename2)
    #         end
    #         check_output_files(output_dir, ideal_dir)
    #     end
    # end
end

addprocs(3)
@everywhere using BioDemuX
@testset "Multi Thread Tests" begin
    @testset "Multi-thread with 1 FASTQ file" begin
        mktempdir() do output_dir
            ideal_dir = "results/demo1_R1"

            files = readdir("FASTQ_files/demo1_R1")
            for file in files
                execute_demultiplexing("FASTQ_files/demo1_R1/$file", "reference_files/demo1.tsv", output_dir)
            end
            check_output_files(output_dir, ideal_dir)
        end
    end

    @testset "Multi-thread with 2 FASTQ files" begin
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
end