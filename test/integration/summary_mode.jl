using BioDemuX
using Test
using Dates

@testset "Summary Mode Tests" begin
    mktempdir() do temp_dir
        bc_file = joinpath(temp_dir, "barcodes.fasta")
        open(bc_file, "w") do io
            write(io, ">BC1\nAAAA\n")
            write(io, ">BC2\nTTTT\n")
        end

        bc_file2 = joinpath(temp_dir, "barcodes2.fasta")
        open(bc_file2, "w") do io
            write(io, ">BC2_1\nCCCC\n")
            write(io, ">BC2_2\nGGGG\n")
        end

        fastq_file = joinpath(temp_dir, "reads.fastq")
        open(fastq_file, "w") do io
            write(io, "@read1\nAAAA\n+\nIIII\n") # Match BC1
            write(io, "@read2\nTTTT\n+\nIIII\n") # Match BC2
            write(io, "@read3\nGGGG\n+\nIIII\n") # Unknown
            write(io, "@read4\nAAAT\n+\nIIII\n") # Ambiguous/Mismatch
        end
        # Create dummy FASTQ files (Dual - Inline)
        # BioDemuX supports inline dual barcodes (both in R1).
        fastq_file1 = joinpath(temp_dir, "reads_R1.fastq")
        fastq_file2 = joinpath(temp_dir, "reads_R2.fastq")
        open(fastq_file1, "w") do io
            write(io, "@read1\nAAAACCCC\n+\nIIIIIIII\n") # Match BC1 (AAAA) + BC2_1 (CCCC)
            write(io, "@read2\nTTTTGGGG\n+\nIIIIIIII\n") # Match BC2 (TTTT) + BC2_2 (GGGG)
        end
        open(fastq_file2, "w") do io
            write(io, "@read1\nNNNN\n+\nIIII\n") # Dummy R2
            write(io, "@read2\nNNNN\n+\nIIII\n") # Dummy R2
        end

        output_dir = joinpath(temp_dir, "output")
        mkpath(output_dir)

        BioDemuX.execute_demultiplexing(
            fastq_file, bc_file, output_dir;
            max_error_rate=0.2,
            summary=true,
            summary_format=:txt
        )

        txt_file = joinpath(output_dir, "summary.txt")
        @test isfile(txt_file)
        content = read(txt_file, String)
        @test occursin("Total Reads: 4", content)
        @test occursin("Matched Reads: 2", content)
        @test occursin("Run Information:", content)
        @test occursin("Barcode File: $bc_file", content)
        @test !occursin("Barcode File 2:", content)

        BioDemuX.execute_demultiplexing(
            fastq_file, bc_file, output_dir;
            max_error_rate=0.2,
            summary=true,
            summary_format=:txt
        )
        content = read(txt_file, String)
        @test occursin("==================================================", content)
        @test count(x -> x == "Run Information:", split(content, "\n")) == 2

        BioDemuX.execute_demultiplexing(
            fastq_file1, fastq_file2, bc_file, output_dir;
            barcode_file2=bc_file2,
            max_error_rate=0.2,
            summary=true,
            summary_format=:json
        )

        json_file = joinpath(output_dir, "summary.json")
        @test isfile(json_file)
        content = read(json_file, String)
        @test occursin("\"total_reads\": 2", content)
        @test occursin("\"matched_reads\": 2", content)
        @test occursin("\"barcode_file2\": \"$bc_file2\"", content)

        BioDemuX.execute_demultiplexing(
            fastq_file1, fastq_file2, bc_file, output_dir;
            barcode_file2=bc_file2,
            max_error_rate=0.2,
            summary=true,
            summary_format=:json
        )
        content = read(json_file, String)
        content = strip(content)
        @test startswith(content, "[")
        @test endswith(content, "]")
        @test count(x -> occursin("\"run_info\"", x), split(content, "\n")) == 2

        BioDemuX.execute_demultiplexing(
            fastq_file1, fastq_file2, bc_file, output_dir;
            barcode_file2=bc_file2,
            max_error_rate=0.2,
            summary=true,
            summary_format=:html
        )

        html_file = joinpath(output_dir, "summary.html")
        @test isfile(html_file)
        content = read(html_file, String)
        @test occursin("Total Reads", content)
        @test occursin("Barcode File 2:", content)
        @test occursin(bc_file2, content)

        BioDemuX.execute_demultiplexing(
            fastq_file1, fastq_file2, bc_file, output_dir;
            barcode_file2=bc_file2,
            max_error_rate=0.2,
            summary=true,
            summary_format=:html
        )
        content = read(html_file, String)
        @test count(x -> occursin("Run Information", x), split(content, "\n")) == 2

        original_stdout = stdout
        (rd, wr) = redirect_stdout()

        try
            BioDemuX.execute_demultiplexing(
                fastq_file1, fastq_file2, bc_file, output_dir;
                barcode_file2=bc_file2,
                max_error_rate=0.2,
                summary=true,
                summary_format=:stdout
            )
        finally
            redirect_stdout(original_stdout)
            close(wr)
        end

        output = read(rd, String)
        @test occursin("BioDemuX Summary Report", output)
        @test occursin("Barcode File 2: $bc_file2", output)
        @test occursin("Total Reads: 2", output)
    end
end
