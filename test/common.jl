using Test
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
