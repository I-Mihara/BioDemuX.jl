using BioDemuX
using Test
using CodecZlib


# Change working directory to the directory of this file (test/)
cd(@__DIR__)

include("common.jl")

@testset "Unit Tests" begin
    include("unit/alignment.jl")
    include("unit/trimming.jl")
    include("unit/hamming.jl")
    include("unit/exact.jl")
end

println("Running tests with $(Threads.nthreads()) threads")

@testset "Integration Tests (Current Process)" begin
    include("integration/single_barcode.jl")
    include("integration/dual_barcode.jl")
    include("integration/summary_mode.jl")

    include("integration/summary_distributions.jl")
    include("integration/hamming_demux.jl")
    include("integration/exact_demux.jl")
end

if Threads.nthreads() == 1
    println("\nRunning multi-threaded tests in a subprocess...")
    @testset "Integration Tests (Subprocess: 4 threads)" begin
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) -t 4 run_subprocess.jl`
        @test success(run(cmd))
    end
else
    println("\nSkipping subprocess tests (already running with multiple threads)")
end