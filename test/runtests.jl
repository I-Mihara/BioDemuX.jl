using BioDemuX
using Test
using CodecZlib

include("common.jl")
include("unit_tests.jl")

println("Running tests with $(Threads.nthreads()) threads")

@testset "Integration Tests (Current Process)" begin
    include("integration_tests.jl")
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