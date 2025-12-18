using PackageCompiler

create_app(".", "biodemux";
    force=true,
    precompile_execution_file="test/runtests.jl",
    incremental=false,
    filter_stdlibs=false
)
