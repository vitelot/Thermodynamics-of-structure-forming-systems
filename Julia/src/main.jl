"""
Montecarlo simulation of structure-forming systems
"""

include("extern.jl");
include("parameters.jl");
include("functions.jl");
include("analysis.jl");

function main()
    B = initializeBox();
    if goParallel
        Results = multiTemperature_parallel(B);
    else
        Results = multiTemperature(B);
    end
    # println("Atoms:$(B.M[1]+B.M[-1]) Magnetization:$(B.M[1]-B.M[-1]) Molecules:$(length(B.molecules))");
    filetext = "_N=$(B.N)_Tmin=$(Tmin)_Tstep=$(Tstep)_Tmax=$(Tmax)_fMA=$(fractionMinNrFreeAtoms)";
    saveResults && CSV.write("results$filetext.csv", Results);
    
    final_results = analyze(Results);
    saveResults && CSV.write("final_results$filetext.csv", final_results);
    
    if doPlot
        if length(unique(final_results.temperature))<2
            println("Less than 2 temperature runs found. No plot.");
            println("$final_results");
        else
            p=plotall(final_results)
            savefig(p,"results$filetext.pdf");
        end    
    end

    final_results
end

Results = main();
