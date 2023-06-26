"""
Montecarlo simulation of structure-forming systems
"""

include("extern.jl");
include("parameters.jl");
include("functions.jl");
include("analysis.jl");

function main()
    B = initializeBox();
    Results = multiTemperature(B);
    # println("Atoms:$(B.M[1]+B.M[-1]) Magnetization:$(B.M[1]-B.M[-1]) Molecules:$(length(B.molecules))");
    
    saveResults && CSV.write("results.csv", Results);
    
    final_results = analyze(Results);
    saveResults && CSV.write("final_results.csv", final_results);
    
    if doPlot
        if length(unique(final_results.temperature))<2
            println("Less than 2 temperature runs found. No plot.");
            println("$final_results");
        else
            p=plotall(final_results)
            savefig(p,"results.pdf");
        end    
    end

    final_results
end

Results = main();
