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
    Results
end

Results = main();
saveResults && CSV.write("results.csv", Results);

final_results = analyze(Results)
saveResults && CSV.write("final_results.csv", final_results);

plotall(final_results)
savefig("results.pdf")
# Results