"""
Montecarlo simulation of structure-forming systems
"""

include("extern.jl");
include("parameters.jl");
include("functions.jl");
include("analysis.jl");

function main()
    B = initializeBox();
    Results = montecarlo(B);
    println("Atoms:$(B.M[1]+B.M[-1]) Magnetization:$(B.M[1]-B.M[-1]) Molecules:$(length(B.molecules))");
    Results
end

Results = main();
saveResults && CSV.write("results.csv", Results);

final_results = analyze(Results)

# plot(Results.average_energy)
# Results