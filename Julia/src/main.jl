"""
Montecarlo simulation of structure-forming systems
"""

include("extern.jl");
include("parameters.jl");
include("functions.jl");

function main()
    B = initializeBox(N);
    Results = montecarlo(B);
    println("Atoms:$(B.M[1]+B.M[-1]) Magnetization:$(B.M[1]-B.M[-1]) Molecules:$(length(B.molecules))");
    Results
end

Results = main();
plot(Results.average_energy)
Results