"""
Montecarlo simulation of structure-forming systems
"""

include("extern.jl");
include("functions.jl");

function main()
    B = initializeBox(N);
    Elist = montecarlo(B);
    println("Atoms:$(B.M[1]+B.M[-1]) Magnetization:$(B.M[1]-B.M[-1]) Molecules:$(length(B.molecules))");
    Elist
end

Elist = main();
plot(Elist)
