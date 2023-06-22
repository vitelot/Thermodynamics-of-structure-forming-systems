"""
Montecarlo simulation of structure-forming systems
"""

include("extern.jl");
include("functions.jl");

function main()
    B = initializeBox(N);
    montecarlo(B);
end

Elist = main();
plot(Elist)
