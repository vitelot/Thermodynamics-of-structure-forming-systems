Real = Float64;

# parameters
N::Int = 20; # total number of atoms

# structures
mutable struct Atom
    id::Int; # identification number
    σ::Int;  # spin (±1, for the moment)
    isfree::Bool; # true if this atom is not combined
end

mutable struct Molecule
    natoms::Int; # number of atoms composing the molecule
    atomlist::Vector{Int}; # id of the atoms composing it
end

mutable struct Box # the simulation box
    atoms::Vector{Atom};
    molecules::Vector{Molecule};
end

