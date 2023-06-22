using Plots;

Real = Float64;

# parameters
N::Int = 200; # total number of atoms
T::Real = 1.0;  # temperature
H::Real = 0.0; # external magnetic field
J::Real = 1.0; # spin-spin coupling
Steps::Int = 2_000_000;
avrgStep = 1000; # used to estimate the average energy

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

