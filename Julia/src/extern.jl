using Plots;

Double = Float64;

DEBUG = 3;

# parameters
N::Int = 200; # total number of atoms
T::Double = 1.0;  # temperature
H::Double = 0.0; # external magnetic field
J::Double = 1.0; # spin-spin coupling
Steps::Int = 2_000_000;
avrgStep = 1000; # used to estimate the average energy
splitProbability = 0.5; # decides if to split molecules or join atoms
spinFlipProbability = 0.5; # decides if to flip atoms or make a molecule move

# structures
mutable struct Atom
    id::Int; # identification number
    σ::Int;  # spin (±1, for the moment)
    # isfree::Bool; # true if this atom is not combined
end

mutable struct Molecule
    # natoms::Int; # number of atoms composing the molecule
    atoms::Vector{Atom}; # the atoms composing the molecule
end

mutable struct Box # the simulation box
    M::Dict{Int,Int}; # magnetisation: number of up and down spin
    atoms::Set{Atom};
    molecules::Set{Molecule};
end

