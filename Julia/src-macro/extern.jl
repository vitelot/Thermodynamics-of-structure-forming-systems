"""
 Structures are defined
"""

using CSV, DataFrames, Plots, Statistics;
# import Base.copy;

###################################
###       Global variables      ###
ProgramVersion = v"0.1.1";
Opt = Dict{String,Any}(); # options read from par.ini

Double = Float64; # I hate writing Float64, I love C
###################################

# structures

mutable struct Box # the simulation box
    ### Static params ###
    N::Int; # total nr of atoms in the box
    T::Double; # temperature
    J::Double; # spin-spin coupling
    H::Double; # external magnetic field
    molEF::Double; # molecule energy formation
    
    ### Dynamical params ###
    nup::Int; # number of up atoms
    ndn::Int; # number of dn atoms
    Nfree::Int; # number of free atoms
    Nmol::Int; # number of molecules
end

