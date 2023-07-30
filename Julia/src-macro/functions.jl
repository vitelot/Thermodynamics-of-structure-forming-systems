function initializeBox()::Box
    Ntot::Int = Opt["Ntot"];
    Tmin::Double = Opt["Tmin"];
    Jcoupling::Double = Opt["Jcoupling"];
    Hfield::Double = Opt["Hfield"];
     
    # zero initial magnetization
    nup::Int = floor(Opt["initialSpinUpFraction"]*Ntot);
    ndn = Ntot-nup;

    Nmol = 0;

    B = Box(
            Ntot,
            Tmin,
            Jcoupling,
            Hfield,
            nup, # random spins
            ndn,
            Ntot, # all atom free 
            Nmol
        );

    return B;
end

"""
    Energy of the system given spin population and parameters
"""
function energy(nup::Int,ndn::Int,J_coupling::Double,H_field::Double)::Double
    n = nup+ndn;
    m = nup-ndn;
    return (n-m*m)*J_coupling/(n-1) - H_field * m;
end

"""
    Energy of the system given the system box
"""
function energy(B::Box)::Double

    return energy(B.nup, B.ndn, B.J, B.H); 
end

"""
Do one montecarlo move and return the delta energy;
Returns 0.0 if the move is not accepted;
The kind of move is proportional to the nr of possible one-move configurations of that kind,
i.e., we minimize the free energy F = E -TS. 
"""
function oneMove!(B::Box)::Double
    fractionMinNrFreeAtoms::Double = Opt["fractionMinNrFreeAtoms"];

    na = atomNumber(B);
    nm = moleculeNumber(B);
    nup = B.nup;
    ndn = B.ndn;

    p::Vector{Double} = [
        # Flip one spin of free particle
        ndn/(nup+1), # Spin flip from down to up
        nup/(ndn+1), # Spin flip from up to down
        # Join two free spins to one molecule
        0.5*nup*(nup-1)/(nm+1), # Join two spins up
        0.5*ndn*(ndn-1)/(nm+1), # Join two spins dn
        0.5*nup*ndn/(nm+1), # Join one spin up and one down
        # Split one molecule to two free particles
        2.0*nm/(nup+2)/(nup+1), # Split two up
        2.0*nm/(ndn+2)/(ndn+1), # Split two dn
        2.0*nm/(nup+1)/(ndn+1) # Split one up and one dn
    ];

    if na <= fractionMinNrFreeAtoms*B.N # may not join atoms to keep the reservoir
        p[3] = p[4] = p[5] = 0.0;
    end

    p = p ./ sum(p);
    for i in 2:length(p)
        p[i] += p[i-1]; # cumul sum
    end
    # println(p);

    r = rand();
    if r < p[1]
        return spinFlipDnToUp!(B);
    elseif r < p[2] 
        return spinFlipUpToDn!(B);
    elseif r < p[3]
        return moleculeJoin2Up!(B);
    elseif r < p[4]
        return moleculeJoin2Dn!(B);   
    elseif r < p[5]
        return moleculeJoinUpDn!(B);
    elseif r < p[6]
        return moleculeSplit2Up!(B);
    elseif r < p[7]
        return moleculeSplit2Dn!(B);
    else
        return moleculeSplit2UpDn!(B);
    end
end

"""
return the magnetisation in the box
"""
function magnetisation(B::Box)::Int
    return B.nup - B.ndn;
end

"""
return the number of molecules in the box
"""
function moleculeNumber(B::Box)::Int
    return B.Nmol;
end

"""
return the number of free atoms in the box
"""
function atomNumber(B::Box)::Int
    return B.Nfree;
end

function multiTemperature(B::Box)::DataFrame
    Tmin::Double  = Opt["Tmin"];
    Tstep::Double = Opt["Tstep"];
    Tmax::Double  = Opt["Tmax"];

    Results  = DataFrame(
                temperature=Double[],
                j_coupling=Double[],
                h_field=Double[],
                step=Int[],
                average_energy=Double[],
                magnetisation=Double[],
                molecules=Int[]
                );
    for temperature = Tmin:Tstep:Tmax
        B.T = temperature;
        @info "Simulating T=$temperature";

        # println(B);
        montecarlo!(B, Results);
        # println(B);
    end

    return Results;
end

# function copy(B::Box)::Box
#     return Box(B.N, B.T, B.J, B.H, copy(B.M), copy(B.atoms), copy(B.molecules));
# end

function multiTemperature_parallel(B::Box)::DataFrame

    Tmin::Double  = Opt["Tmin"];
    Tstep::Double = Opt["Tstep"];
    Tmax::Double  = Opt["Tmax"];
    
    Results  = DataFrame(
                temperature=Double[],
                j_coupling=Double[],
                h_field=Double[],
                step=Int[],
                average_energy=Double[],
                magnetisation=Double[],
                molecules=Int[]
                );
    nthreads = Threads.nthreads();
    if nthreads == 1
        println("You have only one thread selected. If you have more, pls run julia with --threads n");
        println("Going on sequentially: no parallelization\n");
        return multiTemperature(B);
    end

    # VB = Vector{Box}(undef, nthreads);
    VR = Vector{DataFrame}(undef, nthreads);
    for i = 1:nthreads
        # VB[i] = copy(B);
        VR[i] = copy(Results);
    end

    Threads.@threads for temperature = Tmin:Tstep:Tmax
        tid = Threads.threadid();
        Bc = initializeBox();
        Bc.T = temperature;
        @info "Simulating T=$temperature on thread $tid";
        
        # println(Bc);
        montecarlo!(Bc, VR[tid]);
        # println(Bc);
        # Bc = nothing;
    end

    return vcat(VR...);
end

function montecarlo!(B::Box, Results::DataFrame)::Nothing
    thermalisationSteps::Int = Opt["thermalisationSteps"];
    Steps::Int     = Opt["Steps"];
    avrgStep::Int  = Opt["avrgStep"];
    verbosity::Int = Opt["verbosity"];
    onlyAccepted   = Opt["onlyAccepted"];

    wastedmoves = 0;
    for i in 1:thermalisationSteps
        ΔE = oneMove!(B);
        if onlyAccepted
            while ΔE == 0.0
                ΔE = oneMove!(B);
                wastedmoves += 1;
            end
        end
    end
    verbosity >= 1 && @info "$wastedmoves wasted moves during $thermalisationSteps thermalisation steps";
    
    sumEnergy = en = energy(B);
    for i in 1:Steps # do at max Step steps

        ΔE = oneMove!(B);
        if onlyAccepted
            while ΔE == 0.0
                ΔE = oneMove!(B);
            end
        end
        
        en += ΔE;

        sumEnergy += en; 
        # push!(Elist, en);
        if i%avrgStep == 0
            avrg = sumEnergy/i;
            push!(Results, 
                (B.T, B.J, B.H, i, avrg, abs(magnetisation(B)), moleculeNumber(B))
            );
        end
    end
    return;
end

function spinFlipDnToUp!(B::Box)
    nup = B.nup;
    ndn = B.ndn;
    nup_tilde = nup+1;
    ndn_tilde = ndn-1;
    ΔE = energy(nup_tilde,ndn_tilde, B.J,B.H) - energy(nup,ndn, B.J,B.H); 
    if ΔE < 0 || rand()<exp(-ΔE/B.T)
        B.nup = nup_tilde;
        B.ndn = ndn_tilde;
        return ΔE; # spin flip accepted
    end
    return 0.0; # flipping discarded
end

function spinFlipUpToDn!(B::Box)
    nup = B.nup;
    ndn = B.ndn;
    nup_tilde = nup-1;
    ndn_tilde = ndn+1;
    ΔE = energy(nup_tilde,ndn_tilde, B.J,B.H) - energy(nup,ndn, B.J,B.H); 
    if ΔE < 0 || rand()<exp(-ΔE/B.T)
        B.nup = nup_tilde;
        B.ndn = ndn_tilde;
        return ΔE; # spin flip accepted
    end
    return 0.0; # flipping discarded
end

function moleculeJoin2Up!(B::Box)
    
    if B.nup <2 
        @info "No 2 spin up atoms available to join";
        return 0.0;
    end

    nup = B.nup;
    ndn = B.ndn;
    nup_tilde = nup-2;
    ndn_tilde = ndn;
    ΔE = energy(nup_tilde,ndn_tilde, B.J,B.H) - energy(nup,ndn, B.J,B.H); 
    if ΔE < 0 || rand()<exp(-ΔE/B.T)
        B.nup = nup_tilde;
        B.ndn = ndn_tilde;
        B.Nmol += 1;
        B.Nfree -= 2;
        return ΔE; # spin flip accepted
    end
    return 0.0; # flipping discarded

end

function moleculeJoin2Dn!(B::Box)
    
    if B.ndn <2 
        @info "No 2 spin dn atoms available to join";
        return 0.0;
    end

    nup = B.nup;
    ndn = B.ndn;
    nup_tilde = nup;
    ndn_tilde = ndn-2;
    ΔE = energy(nup_tilde,ndn_tilde, B.J,B.H) - energy(nup,ndn, B.J,B.H); 
    if ΔE < 0 || rand()<exp(-ΔE/B.T)
        B.nup = nup_tilde;
        B.ndn = ndn_tilde;
        B.Nmol += 1;
        B.Nfree -= 2;
        return ΔE; # spin flip accepted
    end
    return 0.0; # flipping discarded

end

function moleculeJoinUpDn!(B::Box)
    
    if B.ndn < 1 || B.nup < 1
        @info "No spin dn and up atoms available to join";
        return 0.0;
    end

    nup = B.nup;
    ndn = B.ndn;
    nup_tilde = nup-1;
    ndn_tilde = ndn-1;
    ΔE = energy(nup_tilde,ndn_tilde, B.J,B.H) - energy(nup,ndn, B.J,B.H); 
    if ΔE < 0 || rand()<exp(-ΔE/B.T)
        B.nup = nup_tilde;
        B.ndn = ndn_tilde;
        B.Nmol += 1;
        B.Nfree -= 2;
        return ΔE; # spin flip accepted
    end
    return 0.0; # flipping discarded
end

function moleculeSplit2Up!(B::Box)
    
    if B.Nmol < 1 
        @info "No molecules to split in 2 up atoms";
        return 0.0;
    end

    nup = B.nup;
    ndn = B.ndn;
    nup_tilde = nup+2;
    ndn_tilde = ndn;
    ΔE = energy(nup_tilde,ndn_tilde, B.J,B.H) - energy(nup,ndn, B.J,B.H); 
    if ΔE < 0 || rand()<exp(-ΔE/B.T)
        B.nup = nup_tilde;
        B.ndn = ndn_tilde;
        B.Nmol -= 1;
        B.Nfree += 2;
        return ΔE; # spin flip accepted
    end
    return 0.0; # flipping discarded
end

function moleculeSplit2Dn!(B::Box)
    
    if B.Nmol < 1 
        @info "No molecules to split in 2 up atoms";
        return 0.0;
    end

    nup = B.nup;
    ndn = B.ndn;
    nup_tilde = nup;
    ndn_tilde = ndn+2;
    ΔE = energy(nup_tilde,ndn_tilde, B.J,B.H) - energy(nup,ndn, B.J,B.H); 
    if ΔE < 0 || rand()<exp(-ΔE/B.T)
        B.nup = nup_tilde;
        B.ndn = ndn_tilde;
        B.Nmol -= 1;
        B.Nfree += 2;
        return ΔE; # spin flip accepted
    end
    return 0.0; # flipping discarded
end

function moleculeSplit2UpDn!(B::Box)
    
    if B.Nmol < 1 
        @info "No molecules to split in 2 up atoms";
        return 0.0;
    end

    nup = B.nup;
    ndn = B.ndn;
    nup_tilde = nup+1;
    ndn_tilde = ndn+1;
    ΔE = energy(nup_tilde,ndn_tilde, B.J,B.H) - energy(nup,ndn, B.J,B.H); 
    if ΔE < 0 || rand()<exp(-ΔE/B.T)
        B.nup = nup_tilde;
        B.ndn = ndn_tilde;
        B.Nmol -= 1;
        B.Nfree += 2;
        return ΔE; # spin flip accepted
    end
    return 0.0; # flipping discarded
end
