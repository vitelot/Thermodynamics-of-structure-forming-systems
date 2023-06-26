function initializeBox()::Box
    N = Ntot;
    Atoms = Atom[];
    for i=1:N
        push!(Atoms, Atom(i,rand([-1,1])));
    end
    nup = count(x->x.σ>0, Atoms);
    ndn = count(x->x.σ<0, Atoms);

    B = Box(
            Ntot,
            Tmin,
            Jcoupling,
            Hfield,
            Dict(-1=>ndn, 1=>nup),
            Set(Atoms), 
            Set{Molecule}()
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

    return energy(B.M[1], B.M[-1], B.J, B.H); 
end

function spinFlipDeltaEnergy(B::Box, atom::Atom)

    nup = B.M[1];
    ndn = B.M[-1];
    nup_tilde = nup - (atom.σ + 1)÷2 + (1 - atom.σ)÷2;
    ndn_tilde = ndn - (1 - atom.σ)÷2 + (atom.σ + 1)÷2;
    
    return energy(nup_tilde,ndn_tilde, B.J,B.H) - energy(nup,ndn, B.J,B.H); 
end

"""
Do one montecarlo move and return the delta energy
"""
function oneMove(B::Box)::Double

    if rand() < spinFlipProbability
        return spinFlipMove(B);
    end

    return moleculeMove(B);
end

"""
return the magnetisation in the box
"""
function magnetisation(B::Box)::Int
    return B.M[1] - B.M[-1];
end

"""
return the number of molecules in the box
"""
function moleculeNumber(B::Box)::Int
    return length(B.molecules);
end

function multiTemperature(B::Box)::DataFrame
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
    
    en = energy(B);
    
    for _ in 1:termalisationSteps
        en += oneMove(B);
    end
    
    sumEnergy = en;
    for i in 1:Steps # do at max Step steps
        
        en += oneMove(B);

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

"""
    checks whether a spin flip occurs; does it and returns the ΔE 
"""
function spinFlipMove(B::Box)
    a = rand(B.atoms);
    ΔE = spinFlipDeltaEnergy(B,a);
    if ΔE < 0
        B.M[a.σ] -= 1; # decrease number of old spins
        a.σ = -a.σ;
        B.M[a.σ] += 1; # and increase the new one
        return ΔE; # spin flip accepted
    else
        p = rand();
        if p<exp(-ΔE/B.T)
            B.M[a.σ] -= 1;
            a.σ = -a.σ;
            B.M[a.σ] += 1;
            return ΔE; # spin flip accepted
        end
    end
    return 0.0; # flipping discarded
end

function moleculeMove(B::Box)::Double
    if rand() < moleculeSplitProbability
        return moleculeSplit(B);
    end
    
    return moleculeJoin(B);
end

"""
    checks whether two atoms join; does it and returns the ΔE 
"""
function moleculeJoin(B::Box)::Double
    Atoms = B.atoms;
    if length(Atoms) < 2
        # @info "No atoms available to join";
        return 0.0;
    end

    # select two random different atoms
    a1 = rand(Atoms);
    a2 = rand(Atoms);
    while true
        a2.id == a1.id || break;
        a2 = rand(Atoms);
    end
    
    # current spin population
    nup = B.M[1];
    ndn = B.M[-1];
    # spin population after removal
    nup_tilde = nup - (a1.σ + 1)÷2 - (a2.σ + 1)÷2;
    ndn_tilde = ndn - (1 - a1.σ)÷2 - (1 - a2.σ)÷2;
    # current nr of molecules
    nmol = length(B.molecules);

    # energy is now fast to calculate
    ΔE = energy(nup_tilde,ndn_tilde, B.J,B.H) - energy(nup,ndn, B.J,B.H);

    prob = factorialRatio(nup,nup_tilde)*
           factorialRatio(ndn,ndn_tilde)/
           (2*nmol+2)*
           exp(-ΔE/B.T);

    if rand() < prob
        # accept the joining
        # remove the two atoms from box1
        pop!(Atoms, a1);
        pop!(Atoms, a2);
        mol = Molecule([a1,a2]);
        # add the new molecule into box2
        push!(B.molecules, mol);
        B.M[1] = nup_tilde;
        B.M[-1]= ndn_tilde;
        return ΔE;
    else
        # do not join and do nothing
        return 0.0;
    end
end

"""
    checks whether to split a molecule; 
    if yes, does it, returns the ΔE 
    and assign atoms' spin proportional to spin population 
"""
function moleculeSplit(B)::Double
    nmol = length(B.molecules);
    nmol < 1 && return 0.0;
    # current magnetisation in box1
    nup = B.M[1];
    ndn = B.M[-1];

    ntot = nup+ndn;
    pup = ifelse(ntot>0, nup/ntot, 0.5);
    
    newspin1 = ifelse(rand()<pup, 1, -1);
    newspin2 = ifelse(rand()<pup, 1, -1);
    # vito: actually we know the spin of composing atoms and could use it

    # spin population after addition
    nup_tilde = nup + (newspin1 + 1)÷2 + (newspin2 + 1)÷2;
    ndn_tilde = ndn + (1 - newspin1)÷2 + (1 - newspin2)÷2;

    # energy is now fast to calculate
    ΔE = energy(nup_tilde,ndn_tilde, B.J,B.H) - energy(nup,ndn, B.J,B.H);
    # println("Delta: ",ΔE)

    prob = factorialRatio(nup,nup_tilde)*
           factorialRatio(ndn,ndn_tilde)*
           2*nmol*
           exp(-ΔE/B.T);
    # println("Prob: $prob")

    if rand()<prob
        # accepted
        # println("accepted")
        # select a random molecule and remove it from box2 
        mol = pop!(B.molecules);

        # vito: we stick to the idea of having two atoms only for the moment
        (a1,a2) = mol.atoms;
        a1.σ = newspin1;
        a2.σ = newspin2;
        push!(B.atoms, a1, a2);
        B.M[1] = nup_tilde;
        B.M[-1]= ndn_tilde;
        return ΔE;
    end

    # not accepted
    return 0.0;
end

"""
calculates the ratio n!/m!
"""
function factorialRatio(n::Int, m::Int)::Double
    n==m && return 1.0;
    if n==0 || n==1 
        return 1.0/factorial(m);
    end
    if m==0 || m==1 
        return factorial(n);
    end
    f=1.0;
    if n>m
        for i=m+1:n
            f *= i;
        end
        return f;
    else
        for i = n+1:m
            f /= i
        end
        return f;
    end
end