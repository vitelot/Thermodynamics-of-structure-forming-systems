function initializeBox(N::Int)
    Atoms = Atom[];
    for i=1:N
        push!(Atoms, Atom(i,rand([-1,1])));
    end
    nup = count(x->x.σ>0, Atoms);
    ndn = count(x->x.σ<0, Atoms);

    B = Box(Dict(-1=>ndn, 1=>nup),Set(Atoms), Set{Molecule}());
    
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

    return energy(B.M[1], B.M[-1], J, H); 

    # Es = Eh = 0.0;
    # Atoms = collect(B.atoms);
    # n = length(Atoms);

    # if H > 0.0 # well, it doesn't really matter ... it's so fast anyway
    #     for i in eachindex(Atoms)
    #         si = Atoms[i].σ;
    #         # Eh += si;
    #         for j in 1:i-1
    #             Es += si * Atoms[j].σ;
    #         end
    #     end
    # else
    #     for i in eachindex(Atoms)
    #         si = Atoms[i].σ;
    #         Eh += si;
    #         for j in 1:i-1
    #             Es += si * Atoms[j].σ;
    #         end
    #     end
    # end

    # return -Eh * H - Es * 2.0 * J /(n-1);
end

function spinFlipDeltaEnergy(B::Box, atom::Atom)
    # E = 0.0;
    # n = length(B.atoms);
    # for a in B.atoms
    #     E += a.σ;
    # end 
    # E -= atom.σ;

    # the following is way faster
    nup = B.M[1];
    ndn = B.M[-1];
    nup_tilde = nup - (atom.σ + 1)÷2 + (1 - atom.σ)÷2;
    ndn_tilde = ndn - (1 - atom.σ)÷2 + (atom.σ + 1)÷2;
    
    return energy(nup_tilde,ndn_tilde, J,H) - energy(nup,ndn, J,H); 
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

function montecarlo(B::Box)::Vector{Double}
    # Elist = Double[];
    Avrg  = Double[];

    en = energy(B);
    sumEnergy = en;

    for _ in 1:termalisationSteps
        en += oneMove(B);
    end

    for i in 1:Steps # do at max Step steps
        
        en += oneMove(B);
        
        # energy_check = energy(B);
        # energy_check ≈ en || @info("Energy does not match");

        sumEnergy += en; 
        # push!(Elist, en);
        if i%avrgStep == 0
            avrg = sumEnergy/i;
            # avrg = sum(Elist[end-avrgStep+1:end])/avrgStep;
            push!(Avrg, avrg);
        end
    end
    return Avrg;
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
        if p<exp(-ΔE/T)
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
    ΔE = energy(nup_tilde,ndn_tilde, J,H) - energy(nup,ndn, J,H);

    prob = factorialRatio(nup,nup_tilde)*
           factorialRatio(ndn,ndn_tilde)/
           (2*nmol+2)*
           exp(-ΔE/T);

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
    ΔE = energy(nup_tilde,ndn_tilde, J,H) - energy(nup,ndn, J,H);
    # println("Delta: ",ΔE)

    prob = factorialRatio(nup,nup_tilde)*
           factorialRatio(ndn,ndn_tilde)*
           2*nmol*
           exp(-ΔE/T);
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