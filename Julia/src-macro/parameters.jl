DEBUG = 3;

# parameters
# Ntot::Int         = 200; # total number of atoms
# Hfield::Double    = 0.0; # external magnetic field
# Jcoupling::Double = 0.5; # spin-spin coupling

# Tmin::Double  = 0.0;
# Tstep::Double = 0.002;
# Tmax::Double  = 0.5;

# termalisationSteps::Int = Ntot*Ntot÷2; # run these nr of steps at the beginning discarding the energy
# Steps::Int              = 100*termalisationSteps; # total nr of steps
# avrgStep::Int           = floor(Int, sqrt(Steps)); # used to estimate the average energy

# # moleculeSplitProbability::Double = 0.5; # decides if to split molecules or join atoms
# # spinFlipProbability::Double      = 0.5; # decides if to flip atoms or make a molecule move
# fractionMinNrFreeAtoms = 0.2; # do not accept less than this fraction of free atoms

# saveResults::Bool = true; # saves results into file results.CSV
# doPlot::Bool      = true; # plot and save it on pdf
# goParallel::Bool  = true; # well, it does it, multi-thread

function loadParameters(file::String)
    
    if !isfile(file)
        createIniFile(file)
    end

    for line in eachline(file)
        occursin(r"^#", line) && continue # ignore lines beginning with #
        df = split(line, r"\s+")
        length(df) < 2 && continue # ignore empty lines
        key = df[1] ; val = df[2]
        ####################################################################
        # if(key=="Version")                  Opt[key] = val
        ####################################################################
        if(key=="Ntot")                       Opt[key] = parse(Int, val)
        elseif(key=="Hfield")                 Opt[key] = parse(Double, val)
        elseif(key=="Jcoupling")              Opt[key] = parse(Double, val)
        elseif(key=="initialSpinUpFraction")  Opt[key] = parse(Double, val)
        elseif(key=="fractionMinNrFreeAtoms") Opt[key] = parse(Double, val)
        ####################################################################
        elseif(key=="termalisationSteps")   Opt[key] = parse(Int, val)
        elseif(key=="Steps")                Opt[key] = parse(Int, val)
        elseif(key=="avrgStep")            Opt[key] = parse(Int, val)
        ####################################################################
        elseif(key=="Tmin")                 Opt[key] = parse(Double, val)
        elseif(key=="Tstep")                Opt[key] = parse(Double, val)
        elseif(key=="Tmax")                 Opt[key] = parse(Double, val)
        ####################################################################
        elseif(key=="saveResults")   Opt[key] = parse(Bool, val)
        elseif(key=="doPlot")        Opt[key] = parse(Bool, val)
        elseif(key=="goParallel")    Opt[key] = parse(Bool, val)
        ####################################################################
        else @warn "WARNING: input parameter $key does not exist";
        end
    end
    ts = Opt["termalisationSteps"];
    Opt["termalisationSteps"] = ifelse(ts>0, ts, Opt["Ntot"]*Opt["Ntot"]÷2); # run these nr of steps at the beginning discarding the energy

    st = Opt["Steps"];
    Opt["Steps"]              = ifelse(st>0, st, 100*Opt["termalisationSteps"]); # total nr of steps

    av = Opt["avrgStep"];
    Opt["avrgStep"]           = ifelse(av>0, av, floor(Int, sqrt(Opt["Steps"]))); # used to estimate the average energy

    @info("Parameters loaded, starting the program.")
end

function createIniFile(file::String)

    INI = open(file, "w")
        print(INI,
"""
#key                    value
#############################
Ntot        200     # total nr of particles
Hfield      0.0     # external magnetic field
Jcoupling   1.0     # spin-spin coupling
#############################
initialSpinUpFraction  0.5    # initialize the system with this spin up fraction
fractionMinNrFreeAtoms 0.0    # free atoms cannot drop under this fraction
#############################
termalisationSteps       0    # if <= 0 use the default: Ntot^2/2
Steps                    0    # steps used for each temperature; if <= 0 use the default: 100*termalisationSteps
avrgStep                 0    # steps interval used to calculate averages; if <= 0 use the default: sqrt(Steps)
#############################
Tmin        0.0     # minimum temperature
Tstep       0.2     # temperature increasing step (can be negative in case Tmin>Tmax)
Tmax        2.0     # maximum temperature
#############################
saveResults 0       # save result into file results.csv
doPlot      0       # plot energy, magnetization and nr of molecules as a function of temperature
goParallel  0       # try multi-threading
#############################
"""
)
    close(INI)
    println("Parameter file \"$file\" was missing and a default one was created.\nPlease edit it and rerun.")
    exit()

end
