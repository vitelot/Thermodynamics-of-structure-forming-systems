DEBUG = 3;

# parameters
Ntot::Int         = 200; # total number of atoms
Hfield::Double    = 0.0; # external magnetic field
Jcoupling::Double = 0.5; # spin-spin coupling

Tmin::Double  = 0.1;
Tstep::Double = 0.01;
Tmax::Double  = 3.0;

termalisationSteps::Int = Ntot*Ntot÷2; # run these nr of steps at the beginning discarding the energy
Steps::Int              = 100*termalisationSteps; # total nr of steps
avrgStep::Int           = floor(Int, sqrt(Steps)); # used to estimate the average energy

# moleculeSplitProbability::Double = 0.5; # decides if to split molecules or join atoms
# spinFlipProbability::Double      = 0.5; # decides if to flip atoms or make a molecule move
fractionMinNrFreeAtoms = 0.2; # do not accept less than this fraction of free atoms

saveResults::Bool = true; # saves results into file results.CSV
doPlot::Bool      = true; # plot and save it on pdf
goParallel::Bool  = true; # well, it does it, multi-thread