DEBUG = 3;

# parameters
N::Int = 200; # total number of atoms
T::Double = 1.0;  # temperature
H::Double = 0.0; # external magnetic field
J::Double = 1.0; # spin-spin coupling
Tmin = 0.5;
Tstep = 0.05;
Tmax = 3.0;
termalisationSteps = 1000; # run these nr of steps at the beginning discarding the energy
Steps::Int = 2_000_000;
avrgStep = 1000; # used to estimate the average energy
moleculeSplitProbability = 0.5; # decides if to split molecules or join atoms
spinFlipProbability = 0.5; # decides if to flip atoms or make a molecule move
