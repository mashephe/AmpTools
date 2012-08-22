#!/bin/tcsh -e

# Generate 1,000,000 phase space MC events
generatePhasespace phasespace.gen.root 1000000
plotData phasespace.gen.root plots.phasespace.gen.root

# Apply toy acceptance
toyAcceptance phasespace.gen.root phasespace.acc.root

# Generate physics using gammaKKHelicityAmps with E1 = 95%, M2 = 5%
# Parameters are cfg file, ROOT file, # of events, random seed.
generatePhysics spin1_E1dominated.helicity.cfg physics.spin1_E1dominated.helicity.gen.root 100000 1
plotData physics.spin1_E1dominated.helicity.gen.root plots.spin1_E1dominated.helicity.gen.root

# Generate physics using MultipoleAmps with E1 = 95%, M2 = 5%
# Parameters are cfg file, ROOT file, # of events, random seed.
generatePhysics spin1_E1dominated.multipole.cfg physics.spin1_E1dominated.multipole.gen.root 100000 1
plotData physics.spin1_E1dominated.multipole.gen.root plots.spin1_E1dominated.multipole.gen.root

# Apply toy acceptance
toyAcceptance physics.spin1_E1dominated.helicity.gen.root physics.spin1_E1dominated.helicity.acc.root
toyAcceptance physics.spin1_E1dominated.multipole.gen.root physics.spin1_E1dominated.multipole.acc.root

# Do fit using config file
#
# fit helicity generated with helicity amps
time fitAmplitudesATI spin1_E1dominated.fit.helicity_with_helicity.cfg
#
# fit helicity generated with multipole amps
time fitAmplitudesATI spin1_E1dominated.fit.helicity_with_multipole.cfg
#
# fit multipole generated with helicity amps
time fitAmplitudesATI spin1_E1dominated.fit.multipole_with_helicity.cfg
#
# fit multipole generated with multipole amps
time fitAmplitudesATI spin1_E1dominated.fit.multipole_with_multipole.cfg

echo "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ="
echo "                            FINISHED WITH FITS                                           "
echo "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ="

# Make ROOT file out of fit results
#
# fit helicity generated with helicity amps
time plotResults spin1_E1dominated.fit.helicity_with_helicity.cfg spin1_E1dominated.helicity_with_helicity.fit plots.fitresult.spin1_E1dominated.helicity_with_helicity.root
# fit helicity generated with multipole amps
time plotResults spin1_E1dominated.fit.helicity_with_multipole.cfg spin1_E1dominated.helicity_with_multipole.fit plots.fitresult.spin1_E1dominated.helicity_with_multipole.root
# fit multipole generated with helicity amps
time plotResults spin1_E1dominated.fit.multipole_with_helicity.cfg spin1_E1dominated.multipole_with_helicity.fit plots.fitresult.spin1_E1dominated.multipole_with_helicity.root
# fit multipole generated with multipole amps
time plotResults spin1_E1dominated.fit.multipole_with_multipole.cfg spin1_E1dominated.multipole_with_multipole.fit plots.fitresult.spin1_E1dominated.multipole_with_multipole.root

# Plot the histograms from the ROOT file using ROOT macro
mkdir -p figures/plots.fitresult.spin1_E1dominated.helicity_with_helicity.root
mkdir -p figures/plots.fitresult.spin1_E1dominated.helicity_with_multipole.root
mkdir -p figures/plots.fitresult.spin1_E1dominated.multipole_with_helicity.root
mkdir -p figures/plots.fitresult.spin1_E1dominated.multipole_with_multipole.root
root -b -q 'plotRootFile.C("plots.fitresult.spin1_E1dominated.helicity_with_helicity.root",2)'
root -b -q 'plotRootFile.C("plots.fitresult.spin1_E1dominated.helicity_with_multipole.root",2)'
root -b -q 'plotRootFile.C("plots.fitresult.spin1_E1dominated.multipole_with_helicity.root",2)'
root -b -q 'plotRootFile.C("plots.fitresult.spin1_E1dominated.multipole_with_multipole.root",2)'

echo "Open figures figures/plots.fitresult.spin1_E1dominated.helicity_with_helicity.root/PhotonCosTheta_gen.pdf"
echo "This shows the fit result to the photon angular distribution"
echo "With the cfg file just used, the fit result should give alpha = (-x*x + 6*x - 1) / (3x*x-2*x+3)"
echo "where x is the fraction of M2/E1".
echo "For this example of M2/E1 = 0.05, this should give alpha = -2.4162"
