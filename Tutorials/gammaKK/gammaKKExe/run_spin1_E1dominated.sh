#!/bin/tcsh -e

# Generate physics using gammaKKHelicityAmps with E1 = 95%, M2 = 5%
# Parameters are cfg file, ROOT file, # of events, random seed.
generatePhysics spin1_E1dominated.helicity.cfg physics.spin1_E1dominated.helicity.gen.root 3000000 1
plotData physics.spin1_E1dominated.helicity.gen.root plots.spin1_E1dominated.helicity.gen.root

# Generate physics using MultipoleAmps with E1 = 95%, M2 = 5%
# Parameters are cfg file, ROOT file, # of events, random seed.
generatePhysics spin1_E1dominated.multipole.cfg physics.spin1_E1dominated.multipole.gen.root 3000000 1
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
root -b -q 'plotRootFile.C("plots.fitresult.spin1_E1dominated.helicity_with_helicity.root")'
root -b -q 'plotRootFile.C("plots.fitresult.spin1_E1dominated.helicity_with_multipole.root")'
root -b -q 'plotRootFile.C("plots.fitresult.spin1_E1dominated.multipole_with_helicity.root")'
root -b -q 'plotRootFile.C("plots.fitresult.spin1_E1dominated.multipole_with_multipole.root")'

cd figures/plots.fitresult.spin1_E1dominated.helicity_with_helicity.root
mypdfmerge.pl
mv eachbin.pdf ../../helicity_with_helicity.pdf
cd ../../

cd figures/plots.fitresult.spin1_E1dominated.helicity_with_multipole.root
mypdfmerge.pl
mv eachbin.pdf ../../helicity_with_multipole.pdf
cd ../../

cd figures/plots.fitresult.spin1_E1dominated.multipole_with_helicity.root
mypdfmerge.pl
mv eachbin.pdf ../../multipole_with_helicity.pdf
cd ../../

cd figures/plots.fitresult.spin1_E1dominated.multipole_with_multipole.root
mypdfmerge.pl
mv eachbin.pdf ../../multipole_with_multipole.pdf
cd ../../

open *.pdf
