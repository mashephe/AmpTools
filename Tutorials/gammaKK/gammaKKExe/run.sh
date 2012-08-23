#!/bin/tcsh -e

##########################################################
#
# 2012/08/03
# Author: Kei Moriya (Indiana University)
#
# Example shell script that will run all of the
# tutorial programs for J/psi -> gamma + Ks Ks
#
# The programs are
#
# 1. generatePhaseSpace
# 2. toyAcceptance
# 3. plotData
# 4. generatePhysics
# 5. fitAmplitudesATI
# 6. plotResults
#
# To compile and run the programs, the variables
# - AMPTOOLS
# - MCTOOLKIT
# must be set, and also,
# GAMMAKK must be set to the top directory
# (source the file setup_gammaKK.csh in the top
#  directory to set this).
#
# To compile programs, simply type make,
# or run this script, which has make included.
# This directory is NOT included in the top directory
# Makefile, so this directory cannot be made with the
# top level make.
#
# If the variable VERBOSE is set below, this shell
# will ask for user input at the end of each stage
# so that you can understand which program is running,
# instead of having everything flash by.
#
##########################################################

# Comment out this line to run all programs without pause
# setenv VERBOSE 1

# 0. Make executables
cd $GAMMAKK/gammaKKExe
make
# CHANGE LATER
cd $GAMMAKK/gammaKKExe

#--------------------------------------------------------------------------------------------------
#
# 1. generate 100,000 phase space events
#
$GAMMAKK/gammaKKExe/generatePhaseSpace phasespace.gen.root 100000
if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished generating phase space events";
    echo "Press return key to continue";
    echo "Next step:     apply toy acceptance";
    read inputline > /dev/null
endif

#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#
# 2. apply toy acceptance to phase space
#
$GAMMAKK/gammaKKExe/toyAcceptance phasespace.gen.root phasespace.acc.root
if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished applying toy acceptance to phase space events";
    echo "Press return key to continue";
    echo "Next step:     create plots for phase space MC";
    read inputline > /dev/null
endif
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#
# 3. create plots to compare the flat Dalitz plots
#
mkdir -p ./figures/phasespace.gen.root/
mkdir -p ./figures/phasespace.acc.root/
$GAMMAKK/gammaKKExe/plotData phasespace.gen.root plots.phasespace.gen.root
$GAMMAKK/gammaKKExe/plotData phasespace.acc.root plots.phasespace.acc.root
if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished making plots for phase space events";
    echo "Press return key to continue";
    echo "Next step:     generate physics events";
    read inputline > /dev/null
endif
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#
# 4. generate 10,000 physics events
#
# Use gammaKKHelicityAmp.cc to generate spin-0 resonance
$GAMMAKK/gammaKKExe/generatePhysics helicity0.cfg physics.helicity0.gen.root    10000 1

if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished generating physics events";
    echo "Press return key to continue";
    echo "Next step:     apply toy acceptance to data";
    read inputline > /dev/null
endif
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#
# 5. apply toy acceptance to physics
#
$GAMMAKK/gammaKKExe/toyAcceptance physics.helicity0.gen.root physics.helicity0.acc.root

if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished applying toy acceptance to physics events";
    echo "Press return key to continue";
    echo "Next step:     create plots for data";
    read inputline > /dev/null
endif
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#
# 6. create plots to compare the physics Dalitz plots
#
### helicity amplitudes for spin-0 resonance
mkdir -p ./figures/physics/helicity0.gen.root/
mkdir -p ./figures/physics/helicity0.acc.root/

### TwoPiAngles amplitudes for spin-0 resonance
mkdir -p ./figures/physics/helicity0.gen.root/
mkdir -p ./figures/physics/helicity0.acc.root/

### helicity amplitudes for spin-0 resonance
$GAMMAKK/gammaKKExe/plotData physics.helicity0.gen.root plots.physics.helicity0.gen.root

if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished creating plots for physics events";
    echo "Press return key to continue";
    echo "Next step:     fit the data";
    read inputline > /dev/null
endif
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#
# 7. Fit the data with acceptance and see if we can pull out the
#    original parameters.
#

# helicity amplitudes for spin-0 resonance
$GAMMAKK/gammaKKExe/fitAmplitudesATI helicity0.fit.cfg

if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished doing fits to helicity0 events";
    echo "Press return key to continue";
    echo "Next step:     plot fit results";
    read inputline > /dev/null
endif
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#
# 8. Make plots of the intensity of each amplitude based on
#    the fit results given in the .fit file generated by 7.
#

# helicity amplitudes for spin-0 resonance
$GAMMAKK/gammaKKExe/plotResults helicity0.fit.cfg helicity0.fit plots.fitresult.helicity0.root

if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished creating ROOT files for fit results";
    echo "Press return key to continue";
    echo "Next step:     show the histograms created using ROOT macro";
    read inputline > /dev/null
endif
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#
# 9. Plot the fit results using ROOT macro
#

# helicity amplitudes for spin-0 resonance
mkdir -p ./figures/plots.fitresult.helicity0.root/
root $GAMMAKK/gammaKKExe/plotRootFile.C'("plots.fitresult.helicity0.root",1)'

if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished creating ROOT files for fit results";
#    echo "Press return key to continue";
#    read inputline > /dev/null
endif
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
# End of script.
echo ""
echo ""
echo "Finished running all programs."
echo ""
echo ""
#--------------------------------------------------------------------------------------------------
