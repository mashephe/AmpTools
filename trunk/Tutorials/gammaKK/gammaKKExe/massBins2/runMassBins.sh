#!/bin/tcsh -e

##########################################################
#
# 2012/08/21
# Author: Kei Moriya (Indiana University)
#
# Example shell script that will run an example of
# mass-independent fits for
# J/psi -> gamma + Ks Ks
#
# The script will first divide the data files into
# bins of Ks Ks mass, where the binning is specified
# within the cfg file by the keyword 'binning'.
#
# Once the input ROOT files are divided into mass bins,
# fits are done on each one, with the next fit
# being initialized with the previous fit result.
#
#
# -------------------------------------------------------
#
# Programs used:
#
# 1. generatePhaseSpace
# 2. toyAcceptance
# 3. generatePhysics
# 4. splitByMass
# 5. fitAmplitudesATIMassBins
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
cd $GAMMAKK/gammaKKExe/massBins2

#--------------------------------------------------------------------------------------------------
#
# 1. generate 10,000,000 phase space events
#
# $GAMMAKK/gammaKKExe/generatePhaseSpace phasespace.gen.root 10000000
if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished generating phase space events";
    echo "Press return key to continue";
    read inputline > /dev/null
endif

#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#
# 2. apply toy acceptance to phase space
#
# $GAMMAKK/gammaKKExe/toyAcceptance phasespace.gen.root phasespace.acc.root
if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished applying toy acceptance to phase space events";
    echo "Press return key to continue";
    read inputline > /dev/null
endif
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#
# 3. generate 1,000,000 physics events
#
# Use gammaKKHelicityAmp.cc to generate spin-0 resonance
# $GAMMAKK/gammaKKExe/generatePhysics helicity0.cfg physics.helicity0.gen.root    1000000 1

if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished generating physics events";
    echo "Press return key to continue";
    read inputline > /dev/null
endif
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#
# 4. apply toy acceptance to physics
#
# $GAMMAKK/gammaKKExe/toyAcceptance physics.helicity0.gen.root physics.helicity0.acc.root

if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished applying toy acceptance to physics events";
    echo "Press return key to continue";
    read inputline > /dev/null
endif
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#
# 5. split files in KK mass
#
# data
$GAMMAKK/gammaKKExe/splitByMass helicity0.fit.cfg physics.helicity0.acc.root physics.helicity0.mass_
# acc MC
$GAMMAKK/gammaKKExe/splitByMass helicity0.fit.cfg phasespace.acc.root phasespace.acc.mass_
# gen MC
$GAMMAKK/gammaKKExe/splitByMass helicity0.fit.cfg phasespace.gen.root phasespace.gen.mass_

if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished applying toy acceptance to physics events";
    echo "Press return key to continue";
    read inputline > /dev/null
endif
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#
# 6. Fit the data with acceptance and see if we can pull out the
#    original parameters.
#

# helicity amplitudes for spin-0 resonance
$GAMMAKK/gammaKKExe/fitAmplitudesATISplitMassBins helicity0.fit.cfg

if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished doing fits to helicity0 events";
    echo "Press return key to continue";
    read inputline > /dev/null
endif
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#
# 8. Make plots of the intensity of each amplitude based on
#    the fit results given in the .fit file generated by 7.
#

# helicity amplitudes for spin-0 resonance
#$GAMMAKK/gammaKKExe/plotResults helicity0.fit.cfg helicity0.fit plots.fitresult.helicity0.root

if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished creating ROOT files for fit results";
    echo "Press return key to continue";
    read inputline > /dev/null
endif
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#
# 9. Plot the fit results using ROOT macro
#

# helicity amplitudes for spin-0 resonance
#root $GAMMAKK/gammaKKExe/plotRootFile.C'("plots.fitresult.helicity0.root")'

if( ($?VERBOSE) ) then
    echo "";
    echo "-----------------------------------------------------------";
    echo "Finished creating ROOT files for fit results";
    echo "Press return key to continue";
    read inputline > /dev/null
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