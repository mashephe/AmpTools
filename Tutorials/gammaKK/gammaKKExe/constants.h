#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

/**
 * This is a utility class that contains mathematical constants
 * such as pi, and also physics constants
 * such as the known masses and widths of various particles.
 *
 * The main purpose of this program is that by including this
 * header file, you have access to all constants. This helps in
 * - ensuring consistent use of the same constant
 *   (avoid having pi defined in various files with various
 *    degrees of precision)
 * - (hopefully) no need to look up masses every time you want
 *   to change parameters. Also, consistently have the most
 *   precise masses etc.
 *
 *
 * Usage Example: J/psi -> pi+ pi- pi0
 *
 * ---------------------------------------
 * // ... beginning of code
 * #include "constants.h"
 * // ... other parts of code
 * double parentMass = m_Jpsi;
 * vector<double> daughteMasses;
 * daughterMasses.push_back(m_piPlus);
 * daughterMasses.push_back(m_piMinus);
 * daughterMasses.push_back(m_piZero);
 * // ... other parts of code
 * ---------------------------------------
 *
 */

/////////////////////////////////////////////////////////////////////
//                                                                 //
//                    Mathematical constants                       //
//                                                                 //
/////////////////////////////////////////////////////////////////////

/**
 * Taken from Wikipedia (20 digits)
 * If the variable is called PI, this conflicts with
 * the definition of PI in GPUCustomTypes.h,
 * and will give a compile error.
 */
const double _PI = 3.14159265358979323846;

/////////////////////////////////////////////////////////////////////
//                                                                 //
//                       Physics constants                         //
//                                                                 //
/////////////////////////////////////////////////////////////////////

const double SPEED_OF_LIGHT = 29.9792458; // in cm/ns

/////////////////////////////////////////
//                                     //
//       All masses, widths are        //
//       in GeV, all lifetimes in ns.  //
//       Retrieved from PDG on         //
//       2012/07/06.                   //
//                                     //
/////////////////////////////////////////

// we don't need this, but just for consistency
const double m_photon = 0;

// ===   leptons   ===

const double m_e   = 0.000510998928;
const double m_mu  = 0.1056583715;
const double m_tau = 1.77682;

// ===   mesons    ===
const double m_piPlus  = 0.13957018;
const double m_piMinus = m_piPlus;
const double m_piZero  = 0.1349766;

const double m_eta     = 0.547853;

const double m_rho     = 0.77549; // e+ e- neutral only
const double m_omega   = 0.78265;

const double m_etaPrime = 0.95778;

const double m_phi     = 1.019455;

const double m_KPlus   = 0.493677;
const double m_KMinus  = m_KPlus;
const double m_KZero   = 0.497614;


const double m_etac  = 2.9810;
const double m_Jpsi  = 3.096916;
const double m_chic0 = 3.41475;
const double m_chic1 = 3.51066;
const double m_chic2 = 3.55620;
const double m_hc    = 3.52541;

const double m_etac2S = 3.6389;
const double m_psiPrime = 3.686109;
// ===   baryons   ===

const double m_proton  = 0.938272046;
const double m_neutron = 0.939565379;

const double m_Lambda  = 0.1115683;
const double m_SigmaPlus  = 1.18937;
const double m_SigmaMinus = 1.197449;
const double m_SigmaZero  = 1.192642;

const double G_eta     = 1.30 * pow(10.,-6.); // 1.3 keV
const double G_rho     = 0.1491;              // charged only
const double G_omega   = 0.00849;
const double G_etaPrime = 0.199 * pow(10.,-3.); // 0.199 MeV
const double G_phi     = 0.00426;

const double G_etac = 0.0297;
const double G_Jpsi = 0.0929 * pow(10.,-3.); // 92.9 keV
const double G_chic0 = 0.0104;
const double G_chic1 = 0.86 * pow(10.,-3.);
const double G_chic2 = 1.98 * pow(10.,-3.);
// width of h_c listed only as < 1 MeV
const double G_etac2S = 0.010;
const double G_psiPrime = 304. * pow(10.,-6.); // 304 keV

// all in ns
const double tau_KPlus  = 12.380;
const double tau_KMinus = tau_KPlus;
const double tau_KShort = 0.089564;
const double tau_KLong  = 51.16;

const double tau_Lambda = 0.2632;
const double tau_SigmaPlus = 0.08018;
const double tau_SigmaMinus = 0.1479;

#endif
