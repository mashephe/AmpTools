
#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/NormIntInterface.h"

#include <complex>

using namespace std;

// The goal of this function is to facilitate easier input/output
// checks where one generates MC and then does a fit to check
// the fit result is consistent with generation.  The general problem
// is that the production coefficiencts at generation time have
// a non-trivial relationship with those that come from the fit.
// The easiest way to compare is to compare fit fractions, but
// computing the fit fractions is not possible without the
// normalization integrals, which are usually not computed as
// part of the generation.  This script **ASSUMES** the normalization
// integrals in the .fit output are true integrals of the amplitudes
// configured in the generator.  (This will NOT be true if for example,
// amplitudes have embedded parameters that are floating in the fit
// and they float to values different than those initialized in the
// generation.)  As long as this assumption is valid, then the generated
// fit fractions can be computed with the post-fit integrals and the
// generated production parameters in the config file used for generation.
// The computation of fit fractions from the fit uses standard functions.
// Uncertainty computations and therefore pulls are likely not
// rigorous and should be used as rough indicators of problems.


void compareFitFractions( string genCfgFile, string fitFile ){
  
  FitResults fr( fitFile );
  
  ConfigFileParser cfParser( genCfgFile );
  ConfigurationInfo* cfgInfo = cfParser.getConfigurationInfo();
  
  vector< ReactionInfo* > reactions = cfgInfo->reactionList();
  
  for( vector< ReactionInfo* >::iterator reacInfo = reactions.begin();
      reacInfo != reactions.end(); ++reacInfo ){
    
    string reaction = (*reacInfo)->reactionName();
    
    vector< AmplitudeInfo* > ampList = cfgInfo->amplitudeList( reaction );
    const NormIntInterface* ni = fr.normInt( reaction );
    
    // ** first nested loop over the generated production parameters and
    //    the results of the normalization integrals to pick out all
    //    of the relevant info needed to construct the generated fit
    //    fractions
    
    // need to renormalize the input PDF based using the input parameters
    // and the integrals of the PDF from the fit file
    complex< double > inputTotal( 0, 0 );
    
    map< string, double > nameGenFrac;
    vector< string > ampNameVec;
    
    for( vector< AmplitudeInfo* >::iterator ampInfo = ampList.begin();
        ampInfo != ampList.end(); ++ampInfo ){
      
      string ampName = (**ampInfo).fullName();
      ampNameVec.push_back( ampName );
      
      complex< double > inputProdAmp = atof( (**ampInfo).scale().c_str() ) * (**ampInfo).value();
      
      // this is the unnormalized generated intensity
      nameGenFrac[ampName] = abs( inputProdAmp * conj( inputProdAmp ) *
                                 ni->ampInt( ampName, ampName ) );
      
      // need to compute the "generated total" in order to renormalize
      // the generated production parameters
      for( vector< AmplitudeInfo* >::iterator conjAmpInfo = ampList.begin();
          conjAmpInfo != ampList.end(); ++conjAmpInfo ){
        
        complex< double > inputConjProdAmp = atof( (**conjAmpInfo).scale().c_str() ) *
        (**conjAmpInfo).value();
        
        inputTotal += inputProdAmp * conj( inputConjProdAmp ) *
        ni->ampInt( ampName, (**conjAmpInfo).fullName() );
      }
    }
    
    // do a check to see if this looks like a real number
    if( imag( inputTotal ) > 1E-6 ){
      
      cout << "WARNING:  normalization factor should be real: " << inputTotal << endl;
    }
    
    // ** second loop over the amplitude names -- renormalize to produce
    //    the fit fractions and print things to the screen
    
    pair< double, double > fitTotal = fr.intensity( ampNameVec );
    
    cout << setw(90) << setfill( '*' ) << "\n" << setfill( ' ' );
    cout << "REACTION:  " << reaction << endl << endl;
    
    cout << setw( 40 ) << left << " "
    << setw( 13 ) << right <<  "Generated  "
    << setw( 20 ) << "Fit Result "
    << setw( 12 ) << endl;
    
    cout << setw( 40 ) << left << "Amplitude"
    << setw( 13 ) << right << "Fraction [%]"
    << setw( 20 ) << "Fraction [%]"
    << setw( 12 ) << "Pull" << endl;
    
    cout << setw(90) << setfill( '-' ) << "\n" << setfill( ' ' );
    
    for( vector< string >::iterator name = ampNameVec.begin();
        name != ampNameVec.end(); ++name ) {
      
      // now record the fit fractions
      pair<double,double> fitFrac = fr.intensity( *name );
      fitFrac.first /= fitTotal.first * 0.01;
      fitFrac.second /= fitTotal.first * 0.01;
      
      nameGenFrac[*name] /= abs( inputTotal ) * 0.01;
      
      cout << setw( 40 ) << left <<  *name
      << setw( 10 ) << right << fixed << setprecision( 2 ) << nameGenFrac[*name]
      << setw( 15 ) << fitFrac.first << " +/- " << fitFrac.second
      << setw( 12 ) << ( fitFrac.first - nameGenFrac[*name] ) / fitFrac.second
      << endl;
    }
    
    cout << setw(90) << setfill( '*' ) << "\n" << setfill( ' ' );
    cout << endl;
  }
}
