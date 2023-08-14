#if !defined( LHCONTRIBUTIONMANAGER )
#define LHCONTRIBUTIONMANAGER

#include "IUAmpTools/Neg2LnLikContrib.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/ConfigurationInfo.h"

using namespace std;

class Neg2LnLikContribManager
{
public:
  /** Constructor.
   * Constructs a Neg2LnLikContribManager
   */ 
  Neg2LnLikContribManager();
  
  /** Destructor.
   */
  ~Neg2LnLikContribManager();

  /**
   * This function sets up the Neg2LnLikContribManager based on information
   * provided by a ConfigurationInfo object.  It is intended to be the
   * most user-friendly way of configuring the manager for
   * use.
   *
   * \param[in] configInfo a pointer to a ConfigurationInfo object
   */
  void setupFromConfigurationInfo( const ConfigurationInfo* configInfo );


  void addNeg2LnLikContrib(const string& lhcontName,
                     const vector< string >& args);


  void registerNeg2LnLikContrib( const Neg2LnLikContrib& defaultNeg2LnLikContrib );

    /**
   * This tells an amplitude to use an external pointer to resolve the value
   * of some parameter inside of an amplitude.  The parameter is not a
   * production paramater as discussed above, but instead a parameter that
   * may be floating in the description of the decay, e.g., the mass or width
   * of some state.  This routine is here for convenience, one could also
   * get a list of all the amplitude factors and check each one for the parameter
   * and then use the functionality of the Amplitude class to set the parameter
   * pointer.
   *
   * \param[in] ampName the name of the amplitude that contains a factor that
   * contains the paramter of interest
   *
   * \param[in] parName the name of the parameter
   *
   * \param[in] ampParPtr a pointer to external memory that tells the amplitude
   * where to find the value of the parameter
   *
   * \see setAmpParValue
   * \see Amplitude::setParPtr
   */
  void setParPtr( const string& termName, const string& parName,
                  const double* ampParPtr );

   /**
   * This tells an amplitude to use a particular value for
   * some parameter inside of an amplitude.  The parameter is not a
   * production paramater as discussed above, but instead a parameter that
   * may be floating in the description of the decay, e.g., the mass or width
   * of some state.  This routine is here for convenience, one could also
   * get a list of all the amplitude factors and check each one for the parameter
   * and then use the functionality of the Amplitude class to set the parameter
   * value.
   *
   * \param[in] ampName the name of the amplitude that contains a factor that
   * contains the paramter of interest
   *
   * \param[in] parName the name of the parameter
   *
   * \param[in] ampParValue the value to which to set the parameter
   *
   * \see setAmpParValue
   * \see Amplitude::setParValue
   */
  void setParValue( const string& name, const string& parName,
                               double val );

   bool hasTerm(const string& name){
    return m_mapNameToNeg2LnLikContribs.count(name);
  };
                    
   /**
   * This function will be called whenever a parameter is updated.
   *
   * \see Amplitude::updatePar
   */
  void updatePar( const string& parName ) const;

  private:

  map< string, Neg2LnLikContrib* > m_registeredNeg2LnLikContribs;
  map< string, Neg2LnLikContrib* > m_mapNameToNeg2LnLikContribs;
};

#endif
