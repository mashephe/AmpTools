#if !defined( LHCONTRIBUTION )
#define LHCONTRIBUTION

#include <map>
#include <complex>
#include <vector>
#include <string>
#include <cassert>
#include <set>

#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"
#include "MinuitInterface/MIFunctionContribution.h"

#ifdef GPU_ACCELERATION 
#include "cuda_runtime.h"
#include "GPUManager/CUDA-Complex.cuh"
class GPUManager;
#endif //GPU_ACCELERATION

using std::complex;
using namespace std;

class AmpParameter;
class MinuitMinimizationManager;
class MinuitParameter;

class Neg2LnLikContrib : public MIFunctionContribution
{
  
public:
  
  /**
   * The default constructor.  The user's derived class should contain a
   * default constructor that calls this constructor.
   */
  Neg2LnLikContrib(  ) : m_isDefault(true), MIFunctionContribution(m_minManager) {
    
  }
  
  /**
   * This constructor takes a list of arguments to inititialize a
   * Neg2LnLikContrib and then stores them.  The user's derived class should
   * contain a similar constructor that calls this one.
   */
  Neg2LnLikContrib( const vector< string >& args ) : m_isDefault(false),
  m_args(args), MIFunctionContribution(m_minManager){ }
  
  /**
   * This is the destructor.
   */
  virtual ~Neg2LnLikContrib(){}


  /** 
   * this method delivers the chi2/the likelihood
   */
  double operator()();

  /**
   * will be overridden by user to provide a chi2 or likelihood function
   * using binned data
   * 
   */
  virtual double neg2LnLikelihood();
  
  /**
   * Must be overridden by the user to provide the name of the amplitude.
   * This is necessary to connect ConfigurationInfo to specific class
   * instances when the AmplitudeManager is being setup.
   */
  virtual string name() const = 0;
  
  /**
   * This method returns a string that uniquely identifies the instance
   * of the Neg2LnLikContrib.  It is the Neg2LnLikContrib name and all the arguments
   * concatenated together.  Every Neg2LnLikContrib with the same identifier
   * should behave the same way.
   */
  string identifier() const;
  
  /**
   * This must be overriden by the user and indicates how to convert a list
   * of strings (arguments) into a pointer to a new instance of the
   * users defined Neg2LnLikContrib.
   *
   * The user can avoid writing this method by inheriting from the
   * UserNeg2LnLikContrib class (which derives from this class).
   *
   * \param[in] args a list of string arguments that may, for example, be
   * specified in a configuration file
   *
   *  \see UserNeg2LnLikContrib
   */
  virtual Neg2LnLikContrib* newNeg2LnLikContrib( const vector< string >& args ) const = 0;
  
  /**
   * A function that the user must write that indicates how a Neg2LnLikContrib
   * can duplicate itself.  It returns a pointer to a new instance of the
   * Neg2LnLikContrib that behaves in exactly the same way as the original.
   *
   * The user can avoid writing this method by inheriting from the
   * UserNeg2LnLikContrib class (which derives from this class).
   *
   *  \see UserNeg2LnLikContrib
   */
  virtual Neg2LnLikContrib* clone() const = 0;
  
  /**
   * A function that indicates if this Neg2LnLikContrib was created with the
   * default constructor.
   */
  bool isDefault() const { return ( m_isDefault == true ); }
  
  /**
   * Returns the list of arguments that was passed to the constructor.
   */
  vector<string> arguments() const { return m_args; }
  
  /**
   * The user can override this to do specific one-time tasks that
   * before the fit begins, but after the parameters have been initialized.
   * In general these tasks should go into the this init routine
   * rather than the consructor.
   */
  virtual void init(){}
  
  /**
   * This tells the AmpParameter with the indicated name to obtain its
   * value from some external location in memory.  This functionality is
   * is useful when fitting as the ParameterManager can automatically
   * updated these external values.  The function returns true if the
   * a parameter of the indicated name was found in the list of registered
   * parameters for this Neg2LnLikContrib, false otherwise.  The function calls
   * updatePar after setting the value.
   *
   * \param[in] name the name of the AmpParameter
   * \param[in] ptr the location in memory to find the value of the AmpParameter
   *
   * \see updatePar
   * \see AmplitudeManager::setAmpParPtr
   */
  bool setParPtr( const string& name, const double* ptr ) const;
  
  /**
   * This tells the AmpParameter with the indicated name to set its
   * value to the specified value.  The function returns true if the
   * a parameter of the indicated name was found in the list of registered
   * parameters for this amplitude, false otherwise. Function calls 
   * updatePar after setting the value.
   *
   * \param[in] name the name of the AmpParameter
   * \param[in] val the value to set AmpParameter to
   *
   * \see updatePar
   * \see AmplitudeManager::setAmpParValue
   */
  bool setParValue( const string& name, double val ) const;
  
  /**
   * The user may override this function to recalculate member data in the
   * amplitude class whenever a parameter changes.  For example, the user's
   * amplitude calculation may have expensive integration or other functions
   * that only need to be computed whenever a parameter changes, rather than
   * being computed on the fly for every event.
   *
   * \param[in] par a const reference to the parameter that has been updated
   */
  virtual void updatePar( const AmpParameter& par ) {}
  
  /**
   * \overload
   * 
   * A function used by the framework to signal that a parameter with a
   * particular name has changed.  This function searches the registered
   * parameters to see if this amplitude contains that parameter.  If so
   * it calls the virtual function updatePar with the AmpParameter.  The
   * function returns true of a parameter with the specified name is found.
   *
   * \param[in] name the name of the updated parameter
   */
  
  bool updatePar( const string& name ) const;
  
  /**
   * This is the user-defined function that computes a single real Neg2LnLikContrib
   * as a function of a single coordinate (e.g. sqrt(s)).
   *
   */
  virtual double calcNeg2LnLikContrib( double x ) const ;


  void setMinimizationManager(MinuitMinimizationManager *minManager){
    m_minManager = minManager;
  }
    
protected:
  
  void registerParameter( AmpParameter& par );
  
private:
  
  bool m_isDefault;
  
  vector<string> m_args;
  
  vector< AmpParameter* > m_registeredParams;

  static MinuitMinimizationManager *m_minManager;
  
};


#endif
