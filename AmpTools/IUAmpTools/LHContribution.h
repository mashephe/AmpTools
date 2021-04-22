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
#include "IUAmpTools/BinnedData.h"
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

class LHContribution : public MIFunctionContribution
{
  
public:
  
  /**
   * The default constructor.  The user's derived class should contain a
   * default constructor that calls this constructor.
   */
  LHContribution(  ) : m_isDefault(true), MIFunctionContribution(m_minManager) {
    
  }
  
  /**
   * This constructor takes a list of arguments to inititialize a
   * LHContribution and then stores them.  The user's derived class should
   * contain a similar constructor that calls this one.
   */
  LHContribution( const vector< string >& args ) : m_isDefault(false), 
  m_args(args), MIFunctionContribution(m_minManager){ }
  
  /**
   * This is the destructor.
   */
  virtual ~LHContribution(){}


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
   * of the LHContribution.  It is the LHContribution name and all the arguments
   * concatenated together.  Every LHContribution with the same identifier
   * should behave the same way.
   */
  string identifier() const;
  
  /**
   * This must be overriden by the user and indicates how to convert a list
   * of strings (arguments) into a pointer to a new instance of the
   * users defined LHContribution.
   *
   * The user can avoid writing this method by inheriting from the
   * UserLHContribution class (which derives from this class).
   *
   * \param[in] args a list of string arguments that may, for example, be
   * specified in a configuration file
   *
   *  \see UserLHContribution
   */
  virtual LHContribution* newLHContribution( const vector< string >& args ) const = 0;
  
  /**
   * A function that the user must write that indicates how a LHContribution
   * can duplicate itself.  It returns a pointer to a new instance of the
   * LHContribution that behaves in exactly the same way as the original.
   *
   * The user can avoid writing this method by inheriting from the
   * UserLHContribution class (which derives from this class).
   *
   *  \see UserLHContribution
   */
  virtual LHContribution* clone() const = 0;
  
  /**
   * A function that indicates if this LHContribution was created with the
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
   * parameters for this LHContribution, false otherwise.  The function calls
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
   * This is the user-defined function that computes a single real LHContribution
   * as a function of a single coordinate (e.g. sqrt(s)).
   *
   */
  virtual double calcLHContribution( double x ) const ;


  void setMinimizationManager(MinuitMinimizationManager *minManager){
    m_minManager = minManager;
  }
  
#ifdef GPU_ACCELERATION 
  
  /**
   * If GPU_ACCELERATION flag is set this is the member function that sets
   * the current permutation and then calls the user-defined routine
   * to launch the GPU kernel.  It is the GPU analog of calcAmplitudeAll.
   */
  virtual void calcLHContributionGPU( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
  
  /**
   * The user override this route and use it pass any parameters to a global
   * C function that actually launches the GPU kernel
   */
  virtual void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {
    
    cout << "\nNo GPU function for calculating " << name() << " is defined." << endl;
    assert( false );
  }


  
#endif //GPU_ACCELERATION
  
protected:
  
  /**
   * Any user-defined class derived from Amplitude that has a parameter
   * in the amplitude should register the parameter using this routine.  This
   * shoudl be done for each parameter in the constructor of the user's 
   * Amplitude class.
   *
   * \param[in] par a reference to the AmpParmeter object
   */
  void registerParameter( AmpParameter& par );

  static BinnedData m_data; 
  static AmpParameter parameters[100];
  
private:
  
  bool m_isDefault;
  
  vector<string> m_args;
  
  vector< AmpParameter* > m_registeredParams;

  static MinuitMinimizationManager *m_minManager;

  
};


#endif
