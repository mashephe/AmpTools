#if !(defined USERNEG2LNLIKCONTRIB)
#define USERNEG2LNLIKCONTRIB


#include <string>
#include <vector>
#include "IUAmpTools/Neg2LnLikContrib.h"

using namespace std;

template< class T >
class UserNeg2LnLikContrib : public Neg2LnLikContrib
{
  
public:
  
  /**
   * This is the default constructor.  It should be called in the default
   * constructor of the user's derived class.
   */
  UserNeg2LnLikContrib() : Neg2LnLikContrib() { }
  
  
  /**
   * This constructor takes a list of arguments and stores them.  There should
   * be a corresponding constructor in the user's derived class that calls
   * this constructor.
   */
  UserNeg2LnLikContrib( const vector< string >& args ) : Neg2LnLikContrib( args ) { }
  
  
  /**
   * This is the destructor.
   */
   ~UserNeg2LnLikContrib() {
    
    for( typename map< string, T* >::iterator amp = m_ampInstances.begin();
         amp != m_ampInstances.end(); ++amp ){
      
       delete amp->second;
    }
  }
  
  Neg2LnLikContrib* newNeg2LnLikContrib( const vector< string >& args ) const{
    
    // we need the identifier of the new amplitude first:
    T* newAmp = new T( args );
    
    string ident = newAmp->identifier();
    
    if( m_ampInstances.find( ident ) == m_ampInstances.end() ){
      
      m_ampInstances[ident] = newAmp;
    }
    else{
      
      // already have a functional instance, so delete this one
      delete newAmp;
    }
    
    return m_ampInstances[ident];
  }
  
  
  /**
   * This method can create a clone of an amplitude (of the derived type).
   */
  Neg2LnLikContrib* clone() const {
    
    return ( isDefault() ? new T() : new T( arguments() ) );
   }
  
private:
  
  mutable map< string, T* > m_ampInstances;
  
};

#endif
