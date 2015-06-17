#if !defined(NBODYPHASESPACEFACTORY)
#define NBODYPHASESPACEFACTORY

#include <vector>

#include "TLorentzVector.h"

using namespace std;

class NBodyPhaseSpaceFactory
{
	
 public:
	
  NBodyPhaseSpaceFactory( double parentMass, const vector<double>& childMass );
	
  vector<TLorentzVector> generateDecay() const;
	
 private:
        
  static const double kPi;
	
  double pdk( double a, double b, double c ) const;
  double random( double low, double hi ) const;
	
  double m_parentMass;
  vector<double> m_childMass;    // vector of daughter masses
  int m_Nd;                      // number of decay products
};

#endif
