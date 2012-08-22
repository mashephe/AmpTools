#if !defined(NBODYPHASESPACEFFACTORY)
#define NBODYPHASESPACEFACTORY

#include <vector>

#include "CLHEP/Vector/LorentzVector.h"

using namespace std;
using namespace CLHEP;

class NBodyPhaseSpaceFactory
{
	
 public:
	
  NBodyPhaseSpaceFactory( double parentMass, const vector<double>& childMass );
	
  vector<HepLorentzVector> generateDecay() const;
	
 private:
        
  static const double kPi;
	
  double pdk( double a, double b, double c ) const;
  double random( double low, double hi ) const;
	
  double m_parentMass;
  vector<double> m_childMass;    // vector of daughter masses
  int m_Nd;                      // number of decay products
};

#endif
