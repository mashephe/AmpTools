#if !defined(BINNEDDATA)
#define BINNEDDATA

#include <string>
#include <vector>
#include <cassert>

using namespace std;

struct datapoint{
  double x;
  double y;
  double eyl;
  double eyh;
  double exl;
  double exh;
};

using namespace std;

class BinnedData
{
  
public:
  BinnedData(){ };
  virtual ~BinnedData(){
    m_thedata.clear();
  };
  void addPoint(double x, double y, double exl, double exh, double eyl, double eyh);
  void addPoint(double x, double y, double ey);
  void addPoint(double x, double y, double ex, double ey);
  void addPoint(double x, double y, double ex, double eyl, double eyh);
  datapoint getPoint(unsigned int i);
  unsigned int getN(){
    return m_thedata.size();
  };

protected:

private:
  vector<datapoint> m_thedata;


};

#endif

