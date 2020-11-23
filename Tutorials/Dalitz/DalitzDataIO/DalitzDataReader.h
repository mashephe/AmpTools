#if !defined(DALITZDATAREADER)
#define DALITZDATAREADER

#include <string>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "IUAmpTools/DataReader.h"
#include "IUAmpTools/UserDataReader.h"

using namespace std;

class DalitzDataReader : public UserDataReader< DalitzDataReader >{

public:

  DalitzDataReader() : UserDataReader< DalitzDataReader >() { }

  DalitzDataReader( const vector< string >& args );

  string name() const { return "DalitzDataReader"; }

  virtual Kinematics* getEvent();

  virtual void resetSource();

  virtual unsigned int numEvents() const;

  int eventCounter() const { return m_eventCounter; }


private:

  TFile* m_inFile;
  TTree* m_inTree;
  int m_eventCounter;

  double m_EnP1;
  double m_PxP1;
  double m_PyP1;
  double m_PzP1;

  double m_EnP2;
  double m_PxP2;
  double m_PyP2;
  double m_PzP2;

  double m_EnP3;
  double m_PxP3;
  double m_PyP3;
  double m_PzP3;
  
  double m_weight;

};

#endif
