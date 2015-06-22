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

  float m_EnP1;
  float m_PxP1;
  float m_PyP1;
  float m_PzP1;

  float m_EnP2;
  float m_PxP2;
  float m_PyP2;
  float m_PzP2;

  float m_EnP3;
  float m_PxP3;
  float m_PyP3;
  float m_PzP3;
  
  float m_weight;

};

#endif
