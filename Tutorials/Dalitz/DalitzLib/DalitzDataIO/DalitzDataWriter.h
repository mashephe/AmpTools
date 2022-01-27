#if !defined(DALITZDATAWRITER)
#define DALITZDATAWRITER

#include "IUAmpTools/Kinematics.h"

#include "TTree.h"
#include "TFile.h"

class DalitzDataWriter{

public:
		
  DalitzDataWriter( const string& outFile );
  ~DalitzDataWriter();

  void writeEvent( const Kinematics& kin );

  int eventCounter() const { return m_eventCounter; }

private:

  TFile* m_outFile;
  TTree* m_outTree;
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
  
  double m_s12;
  double m_s23;

};

#endif
