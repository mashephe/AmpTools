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
  
  float m_s12;
  float m_s23;

};

#endif
