#if !defined(GAMMAKKDATAREADER)
#define GAMMAKKDATAREADER

#include <string>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "IUAmpTools/DataReader.h"
#include "IUAmpTools/UserDataReader.h"

using namespace std;

/**
 * This class inherits from IUAmpTools/DataReader
 * and is for reading in event 4-vectors,
 * in conjunction with gammaKKDataWriter,
 * which writes out events.
 *
 * This is accomplished by setting the ROOT
 * file name and Tree name in the constructor,
 * and reading in 3 4-vectors for each event.
 * The 4-vectors are formatted in E, x, y, z,
 * for the particles 1, 2, 3.
 *
 * The variable m_eventCounter keeps track of
 * where within the set of events we are, and
 * ensures that we do not try to get more
 * events (NULL will be returned).
 *
 * By calling the getEvent() method, the
 * (m_eventCounter)th event will be returned
 * in the form of a IUAmpTools/Kinematics pointer.
 * This object will have as a member a vector of
 *  4-vectors that were pushed in.
 *
 */

class gammaKKDataReader : public UserDataReader<gammaKKDataReader>{

public:
 gammaKKDataReader() : UserDataReader<gammaKKDataReader>(){}

  gammaKKDataReader(const vector<string> &args);

  string name() const {return "gammaKKDataReader";}
  
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

  // Optional variables to specify the range of
  // M23 that should be read in.
  bool  m_rangeSpecified;
  float m_min;
  float m_max;
  // Number of events that will be read in when
  // using the constructor with 3 arguments
  int   m_numEvents;

};

#endif
