#if !defined(GAMMAKKDATAWRITER)
#define GAMMAKKDATAWRITER

#include "IUAmpTools/Kinematics.h"

#include "TTree.h"
#include "TFile.h"

/**
 * This class is for writing out event 4-vectors,
 * in conjunction with gammaKKDataWriter,
 * which reads in events.
 *
 * The constructor takes in the outfile name desired,
 * and creates a ROOT file with that name, which contains
 * a Tree named "nt". Within this tree, branches for
 * E, x, y, z are created for 3 particles, which are
 * filled for each event. The variable m_eventCounter
 * is initially set to 0.
 * Any other variable that we wish to write out can be
 * set here, and should be implemented in the writeEvent()
 * function.
 *
 * To fill events, we rely on the AmpTools/Kinematics class,
 * which contains a vector of 4-vectors for the event we want.
 * The event counter is incremented. This variable can be
 * examined via the function eventCounter(), but no checks
 * are put in to inspect its value.
 *
 */

class gammaKKDataWriter{

public:
		
  gammaKKDataWriter( const string& outFile );
  ~gammaKKDataWriter();

  void writeEvent( const Kinematics& kin );

  void write();

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

  // Dalitz variables (cyclic)
  float m_s12;
  float m_s23;
  float m_s31;

};

#endif
