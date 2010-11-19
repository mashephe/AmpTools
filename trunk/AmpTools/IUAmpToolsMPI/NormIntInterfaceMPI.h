#if !defined(NORMINTINTERFACEMPI)
#define NORMINTINTERFACEMPI

#include "IUAmpTools/NormIntInterface.h"

class NormIntInterfaceMPI : public NormIntInterface
{

 public:

  enum IntType { kNormInt, kAmpInt };

  NormIntInterfaceMPI( const string& normIntFile );
  NormIntInterfaceMPI( DataReader* genMCData, DataReader* accMCData, 
		       const AmplitudeManager& ampManager );

  ~NormIntInterfaceMPI();

  // Note: for managing cases where NI must be recomputed at each 
  // iteration, this class must override the normInt method
  // of NormIntInterface to ensure integrals are properly
  // recomputed

 private:

  void setupMPI();
  void sumIntegrals( IntType type );

  bool m_mpiSetup;

  int m_rank;
  int m_numProc;
  bool m_isMaster;
};

#endif
