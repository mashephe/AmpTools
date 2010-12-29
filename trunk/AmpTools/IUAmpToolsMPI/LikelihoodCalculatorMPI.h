#if !defined(LIKELIHOODCALCULATORMPI)
#define LIKELIHOODCALCULATORMPI

#include "IUAmpToolsMPI/MPITag.h"
#include "IUAmpToolsMPI/ParameterManagerMPI.h"
#include "IUAmpTools/LikelihoodCalculator.h"

class MISubject;

class LikelihoodCalculatorMPI : public LikelihoodCalculator
{

 public:

  enum { kFirstId = MPITag::kMaxTags };

  // the master node should use this constructor
  LikelihoodCalculatorMPI( const AmplitudeManager& ampManager,
                          const NormIntInterface& normInt,
                          DataReader& dataReader,
                          ParameterManagerMPI& parManager );
  
  ~LikelihoodCalculatorMPI();

  // override this operator which is where the likelihood gets
  // calculated -- this will be called on the master node by
  // the MinuitMinimiationManager
  double operator()();

  // also override the update method which does precalculation in the
  // parent class
  void update( const MISubject* );

  // the following functions are used by the LikelihoodManager to trigger
  // portions of the likeihood calculation -- they should only be called
  // on the worker nodes
  void updateParameters();
  void computeLikelihood();

 private:

  static int m_idCounter;

  void setupMPI();

  ParameterManagerMPI& m_parManager;
  bool m_functionEvaluated;
  int m_thisId;

  int m_rank;
  int m_numProc;
  bool m_isMaster;
};

#endif
