#ifndef DATAREADERMPI
#define DATAREADERMPI

#include "mpi.h"

#include "IUAmpTools/Kinematics.h"
#include "IUAmpToolsMPI/MPITag.h"

/**
 * This is a helper struct that is used to pass kinematic data between
 * different processes using MPI.
 */

struct KinStruct {
  
  int nPart;
  float weight;
  float e[Kinematics::kMaxParticles];
  float px[Kinematics::kMaxParticles];
  float py[Kinematics::kMaxParticles];
  float pz[Kinematics::kMaxParticles];
};

/**
 * This is a template that can be used to turn a user-defined DataReader
 * object into a parallel data reader using MPI.  Use of this template will
 * cause the data to be distributed evenly amongst the worker nodes.  Each
 * of the instances of this class on the worker nodes will behave as if
 * they have only a subset of the data.  The instance of it on the master
 * node will behave as if it has all of the data.  Be sure that the getEvent, 
 * resetSource, and numEvents methods in the user-defined class are declared
 * virtual.
 *
 * \ingroup IUAmpToolsMPI
 */

class DataReader;

template < class T >
class DataReaderMPI : public T
{
  
public:

  /**
   * This is the default constructor.
   */

  DataReaderMPI() : T() { m_isDefault = true; m_isMaster = true; }

  /**
   * This is the constructor for the templated class, which takes as
   * input a vector of arguments (as in the base class).
   */

  DataReaderMPI( const vector<string>& args );
  
  virtual ~DataReaderMPI();
  
  // override these functions with parallelized versions -- it is very 
  // important that the methods in the base class be declared virtual!
  
  Kinematics* getEvent();

  void resetSource();
  
  unsigned int numEvents() const;

  /**
   * This method can create a new data reader (of the derived type).
   */
  virtual DataReader* newDataReader( const vector< string >& args ) const{
    return new DataReaderMPI<T>( args );
  }


  /**
   * This method can create a clone of a data reader (of the derived type).
   */
  virtual DataReader* clone() const{
    return ( isDefault() ? new DataReaderMPI<T>() : new DataReaderMPI<T>( arguments() ) );
  }

  /**
   * Returns the list of arguments that was passed to the constructor.
   */
  virtual vector<string> arguments() const { return m_args; }

  /**
   * Returns true if this instance was created using the default constructor
   * and returns false otherwise.
   */
  virtual bool isDefault() const { return ( m_isDefault == true ); }

  
private:

          bool m_isDefault;
          vector<string> m_args;
  
  // some helper functions:
  
  void defineMPIType();
  void distributeData();
  void receiveData();
  
  Kinematics* createKin( KinStruct* kinStruct );
  void fillStruct( KinStruct* kinStruct, Kinematics* kin );
  
  MPI_Datatype MPI_KinStruct;
  
  int m_rank;
  int m_numProc;
  bool m_isMaster;
  
  vector<Kinematics*> m_ptrCache;
  vector<Kinematics*>::iterator m_ptrItr;
  
  unsigned int m_numEvents;
  
};


template< class T >
DataReaderMPI<T>::DataReaderMPI( const vector< string >& args ) : 
  T( args ),
  m_ptrCache( 0 ),
  m_ptrItr( m_ptrCache.begin() ),
  m_isDefault(false),
  m_args(args)
{
   
  MPI_Comm_rank( MPI_COMM_WORLD, &m_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &m_numProc );
  
  m_isMaster = ( m_rank == 0 );
  
  defineMPIType();
  
  if( m_isMaster ) distributeData();
  else receiveData();  
}

template< class T >
DataReaderMPI<T>::~DataReaderMPI(){
  
  // clean up data cache on worker nodes
  
  if( !m_isMaster && !m_isDefault ){
     
    for( vector<Kinematics*>::iterator ptrItr = m_ptrCache.begin();
        ptrItr != m_ptrCache.end();
        ++ptrItr ){
      
      delete *ptrItr;
    }
    
    m_ptrCache.clear();
  }
}

template< class T >
Kinematics* DataReaderMPI<T>::getEvent()
{
  
  if( m_isMaster ) return T::getEvent();
  
  if( m_ptrItr != m_ptrCache.end() ){
    
    // the standard behavior for DataReaders is that the calling class
    // takes ownership of the memory returned by getEvent -- this means
    // that this memory will be deleted after the call.  We don't want
    // our cache of data on the workers to be deleted since this 
    // doesn't allow for a reset and re-read of the data.  Return a new
    // Kinematics object then.
    
    return new Kinematics(**m_ptrItr++);
  }
  else{
    
    return NULL;
  }
}

template< class T >
void DataReaderMPI<T>::resetSource()
{
  
  if( m_isMaster ){
    
    //    cout << "Resetting master data source " << m_rank << endl;
    
    // on the master, we want the data reader to behave
    // like a normal instance of the non-MPI version so use the 
    // user defined source reset method
    
    T::resetSource();
  }
  else{
    
    //    cout << "Resetting data source on process with rank " << m_rank << endl;
    
    // on the workers the "source" is the local cache of Kinematics
    // objects so this should be reset
    
    m_ptrItr = m_ptrCache.begin();
  }
}

template< class T >
void DataReaderMPI<T>::distributeData()
{
  
  assert( m_numProc > 1 );
  
  MPI_Status status;
  
  int totalEvents = T::numEvents();
  int stepSize = totalEvents / ( m_numProc - 1 );
  int remainder = totalEvents % ( m_numProc - 1 );
  
  KinStruct* kinArray = new KinStruct[stepSize+1];
  
  for( int i = 1; i < m_numProc; ++i ){
    
    int nEvents = ( i > remainder ? stepSize : stepSize + 1 );
    
    cout << "Sending process " << i << " " << nEvents << " events" << endl;    
    flush( cout );
    
    MPI_Send( &nEvents, 1, MPI_INT, i, MPITag::kIntSend, MPI_COMM_WORLD );
    for( int j = 0; j < nEvents; ++j ){
      
      Kinematics* event = T::getEvent();
      
      // the pointer should never be null as long as we can count OK
      assert( event != NULL );
      fillStruct( &(kinArray[j]), event );
      delete event;
    }
    
    MPI_Send( kinArray, nEvents, MPI_KinStruct, i, MPITag::kDataSend,
              MPI_COMM_WORLD );
    
    //    cout << "Waiting for acknowledge from " << i << endl;
    //    flush( cout );
    
    // wait for acknowledgment before continuing
    MPI_Recv( &nEvents, 1, MPI_INT, i, MPITag::kAcknowledge, 
             MPI_COMM_WORLD, &status );
    
    //    cout << "Send to process " << i << " finished." << endl;
    //    flush( cout );
  }
  
  delete[] kinArray;
}

template< class T >
void DataReaderMPI<T>::receiveData()
{
  
  int nEvents;
  
  MPI_Status status;
  
  MPI_Recv( &nEvents, 1, MPI_INT, 0, MPITag::kIntSend, 
           MPI_COMM_WORLD, &status );
  
  //  cout << "Process " << m_rank << " waiting for " 
  //       << nEvents << " events." << endl;
  
  KinStruct* kinArray = new KinStruct[nEvents];
  
  MPI_Recv( kinArray, nEvents, MPI_KinStruct, 0, MPITag::kDataSend,
            MPI_COMM_WORLD, &status );
  
  for( int i = 0; i < nEvents; ++i ){
    
    m_ptrCache.push_back( createKin( &(kinArray[i]) ) );
  }
  
  cout << "Process " << m_rank << " received "
       << m_ptrCache.size() << " events." << endl;
  flush( cout );
  
  // adjust the pointer iterator to point to the beginning of the cache
  resetSource();
  m_numEvents = static_cast< unsigned int >( m_ptrCache.size() );
  
  // send acknowledgment
  MPI_Send( &nEvents, 1, MPI_INT, 0, MPITag::kAcknowledge, MPI_COMM_WORLD );
  
  delete[] kinArray;
}

template< class T >
unsigned int DataReaderMPI<T>::numEvents() const
{
  if( m_isMaster ){
    
    return T::numEvents();
  }
  else{
    
    return m_numEvents;
  }
}

template< class T >
void DataReaderMPI<T>::fillStruct( KinStruct* kinStruct, Kinematics* kin )
{
  
  const vector<TLorentzVector>& partList = kin->particleList();
  kinStruct->nPart = partList.size();
  kinStruct->weight = kin->weight();
  
  assert( partList.size() <= Kinematics::kMaxParticles );
  
  for( vector<TLorentzVector>::const_iterator vec = partList.begin();
      vec != partList.end(); ++vec ){
    
    int i = ( vec - partList.begin() );
    
    kinStruct->e[i] = (*vec).E();
    kinStruct->px[i] = (*vec).Px();
    kinStruct->py[i] = (*vec).Py();
    kinStruct->pz[i] = (*vec).Pz();
  }
}

template< class T >
Kinematics* DataReaderMPI<T>::createKin( KinStruct* kinStruct )
{
  
  vector<TLorentzVector> partList;
  
  for( int i = 0; i < kinStruct->nPart; ++i ){
    
    partList.push_back( TLorentzVector( kinStruct->px[i],
                                          kinStruct->py[i],
                                          kinStruct->pz[i],
                                          kinStruct->e[i] ) );
  }
  
  return new Kinematics( partList, kinStruct->weight );
}

template< class T > 
void DataReaderMPI<T>::defineMPIType()
{
  
  KinStruct kinStruct;
  
  // arrays used to define info about the six elements in the struct
  int length[6];
  MPI_Aint loc[6];
  MPI_Datatype type[6];
  
  MPI_Aint baseAddress;
  MPI_Address( &kinStruct, &baseAddress );
  
  length[0] = 1;
  MPI_Address( &kinStruct.nPart, &loc[0] );
  loc[0] -= baseAddress;
  type[0] = MPI_INT;
  
  length[1] = 1;
  MPI_Address( &kinStruct.weight, &loc[1] );
  loc[1] -= baseAddress;
  type[1] = MPI_FLOAT;
  
  length[2] = Kinematics::kMaxParticles;
  MPI_Address( &kinStruct.e, &loc[2] );
  loc[2] -= baseAddress;
  type[2] = MPI_FLOAT;
  
  length[3] = Kinematics::kMaxParticles;
  MPI_Address( &kinStruct.px, &loc[3] );
  loc[3] -= baseAddress;
  type[3] = MPI_FLOAT;
  
  length[4] = Kinematics::kMaxParticles;
  MPI_Address( &kinStruct.py, &loc[4] );
  loc[4] -= baseAddress;
  type[4] = MPI_FLOAT;
  
  length[5] = Kinematics::kMaxParticles;
  MPI_Address( &kinStruct.pz, &loc[5] );
  loc[5] -= baseAddress;
  type[5] = MPI_FLOAT;
  
  MPI_Type_struct( 6, length, loc, type, &MPI_KinStruct );
  MPI_Type_commit( &MPI_KinStruct );
}

#endif
