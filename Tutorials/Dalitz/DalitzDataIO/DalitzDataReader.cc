
#include <vector>
#include <cassert>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "DalitzDataIO/DalitzDataReader.h"


DalitzDataReader::DalitzDataReader( const string& inFileName,
                                    const string& inTreeName ) :
                                    m_eventCounter( 0 ){

  TH1::AddDirectory( kFALSE );

  m_inFile = new TFile( inFileName.c_str() );
  m_inTree = static_cast<TTree*>( m_inFile->Get( inTreeName.c_str() ) );

  m_inTree->SetBranchAddress( "EnP1", &m_EnP1 );
  m_inTree->SetBranchAddress( "PxP1", &m_PxP1 );
  m_inTree->SetBranchAddress( "PyP1", &m_PyP1 );
  m_inTree->SetBranchAddress( "PzP1", &m_PzP1 );

  m_inTree->SetBranchAddress( "EnP2", &m_EnP2 );
  m_inTree->SetBranchAddress( "PxP2", &m_PxP2 );
  m_inTree->SetBranchAddress( "PyP2", &m_PyP2 );
  m_inTree->SetBranchAddress( "PzP2", &m_PzP2 );

  m_inTree->SetBranchAddress( "EnP3", &m_EnP3 );
  m_inTree->SetBranchAddress( "PxP3", &m_PxP3 );
  m_inTree->SetBranchAddress( "PyP3", &m_PyP3 );
  m_inTree->SetBranchAddress( "PzP3", &m_PzP3 );

}


void
DalitzDataReader::resetSource(){

  m_eventCounter = 0;

}


Kinematics*
DalitzDataReader::getEvent(){

  if( m_eventCounter < static_cast<int>( m_inTree->GetEntries() ) ){

    m_inTree->GetEntry( m_eventCounter++ );

    vector< HepLorentzVector > particleList;
    particleList.push_back( HepLorentzVector( m_PxP1, m_PyP1, m_PzP1, m_EnP1 ) );
    particleList.push_back( HepLorentzVector( m_PxP2, m_PyP2, m_PzP2, m_EnP2 ) );
    particleList.push_back( HepLorentzVector( m_PxP3, m_PyP3, m_PzP3, m_EnP3 ) );

    return new Kinematics( particleList );

  }

  else{
    return NULL;
  }

}


unsigned int
DalitzDataReader::numEvents() const{
  return static_cast< unsigned int >( m_inTree->GetEntries() );
}
