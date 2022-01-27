
#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "DalitzDataIO/DalitzDataReader.h"
#include "TSystem.h"

using namespace std;

DalitzDataReader::DalitzDataReader( const vector< string >& args ) :
UserDataReader< DalitzDataReader >(args),
m_eventCounter( 0 ){
  
  assert(args.size() == 1);
  string inFileName(args[0]);
  string inTreeName("nt");

  TH1::AddDirectory( kFALSE );

  gSystem->Load( "libTree" );

  ifstream fileexists( inFileName.c_str() );
  if (fileexists){
    m_inFile = new TFile( inFileName.c_str() );
    m_inTree = static_cast<TTree*>( m_inFile->Get( inTreeName.c_str() ) );
  }
  else{
    cout << "DalitzDataReader WARNING:  Cannot find file... " << inFileName << endl;
    m_inFile = NULL;
    m_inTree = NULL;
  }

  if (m_inTree){

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

    m_inTree->SetBranchAddress( "weight", &m_weight );
  }

}


void
DalitzDataReader::resetSource(){

  m_eventCounter = 0;

}


Kinematics*
DalitzDataReader::getEvent(){

  if( m_eventCounter < numEvents() ){

    m_inTree->GetEntry( m_eventCounter++ );

    vector< TLorentzVector > particleList;
    particleList.push_back( TLorentzVector( m_PxP1, m_PyP1, m_PzP1, m_EnP1 ) );
    particleList.push_back( TLorentzVector( m_PxP2, m_PyP2, m_PzP2, m_EnP2 ) );
    particleList.push_back( TLorentzVector( m_PxP3, m_PyP3, m_PzP3, m_EnP3 ) );

    return new Kinematics( particleList, m_weight );

  }

  else{
    return NULL;
  }

}


unsigned int
DalitzDataReader::numEvents() const{
  if (!m_inTree) return 0;
  return static_cast< unsigned int >( m_inTree->GetEntries() );
}

