
#include <vector>
#include <cassert>
#include <iostream>

#include "gammaKKDataIO/gammaKKDataWriter.h"
#include "TLorentzVector.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

using namespace std;

gammaKKDataWriter::gammaKKDataWriter( const string& outFile ){

  TH1::AddDirectory( kFALSE );

  m_outFile = new TFile( outFile.c_str(), "recreate" );
  m_outTree = new TTree( "nt", "nt" );

  m_outTree->Branch( "EnP1", &m_EnP1, "EnP1/F" );
  m_outTree->Branch( "PxP1", &m_PxP1, "PxP1/F" );
  m_outTree->Branch( "PyP1", &m_PyP1, "PyP1/F" );
  m_outTree->Branch( "PzP1", &m_PzP1, "PzP1/F" );

  m_outTree->Branch( "EnP2", &m_EnP2, "EnP2/F" );
  m_outTree->Branch( "PxP2", &m_PxP2, "PxP2/F" );
  m_outTree->Branch( "PyP2", &m_PyP2, "PyP2/F" );
  m_outTree->Branch( "PzP2", &m_PzP2, "PzP2/F" );

  m_outTree->Branch( "EnP3", &m_EnP3, "EnP3/F" );
  m_outTree->Branch( "PxP3", &m_PxP3, "PxP3/F" );
  m_outTree->Branch( "PyP3", &m_PyP3, "PyP3/F" );
  m_outTree->Branch( "PzP3", &m_PzP3, "PzP3/F" );

  m_outTree->Branch( "s12", &m_s12, "s12/F" );
  m_outTree->Branch( "s23", &m_s23, "s23/F" );
  m_outTree->Branch( "s31", &m_s31, "s31/F" );

  m_eventCounter = 0;

}

void gammaKKDataWriter::write(){
  // cout << "writing " << m_outFile->GetName() << endl;
  m_outFile->cd();
  m_outTree->Write();
  m_outFile->Close();

}

gammaKKDataWriter::~gammaKKDataWriter(){
  // cout << "deleting " << m_outFile->GetName() << endl;
  m_outFile->cd();
  m_outTree->Write();
  m_outFile->Close();

}


void
gammaKKDataWriter::writeEvent( const Kinematics& kin ){

  vector< TLorentzVector > particleList = kin.particleList();

  m_EnP1 = particleList[0].E();
  m_PxP1 = particleList[0].Px();
  m_PyP1 = particleList[0].Py();
  m_PzP1 = particleList[0].Pz();

  m_EnP2 = particleList[1].E();
  m_PxP2 = particleList[1].Px();
  m_PyP2 = particleList[1].Py();
  m_PzP2 = particleList[1].Pz();

  m_EnP3 = particleList[2].E();
  m_PxP3 = particleList[2].Px();
  m_PyP3 = particleList[2].Py();
  m_PzP3 = particleList[2].Pz();

  m_s12 = (particleList[0]+particleList[1]).M2();
  m_s23 = (particleList[1]+particleList[2]).M2();
  m_s31 = (particleList[2]+particleList[0]).M2();

  m_outTree->Fill();
    
  m_eventCounter++;

}
