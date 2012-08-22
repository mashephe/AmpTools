
#include <vector>
#include <cassert>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "gammaKKDataIO/gammaKKDataReader.h"

gammaKKDataReader::gammaKKDataReader(const vector<string> &args) :
  UserDataReader<gammaKKDataReader>(args), m_eventCounter( 0 ){

  // The number of arguments should be
  // 1 ... file name only, read all data
  // 2 ... file name, and also min and max of M23 (very specific for this example)
  assert(args.size()==1 || args.size()==3);
  m_rangeSpecified = false;

  TH1::AddDirectory( kFALSE );

  m_inFile = new TFile(args[0].c_str());
  m_inTree = static_cast<TTree*>( m_inFile->Get("nt") );
  m_numEvents = m_inTree->GetEntries();
  

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

  if(args.size()==3){
    
    m_rangeSpecified = true;
    m_min = atof(args[1].c_str());
    m_max = atof(args[2].c_str());

    cout << "*******************************************************************" << endl;
    cout << "Constructor for gammaKKReader" << endl;
    cout << "Mass ranges: min         = " << m_min << ", max = " << m_max << endl;
    cout << "Number of events in file = " << m_inTree->GetEntries() << endl;

    m_numEvents = 0;
    while(m_eventCounter < static_cast<int>( m_inTree->GetEntries() )){
      m_inTree->GetEntry(m_eventCounter++);
      
      vector< HepLorentzVector > particleList;
      particleList.push_back( HepLorentzVector( m_PxP1, m_PyP1, m_PzP1, m_EnP1 ) );
      particleList.push_back( HepLorentzVector( m_PxP2, m_PyP2, m_PzP2, m_EnP2 ) );
      particleList.push_back( HepLorentzVector( m_PxP3, m_PyP3, m_PzP3, m_EnP3 ) );
      
      
      float m23 = (particleList[1] + particleList[2]).m();
      if(m_min <= m23 && m23 < m_max){
	m_numEvents++;
      }
    }
    cout << "Number of events kept    = " << m_numEvents << endl;
    cout << "*******************************************************************" << endl;
    resetSource();


	/*
    cout << "Start of creating cutTree.........................................." << endl;
    m_rangeSpecified = true;
    m_min = atof(args[1].c_str());
    m_max = atof(args[2].c_str());
    TString stringMass("");
    TString stringCut("");
    stringMass +=
      "sqrt((EnP2+EnP3)*(EnP2+EnP3)-(PxP2+PxP3)*(PxP2+PxP3)-(PyP2+PyP3)*(PyP2+PyP3)-(PzP2+PzP3)*(PzP2+PzP3))";
    stringCut += stringMass;  stringCut += ">";  stringCut += args[1]; stringCut += "&&";
    stringCut += stringMass;  stringCut += "<";  stringCut += args[2];
    TFile* cutFile = new TFile("cutfile.root","recreate");
    TTree* cutTree = m_inTree->CopyTree(stringCut);
    // delete m_inTree;
    m_inTree = cutTree;
    cout << "End of creating cutTree.........................................." << endl;
	*/
  }
}


void
gammaKKDataReader::resetSource(){

  m_eventCounter = 0;

}


Kinematics*
gammaKKDataReader::getEvent(){

  // If no range is specified for M23, then read all events in
  if(m_rangeSpecified == false){
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
    // If a range is specified for M23, then read in only events that
    // have M23 within that range
  }else{

    // cout << "there are " << m_inTree->GetEntries() << " events in tree" << endl;
    // Get each event, and read it in only if M23 is within the specified range
    // (We will save the value of m_eventCounter for each event that is accepted,
    //  and start looking at event m_eventCounter+1 if m_eventCounter is not 
    //  the first event)
    // cout << "-------------------------------------------------------------------------" << endl;
    while(m_eventCounter < static_cast<int>( m_inTree->GetEntries() )){

      m_inTree->GetEntry(m_eventCounter++);

      vector< HepLorentzVector > particleList;
      particleList.push_back( HepLorentzVector( m_PxP1, m_PyP1, m_PzP1, m_EnP1 ) );
      particleList.push_back( HepLorentzVector( m_PxP2, m_PyP2, m_PzP2, m_EnP2 ) );
      particleList.push_back( HepLorentzVector( m_PxP3, m_PyP3, m_PzP3, m_EnP3 ) );


      float m23 = (particleList[1] + particleList[2]).m();
      if(m_min <= m23 && m23 < m_max){
	// cout << "kept event, eventCounter = " << m_eventCounter << endl;
	// cout << "event counter: " << m_eventCounter << endl;
	return new Kinematics( particleList );
      }
      else{
	// cout << "rejected event, eventCounter = " << m_eventCounter << endl;
      }
    }
    // cout << "Got to here" << endl;
    return NULL;
  } // end of m_rangeSpecified is true

  // cout << "We are here" << endl;
  // If for some reason we get here, return null
  return NULL;

}


unsigned int
gammaKKDataReader::numEvents() const{
  return m_numEvents;
}
