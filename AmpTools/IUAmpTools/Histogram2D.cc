//******************************************************************************
// This file is part of AmpTools, a package for performing Amplitude Analysis
//
// Copyright Trustees of Indiana University 2010, all rights reserved
//
// This software written by Matthew Shepherd, Ryan Mitchell, and
//                  Hrayr Matevosyan at Indiana University, Bloomington
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice and author attribution, this list of conditions and the
//    following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice and author attribution, this list of conditions and the
//    following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// 3. Neither the name of the University nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// Creation of derivative forms of this software for commercial
// utilization may be subject to restriction; written permission may be
// obtained from the Trustees of Indiana University.
//
// INDIANA UNIVERSITY AND THE AUTHORS MAKE NO REPRESENTATIONS OR WARRANTIES,
// EXPRESS OR IMPLIED.  By way of example, but not limitation, INDIANA
// UNIVERSITY MAKES NO REPRESENTATIONS OR WARRANTIES OF MERCANTABILITY OR
// FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THIS SOFTWARE OR
// DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS,
// OR OTHER RIGHTS.  Neither Indiana University nor the authors shall be
// held liable for any liability with respect to any claim by the user or
// any other party arising from use of the program.
//******************************************************************************

#include <iostream>
#include <math.h>
#include <assert.h>

#include "IUAmpTools/Histogram2D.h"

using namespace std;

Histogram2D::Histogram2D() :
m_nBinsY( 0 ),
m_yLow( 0 ),
m_yHigh( 0 )
{
  m_dimensions=2;
}

Histogram2D::Histogram2D( HistStruct& hist )
{
  
  m_nBinsX = hist.nBinsX;
  m_nBinsY = hist.nBinsY;
  m_xLow = hist.xLow;
  m_xHigh = hist.xHigh;
  m_nBinsY = hist.nBinsY;
  m_yLow = hist.yLow;
  m_yHigh = hist.yHigh;
  m_nBins = hist.nBins;
  m_entries = hist.entries;
  m_dimensions = 2;
  m_binContents.resize(m_nBins);
  m_sumWeightSq.resize(m_nBins);
  
  for( int i = 0; i < m_nBins; ++i ){
    
    m_binContents[i] = hist.contents[i];
    m_sumWeightSq[i] = hist.sumW2[i];
  }
  
}

Histogram2D::Histogram2D( int nBinsX, double xLow, double xHigh,
                          int nBinsY, double yLow, double yHigh,
                          string name, string title ) :
Histogram( name, title )
{
  m_xLow=xLow;
  m_xHigh=xHigh;
  m_yLow=yLow;
  m_yHigh=yHigh;
  m_nBinsX=nBinsX;
  m_nBinsY=nBinsY;
  m_nBins=nBinsX*nBinsY;
  m_dimensions = 2;
  m_entries=0;
  m_binContents.resize( nBinsX*nBinsY, 0 );
  m_sumWeightSq.resize( nBinsX*nBinsY, 0 );
  m_binSizeX = ( xHigh - xLow ) / nBinsX;
  m_binSizeY = ( yHigh - yLow ) / nBinsY;
}

void
Histogram2D::fill(vector < double > values, double weight ){
  
  assert (values.size()==2);
  double valueX=values[0];
  double valueY=values[1];

  if( ( valueX < m_xHigh ) && ( valueX >= m_xLow ) && ( valueY < m_yHigh ) && ( valueY >= m_yLow ) ){

    int ibinX=(int)((valueX - m_xLow)/m_binSizeX);
    int ibinY=(int)((valueY - m_yLow)/m_binSizeY);
    int ibin=ibinY*m_nBinsX+ibinX; //a-la ROOT

    m_binContents[ibin] += weight;
    m_sumWeightSq[ibin] += weight*weight;
  }
  // count overflows and underflows in the number of entries
  m_entries += weight;
}




TH1* Histogram2D::toRoot( void ) const {
  
  TH2F *plot2D=new TH2F( name().c_str(), title().c_str(),
                         m_nBinsX, m_xLow, m_xHigh,m_nBinsY, m_yLow, m_yHigh );
  int ibinX,ibinY;
  for( unsigned int i = 0; i < m_binContents.size(); ++i ){
    //  int ibin=ibinY*m_nBinsX+ibinX; I used this. So:
    ibinX=i % m_nBinsX;
    ibinY=i / m_nBinsX;
    plot2D->SetBinContent( ibinX+1,ibinY+1, m_binContents[i] );       //+1 since ROOT considers the bin 0-th as underflow.
    plot2D->SetBinError( ibinX+1,ibinY+1, sqrt( m_sumWeightSq[i] ) );
  }
		return plot2D;
}

HistStruct
Histogram2D::toStruct( void ) const {
  
  if( m_nBins > MAXBINS ){
    
    cout << "Too many bins in histogram -- increase MAXBINS" << endl;
    assert( false );
  }
  
  HistStruct hist2D;
  hist2D.xLow = m_xLow;
  hist2D.xHigh = m_xHigh;
  hist2D.nBinsX = m_nBinsX;
  
  hist2D.yLow = m_yLow;
  hist2D.yHigh = m_yHigh;
  hist2D.nBinsY = m_nBinsY;
  
  hist2D.nBins = m_nBins;
  
  hist2D.entries = m_entries;
  
  for( int i = 0; i < m_nBins; ++i ){
    
    hist2D.contents[i] = static_cast< float >( m_binContents[i] );
    hist2D.sumW2[i] = static_cast< float >( m_sumWeightSq[i] );
  }
  
  return hist2D;
}


Histogram* Histogram2D::Clone() const{
  return new Histogram2D(*this);
}







