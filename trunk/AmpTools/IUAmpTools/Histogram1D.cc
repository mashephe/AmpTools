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

#include "IUAmpTools/Histogram.h"
#include "IUAmpTools/Histogram1D.h"

using namespace std;

Histogram1D::Histogram1D() 
{
	m_dimensions=1;
}

Histogram1D::Histogram1D( HistStruct& hist )
{
	m_binContents.resize(hist.nBins);
	m_binContents.clear();
	m_entries = 0;
	m_dimensions=1;
	m_nBinsX = hist.nBins;
	m_nBins = hist.nBins;
	m_xLow = hist.xLow;
	m_xHigh = hist.xHigh;
    
    for( int i = 0; i < m_nBins; ++i ){
        
        m_binContents[i] = hist.contents[i];
    }    
  
    m_entries = hist.entries;
}

Histogram1D::Histogram1D( int nBins, double xLow, double xHigh ) 
{
	m_nBinsX = nBins;
	m_nBins = nBins;
	m_xLow= xLow;
	m_xHigh= xHigh;
	m_entries= 0;
	m_dimensions= 1;
	m_binContents.resize(nBins);
	m_binContents.clear();
	m_binSizeX = ( xHigh - xLow ) / nBins;
}

void
Histogram1D::fill( vector < double > values, double weight ){
        
	assert (values.size()==1);	//This is a 1D histogram!
	
	double value=values.at(0);
	
	if( ( value < m_xHigh ) && ( value >= m_xLow ) ){

        m_binContents[(int)((value - m_xLow)/m_binSizeX)] += weight;
    }
    // count overflows and underflows in the number of entries
    m_entries += weight;
}




TH1* Histogram1D::toRoot( void ) const {
 
    TH1F *plot= new TH1F( "hist", "Histogram", m_nBins, m_xLow, m_xHigh );
    
    for( unsigned int i = 0; i < m_binContents.size(); ++i ){
        
        plot->SetBinContent( i+1, m_binContents[i] );
        plot->SetBinError( i+1, sqrt( m_binContents[i] ) );
    }
      
	return (TH1*)plot;
}

HistStruct
Histogram1D::toStruct( void ) const {
 
    if( m_nBins > MAXBINS ){
        
        cout << "Too many bins in histogram -- increase MAXBINS" << endl;
        assert( false );
    }
    
    HistStruct hist;
    hist.xLow = m_xLow;
    hist.xHigh = m_xHigh;
    hist.nBins = m_nBins;
    hist.nBinsX = m_nBins;
    hist.entries = m_entries;
    
    for( int i = 0; i < m_nBins; ++i ){
        
        hist.contents[i] = static_cast< float >( m_binContents[i] );
    }
    
    return hist;
}

Histogram* Histogram1D::Clone() const{
	Histogram* m_histo;
	m_histo=(Histogram*)(new Histogram1D(*this));
	return m_histo;
}