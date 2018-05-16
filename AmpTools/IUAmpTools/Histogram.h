#if !defined(HISTOGRAM)
#define HISTOGRAM

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

#include <vector>

#include "TH1F.h"

using namespace std;

#define MAXBINS 200*200

struct HistStruct {
	int nBins,nBinsX,nBinsY;
	float xLow,xHigh;
	float yLow,yHigh;
	float entries;
	float contents[MAXBINS];
  float sumW2[MAXBINS];
};


class Histogram 
{
    
public:
  
  Histogram( string name = "Histogram", string title = "hist" );
  
	virtual ~Histogram(){};
	/*Pure virtual methods*/
	virtual void fill( vector< double > value, double weight = 1.0 )= 0;
	
	virtual TH1* toRoot() const=0;
	virtual HistStruct toStruct() const=0;
	virtual Histogram* Clone() const=0;
	
	void normalize( double scaleFactor );
	double entries(){ return( m_entries ); }
	void clear();
	void operator+=( HistStruct& hStruct );
  
  string title() const { return m_title; }
  string name()  const { return m_name;  }
  
	bool empty() const { return( m_entries == 0 ); }
  
protected:
  
  int m_nBins,m_nBinsX;
  double m_xLow;
  double m_xHigh;
  
  double m_entries;
  vector< double > m_binContents;  // this is the sum of the weights
  vector< double > m_sumWeightSq;
  double m_binSizeX;
	
  int m_dimensions;

private:
  
  string m_title;
  string m_name;
};

#endif
