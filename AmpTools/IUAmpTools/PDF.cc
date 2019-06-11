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

#include <map>
#include <string>
#include <cassert>
#include <iostream>
#include <sstream>

#include "IUAmpTools/PDF.h"
#include "IUAmpTools/Kinematics.h"

#ifdef VTRACE
#include "vt_user.h"
#endif

void
PDF::calcPDFAll( GDouble* pdData, GDouble* pdAmps, int iNEvents, int iNParticles,
                 GDouble* pdUserVars ) const
{
  
#ifdef VTRACE
  string info = name();
  info += "::calcPDFAll";
  VT_TRACER( info.c_str() );
#endif
  
  unsigned int numVars = numUserVars();
  
  GDouble** pKin = new GDouble*[iNParticles];
  
  int i, iEvent;
  for( iEvent=0; iEvent<iNEvents; iEvent++ ){
    
    for( i = 0; i < iNParticles; i++ ){
      
      pKin[i] = &(pdData[4*iNParticles*iEvent]);
    }

    unsigned int userIndex = iNEvents*numVars + iEvent*numVars;
    
    if( numVars != 0 ){

      pdAmps[iEvent] = calcPDF( pKin, &(pdUserVars[userIndex]) );
    }
    else{
      
      pdAmps[iEvent] = calcPDF( pKin );
    }
  }
  
  delete[] pKin;
}

GDouble
PDF::calcPDF( GDouble** pKin, GDouble* userVars ) const {
  
  cout
  << "***********************************************************\n"
  << "ERROR in the construction of the class that defines\n"
  << "the PDF named " << name() << ".\n"
  << "One of the following two cases result in this error.\n\n"
  << "(1) The numUserVars() method of the class indicates that\n"
  << "    at least one user-defined variable will be calculated,\n"
  << "    but the calcPDF method hasn't been defined such\n"
  << "    that it can accept a pointer to the user-defined data\n"
  << "    block.  Please define the function:\n"
  << "      " << name() << "::\n"
  << "         calcPDF( GDouble** pKin, GDouble* userVars )\n\n"
  << "(2) No calcPDF function (with or without user data\n"
  << "    data pointer is defined in the class.\n"
  << "***********************************************************\n"
  << endl;
  
  assert( false );
}

GDouble
PDF::calcPDF( GDouble** pKin ) const {
  
  // It is possible to end up here if the user has
  // defined calcPDF such that it takes two
  // arguments and the number of user variables
  // to calculate is zero. (This is the else clause
  // in the next method.)  In this case try to call
  // the user's calcAmplitude function by passing
  // in a NULL pointer to the user data block.
  // If that isn't defined either then the error
  // above will print and the program will exit.
  
  return calcPDF( pKin, NULL );
}

GDouble
PDF::calcPDF( const Kinematics* pKin, GDouble* userVars ) const {
  
#ifdef VTRACE
  string info = name();
  info += "::calcPDF";
  VT_TRACER( info.c_str() );
#endif
  
  vector<TLorentzVector> particleList = pKin->particleList();
  GDouble** pData = new GDouble*[particleList.size()];
  
  for (int i = 0; i < particleList.size(); i++){
    pData[i] = new GDouble[4];
    pData[i][0] = particleList[i].E();
    pData[i][1] = particleList[i].Px();
    pData[i][2] = particleList[i].Py();
    pData[i][3] = particleList[i].Pz();
  }
  
  GDouble value;
  
  if( userVars != NULL ){
  
    value = calcPDF( pData, userVars );
  }
  else{
    
    value = calcPDF( pData );
  }
  
  for (int i = 0; i < particleList.size(); i++){
    delete[] pData[i];
  }
  delete[] pData;
  
  return value;
  
}

#ifdef GPU_ACCELERATION
void
PDF::calcPDFGPU( dim3 dimGrid, dim3 dimBlock, GPU_PDF_PROTO ) const {

#ifdef VTRACE
  string info = name();
  info += "::calcPDFGPU";
  VT_TRACER( info.c_str() );
#endif
  
  launchGPUKernel( dimGrid, dimBlock, GPU_PDF_ARGS );
}
#endif

