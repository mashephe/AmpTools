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

#include "IUAmpTools/report.h"
#include "GPUManager/GPUCustomTypes.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <string.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace std;

static const char* kModule = "report";

static ofstream nullStream;
static bool initReportSplash = false;
static ReportLevel currentLevel = INFO;

static ReportLevel currentLevelMPI = WARNING;
static int mpiRank_AT = 0;

ostream& report( ReportLevel level ){

  if( !initReportSplash ) initReport();

#ifdef USE_MPI
  if( mpiRank_AT > 0 && level < currentLevelMPI )
    return( nullStream );
#endif
    
  if( level >= currentLevel ) return( cout );
  else return( nullStream );
}

ostream& report( ReportLevel level,
                const char* module ){

  return report( level, string( module ) ) ;
}

ostream& report( ReportLevel level,
                 const string& module)
{

  static string lastModule( "" ) ;
  
  if( module != lastModule && level >= currentLevel ) {
    
    report( level ) << endl << "[ " << module;
    
#ifdef USE_MPI
    report( level ) << " (MPI rank: " << setw(4) << mpiRank_AT << ")";
#endif
    
    report( level ) << " ]:\n\n";
    lastModule = module ;
  }
  
  return( report( level ) << "\t" );
}

void initReport(){

  nullStream.open( "/dev/null", ios::out );
  
  if( const char* envLevel = getenv( "AMPTOOLS_REPORT_LEVEL" ) ){
    
    if( strcmp( envLevel, "DEBUG" ) == 0 ) currentLevel = DEBUG;
    else if( strcmp( envLevel, "INFO" ) == 0 ) currentLevel = INFO;
    else if( strcmp( envLevel, "NOTICE" ) == 0 ) currentLevel = NOTICE;
    else if( strcmp( envLevel, "WARNING" ) == 0 ) currentLevel = WARNING;
    else if( strcmp( envLevel, "ERROR" ) == 0 ) currentLevel = ERROR;
    else report( WARNING, kModule ) << "Unknown report level:  " << envLevel << endl;
  }
  
#ifdef USE_MPI
  const char* mpiSupport = "YES";
  MPI_Comm_rank( MPI_COMM_WORLD, &mpiRank_AT );
  if( const char* envLevel = getenv( "AMPTOOLS_REPORT_LEVEL_MPI" ) ){
    
    if( strcmp( envLevel, "DEBUG" ) == 0 ) currentLevelMPI = DEBUG;
    else if( strcmp( envLevel, "INFO" ) == 0 ) currentLevelMPI = INFO;
    else if( strcmp( envLevel, "NOTICE" ) == 0 ) currentLevelMPI = NOTICE;
    else if( strcmp( envLevel, "WARNING" ) == 0 ) currentLevelMPI = WARNING;
    else if( strcmp( envLevel, "ERROR" ) == 0 ) currentLevelMPI = ERROR;
    else report( WARNING, kModule ) << "Unknown MPI report level:  " << envLevel << endl;
  }
  if( mpiRank_AT > 0 ) return;
#else
  const char* mpiSupport = "NO";
#endif

  const char* gpuSupport =
#ifdef GPU_ACCELERATION
  GPU_ARCH;
#else
  "NO";
#endif

#ifndef __ACLIC__
  
  if( getenv("AMPTOOLS_DISABLE_SPLASH") != NULL ){
    
    cout << "   "; printLine();
    cout << setw(22) << right << "|        ^         " << setw(46) << right <<  "|" << endl;
    cout << setw(22) << right << "|       / \\        " << setw(15) << right << "Version:  " << setw(25) << left << ATVERSION << setw(6) << right << "|" << endl;
    cout << setw(22) << right << "|      /---\\       " << setw(46) << right << "|" << endl;
    cout << setw(22) << right << "|     /     \\      " << setw(15) << right << "GDouble:  " << sizeof(GDouble) << setw(19) << left << " bytes" << setw(11) << right << "|" << endl;
    cout << setw(22) << right << "|    /       \\ MP  " << setw(15) << right << "MPI:  " << setw(20) << left << mpiSupport << setw(11) << right << "|" << endl;
    cout << setw(22) << right << "|     -------      " << setw(15) << right << "GPU:  " << setw(20) << left << gpuSupport << setw(11) << right << "|" << endl;
    cout << setw(22) << right << "|        |         " << setw(46) << right << "|" << endl;
    cout << setw(22) << right << "|        |         " << setw(46) << right << "doi.org/10.5281/zenodo.5039377          |" << endl;
    cout << setw(22) << right << "|        | OOLS    " << setw(46) << right << "|" << endl;
    cout << "   ";  printLine();
  }
  
#endif
  
  initReportSplash = true;
}

void printLine(){
  
  for( int i = 0; i < 65; ++i ) cout << "=";
  cout << endl;
}

