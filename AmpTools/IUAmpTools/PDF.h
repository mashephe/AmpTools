#if !defined( PDF_H )
#define PDF_H

//******************************************************************************
// This file is part of AmpTools, a package for performing PDF Analysis
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
#include <complex>
#include <vector>
#include <string>
#include <cassert>

#include "IUAmpTools/Term.h"
#include "IUAmpTools/Kinematics.h"
#include "GPUManager/GPUCustomTypes.h"

#ifdef GPU_ACCELERATION
#include "cuda_runtime.h"
#include "GPUManager/CUDA-Complex.cuh"
class GPUManager;
#endif //GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;
class AmpParameter;

/**
 * This class represents a user defined PDF.  In its most abstract
 * sense it is a mechanism to turn a set of four vectors describing an
 * event into a probability.
 *
 * \ingroup IUAmpTools
 */

class PDF : public Term
{
  
public:
  
  /**
   * The default constructor.  The user's derived class should contain a
   * default constructor that calls this constructor.
   */
  PDF() : Term() { }
  
  /**
   * This constructor takes a list of arguments to inititialize an
   * amplitude and then stores them.  The user's derived class should
   * contain a similar constructor that calls this one.
   */
  PDF( const vector< string >& args ) : Term( args ) { }
  
  /**
   * This is the destructor.
   */
  virtual ~PDF(){}
  
  /**
   * This must be overriden by the user and indicates how to convert a list
   * of strings (arguments) into a pointer to a new instance of the
   * users defined PDF.
   *
   * The user can avoid writing this method by inheriting from the
   * UserPDF class (which derives from this class).
   *
   * \param[in] args a list of string arguments that may, for example, be
   * specified in a configuration file
   *
   *  \see UserPDF
   */
  virtual PDF* newPDF( const vector< string >& args ) const = 0;
  
  /**
   * A function that the user must write that indicates how an PDF
   * can duplicate itself.  It returns a pointer to a new instance of the
   * PDF that behaves in exactly the same way as the original.
   *
   * The user can avoid writing this method by inheriting from the
   * UserPDF class (which derives from this class).
   *
   *  \see UserPDF
   */
  virtual PDF* clone() const = 0;

  // speed may be enhanced if the two functions below are combined
  // as this avoids an extra function call, but puts more complicated
  // calcPDFAll loops in user code
  
  /**
   * This loops over all events and calculates the PDF for each event.
   * It does so by calling the user-defined calcPDF routine for each
   * permutation of particles.  The function is virtual since, in principle,
   * the user may choose to override it and perform the loop and computation
   * directly and more efficiently than using a separate function call
   * for each event.
   *
   * \param[in] pdData a pointer to the array of data.  This is a long list
   * of GDoubles E, px, py, pz repeated sequentially for each particle in
   * the event and then the series repeated for each event in the data set.
   *
   * \param[out] pdAmps a pointer to a block of memory where the calculated
   * PDFs should go.  The values are stored as real and imaginary parts
   * repeated for each permutation and then for each event.
   * Total length of memory needed is 2 * sizeof( GDDouble ) * N_evts * N_perms
   *
   * \param[in] iNEvents the number of events in the data set
   *
   * \see calcPDFAll
   * \see PDFManager::addAmpPermutation
   */
  virtual void calcPDFAll( GDouble* pdData, GDouble* pdAmps, int iNEvents,
                           int iNParticles, GDouble* pdUserVars = 0 ) const;
  
  /**
   * This is the user-defined function that computes a single PDF
   * for a set of four-vectos that describe the event kinematics.  As discussed
   * above this function should be factorized as much as possible.
   *
   * \param[in] pKin a pointer to a single event.  pKin[0][0-3] define E, px,
   * py, pz for the first particle, pKin[1][0-3] for the second, and so on
   *
   */
  virtual GDouble calcPDF( GDouble** pKin ) const;
  
  
  /**
   * This is the user-defined function that computes a single PDF
   * for a set of four-vectos that describe the event kinematics.
   * The user must override this function in his or her PDF class.
   *
   * For the user to utilize user-defiend data in the amplitude calculation,
   * this function must be overridden by the derived class.  Either this function
   * or the function above must be defined for any Amplitude class.
   *
   * \param[in] pKin a pointer to a single event.  pKin[0][0-3] define E, px,
   * py, pz for the first particle, pKin[1][0-3] for the second, and so on
   *
   * \param[in] userVars is an optional pointer to the user data block associated
   * with this event.  It can be used to store intermediate portions of
   * the calculation in the case that calcPDF must be called multiple times
   * during the course of a fit.  The userVars memory block is filled in
   * calcUserVars.
   */
  virtual GDouble calcPDF( GDouble** pKin,  GDouble* userVars ) const;
  
  /**
   * \overload
   *
   * This is a user-friendly interface to calcPDF used for diagnostics.
   * It takes a pointer to a Kinematics object, converts it into a GDouble** format,
   * then calls the standard calcPDF method.
   *
   * \param[in] pKin a pointer to a Kinematics object
   *
   * \see calcPDF
   */
  
  GDouble calcPDF( const Kinematics* pKin, GDouble* userVars = 0 ) const;
  
  
#ifdef GPU_ACCELERATION
  
  /**
   * If GPU_ACCELERATION flag is set this is the member function that sets
   * the current permutation and then calls the user-defined routine
   * to launch the GPU kernel.  It is the GPU analog of calcPDFAll.
   */
  virtual void calcPDFGPU( dim3 dimGrid, dim3 dimBlock, GPU_PDF_PROTO ) const;
  
  /**
   * The user override this route and use it pass any parameters to a global
   * C function that actually launches the GPU kernel
   */
  virtual void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_PDF_PROTO ) const {
    
    cout << "\nNo GPU function for calculating " << name() << " is defined." << endl;
    assert( false );
  }
  
#endif //GPU_ACCELERATION

};


#endif
