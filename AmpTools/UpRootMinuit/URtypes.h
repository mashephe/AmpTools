/* @(#)root/base:$Name:  $:$Id: URtypes.h,v 1.1.1.1 2006/11/15 18:01:50 mashephe Exp $ */

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef UPROOT_Rtypes
#define UPROOT_Rtypes

#include "UpRootMinuit/URConfig.h"

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Rtypes                                                               //
//                                                                      //
// Basic types used by ROOT.                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <cstddef>


//---- types -------------------------------------------------------------------

typedef char           Char_urt;      //Signed Character 1 byte (char)
typedef unsigned char  UChar_urt;     //Unsigned Character 1 byte (unsigned char)
typedef short          Short_urt;     //Signed Short integer 2 bytes (short)
typedef unsigned short UShort_urt;    //Unsigned Short integer 2 bytes (unsigned short)
#ifdef R__INT16
typedef long           Int_urt;       //Signed integer 4 bytes
typedef unsigned long  UInt_urt;      //Unsigned integer 4 bytes
#else
typedef int            Int_urt;       //Signed integer 4 bytes (int)
typedef unsigned int   UInt_urt;      //Unsigned integer 4 bytes (unsigned int)
#endif
#ifdef R__B64    // Note: Long_urt and ULong_urt are currently not portable types
typedef int            Seek_urt;      //File pointer (int)
typedef long           Long_urt;      //Signed long integer 8 bytes (long)
typedef unsigned long  ULong_urt;     //Unsigned long integer 8 bytes (unsigned long)
#else
typedef int            Seek_urt;      //File pointer (int)
typedef long           Long_urt;      //Signed long integer 4 bytes (long)
typedef unsigned long  ULong_urt;     //Unsigned long integer 4 bytes (unsigned long)
#endif
typedef float          Float_urt;     //Float 4 bytes (float)
typedef double         Double_urt;    //Float 8 bytes (double)
typedef char           Text_urt;      //General string (char)
typedef bool           Bool_urt;      //Boolean (0=false, 1=true) (bool)
typedef unsigned char  Byte_urt;      //Byte (8 bits) (unsigned char)
typedef short          Version_urt;   //Class version identifier (short)
typedef const char     Option_urt;    //Option string (const char)
typedef int            Ssiz_urt;      //String size (int)
typedef float          Real_urt;      //TVector and TMatrix element type (float)
#if defined(R__WIN32) && !defined(__CINT__)
typedef __int64            Long64_urt;  //Portable signed long integer 8 bytes
typedef unsigned __int64   ULong64_urt; //Portable unsigned long integer 8 bytes
#else
typedef long long          Long64_urt;  //Portable signed long integer 8 bytes
typedef unsigned long long ULong64_urt; //Portable unsigned long integer 8 bytes
#endif

typedef void         (*VoidFuncPtr_urt)();  //pointer to void function


//---- constants ---------------------------------------------------------------

#ifndef NULL
#define NULL 0
#endif

const Bool_urt kurTRUE   = 1;
const Bool_urt kurFALSE  = 0;

const Int_urt     kurMaxUShort   = 65534;
const Int_urt     kurMaxShort    = kurMaxUShort >> 1;
const UInt_urt    kurMaxUInt     = ~0;
const Int_urt     kurMaxInt      = Int_urt(kurMaxUInt >> 1);
const ULong64_urt kurMaxULong64  = ~R__LL(0); 
const Long64_urt  kurMaxLong64   = Long64_urt(kurMaxULong64 >> 1);

const std::size_t    kurBitsPerByte = 8;
const Ssiz_urt    kurNPOS        = ~(Ssiz_urt)0;

#endif
