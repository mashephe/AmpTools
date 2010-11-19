#if !defined(MINUITINTERFACE_MIPOINTERLISTITERATOR_H)
#define MINUITINTERFACE_MIPOINTERLISTITERATOR_H

// This file is a part of MinuitInterface - a front end for the Minuit minimization
//       package (Minuit itself was authored by Fred James, of CERN)
// 
// 
// Copyright Cornell University 1993, 1996, All Rights Reserved.
// 
// This software written by Lawrence Gibbons, Cornell University.
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
// obtained from Cornell University.
// 
// CORNELL MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  By way
// of example, but not limitation, CORNELL MAKES NO REPRESENTATIONS OR
// WARRANTIES OF MERCANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
// THE USE OF THIS SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS,
// COPYRIGHTS, TRADEMARKS, OR OTHER RIGHTS.  Cornell University shall not be
// held liable for any liability with respect to any claim by the user or any
// other party arising from use of the program.
//

typedef long relative_address;

template < class T_iter, class T_return, class T_valueType >
class MIPointerListIterator
{

   // ---------- private typedefs ---------------
   typedef MIPointerListIterator< T_iter, T_return, T_valueType > self;
   
public:
      // ---------- constants, enums and typedefs --------------
   typedef T_valueType      value_type;
   typedef relative_address difference_type;
   typedef T_return&        reference;
   typedef T_return*        pointer;
   
#if defined(NO_ITERATOR_TRAITS_BUG)
   typedef  typename iterator_traits<T_iter>::iterator_category iterator_category;
#endif
   
   // ---------- Constructors and destructor ----------------
   // default constructor
   MIPointerListIterator();
   // constructor based on the STL implemented iterator for the
   // list of T* (takes double-dereferencing to get to T)
   MIPointerListIterator( const T_iter& anIterator );
   
   // copy constructor
   MIPointerListIterator( const self& rhs );
   
   // assignment operator
   const self& operator=( const self& rhs );
   
   virtual ~MIPointerListIterator();
   
   // ---------- member functions ---------------------------
   self& operator++();
   self  operator++(int);
   
   self& operator--();
   self  operator--(int);
   
   self& operator+=( relative_address n);
   self& operator-=( relative_address n);
   
   // dereferencing operators
   T_return& operator*();
   T_return* operator->();
   
   // ---------- const member functions ---------------------
   
   // comparison operators
   bool operator==( const self& anMIPointerListIterator ) const;
   bool operator!=( const self& anMIPointerListIterator ) const;
   
   // iterator addition/subtraction
   self  operator+(  relative_address n) const;
   self  operator-(  relative_address n) const;
   
protected:
   // --------- not intended for public use -----------
   const T_iter& bareIterator() const;
   
private:
   // ---------- data members -------------------------------
   T_iter m_iterator;
   
};

// inline function definitions

//
// Package:     <MinuitInterface>
// Module:      MIPointerListIterator
// 
// Description: A smart pointer that implements an iterator over a
//              container of T*, eliminating the need for users to
//              double dereference.
//
// Implementation:
//     This just holds an STL list iterator, with the dereferencing
//     operators performing the first level of iteration, returning
//     either a T* (->) or a T& (*).
//
//
// constructors and destructor
//
template<class T_iter, class T_return, class T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>::MIPointerListIterator()
{}

template<class T_iter, class T_return, class T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>::MIPointerListIterator( const T_iter& anIterator ) :
m_iterator( anIterator )
{}

template<class T_iter, class T_return, class T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>::MIPointerListIterator( const MIPointerListIterator<T_iter, T_return, T_valueType>& rhs ) :
m_iterator( rhs.m_iterator )
{}

template<class T_iter, class T_return, class T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>::~MIPointerListIterator()
{
}

//
// assignment operators
//
template<class T_iter, class T_return, class T_valueType>
const MIPointerListIterator<T_iter, T_return, T_valueType>& 
MIPointerListIterator<T_iter, T_return, T_valueType>::operator=(const MIPointerListIterator<T_iter, T_return, T_valueType>& rhs)
{
   if( this != &rhs ) 
   {
      
      m_iterator = rhs.m_iterator;
   }
   
   return *this;
}

//
// member functions
//
template<class T_iter, class T_return, class T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>&
MIPointerListIterator<T_iter, T_return, T_valueType>::operator++()        // prefix
{
   ++m_iterator;
   return *this;
}

template<class T_iter, class T_return, class T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>::operator++(int)     // postfix
{
   MIPointerListIterator<T_iter, T_return, T_valueType> before( *this );
   ++(*this); // use prefix operator
   return( before );
}


template<class T_iter, class T_return, class T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>&
MIPointerListIterator<T_iter, T_return, T_valueType>::operator--()        // prefix
{
   --m_iterator;
   return *this;
}


template<class T_iter, class T_return, class T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>::operator--(int)     // postfix
{
   MIPointerListIterator<T_iter, T_return, T_valueType> before( *this );
   --(*this); // use prefix operator
   return ( before );
}


template<class T_iter, class T_return, class T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>&
MIPointerListIterator<T_iter, T_return, T_valueType>::operator+=( relative_address n )
{
   for ( relative_address i = 0; i != n; ++i, ++*this ) {}
   return *this;
}


template<class T_iter, class T_return, class T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>&
MIPointerListIterator<T_iter, T_return, T_valueType>::operator-=( relative_address n )
{
   for ( relative_address i = 0; i != n; ++i, --*this ) {}
   return *this;
}


//
// const member functions
//

template<class T_iter, class T_return, class T_valueType>
T_return&
MIPointerListIterator<T_iter, T_return, T_valueType>::operator*()
{
   return( **m_iterator );
}


template<class T_iter, class T_return, class T_valueType>
T_return*
MIPointerListIterator<T_iter, T_return, T_valueType>::operator->()
{
   return( *m_iterator );
}


template<class T_iter, class T_return, class T_valueType>
bool
MIPointerListIterator<T_iter, T_return, T_valueType>::operator==( const self& anMIPointerListIterator ) const
{
   return ( m_iterator == anMIPointerListIterator.m_iterator );
}


template<class T_iter, class T_return, class T_valueType>
bool
MIPointerListIterator<T_iter, T_return, T_valueType>::operator!=( const self& anItr ) const
{
   return (m_iterator != anItr.m_iterator); 
}


template<class T_iter, class T_return, class T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>::operator+( relative_address n ) const
{
   MIPointerListIterator<T_iter, T_return, T_valueType> tempItr( *this );
   for( relative_address i = 0; i != n; ++i, ++tempItr ) {}
   return tempItr;
}


template<class T_iter, class T_return, class T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>
MIPointerListIterator<T_iter, T_return, T_valueType>::operator-( relative_address n ) const
{
   MIPointerListIterator<T_iter, T_return, T_valueType> tempItr( *this );
   for( relative_address i = 0; i != n; ++i, --tempItr ) {}
   return tempItr;
}


template<class T_iter, class T_return, class T_valueType>
const T_iter&
MIPointerListIterator<T_iter, T_return, T_valueType>::bareIterator() const
{
   return( m_iterator );
}

//
// static member functions
//

#endif /* MINUITINTERFACE_MIPOINTERLISTITERATOR_H */
