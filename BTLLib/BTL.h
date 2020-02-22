//
// BTL.h
//
// It is included by most of the library and is used to control the behaviour
// of the library as a whole, such as the floating point precision and whether
// the library is in debug mode.
//
// This code is part of the Bioinformatics Template Library (BTL).
//
// Copyright (C) 1997-1999 Birkbeck College, Malet Street, London, U.K.
// Copyright (C) 2005 University College, Gower Street, London, U.K.
//
// This library is free software; you can redistribute it and/or modify it under the 
// terms of the GNU Library General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option) 
// any later version.  This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public License for 
// more details.You should have received a copy of the GNU Library General Public 
// License along with this library; if not, write to the Free Software Foundation,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
///////////////////////////////////////////////////////////////////////////
#if !defined(BTL_BTL_H)
#define BTL_BTL_H 1

// The floating point precision of the class library as a whole is controlled
// by this typedef statement.

typedef double BTL_REAL;
typedef unsigned int BTL_INT;

// If the DEBUG symbol is defined then the library will carry out extra error
// checking such as some array bounds checking. This checking will obviously
// make the library more robust but will also increase excecution times. Simply
// comment this out if the extra checking is not required.

#define DEBUG_VERSION

// These macros can be used for simple debugging and run time error messages.
// They simply print the line and file where it is called and the value of the
// given variable. In the case of FATAL_ERROR and WARNING this variable should
// be a message giving some information concerning the problem. The FATAL_ERROR
// macro causes the immediate termination of the program and should only be used
// for serious problems.
//

#define BUGLINE  "\nin File: " << __FILE__ << " at line: " << __LINE__ << endl

#define DEBUG(x) {       cerr << "DEBUG    : " << (x) << BUGLINE; }
#define FATAL_ERROR(x) { cerr << "!!ERROR  : " << (x) << BUGLINE; exit(1); }
#define WARNING(x) {     cerr << "!WARNING : " << (x) << BUGLINE; }

// This define is used to control the usage of the keyword BTL_TYPENAME 
//  which is not implemented is some compilers. 

#define BTL_TYPENAME typename

// Control the use of namespaces

#define _BTL_BEGIN_NAMESPACE namespace btl {
#define _BTL_END_NAMESPACE }
//#define _BTL_BEGIN_NAMESPACE 
//#define _BTL_END_NAMESPACE 


#endif
