/*---------------------------------------------------------------------------*/
/*                                                                           */
/*               M I L A N - A Library of Matching Algorithms                */
/*                                                                           */
/*               Mahantesh Halappanavar <mhalappa@cs.odu.edu>                */
/*                    Department of Computer Science                         */
/*                       Old Dominion University                             */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Copyright (C) 2004 Mahantesh Halappanavar                                 */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU General Public License               */
/* as published by the Free Software Foundation; either version 2            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 59 Temple Place-Suite 330,Boston,MA 02111-1307,USA.     */
/*                                                                           */
/*---------------------------------------------------------------------------*/

#ifndef _primitiveDataType_Definition_
#define _primitiveDataType_Definition_

#include "preProcessorDirectives.h"

using namespace std;

//Regular integer:
#ifndef INTEGER_H
#define INTEGER_H
	typedef int MilanInt;
#endif

//Regular long Integer:
#ifndef LONG_INT_H
#define LONG_INT_H
	#ifdef BIT64
	typedef long int MilanLongInt;
	#else
	typedef int MilanLongInt;
	#endif
#endif

//Regular boolean
#ifndef BOOL_H
#define BOOL_H
	typedef bool MilanBool;
#endif

//Regular double and the Absolute Function:
#ifndef REAL_H
#define REAL_H
	typedef double MilanReal;
	inline MilanReal MilanAbs(MilanReal value)
	{
	  return fabs(value);
	}
#endif

//// Define the limits:
#ifndef LIMITS_H
#define LIMITS_H

//Integer Maximum and Minimum:
#define MilanIntMax INT_MAX
#define MilanIntMin INT_MIN

#ifdef BIT64
	#define MilanLongIntMax LONG_MAX
	#define MilanLongIntMin -LONG_MAX
#else
	#define MilanLongIntMax INT_MAX
	#define MilanLongIntMin -INT_MAX
#endif

//Double Maximum and Minimum:
//Note: You can alternative use INFINITY defined in math.h
//It has been my experience that this is not very portable.
//Therefore I have adopted for LDBL_MAX and LDBL_MIN as +/- infinity.

//Largest positive number: LDBL_MAX = +infinity
//Smallest positive number: LDBL_MIN 
//Smallest negative number: -LDBL_MAX = -infinity
//Largest negative number: -LDBL_MIN  (just next to zero on the other side?)

// +INFINITY
const double PLUS_INFINITY = numeric_limits<double>::infinity();
//if(numeric_limits<float>::has_infinity)
// PLUS_INFINITY=numeric_limits<float>::infinity();
//else cerr<<"infinity for float isn’t supported";

const double MINUS_INFINITY = -PLUS_INFINITY;  


//#define MilanRealMax LDBL_MAX
#define MilanRealMax PLUS_INFINITY

// -INFINITY
//Instead of assigning smallest possible positive number, assign smallest negative number
//although we only consider postive weights, just for correctness of understand.
//#define MilanRealMin -LDBL_MAX  
//#define MilanRealMin LDBL_MIN  
#define MilanRealMin MINUS_INFINITY  

//const double PLUS_INFINITY = LDBL_MAX;   //deprecated



#endif

#endif
