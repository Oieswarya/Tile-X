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

#ifndef _sorting_algorithms_
#define _sorting_algorithms_

#include "defs.h"
//#include "primitiveDataTypeDefinitions.h"

typedef long MilanLongInt;
typedef bool MilanBool;
typedef double MilanReal;

using namespace std;
/* Class Definitions */
class indexWeight;

class indexWeight
{
private:
	MilanLongInt mIndex; //The mIndex
	MilanReal mWeight ;   //The weight based on which to sort

	//indexWeight( const indexWeight& src);  
	//indexWeight& operator=(const indexWeight& rhs);

public:
	//Default Constructor:
	indexWeight();

	//Copy Constructor:
	//indexWeight( indexWeight& iw) 
	//{mIndex = iw->getIndex(); mWeight  = iw->getWeight(); }
	
	//Build the compressed Column format from given vectors and variables:
	indexWeight(MilanLongInt ind, MilanReal wt);
	
	//Matrix Destructor
	~indexWeight();
	
	/* ****************************************** */
	// Access and Information Functions
	
	//Get the Value:
	MilanLongInt getIndex() { return mIndex; }
	MilanReal getWeight() { return mWeight ; }
		
	//Set the Value:
	void setIndex(MilanLongInt i) { mIndex = i; }
	void setWeight(MilanReal w) { mWeight  = w; }

	//Display Values
	void displayIndexWeight() const;
	
	//Overloaded operator << to print matrix elements to the standard output
	// in the format: <i> <j> <value>
	//friend ostream& operator<<(stream &os, const compressedColumn &Cc);
};  //End of Class declaration.

/* Function Declarations */
///////////////////////////////////////////////////////////////////////////////////////////
//Quick Sort Algorithm:
void QuickSort(vector<indexWeight>::iterator SortThis, MilanLongInt ArraySize);
//Recursive Call function for Quick Sort:
void QuickSortRecur(vector<indexWeight>::iterator SortThis, MilanLongInt LHS, MilanLongInt RHS);

///////////////////////////////////////////////////////////////////////////////////////////
void QuickSort3Way(vector<indexWeight>::iterator SortThis, MilanLongInt ArraySize, MilanBool descending);
void QuickSortRecur3Way(vector<indexWeight>::iterator SortThis, MilanLongInt LHS, MilanLongInt RHS, MilanBool descending);

//Inline Function to Swap two nodes
inline void SwapTwoNodes(vector<indexWeight>::iterator SortThis, MilanLongInt first, MilanLongInt second);

#endif
