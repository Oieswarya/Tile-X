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

/* QUICK SORT PART: 
Adapted from: Michael Lamont's WebSite: http://linux.wku.edu/~lamonml/ 
Source: http://linux.wku.edu/~lamonml/algor/sort/sort.html
*/

#include "utilitySortingAlgorithms.h"

using namespace std;
/* Class Member Definitions */
indexWeight::indexWeight()
{
	mIndex = 0;
	mWeight = 0.0f;
}

//Build the compressed Column format from given vectors and variables:
indexWeight::indexWeight(MilanLongInt ind, MilanReal wt)
{
	mIndex = ind;
	mWeight = wt;
}

//Matrix Destructor
indexWeight::~indexWeight()
{
	//Do not have to do anything
}

/* ****************************************** */
// Access and Information Functions

//Display Values
void indexWeight::displayIndexWeight() const
{
	cout<<mIndex<<"\t"<<mWeight<<"\n";
}

/* Function Definitions 
Adapted from: Michael Lamont's WebSite: http://linux.wku.edu/~lamonml/ 
Source: http://linux.wku.edu/~lamonml/algor/sort/sort.html
*/

//Sorts a given list in DESCENDING order.
void QuickSort(vector<indexWeight>::iterator SortThis, MilanLongInt ArraySize)
{
	QuickSortRecur(SortThis, 0, ArraySize - 1);
}

void QuickSortRecur(vector<indexWeight>::iterator SortThis, MilanLongInt LHS, MilanLongInt RHS)
{
  MilanLongInt pivoted, LHS_old, RHS_old, switchHere;
  MilanReal pivot;

  LHS_old = LHS;
  RHS_old = RHS;
  pivot = SortThis[LHS].getWeight();
  pivoted = SortThis[LHS].getIndex();
  while (LHS < RHS)
  {
    while ((SortThis[RHS].getWeight() <= pivot) && (LHS < RHS))
      RHS--;
    if (LHS != RHS)
    {
      SortThis[LHS].setWeight(SortThis[RHS].getWeight());
	  SortThis[LHS].setIndex(SortThis[RHS].getIndex());
      LHS++;
    }
    while ((SortThis[LHS].getWeight() >= pivot) && (LHS < RHS))
      LHS++;
    if (LHS != RHS)
    {
      SortThis[RHS].setWeight(SortThis[LHS].getWeight());
	  SortThis[RHS].setIndex(SortThis[LHS].getIndex());
      RHS--;
    }
  }
  SortThis[LHS].setWeight(pivot);
  SortThis[LHS].setIndex(pivoted);

  //Recursion:	
  switchHere = LHS;
  LHS = LHS_old;
  RHS = RHS_old;
  if (LHS < switchHere)
    QuickSortRecur(SortThis, LHS, switchHere-1);
  if (RHS > switchHere)
    QuickSortRecur(SortThis, switchHere+1, RHS);
}

////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// Quick Sort with Three-way Partitioning ////////////////////////////////////
//////// Source: http://www.cs.princeton.edu/~rs/talks/QuicksortIsOptimal.pdf //////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//Sorts a given list in DESCENDING order.
void QuickSort3Way(vector<indexWeight>::iterator SortThis, MilanLongInt ArraySize, 
                                                                          MilanBool descending)
{
	//cout<<"\n Within Sort Function \n";
	QuickSortRecur3Way(SortThis, 0, ArraySize - 1, descending);
	//cout<<"\n Returning back from Sorting \n";
}

void QuickSortRecur3Way(vector<indexWeight>::iterator SortThis, MilanLongInt LHS, 
                                                      MilanLongInt RHS, MilanBool descending)
{
	//Termination Condition:
	if ( RHS <= LHS )
		return;
	//Continue if not reached termination:
	//cout<<" Sorting between "<<LHS<<" - "<<RHS<<endl;
	MilanLongInt i = LHS - 1, j = RHS, p = LHS - 1, q = RHS, k;	
	MilanReal v = SortThis[RHS].getWeight();
	for ( ; ; )
	{
		//cout<<" Value of i: "<<i;
		if ( descending )
			while ( SortThis[++i].getWeight() > v ) ;   //Descending
		else
			while ( SortThis[++i].getWeight() < v ) ; //Ascending
		
		//cout<<" after: "<<i<<endl;
		//cout<<" Value of j: "<<j;
		if ( descending )
		{
			while ( v > SortThis[--j].getWeight() )    //Descending
				if ( j == LHS )
					break;
		}
		else
		{
			while ( v < SortThis[--j].getWeight() )  //Ascending
				if ( j == LHS )
					break;
		}
		//cout<<" after: "<<j<<endl;
		if ( i >= j )
			break;
		//Swap the Nodes A[i] with A[j]:
		SwapTwoNodes(SortThis, i, j);
		//cout<<" Value of p: "<<p;
		if ( SortThis[i].getWeight() == v )
		{
			p++;
			SwapTwoNodes(SortThis, p, i);
		}
		//cout<<" after: "<<p<<endl;
		//cout<<" Value of q: "<<q;
		if ( v == SortThis[j].getWeight() )
		{
			q--;
			SwapTwoNodes(SortThis, j, q);
		}
		//cout<<" after: "<<q<<endl;
	} // End of for loop
	SwapTwoNodes(SortThis, i, RHS);
	j = i - 1;
	i++;
	for ( k=LHS; k < p; k++, j-- )
		SwapTwoNodes(SortThis, k, j);
	for ( k=RHS-1; k > q; k--, i++ )
		SwapTwoNodes(SortThis, i, k);

	QuickSortRecur3Way(SortThis, LHS, j, descending);
	QuickSortRecur3Way(SortThis, i, RHS, descending);
}

inline void SwapTwoNodes(vector<indexWeight>::iterator SortThis, MilanLongInt first, MilanLongInt second)
{
	//cout<<" Swapping "<<first<<" - "<<second<<endl;
	MilanReal tempWt;
	MilanLongInt tempInd;
	tempWt = SortThis[first].getWeight();
	tempInd = SortThis[first].getIndex(); 
	SortThis[first].setWeight(SortThis[second].getWeight());
	SortThis[first].setIndex(SortThis[second].getIndex());
	SortThis[second].setWeight(tempWt);
	SortThis[second].setIndex(tempInd);
}
