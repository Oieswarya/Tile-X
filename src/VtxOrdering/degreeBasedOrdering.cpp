/*---------------------------------------------------------------------------*/
/*                                                                           */
/*                          Mahantesh Halappanavar                           */
/*                        High Performance Computing                         */
/*                Pacific Northwest National Lab, Richland, WA               */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Copyright (C) 2010 Mahantesh Halappanavar                                 */
/*                                                                           */
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

#include "defs.h"
#include "vtxOrdering.h"
#include "utilitySortingAlgorithms.h"

using namespace std;

// Perform ordering based on degree: the ordering for duplicates will be retained from the sorting algorithm
// Returns the old2New Pointer such that the location will indicate the new id of an old vertex
//        old id is implicit in the array location
void degreeBasedOrdering(graph *G, long *pOrder, int nThreads) {
    printf("Within degreeBasedOrdering() \n");
    
    double time1=0, time2=0, totalTime=0;
    long    NV        = G->numVertices;
    long    NE        = G->numEdges;
    long    *vtxPtr   = G->edgeListPtrs;
    edge    *vtxInd   = G->edgeList;
    printf("Vertices:%ld  Edges:%ld\n", NV, NE);
    
    //STEP 1: Sort the vertices in order of their degree
    //Compute the degree of each vertex:
    time1 = omp_get_wtime();
    //Vector to sort the newly discovered vertices
    //To make it Minimization problem, substract the degree from NR, the max degree possible.
    vector<indexWeight> sortedWeightList; //Data structure to link index to weights
    MilanBool descending = false; //To arrange the weights in the non-decreasing order of degree
    try {
        sortedWeightList.reserve(NV);
    } catch ( length_error ) {
        cerr<<"Within Function: degreeBasedOrdering() \n";
        cerr<<"Error: Not enough memory to allocate for sortedWeightList \n";
        exit(1);
    }
    for (long i=0; i<NV; i++) {
        indexWeight newNode( i, double(vtxPtr[i+1] - vtxPtr[i]) ); //Add for sorting
        sortedWeightList.push_back(newNode);
    }
    //Now sort the queue:
    QuickSort3Way(sortedWeightList.begin(), NV,  descending);
    time2 = omp_get_wtime();
    printf("Time for sorting: %lf\n",NV, time2-time1);
    totalTime += time2-time1;
    
    //STEP 2: Reorder the vertices based on sorting
    time1 = omp_get_wtime();
    long *R = (long *) malloc (NV * sizeof(long)); assert(R != 0);
    //Add vertices in a sorted order
    //printf("***********************************************\n");
    for (long i=0; i < NV; i++ ) { //O(n)
        R[i] = sortedWeightList[i].getIndex(); //Add it to the ranked list
        //printf("%ld -- %lf \n", sortedWeightList[i].getIndex(), sortedWeightList[i].getWeight());
    }
    //printf("***********************************************\n");
    sortedWeightList.clear(); //Clear up the space for next frontier
    
    //STEP 3: Build the old2NewMap
    for (long i=0; i<NV; i++) {
        pOrder[R[i]] = i; //old2New index mapping
    }
    time2 = omp_get_wtime();
    printf("Time for reordering: %lf\n", NV, time2-time1);
    totalTime += time2-time1;
    
    printf("***********************************************\n");
    printf("Total Time  : %lf sec\n", totalTime);
    printf("***********************************************\n");
    
    //Cleanup:
    free(R);
    
} //End of degreeBasedOrdering()
