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
#include "input_output.h"
#include "basic_comm.h"
#include "basic_util.h"
#include "utilityClusteringFunctions.h"
#include "color_comm.h"
#include "sync_comm.h"
#include "vtxOrdering.h"
#include <dataStructureHeap.h>

using namespace std;

// Perform ordering based on the output of Grappolo
// Returns the old2New Pointer such that the location will indicate the new id of an old vertex
//        old id is implicit in the array location
void clusteringBasedOrderingSimple(graph *G, long *pOrder, int nThreads, int coloring) {
    long NV = G->numVertices;
    long NS = G->sVertices;
    long NT = NV - NS;
    
    double time1, time2;
    //Retain the original copy of the graph:
    graph* G_original = (graph *) malloc (sizeof(graph)); //The original version of the graph
    time1 = omp_get_wtime();
    duplicateGivenGraph(G, G_original);
    time2 = omp_get_wtime();
    printf("Time to duplicate : %lf\n", time2-time1);
    
    //STEP 1: Perform clustering using Grappolo:
    long *C = (long *) malloc (NV * sizeof(long)); assert(C != 0);
#pragma omp parallel for
    for (long i=0; i<NV; i++) {
        C[i] = -1;
    }
    int threadsOpt = 0;
    if(coloring > 0){
        runMultiPhaseColoring(G_original, C, coloring, 32, 0, 10000, 1000, 0.0000001, 0.0000001, nThreads, threadsOpt);
    } else {
        runMultiPhaseBasic(G_original, C, 0, 10000, 1000, 0.0000001, 0.0000001, nThreads, threadsOpt);
    }
    
    //STEP 2: Now build a compact representation of clusters:
    //        Create a CSR-like datastructure for communities in C
    //Compute number of communities:
    //Assume zero is a valid community id
    long nC=-1;
    bool isZero = false;
    bool isNegative = false;
    for(long i = 0; i < NV; i++) {
        if(C[i] == 0)
            isZero = true; //Check if zero is a valid community
        if(C[i] < 0)
            isNegative = true; //Check if zero is a valid community
        if (C[i] > nC) {
            nC = C[i];
        }
    }
    printf("Largest community id observed: %ld\n", nC);
    if(isZero) {
        printf("Zero is a valid community id\n");
        nC++;
    }
    if(isNegative) {
        printf("Some vertices have not been assigned communities\n");
        nC++; //Place to store all the unassigned vertices
    }
    assert(nC>0);
    printf("Number of unique communities in C= %d\n", nC);
    
    long * commPtr = (long *) malloc ((nC+1) * sizeof(long)); assert(commPtr != 0);
    long * commIndex = (long *) malloc (NV * sizeof(long)); assert(commIndex != 0);
    long * commAdded = (long *) malloc (nC * sizeof(long)); assert(commAdded != 0);
    
    // Initialization
#pragma omp parallel for
    for (long i=0; i<NV; i++) { //initialize
        commIndex[i] = -1;
    }
#pragma omp parallel for
    for(long i = 0; i < nC; i++) {
        commPtr[i] = 0;
        commAdded[i] = 0;
    }
    commPtr[nC] = 0;
    // Count the size of each community
#pragma omp parallel for
    for(long i = 0; i < NV; i++) {
        if(C[i] < 0) { //A negative value
            __sync_fetch_and_add(&commPtr[nC],1); //Unassigned vertices
        } else { //A positive value
            if(isZero)
                __sync_fetch_and_add(&commPtr[C[i]+1],1); //Zero-based indexing
            else
                __sync_fetch_and_add(&commPtr[C[i]],1); //One-based indexing
        }
    }//End of for(i)
    //Prefix sum:
    for(long i=0; i<nC; i++) {
        commPtr[i+1] += commPtr[i];
    }
    //STEP 3: Group vertices with the same community in an order
#pragma omp parallel for
    for (long i=0; i<NV; i++) {
        long tc = (long)C[i];
        if(tc < 0) { //A negative value
            tc = nC-1;
            long Where = commPtr[tc] + __sync_fetch_and_add(&(commAdded[tc]), 1);
            assert(Where < NV);
            commIndex[Where] = i; //The vertex id
        } else {
            if(!isZero)
                tc--; //Convert to zero based index
            long Where = commPtr[tc] + __sync_fetch_and_add(&(commAdded[tc]), 1);
            assert(Where < NV);
            commIndex[Where] = i; //The vertex id
        }
    }
    
    //STEP 4: Create the old2New map:
#pragma omp parallel for
    for (long i=0; i<NV; i++) {
        pOrder[commIndex[i]] = i;
        //printf("(%ld) --> (%ld)\n", commIndex[i]+1, i+1);
    }
    
    //Cleanup:
    free(C);
    free(commIndex);
    free(commPtr);
    free(commAdded);
} //End of clusteringBasedOrderingSimple()


// Perform ordering combining ideas from clustering (Grappolo) and RCM
// Returns the old2New Pointer such that the location will indicate the new id of an old vertex
//        old id is implicit in the array location
void clusteringBasedOrderingWithRcm(graph *G, long *pOrder, int nThreads, int coloring) {
    long NV = G->numVertices;
    long NS = G->sVertices;
    long NT = NV - NS;
    
    double time1, time2;
    //Retain the original copy of the graph:
    graph* G_original = (graph *) malloc (sizeof(graph)); //The original version of the graph
    time1 = omp_get_wtime();
    duplicateGivenGraph(G, G_original);
    time2 = omp_get_wtime();
    printf("Time to duplicate : %lf\n", time2-time1);
    
    //STEP 1: Perform clustering using Grappolo:
    long *C = (long *) malloc (NV * sizeof(long)); assert(C != 0);
#pragma omp parallel for
    for (long i=0; i<NV; i++) {
        C[i] = -1;
    }
    int threadsOpt = 0;
    if(coloring > 0){
        runMultiPhaseColoring(G_original, C, coloring, 32, 0, 10000, 1000, 0.0000001, 0.0000001, nThreads, threadsOpt);
    } else {
        runMultiPhaseBasic(G_original, C, 0, 10000, 1000, 0.0000001, 0.0000001, nThreads, threadsOpt);
    }
    
    //Step 2: Build the coarsened graph to enable RCM on it
    //Build a new graph from clusters:
    //Renumber the clusters contiguiously
    long numClusters = renumberClustersContiguously(C, G->numVertices);
    printf("Number of unique clusters: %ld\n", numClusters);
    graph *Gnew = (graph *) malloc (sizeof(graph)); assert(Gnew != 0);
    //double tmpTime =  buildNextLevelGraphOpt(G, Gnew, C, numClusters, nThreads);
    double tmpTime = omp_get_wtime();
    buildCoarsenedGraph(G, Gnew, C, numClusters, nThreads);
    tmpTime = omp_get_wtime() - tmpTime;
    printf("Time to build coarsened graph: %lf\n", tmpTime);
    displayGraphCharacteristics(Gnew);
    displayGraph(Gnew);
    
    char outGorder[256];
    sprintf(outGorder,"TmpGraph.edges"); //Output file name
    writeGraphGorderEdgeList(Gnew, outGorder);
    
    //Step 3: Perform RCM:
    // Datastructures to store old2New map
    long NVnew = Gnew->numVertices;
    long *old2NewMap = (long *) malloc (NVnew * sizeof(long)); assert(old2NewMap != 0);
    //Initialize the Vectors:
#pragma omp parallel for
    for (long i=0; i<NVnew; i++) {
        old2NewMap[i] = -1; //Initialize the rank as -1
    }
    //Call the RCM algorithm:
    algoReverseCuthillMcKeeStrictGraph(Gnew, old2NewMap, nThreads);
    //algoReverseCuthillMcKeeGraph(Gnew, old2NewMap, nThreads);
    
    //Clean up:
    if(Gnew != 0) {
        free(Gnew->edgeListPtrs);
        free(Gnew->edgeList);
        free(Gnew);
    }
    
    //Step 4: Renumber the cluster ids to match the RCM order:
    for (long i=0; i<numClusters; i++) {
        if (C[i] >= 0) //Ignore negative community numbers
            C[i] = old2NewMap[C[i]];
    }
    printf("Reordered cluster ids based on RCM order\n");
    
    //STEP 5: Now build a compact representation of clusters:
    //        Create a CSR-like datastructure for communities in C
    //First Compute number of communities --- Assume zero is a valid community id
    long nC=-1;
    bool isZero = false;
    bool isNegative = false;
    for(long i = 0; i < NV; i++) {
        if(C[i] == 0)
            isZero = true; //Check if zero is a valid community
        if(C[i] < 0)
            isNegative = true; //Check if zero is a valid community
        if (C[i] > nC) {
            nC = C[i];
        }
    }
    printf("Largest community id observed: %ld\n", nC);
    if(isZero) {
        printf("Zero is a valid community id\n");
        nC++;
    }
    if(isNegative) {
        printf("Some vertices have not been assigned communities\n");
        nC++; //Place to store all the unassigned vertices
    }
    assert(nC>0);
    printf("Number of unique communities in C= %d\n", nC);
    if(nC != NVnew) {
        printf("Number of communities from buildNextLevelGraphOpt is different (%ld - %ld)\n", nC, NVnew);
    }
    
    long * commPtr = (long *) malloc ((nC+1) * sizeof(long)); assert(commPtr != 0);
    long * commIndex = (long *) malloc (NV * sizeof(long)); assert(commIndex != 0);
    long * commAdded = (long *) malloc (nC * sizeof(long)); assert(commAdded != 0);
    
    // Initialization
#pragma omp parallel for
    for (long i=0; i<NV; i++) { //initialize
        commIndex[i] = -1;
    }
#pragma omp parallel for
    for(long i = 0; i < nC; i++) {
        commPtr[i] = 0;
        commAdded[i] = 0;
    }
    commPtr[nC] = 0;
    // Count the size of each community
#pragma omp parallel for
    for(long i = 0; i < NV; i++) {
        if(C[i] < 0) { //A negative value
            __sync_fetch_and_add(&commPtr[nC],1); //Unassigned vertices
        } else { //A positive value
            if(isZero)
                __sync_fetch_and_add(&commPtr[C[i]+1],1); //Zero-based indexing
            else
                __sync_fetch_and_add(&commPtr[C[i]],1); //One-based indexing
        }
    }//End of for(i)
    //Prefix sum:
    for(long i=0; i<nC; i++) {
        commPtr[i+1] += commPtr[i];
    }
    //STEP 3: Group vertices with the same community in an order
#pragma omp parallel for
    for (long i=0; i<NV; i++) {
        long tc = (long)C[i];
        if(tc < 0) { //A negative value
            tc = nC-1;
            long Where = commPtr[tc] + __sync_fetch_and_add(&(commAdded[tc]), 1);
            assert(Where < NV);
            commIndex[Where] = i; //The vertex id
        } else {
            if(!isZero)
                tc--; //Convert to zero based index
            long Where = commPtr[tc] + __sync_fetch_and_add(&(commAdded[tc]), 1);
            assert(Where < NV);
            commIndex[Where] = i; //The vertex id
        }
    }
    
    //STEP 4: Create the old2New map:
#pragma omp parallel for
    for (long i=0; i<NV; i++) {
        pOrder[commIndex[i]] = i;
        //printf("(%ld) --> (%ld)\n", commIndex[i]+1, i+1);
    }
    
    //Cleanup:
    free(C);
    free(commIndex);
    free(commPtr);
    free(commAdded);
} //End of clusteringBasedOrderingWithRcm()
