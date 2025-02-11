// ***********************************************************************
//
//            Grappolo: A C++ library for graph clustering
//               Mahantesh Halappanavar (hala@pnnl.gov)
//               Pacific Northwest National Laboratory
//
// ***********************************************************************
//
//       Copyright (2014) Battelle Memorial Institute
//                      All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************

#include "defs.h"

using namespace std;

//Half-approximate matching algorithm
//Returns the matching information in Mate:
void algoHalfApproxMatchingSuitor( graph* G, long *Mate )
{
    int nthreads;
#pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }
    printf("algoEdgeApproxDominatingEdgesSuitorNew(): Number of threads: %d\n", nthreads);
    
    double time1=0, time2=0;
    //Get the iterators for the graph:
    long NVer     = G->numVertices;
    long NEdge    = G->numEdges;       //Returns the correct number of edges (not twice)
    long *verPtr  = G->edgeListPtrs;   //Vertex Pointer: pointers to endV
    edge *verInd  = G->edgeList;       //Vertex Index: destination id of an edge (src -> dest)
    printf("|V|= %ld  |E|= %ld\n", NVer, NEdge);
    
    //Allocate Data Structures:
    time1 = omp_get_wtime();
    omp_lock_t *nlocks = (omp_lock_t *) malloc (NVer * sizeof(omp_lock_t));
    assert( nlocks != 0);
    long *s            = Mate;
    double *ws         = (double *) malloc (NVer * sizeof(double)); assert(ws != 0);
    //Initialize the Vectors:
#pragma omp parallel for
    for (long i=0; i<NVer; i++) {
        ws[i]= 0.0f;               //Set current weight of best suitor to zero
        omp_init_lock(&nlocks[i]); //Initialize locks
    }
    time1  = omp_get_wtime() - time1;
    
    time2 = omp_get_wtime();
#pragma omp parallel for
    for(long x=0; x<NVer; x++) {  // Loop over vertices
        long i = x;
        int done = 0; //FALSE;
        while (!done) {
            double heaviest = ws[i];
            long partner     =  s[i];
            long next_vertex;
            //printf("Processing node %d \n",i);
            long adj1 = verPtr[i];
            long adj2 = verPtr[i+1];
            for(long j=adj1; j<adj2; j++) { // Loop over neighbors of vertex i
                long y = verInd[j].tail;    // y is the current neighbor of i
                double edgeWeight = verInd[j].weight; //weight of y
                if( (edgeWeight < heaviest)||(edgeWeight < ws[y]) )
                    continue;
                if( (edgeWeight == heaviest)&&(y < partner) )
                    continue;
                if( (edgeWeight == ws[y])&&(i < s[y]) )
                    continue;
                // Check if w(i,y) is the best so far, and if it is a better option for y
                heaviest = edgeWeight;      // Store the best so far
                partner = y;
            } // loop over neighbors
            done = 1; //TRUE;
            if ( heaviest > 0 ) {
                omp_set_lock(&nlocks[partner]);    // Locking partner
                if( (heaviest > ws[partner])||( (heaviest == ws[partner])&&(i>s[partner]) ) ) {
                    if (s[partner] >= 0 ) {
                        next_vertex = s[partner];
                        done = 0; //FALSE;
                    }
                    s[partner]  = i;         // i is now the current suitor of s
                    ws[partner] = heaviest;  // The weight of the edge (i,partner)
                }
                else {   // can no longer use the result for node i, must find new partner
                    done = 0; // FALSE;
                    next_vertex = i;
                }
                omp_unset_lock(&nlocks[partner]); // Unlocking partner
            }
            if( !done ) { // Continue with the next vertex
                i = next_vertex;
            }
        } // while not done
    } // loop over vertices
    time2  = omp_get_wtime() - time2;
    
    printf("***********************************************\n");
    printf("Time for Initialization    : %lf sec\n", time1);
    printf("Time for Matching          : %lf sec\n", time2);
    printf("Total Time                 : %lf sec\n", time1+time2);
    printf("***********************************************\n");
    
    //Clean Up:
    free(ws);
#pragma omp parallel for
    for (long i=0; i<NVer; i++) {
        omp_destroy_lock(&nlocks[i]); //Initialize locks
    }
    //free(nlocks);
} //End of algoHalfApproxMatchingSuitor()


