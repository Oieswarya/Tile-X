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
#include "input_output.h"
#include "basic_comm.h"
#include "basic_util.h"
#include "utilityClusteringFunctions.h"
#include "color_comm.h"
#include "sync_comm.h"
#include "utilityGraphPartitioner.h"

using namespace std;
//#define USEHDF5

int main(int argc, char** argv) {
    
    //Parse Input parameters:
    clustering_parameters opts;
    if (!opts.parse(argc, argv)) {
        return -1;
    }
    int nT = 1; //Default is one thread
#pragma omp parallel
    {
        nT = omp_get_num_threads();
    }
    if (nT <= 1) {
        printf("The number of threads should be greater than one.\n");
        return 0;
    }
    double time1, time2;
    graph* G = (graph *) malloc (sizeof(graph));
    
    /* Step 2: Parse the graph in Matrix Market format */
    int fType = opts.ftype; //File type
    char *inFile = (char*) opts.inFile;
    bool isSym = true; //Assume symmetric by default
    switch (fType) {
        case 1: parse_MatrixMarket_Sym_AsGraph(G, inFile); break;
        case 2: parse_Dimacs9FormatDirectedNewD(G, inFile); break;
        case 3: parse_PajekFormat(G, inFile); break;
        case 4: parse_PajekFormatUndirected(G, inFile); break;
        case 5: loadMetisFileFormat(G, inFile); break;
        case 6: parse_UndirectedEdgeList(G, inFile); break;
            //parse_UndirectedEdgeListDarpaHive(G, inFile); break;
        case 7: /*parse_DirectedEdgeList(G, inFile);*/ printf("This routine is under development.\n"); exit(1); break;
        case 8: parse_SNAP(G, inFile); break;
        case 9: parse_EdgeListBinaryNew(G,inFile); break;
        case 10:
#ifdef USEHDF5
            //parse_EdgeListCompressedHDF5(G,inFile);
            parse_EdgeListCompressedHDF5NoDuplicates(G,inFile);
#endif
            break;
        case 11: parse_UndirectedEdgeListFromJason(G, inFile); break;
        case 12: parse_UndirectedEdgeListWeighted(G, inFile); break; // for John F's graphs
        case 13: parse_UndirectedEdgeListDarpaHive(G, inFile); break;
        case 14: parse_EdgeListFromGorder(G, inFile); break;
        default:  cout<<"A valid file type has not been specified"<<endl; exit(1);
    }
    
    displayGraphCharacteristics(G);
    int threadsOpt = 0;
    if(opts.threadsOpt)
        threadsOpt = 1;
    threadsOpt =1;
    
    long NV = G->numVertices;
    long *C = (long *) malloc (NV * sizeof(long)); assert(C != 0);
    long *old2New   = (long *) malloc (NV * sizeof(long)); assert(old2New != 0);
    double timePartition=0, timeRebuild=0, tmpTime=0;
    
    //Set the parameter for #partitions
    //int myVec[7]={4,8,16,32,64,128,256};
    int myVec[1]={32};
    
    for (int i=0; i<2; i++) {
        //Initialize:
        for (long j=0; j<NV; j++) {
            C[j] = -1;
            old2New[j]   = -1;
        }
        //Step 1: Call the graph partitioner with a given number of partitions
        tmpTime = omp_get_wtime();
        MetisGraphPartitioner( G, C, myVec[i]);
        tmpTime = omp_get_wtime() - tmpTime;
        timePartition += tmpTime;
        
        
        //STEP 2: Now build a compact representation of partitions:
        //        Create a CSR-like datastructure for partitions
        //Compute number of communities:
        //Assume zero is a valid community id
        long nC=-1;
        bool isZero = false;
        bool isNegative = false;
        for(long i = 0; i < NV; i++) {
            if(C[i] == 0)
                isZero = true; //Check if zero is a valid partition
            if(C[i] < 0)
                isNegative = true; //Check if zero is a valid partition
            if (C[i] > nC) {
                nC = C[i];
            }
        }
        printf("Largest partition id observed: %ld\n", nC);
        if(isZero) {
            printf("Zero is a valid community id\n");
            nC++;
        }
        if(isNegative) {
            printf("Some vertices have not been assigned communities\n");
            nC++; //Place to store all the unassigned vertices
        }
        assert(nC>0);
        printf("Number of unique partitions in C= %d\n", nC);
        
        tmpTime = omp_get_wtime();
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
            old2New[commIndex[i]] = i;
            //printf("(%ld) --> (%ld)\n", commIndex[i]+1, i+1);
        }
        tmpTime = omp_get_wtime() - tmpTime;
        timeRebuild += tmpTime;
        
        //Step 3: Write the reordered graph to a file:
        char outFile[256];
        sprintf(outFile,"%s_METIS_%d.edges", opts.inFile, myVec[i]);
        printf("Processing with %d partitions; will be stored in file: %s\n", nC, outFile);
        writeGraphGorderEdgeListReordered(G, outFile, old2New);
    }
    printf("=====================================================\n");
    printf("Time for partitioning: %lf\n", timePartition);
    printf("Time for rebuilding  : %lf\n", timeRebuild);
    printf("Total time           : %lf\n", timePartition+timeRebuild);
    printf("=====================================================\n");
    
    //Cleanup:
    if(G != 0) {
        free(G->edgeListPtrs);
        free(G->edgeList);
        free(G);
    }
    free(C);
    free(old2New);
    
    return 0;
}//End of main()
