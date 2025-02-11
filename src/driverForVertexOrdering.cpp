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
#include "vtxOrdering.h"

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
    
    // File Loading
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
        case 7: printf("Not available\n"); exit(1); break;                              //parse_DirectedEdgeList(G, inFile); break;
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
    //displayGraph(G);
    
    //Print statistics only:
    if(opts.vtxOrderingAlgorithm == 0) {
        computeMinLAScores(G);
        computeMinLAScoresNew(G);
    }
    
    //RCM: Print out the file
    if(opts.vtxOrderingAlgorithm == 1) {
        computeMinLAScores(G);
        
        // Datastructures to store old2New map
        long NV = G->numVertices;
        long *old2NewMap = (long *) malloc (NV * sizeof(long)); assert(old2NewMap != 0);
        //Initialize the Vectors:
#pragma omp parallel for
        for (long i=0; i<NV; i++) {
            old2NewMap[i] = -1; //Initialize the rank as -1
        }
        
        //Call the RCM algorithm:
        //algoReverseCuthillMcKeeStrictGraph( G, old2NewMap, nT);
        algoReverseCuthillMcKeeGraph( G, old2NewMap, nT);
        
        //Now output the graph in Gorder format:
        char outFileRcm[256];
        sprintf(outFileRcm,"%s_Rcm.edges", opts.inFile);
        writeGraphGorderEdgeListReordered(G, outFileRcm, old2NewMap);
        free(old2NewMap);
    } //End of RCM
    
    //Grappolo-based reordering:
    if(opts.vtxOrderingAlgorithm == 2) {
        computeMinLAScores(G);
        
        // Datastructures to store old2New map
        long NV = G->numVertices;
        long *old2NewMap = (long *) malloc (NV * sizeof(long)); assert(old2NewMap != 0);
        //Initialize the Vectors:
#pragma omp parallel for
        for (long i=0; i<NV; i++) {
            old2NewMap[i] = -1; //Initialize the rank as -1
        }
        clusteringBasedOrderingSimple(G, old2NewMap, nT, opts.coloring);
        
        //Now output the graph in Gorder format:
        char outFileGrappolo[256];
        sprintf(outFileGrappolo,"%s_Grappolo.edges", opts.inFile);
        writeGraphGorderEdgeListReordered(G, outFileGrappolo, old2NewMap);
        free(old2NewMap);
    } //End of Grappolo
    
    //Write file in Gorder format:
    if(opts.vtxOrderingAlgorithm == 3) {
        computeMinLAScores(G);
        
        char outGorder[256];
        sprintf(outGorder,"%s_Gorder.edges", opts.inFile); //Output file name
        writeGraphGorderEdgeList(G, outGorder);
    }
    
    //RCM-Strict: Print out the file
    if(opts.vtxOrderingAlgorithm == 4) {
        computeMinLAScores(G);
        
        // Datastructures to store old2New map
        long NV = G->numVertices;
        long *old2NewMap = (long *) malloc (NV * sizeof(long)); assert(old2NewMap != 0);
        //Initialize the Vectors:
#pragma omp parallel for
        for (long i=0; i<NV; i++) {
            old2NewMap[i] = -1; //Initialize the rank as -1
        }
        
        //Call the RCM algorithm:
        //algoReverseCuthillMcKeeStrictGraph( G, old2NewMap, nT);
        algoReverseCuthillMcKeeStrictGraph( G, old2NewMap, nT);
        
        //Now output the graph in Gorder format:
        char outFileRcm[256];
        sprintf(outFileRcm,"%s_RcmStrict.edges", opts.inFile);
        writeGraphGorderEdgeListReordered(G, outFileRcm, old2NewMap);
        free(old2NewMap);
    } //End of RCM
    
    //Grappolo-ML reordering:
    if(opts.vtxOrderingAlgorithm == 5) {
        computeMinLAScores(G);
        
        // Datastructures to store old2New map
        long NV = G->numVertices;
        long *old2NewMap = (long *) malloc (NV * sizeof(long)); assert(old2NewMap != 0);
        //Initialize the Vectors:
#pragma omp parallel for
        for (long i=0; i<NV; i++) {
            old2NewMap[i] = -1; //Initialize the rank as -1
        }
        clusteringBasedOrderingWithRcm(G, old2NewMap, nT, opts.coloring);
        
        //Now output the graph in Gorder format:
        char outFileGrappolo[256];
        sprintf(outFileGrappolo,"%s_GrappoloML.edges", opts.inFile);
        writeGraphGorderEdgeListReordered(G, outFileGrappolo, old2NewMap);
        free(old2NewMap);
    } //End of Grappolo-ML
    
    //Degree-based: Print out the file
    if(opts.vtxOrderingAlgorithm == 6) {
        computeMinLAScores(G);
        
        // Datastructures to store old2New map
        long NV = G->numVertices;
        long *old2NewMap = (long *) malloc (NV * sizeof(long)); assert(old2NewMap != 0);
        //Initialize the Vectors:
#pragma omp parallel for
        for (long i=0; i<NV; i++) {
            old2NewMap[i] = -1; //Initialize the rank as -1
        }
        
        //Call the Degree-based algorithm:
        degreeBasedOrdering(G, old2NewMap, nT);
        
        //Now output the graph in Gorder format:
        char outFileDegree[256];
        sprintf(outFileDegree,"%s_Degree.edges", opts.inFile);
        writeGraphGorderEdgeListReordered(G, outFileDegree, old2NewMap);
        free(old2NewMap);
    } //End of Degree
    
    //Cleanup:
    if(G != 0) {
        free(G->edgeListPtrs);
        free(G->edgeList);
        free(G);
    }
    return 0;
}//End of main()
