# Tile-X: A Vertex Reordering Approach for Scalable Long Read Assembly

# Oieswarya Bhowmik, Ananth Kalyanaraman

#      (oieswarya.bhowmik@wsu.edu, tazin.rahman@wsu.edu, ananth@wsu.edu)

# Washington State University

#

# **************************************************************************************************

# Copyright (c) 2025. Washington State University ("WSU"). All Rights Reserved.
# Permission to use, copy, modify, and distribute this software and its documentation
# for educational, research, and not-for-profit purposes, without fee, is hereby
# granted, provided that the above copyright notice, this paragraph and the following
# two paragraphs appear in all copies, modifications, and distributions. For
# commercial licensing opportunities, please contact The Office of Commercialization,
# WSU, 280/286 Lighty, PB Box 641060, Pullman, WA 99164, (509) 335-5526,
# commercialization@wsu.edu<mailto:commercialization@wsu.edu>, https://commercialization.wsu.edu/

# IN NO EVENT SHALL WSU BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
# OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
# THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF WSU HAS BEEN ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# WSU SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND
# ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". WSU HAS NO
# OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
# **************************************************************************************************







#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <thread>
#include <mutex>
#include <stack>
#include <atomic>
#include <cmath>
#include <climits>
#include <deque>
#include <utility> // For std::pair

std::mutex mtx;
std::atomic<int> nodes_processed(0);

struct Edge {
    int neighbor;
    int weight;
};

using Graph = std::unordered_map<int, std::vector<Edge>>;

// Custom hash function for a pair of integers
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

std::pair<int, int> create_edge(int a, int b) {
    return std::minmax(a, b);
}

Graph load_graph(const std::string &filename, std::unordered_set<int> &unique_nodes) {
    Graph graph;
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    int node_a, node_b, weight;
    std::unordered_set<std::pair<int, int>, pair_hash> edges;
    int total_edges = 0;

    while (infile >> node_a >> node_b >> weight) {
        auto edge = create_edge(node_a, node_b);
        if (edges.find(edge) == edges.end()) {
            edges.insert(edge);
            graph[node_a].push_back({node_b, weight});
            graph[node_b].push_back({node_a, weight});
            unique_nodes.insert(node_a);
            unique_nodes.insert(node_b);
            total_edges++;
        }
    }

    std::cout << "Total unique edges loaded: " << total_edges << std::endl;
    std::cout << "Total nodes after loading: " << unique_nodes.size() << std::endl;

    return graph;
}

void save_paths(const std::vector<std::vector<int>> &all_paths, const std::string &filename) {
    std::ofstream outfile(filename);
    if (!outfile) {
        std::cerr << "Error: Unable to open output file " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
    for (const auto &path : all_paths) {
        for (size_t i = 0; i < path.size(); ++i) {
            if (i > 0) outfile << ",";
            outfile << path[i];
        }
        outfile << "\n";
    }
    std::cout << "Paths saved to " << filename << std::endl;
}

void save_unvisited_nodes(const std::vector<int> &unvisited_nodes, const std::string &filename) {
    std::ofstream outfile(filename);
    if (!outfile) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < unvisited_nodes.size(); ++i) {
        if (i > 0) outfile << ",";
        outfile << unvisited_nodes[i];
    }
    outfile << "\n";
    std::cout << "Unvisited nodes saved to " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <graph_file> <output_file> <unvisited_nodes_file>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string graph_file = argv[1];
    std::string output_file = argv[2];
    std::string unvisited_nodes_file = argv[3];

    std::unordered_set<int> unique_graph_nodes;
    std::unordered_set<int> unique_path_nodes;

    Graph graph = load_graph(graph_file, unique_graph_nodes);
    std::cout << "Graph loaded successfully." << std::endl;
    std::cout << "Total unique nodes in the graph: " << unique_graph_nodes.size() << std::endl;

    std::vector<std::vector<int>> all_paths;

    // Placeholder for the actual DFS function
    // multithreaded_dfs(graph, unique_graph_nodes, all_paths, unique_path_nodes);

    save_paths(all_paths, output_file);

    // Placeholder for finding unvisited nodes
    std::vector<int> unvisited_nodes;
    save_unvisited_nodes(unvisited_nodes, unvisited_nodes_file);
    
    std::cout << "Traversal complete." << std::endl;
    return 0;
}
