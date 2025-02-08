# Tile-X: A Vertex Reordering Approach for Scalable Long Read Assembly

Tile-X is a graph-based approach for optimizing long-read genome assembly by reordering sequencing reads prior to assembly. By leveraging vertex reordering techniques, Tile-X enhances parallelism, reduces memory usage, and improves the contiguity of assembled genomes while maintaining high accuracy.

## Features

1. Graph-Theoretic Read Reordering: Computes an overlap graph and applies vertex reordering techniques to improve assembly efficiency.
2. Multiple Reordering Strategies: Implements standard reordering heuristics like Reverse Cuthill-McKee (RCM), Metis, and Grappolo, as well as a novel Farthest Neighbor (Tile-Far) heuristic for sparsified assembly.
3. Scalability: Reduces computational overhead and enables efficient assembly of large genomes.







**Tips:**
1. On some clusters, you may need to load specific modules before installing dependencies and and then also while running Maptcha.
2. Ensure that you have the appropriate permissions to execute the job script.

Maptcha utilizes the following tools:

- **JEM-Mapper**: [JEM-Mapper GitHub Repository](https://github.com/TazinRahman1105050/JEM-Mapper)
- **Hifiasm**: [Hifiasm GitHub Repository](https://github.com/chhylp123/hifiasm)
