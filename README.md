# Tile-X: A Vertex Reordering Approach for Scalable Long Read Assembly

Tile-X is a graph-based approach for optimizing long-read genome assembly by reordering sequencing reads prior to assembly. By leveraging vertex reordering techniques, Tile-X enhances parallelism, reduces memory usage, and improves the contiguity of assembled genomes while maintaining high accuracy.

## Features

1. Graph-Theoretic Read Reordering: Computes an overlap graph and applies vertex reordering techniques to improve assembly efficiency.
2. Multiple Reordering Strategies: Implements standard reordering heuristics like Reverse Cuthill-McKee (RCM), Metis, and Grappolo, as well as a novel Farthest Neighbor (Tile-Far) heuristic for sparsified assembly.
3. Scalability: Reduces computational overhead and enables efficient assembly of large genomes.

### Step-by-Step Guide

1. **Clone the Tile-X Repository:**

   ```bash
   git clone https://github.com/Oieswarya/Tile-X.git
   cd Tile-X

2. **Compile the source files and setup directories:**

   ```bash
   make all

3. **Check if Tile-X is properly installed:**

   ```bash
   ./tileX.sh -h


### Usage
Run the tileX.sh script from the root directory:

```bash
./tileX.sh -lr path/to/longreads.fa [options]


-lr,--longreads    Path to the long reads input file
Options:
-o, --output       Output directory (default: $HOME/Maptcha/Output/)
-t, --threads      Number of threads to use (default: 16)
-n, --nodes        Number of nodes to use (default: 2)
-p, --processes    Number of processes per node (default: 2)
-tile, --module    Tile-X module to use (default: Tile-Far)
                     Options: Tile-Far (default), Tile-RCM (rcm), Tile-Metis (met), Tile-Grappolo (grap)
-h, --help         Show this help message
```

Note:
This code has been tested on high-performance cluster (HPC) systems with MPI and OpenMP compatibility and has been tested for both PBS and SLURM job scheduling systems.


### For a quick test, you can use the provided test input. Navigate within the Tile-X repository and run the `tileX.sh` script. 

```bash
~/Tile-X/tileX.sh ~/Tile-X/TestInput/CoxiellaBurnetii_longreads.fa
```

The final scaffolds will be located here: `~/Tile-X/Output/Final/finalAssembly.fa`, within the Output folder of the Maptcha directory.




**Tips:**
1. On some clusters, you may need to load specific modules before installing dependencies and and then also while running Maptcha.
2. Ensure that you have the appropriate permissions to execute the job script.

Tile-X utilizes the following tools:

- **JEM-Mapper**: [JEM-Mapper GitHub Repository](https://github.com/TazinRahman1105050/JEM-Mapper)
- **Hifiasm**: [Hifiasm GitHub Repository](https://github.com/chhylp123/hifiasm)
