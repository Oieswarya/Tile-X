# Variables
SRC_DIR = src
BIN_DIR = bin
TILEX_SCRIPT = tileX.sh
HIFIASM_DIR = Hifiasm
HIFIASM_BIN = $(HIFIASM_DIR)/hifiasm
#JEM_MAPPER_DIR = JEM-Mapper

# Compiler and flags
CXX = g++
CXXFLAGS = -fopenmp -O3

# Source files
GRAPH_CONSTR_WH_SRC = $(SRC_DIR)/GraphConstrWH.cpp
GRAPH_LRID_SRC = $(SRC_DIR)/graphLRID.cpp

# Binaries (all executables will go inside bin/)
GRAPH_CONSTR_WH_BIN = $(BIN_DIR)/GraphConstrWH
GRAPH_LRID_BIN = $(BIN_DIR)/graphLRID

# Default target
all: setup compile hifiasm jem_mapper reordering permissions

# Setup directories
setup:
	mkdir -p $(BIN_DIR)

# Compile source files
compile: $(GRAPH_CONSTR_WH_BIN) $(GRAPH_LRID_BIN)

$(GRAPH_CONSTR_WH_BIN): $(GRAPH_CONSTR_WH_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(GRAPH_LRID_BIN): $(GRAPH_LRID_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

# Compile Hifiasm
hifiasm:
	cd $(HIFIASM_DIR) && make

# Compile JEM-Mapper
#jem_mapper:
	#cd $(JEM_MAPPER_DIR) && make ksize=15

# Compile Reordering executables (now in src/ and placed inside bin/)
reordering:
	cd $(SRC_DIR) && make clean && make CC=gcc CPP=g++

	mv $(SRC_DIR)/bin/* $(BIN_DIR)/
	@echo "Reordering executables are now in $(BIN_DIR)"

# Ensure correct permissions and format
permissions:
	chmod +x $(HIFIASM_BIN)
	dos2unix $(TILEX_SCRIPT)
	chmod +x $(TILEX_SCRIPT)

# Clean up binaries and output
clean:
	rm -rf $(BIN_DIR)
	cd $(SRC_DIR) && make clean

# Install necessary Python dependencies (if needed)
install-dependencies:
	pip install -r requirements.txt

# Help
help:
	@echo "Usage:"
	@echo "  make all                  - Setup directories, compile binaries, compile Hifiasm, compile JEM-Mapper, compile Reordering, set permissions"
	@echo "  make setup                - Setup directories"
	@echo "  make compile              - Compile source files"
	@echo "  make hifiasm              - Compile Hifiasm"
	@echo "  make jem_mapper           - Compile JEM-Mapper with ksize=15"
	@echo "  make reordering           - Compile Reordering executables"
	@echo "  make permissions          - Set permissions and format"
	@echo "  make clean                - Clean up binaries and output"
	@echo "  make install-dependencies - Install necessary Python dependencies"
	@echo "  make help                 - Show this help message"

.PHONY: all setup compile hifiasm jem_mapper reordering permissions clean install-dependencies help
