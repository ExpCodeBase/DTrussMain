# README #

This is the codebase of the Dynamic D-truss Search. 

## File Organization ##

* **ddecomp/ddecomp.**: the implementation of D-truss decomposition
* **ddecomp/dsample.cc**: perform the decomposition and initialize the index
* **dorder/dorder.cc**: the implementation of D-truss maintenance
* **dorder/dgraph.h**: the supportive functions for D-truss management
* **dorder/dtest.cc**: perform the maintenance based on the index

## How to Use the Code? ##

* First, running the code under the path `ddecomp` to construct an index for the graph before update.
* Then, using this index as input, you can construct the D-Index to handle update edges, the detailed commands are as follows.

### Command Lines ###

* Perform the decomposition and initialize the index, under the path `./ddecomp/`:

  `./dsample <DATA_PATH> <INDEX_PATH>`

* Perform the maintenance based on the D-Index, under the path `./dorder/`:
  
  unit delete:`./dm udelete <OLD_INDEX_PATH> <UPDATE_EDGE_PATH> <GROUND_TRUTH_PATH> <UPDATED_INDEX_PATH> <HELPER_FINDEX_PATH>`

  batch delete:`./dm bdelete <OLD_INDEX_PATH> <UPDATE_EDGE_PATH> <GROUND_TRUTH_PATH> <UPDATED_INDEX_PATH> <HELPER_FINDEX_PATH>`

  unit insert:`./dm uinsert <OLD_INDEX_PATH> <UPDATE_EDGE_PATH> <GROUND_TRUTH_PATH> <UPDATED_INDEX_PATH> <HELPER_FINDEX_PATH>`

  batch insert:`./dm binsert <OLD_INDEX_PATH> <UPDATE_EDGE_PATH> <GROUND_TRUTH_PATH> <UPDATED_INDEX_PATH> <HELPER_FINDEX_PATH>`

