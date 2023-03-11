#include <unistd.h>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "defs.h"
#include "dgraph.h"
#include "dorder.h"


// a sample program
int main(int argc, char** argv) {
  ASSERT(7 == argc);
  const std::string op = argv[1];
  const std::string old_index_file = argv[2];
  const std::string update_file = argv[3];
  const std::string ground_truth_file = argv[4];
  const std::string final_file = argv[5];
  const std::string f_index_file = argv[6];
  printf("*****************************************************************\n");
  printf("old index file: %s\n", old_index_file.c_str());
  printf("update file: %s\n", update_file.c_str());
  printf("ground truth file: %s\n", ground_truth_file.c_str());
  printf("final file: %s\n", final_file.c_str());
  printf("f index file: %s\n", f_index_file.c_str());
  printf("*****************************************************************\n");
  // read the header
  std::ifstream infile(old_index_file, std::ios::binary);
  uint32_t n = -1, m = -1;  // # of vertices and edges
  infile.read(reinterpret_cast<char*>(&n), sizeof n)
        .read(reinterpret_cast<char*>(&m), sizeof m);
  infile.close();
  // read the graph and the index
  // We set the param "l" to "m * 2" here, where "l" is the maximum number of
  // edges that a graph can hold. We impose this constraint only for ease of
  // implementation. TODO: remove this constraint in the future.
  
  // apply the updates
  const auto beg = std::chrono::steady_clock::now();
  std::vector<std::vector<truss_maint::ArrayEntry>> adj_out;
  std::vector<uint32_t> fs_; uint32_t maxf= 0;
  truss_maint::mainFlow(final_file, adj_out, fs_, maxf);
  truss_maint::mainDOrder(n, m, adj_out, fs_, maxf, old_index_file, update_file, ground_truth_file, final_file, f_index_file, op);
  const auto end = std::chrono::steady_clock::now();
  const auto dif = end - beg;
  printf("Applying the updates costs \x1b[1;31m%f\x1b[0m ms.\n",
         std::chrono::duration<double, std::milli>(dif).count());
  // verify the result
  
  printf("*****************************************************************\n");
}
