#include "ddecom.h"
#include <chrono>
#include <cstdio>

int main(int /* argc */, char** argv) {
  // read the graph and truss-decompose it
  const auto beg = std::chrono::steady_clock::now();
  truss_maint::decomp::Decomp index(argv[1]);
  const auto end = std::chrono::steady_clock::now();
  const auto dif = end - beg;
  printf("Re-decomposing costs \x1b[1;31m%f\x1b[0m ms.\n",
         std::chrono::duration<double, std::milli>(dif).count());
  index.DWriteToFile(argv[2]);
}
