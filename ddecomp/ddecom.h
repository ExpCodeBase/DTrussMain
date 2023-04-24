#ifndef TRUSS_MAINT_DECOMP_DECOMP_H_
#define TRUSS_MAINT_DECOMP_DECOMP_H_

#include <cstdint>
#include <string>
#include <vector>

namespace truss_maint {
namespace decomp {
// class Decomp is to truss-decompose a graph; as a byproduct,
// it produces a truss-decomposition order



class Decomp final {
 public:

  Decomp(const std::string& file_name);
  Decomp(const Decomp&) = delete;
  Decomp& operator=(const Decomp&) = delete;
  // write the results to disk
  void cWriteToFile(const std::string& file_name) const;
  void fWriteToFile(const std::string& file_name) const;
  void DWriteToFile(const std::string& file_name) const;
  // adjacency array entry type
  typedef struct final {
    std::uint32_t vid;
    std::uint32_t eid;
  } ArrayEntry;
  std::vector<uint32_t> intersection(std::vector<Decomp::ArrayEntry>& nums1, std::vector<Decomp::ArrayEntry>& nums2);
  std::vector<uint32_t> intersectionQuali(std::vector<Decomp::ArrayEntry>& nums1, std::vector<Decomp::ArrayEntry>& nums2, std::vector<bool>& qualify);
  std::vector<std::pair<uint32_t, uint32_t>> intersecedge(std::vector<Decomp::ArrayEntry>& nums1, std::vector<Decomp::ArrayEntry>& nums2);
  std::vector<std::pair<uint32_t, uint32_t>> intersecedgeQuali(std::vector<Decomp::ArrayEntry>& nums1, std::vector<Decomp::ArrayEntry>& nums2, std::vector<bool>& qualify);
  void cycleDecomp(std::vector<std::vector<Decomp::ArrayEntry>> adj_in,
                        std::vector<std::vector<Decomp::ArrayEntry>> adj_out,
                        std::vector<uint32_t> verts,
                        std::vector<std::uint32_t>& cs_,
                        std::vector<std::uint32_t>& crem_,
                        std::vector<std::uint32_t>& cts_,
                        std::vector<std::uint32_t>& cord_,
                        std::vector<bool>& qualify);
  void flowDecomp(std::vector<std::vector<Decomp::ArrayEntry>> adj_in,
                        std::vector<std::vector<Decomp::ArrayEntry>> adj_out,
                        std::vector<uint32_t> verts,
                        std::vector<std::uint32_t>& fs_,
                        std::vector<std::uint32_t>& frem_,
                        std::vector<std::uint32_t>& fts_,
                        std::vector<std::uint32_t>& ford_,
                        std::uint32_t& maxf);                    
 private:
  
  // data members
  std::uint32_t n_;  // the # of vertices
  std::uint32_t m_;  // the # of edges
  // the adjacency array representation
  std::vector<std::vector<ArrayEntry>> adj_;
  std::vector<std::vector<ArrayEntry>> adj_in;
  std::vector<std::vector<ArrayEntry>> adj_out;
  // the cycle support and flow support
  std::vector<std::uint32_t> cs_;
  std::vector<std::uint32_t> fs_;
  // the D-trussness
  std::vector<std::vector<std::uint32_t>> D_;
  // the truss numbers
  std::vector<std::uint32_t> k_;
  // the remaining supports
  std::vector<std::uint32_t> rem_;
  std::vector<std::uint32_t> crem_;
  std::vector<std::uint32_t> frem_;
  std::vector<std::vector<std::uint32_t>> Drem_;
  // the triangle supports
  std::vector<std::uint32_t> ts_;
  std::vector<std::uint32_t> cts_;
  std::vector<std::uint32_t> fts_;
  std::vector<std::vector<std::uint32_t>> Dts_;
  // the edge peeling order
  std::vector<std::uint32_t> ord_;
  std::vector<std::uint32_t> cord_;
  std::vector<std::uint32_t> ford_;
  std::vector<std::vector<std::uint32_t>> Dord_;
  // the set of edges
  std::vector<std::pair<std::uint32_t, std::uint32_t>> edges_;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> edges_out;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> edges_in;
};


}  // namespace decomp
}  // namespace truss_maint

#endif
