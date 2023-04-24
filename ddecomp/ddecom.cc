#include "ddecom.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <utility>
#include <unordered_set>

#define ASSERT(truth) \
    if (!(truth)) { \
      std::cerr << "\x1b[1;31mASSERT\x1b[0m: " \
                << "LINE " << __LINE__ \
                << ", " << __FILE__ \
                << std::endl; \
      std::exit(EXIT_FAILURE); \
    } else

#define ASSERT_MSG(truth, msg) \
    if (!(truth)) { \
      std::cerr << "\x1b[1;31mASSERT\x1b[0m: " \
                << "LINE " << __LINE__ \
                << ", " << __FILE__ << '\n' \
                << "\x1b[1;32mINFO\x1b[0m: " << msg \
                << std::endl; \
      std::exit(EXIT_FAILURE); \
    } else

namespace truss_maint {
namespace decomp {
// for convenience
using std::uint32_t;

// truss decomposition and the corresponding order
Decomp::Decomp(const std::string& file_name) {

  std::ifstream infile(file_name, std::ios::in);
  ASSERT_MSG(infile.is_open(), "cannot open the file");
  // read the size of the graph
  infile >> n_ >> m_;
  ASSERT_MSG(!infile.eof(), "invalid graph file");
  // read the edges
  while (true) {
    uint32_t v1, v2;
    infile >> v1 >> v2;
    if (infile.eof()) break;
    edges_.push_back({v1, v2});
  }
  infile.close();

  ASSERT_MSG(edges_.size() == m_, "invalid graph file (# of edges are not consistent)");
  for (const auto edge : edges_) {
    ASSERT_MSG(edge.first != edge.second, "self-loop exist in the graph");
    ASSERT_MSG(edge.first < n_ && edge.second < n_, "invalid vertex ID");
  }

  std::sort(edges_.begin(), edges_.end());
 
  ASSERT_MSG(std::unique(edges_.begin(), edges_.end()) == edges_.end(),
             "duplicate edges exist in the graph");


  // initialize adjacency arrays
  adj_out.resize(n_);
  adj_in.resize(n_);
  for (uint32_t eid = 0; eid < m_; ++eid) {
    const uint32_t v1 = edges_[eid].first;
    const uint32_t v2 = edges_[eid].second;
    adj_out[v1].push_back({v2, eid});
    adj_in[v2].push_back({v1, eid});
  }

  for (uint32_t vid = 0; vid < n_; ++vid) {
    adj_out[vid].shrink_to_fit();
    adj_in[vid].shrink_to_fit();
    std::sort(adj_out[vid].begin(), adj_out[vid].end(),
              [](const ArrayEntry& ae1, const ArrayEntry& ae2) {
                return ae1.vid < ae2.vid;
              });
    std::sort(adj_in[vid].begin(), adj_in[vid].end(),
              [](const ArrayEntry& ae1, const ArrayEntry& ae2) {
                return ae1.vid < ae2.vid;
              });
  }

  

  // D-truss decomposition

  std::vector<uint32_t> verts(n_);
  std::iota(verts.begin(), verts.end(), 0);
  
  uint32_t maxf= 0;
  fs_.resize(m_, 0); ford_.resize(m_); frem_.resize(m_, 0); fts_.resize(m_, 0);
  flowDecomp(adj_in, adj_out, verts, fs_, frem_, fts_, ford_, maxf);
  
  std::vector<bool> qualify(m_, true);
  D_.resize(maxf+1); Dord_.resize(maxf+1); Drem_.resize(maxf+1); Dts_.resize(maxf+1);
  

  for(uint32_t i = 0; i <= maxf; i++){

    for(uint32_t eid =0; eid< m_; eid++){
      if(fs_[eid] < i) qualify[eid] = false;      
    }

    cs_.resize(m_, 0); cord_.resize(m_); crem_.resize(m_, 0); cts_.resize(m_, 0);
    
    // 1. count cycle supports
    for (const uint32_t u : verts) {    
      for (const auto ae : adj_out[u]) {
        const uint32_t v = ae.vid;
        const uint32_t e = ae.eid;
        if(!qualify[e]) continue;
        std::vector<uint32_t> W_;
        W_ = intersectionQuali(adj_in[u], adj_out[v], qualify); 
        cs_[e] += W_.size();
      }
    }
    uint32_t maxc = *max_element(cs_.cbegin(), cs_.cend());
    if(!maxc){
      D_[i] = cs_; Drem_[i] = crem_; Dts_[i] = cts_; Dord_[i] = cord_;
      break;
    }
    cycleDecomp(adj_in, adj_out, verts, cs_, crem_, cts_, cord_, qualify);
    D_[i].resize(m_, 0); Dord_[i].resize(m_); Drem_[i].resize(m_, 0); Dts_[i].resize(m_, 0);
    D_[i] = cs_; Drem_[i] = crem_; Dts_[i] = cts_; Dord_[i] = cord_;
    cs_.clear(); cord_.clear(); crem_.clear(); cts_.clear();
    
  }
  decltype(adj_in)().swap(adj_in);
  decltype(adj_out)().swap(adj_out);
}

void Decomp::cycleDecomp(std::vector<std::vector<Decomp::ArrayEntry>> adj_in,
                        std::vector<std::vector<Decomp::ArrayEntry>> adj_out,
                        std::vector<uint32_t> verts,
                        std::vector<std::uint32_t>& cs_,
                        std::vector<std::uint32_t>& crem_,
                        std::vector<std::uint32_t>& cts_,
                        std::vector<std::uint32_t>& cord_,
                        std::vector<bool>& qualify){
  

  // 2. decomposition
  const uint32_t maxc = *max_element(cs_.cbegin(), cs_.cend()); 
  std::vector<uint32_t> cbin(maxc + 1, 0);
  // 2.1. cycle decompose: build cbin and cord_, sort the edges according to their supports
  for (uint32_t eid = 0; eid < m_; ++eid) {
    if(!qualify[eid]) continue;
    ++cbin[cs_[eid]];
  }
  for (uint32_t i = 0, start = 0; i <= maxc; ++i) {
    start += cbin[i];
    cbin[i] = start - cbin[i];
  }
  std::vector<uint32_t> cpos(m_);
  for (uint32_t eid = 0; eid < m_; ++eid) {
    if(!qualify[eid]) continue;
    cpos[eid] = cbin[cs_[eid]];
    cord_[cpos[eid]] = eid;
    ++cbin[cs_[eid]];
  }
  std::rotate(cbin.rbegin(), cbin.rbegin() + 1, cbin.rend());
  cbin[0] = 0;
  // 2.2. cycle decompose: build crem_ and cts_ via peeling
  std::vector<bool> cremoved(m_, false);
  uint32_t c = 0;
 
  for (uint32_t i = 0; i < m_; ++i) { 
    const uint32_t eid = cord_[i];
    if(!qualify[eid]) continue;
    c = std::max(c, cs_[cord_[i]]);
    ++cbin[cs_[eid]];
    cremoved[eid] = true;
    // find triangles containing the edge with ID eid
    std::vector<std::pair<uint32_t, uint32_t>> ctris; {
      const uint32_t v1 = edges_[eid].first;
      const uint32_t v2 = edges_[eid].second;
      ctris = intersecedgeQuali(adj_in[v1], adj_out[v2], qualify); 
    }

    for (const auto tri : ctris) {
      const uint32_t e1 = tri.first;
      const uint32_t e2 = tri.second;
      if(!qualify[e1] || !qualify[e2]) continue;
      if (cs_[e1] >= c && cs_[e2] >= c) ++cts_[eid];
      if (cremoved[e1] || cremoved[e2]) continue;
      ++crem_[eid];
      for (const uint32_t e : {e1, e2}) {
        if (cs_[e] > c) {
          const uint32_t pe3 = cbin[cs_[e]];
          const uint32_t pe = cpos[e];
          if (pe3 != pe) {
            const uint32_t e3 = cord_[pe3];
            cord_[pe] = e3;
            cpos[e3] = pe;
            cord_[pe3] = e;
            cpos[e] = pe3;
          }
          ++cbin[cs_[e]]; 
          --cs_[e];
        }
      }
    }
  }
}


// similar to cycleDecomp
void Decomp::flowDecomp(std::vector<std::vector<Decomp::ArrayEntry>> adj_in,
                        std::vector<std::vector<Decomp::ArrayEntry>> adj_out,
                        std::vector<uint32_t> verts,
                        std::vector<std::uint32_t>& fs_,
                        std::vector<std::uint32_t>& frem_,
                        std::vector<std::uint32_t>& fts_,
                        std::vector<std::uint32_t>& ford_,
                        std::uint32_t& maxf){

}



void Decomp::DWriteToFile(const std::string& file_name) const {
  std::ofstream outfile(file_name, std::ios::binary);
  outfile.write(reinterpret_cast<const char*>(&n_), sizeof n_)
         .write(reinterpret_cast<const char*>(&m_), sizeof m_);

  for (uint32_t i =0; i< Dord_.size(); i++) {

    for (const uint32_t e : Dord_[i]){
      const uint32_t buf[] = {edges_[e].first, edges_[e].second, D_[i][e], Drem_[i][e], Dts_[i][e]}; 
      outfile.write(reinterpret_cast<const char*>(buf), sizeof buf);
    }
    
  }
  outfile.close();
}

std::vector<uint32_t> Decomp::intersection(std::vector<Decomp::ArrayEntry>& nums1, std::vector<Decomp::ArrayEntry>& nums2) {
    if (nums1.empty() || nums2.empty()){
        return std::vector<uint32_t>();
    }
    std::vector<uint32_t> a1, a2;
    for(auto e: nums1){
      a1.push_back(e.vid);
    }
    for(auto e: nums2){
      a2.push_back(e.vid);
    }
    std::unordered_set<uint32_t> set{a1.cbegin(), a1.cend()};
    std::vector<uint32_t> intersections;
    for (auto n: a2){
        if (set.erase(n) > 0){ // if n exists in set, then 1 is returned and n is erased; otherwise, 0.
            intersections.push_back(n);
        } 
    }
    return intersections;
}

std::vector<uint32_t> Decomp::intersectionQuali(std::vector<Decomp::ArrayEntry>& nums1, std::vector<Decomp::ArrayEntry>& nums2, std::vector<bool>& qualify) {
    if (nums1.empty() || nums2.empty()){
        return std::vector<uint32_t>();
    }
    std::vector<uint32_t> a1, a2;
    for(auto e: nums1){
      if(qualify[e.eid]) a1.push_back(e.vid);
    }
    for(auto e: nums2){
      if(qualify[e.eid]) a2.push_back(e.vid);
    }
    std::unordered_set<uint32_t> set{a1.cbegin(), a1.cend()};
    std::vector<uint32_t> intersections;
    for (auto n: a2){
        if (set.erase(n) > 0){ // if n exists in set, then 1 is returned and n is erased; otherwise, 0.
            intersections.push_back(n);
        } 
    }
    return intersections;
}

std::vector<std::pair<uint32_t, uint32_t>> Decomp::intersecedge(std::vector<Decomp::ArrayEntry>& nums1, std::vector<Decomp::ArrayEntry>& nums2) {
    
    if (nums1.empty() || nums2.empty()){
        return std::vector<std::pair<uint32_t, uint32_t>>();
    }

    std::vector<uint32_t> a1, a2;
    for(auto e: nums1){
      a1.push_back(e.vid);
    }
    for(auto e: nums2){
      a2.push_back(e.vid);
    }

    std::unordered_set<uint32_t> set{a1.cbegin(), a1.cend()};
    std::vector<std::pair<uint32_t, uint32_t>> intersecedges;
    

    for (auto n: a2){
        if (set.erase(n) > 0){ // if n exists in set, then 1 is returned and n is erased; otherwise, 0.
            std::pair<uint32_t, uint32_t> arr;
            uint32_t ar1, ar2;
            for(auto e: nums1){
              if(n == e.vid) ar1 = e.eid;
            }
            for(auto e: nums2){
              if(n == e.vid) ar2 = e.eid;
            }
            arr = std::make_pair(ar1, ar2);
            intersecedges.push_back(arr);
        } 
    }
    return intersecedges;
}

std::vector<std::pair<uint32_t, uint32_t>> Decomp::intersecedgeQuali(std::vector<Decomp::ArrayEntry>& nums1, std::vector<Decomp::ArrayEntry>& nums2, std::vector<bool>& qualify) {
    
    if (nums1.empty() || nums2.empty()){
        return std::vector<std::pair<uint32_t, uint32_t>>();
    }

    std::vector<uint32_t> a1, a2;
    for(auto e: nums1){
      if(qualify[e.eid]) a1.push_back(e.vid);
    }
    for(auto e: nums2){
      if(qualify[e.eid]) a2.push_back(e.vid);
    }

    std::unordered_set<uint32_t> set{a1.cbegin(), a1.cend()};
    std::vector<std::pair<uint32_t, uint32_t>> intersecedges;
    

    for (auto n: a2){
        if (set.erase(n) > 0){ // if n exists in set, then 1 is returned and n is erased; otherwise, 0.
            std::pair<uint32_t, uint32_t> arr;
            uint32_t ar1, ar2;
            for(auto e: nums1){
              if(n == e.vid) ar1 = e.eid;
            }
            for(auto e: nums2){
              if(n == e.vid) ar2 = e.eid;
            }
            arr = std::make_pair(ar1, ar2);
            intersecedges.push_back(arr);
        } 
    }
    return intersecedges;
}


}  // namespace decomp
}  // namespace truss_maint
