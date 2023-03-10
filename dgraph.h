#ifndef TRUSS_MAINT_GRAPH_H_
#define TRUSS_MAINT_GRAPH_H_

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <utility>
#include <vector>
#include <unordered_set>

#include "defs.h"

namespace truss_maint {
using std::int32_t;
using std::uint32_t;
// edge type
typedef std::pair<uint32_t, uint32_t> EdgT;

typedef struct final {
    uint32_t vid;
    uint32_t eid;
  } ArrayEntry;

// graph class
class Graph final {
 public:
  
  // construct a graph with only n isolated vertices
  Graph(const uint32_t n, const uint32_t l, const std::vector<int32_t>& k)
      : l_(l), n_(n), m_(0), k_(k) {
    ASSERT_MSG(0 < n_ && n_ < (static_cast<uint32_t>(1) << 29),
               "invalid argument");
    ASSERT_MSG(0 < l_ && l_ < (static_cast<uint32_t>(1) << 29),
               "invalid argument");
    // available edges
    free_edges_.resize(l_);
    std::iota(free_edges_.rbegin(), free_edges_.rend(), 0);
    free_.resize(l_, true);
    // adjacency array
    adj_.resize(n_);
    adj_in.resize(n_);
    adj_out.resize(n_);
    // edge information
    edge_info_.resize(l_, {UINT32_MAX, UINT32_MAX});
  }
  ~Graph() {}




  // get the endpoints of the edge with ID eid
  inline EdgT Get(const uint32_t eid) const {
    ASSERT_MSG(UINT32_MAX != edge_info_.at(eid).first, "invalid edge ID");
    return edge_info_[eid];
  }

  // get the ID of the edge with endpoints v1 and v2
  uint32_t Get(uint32_t v1, uint32_t v2) const {
    for (const auto ae : adj_out[v1]) {
      if (ae.vid == v2) return ae.eid;
    }
    ASSERT(false);
    return UINT32_MAX;
  }

  // get the triangles containing the edge with ID eid
  std::vector<std::pair<uint32_t, uint32_t>>
  GetTriangles(const uint32_t eid) const {
    ASSERT_MSG(UINT32_MAX != edge_info_.at(eid).first, "invalid edge ID");
    const uint32_t v1 = edge_info_[eid].first;
    const uint32_t v2 = edge_info_[eid].second;
    ASSERT(v1 < n_ && v2 < n_);
    std::vector<std::pair<uint32_t, uint32_t>> triangles;
    // find common neighbors
    size_t p1 = 0, p2 = 0;
    while (p1 < adj_[v1].size() && p2 < adj_[v2].size()) {
      if (adj_[v1][p1].vid == adj_[v2][p2].vid) {
        triangles.push_back({adj_[v1][p1].eid, adj_[v2][p2].eid});
        ++p1; ++p2;
      } else if (adj_[v1][p1].vid < adj_[v2][p2].vid) {
        ++p1;
      } else {
        ++p2;
      }
    }
    return triangles;
  }
  // get the triangles which contain the edge with ID eid
  std::vector<std::pair<uint32_t, uint32_t>>
  GetTriangles(const uint32_t eid, const int32_t k) const {
    ASSERT_MSG(UINT32_MAX != edge_info_[eid].first, "invalid edge ID");
    const uint32_t v1 = edge_info_[eid].first;
    const uint32_t v2 = edge_info_[eid].second;
    ASSERT(v1 < n_ && v2 < n_);
    std::vector<std::pair<uint32_t, uint32_t>> triangles;
    // find common neighbors
    size_t p1 = 0, p2 = 0;
    while (p1 < adj_[v1].size() && p2 < adj_[v2].size()) {
      if (adj_[v1][p1].vid == adj_[v2][p2].vid) {
        // the truss numbers of the other two edges should be at least k
        if (k_[adj_[v1][p1].eid] >= k && k_[adj_[v2][p2].eid] >= k) {
          triangles.push_back({adj_[v1][p1].eid, adj_[v2][p2].eid});
        }
        ++p1; ++p2;
      } else if (adj_[v1][p1].vid < adj_[v2][p2].vid) {
        ++p1;
      } else {
        ++p2;
      }
    }
    return triangles;
  }
  // insert an edge (v1, v2) and return its edge ID;


  // Counting cycles: two overloads, with or without k
  std::vector<std::pair<uint32_t, uint32_t>> 
  GetCycles(const uint32_t eid) const {
    ASSERT_MSG(UINT32_MAX != edge_info_.at(eid).first, "invalid edge ID");
    const uint32_t v1 = edge_info_[eid].first;
    const uint32_t v2 = edge_info_[eid].second;
    ASSERT(v1 < n_ && v2 < n_);

    if (adj_in[v1].empty() || adj_out[v2].empty()){
        return std::vector<std::pair<uint32_t, uint32_t>>();
    }
    std::vector<uint32_t> a1, a2;
    for(auto e: adj_in[v1]){
      a1.push_back(e.vid);
    }
    for(auto e: adj_out[v2]){
      a2.push_back(e.vid);
    }
    std::unordered_set<uint32_t> set{a1.cbegin(), a1.cend()};
    std::vector<std::pair<uint32_t, uint32_t>> triangles;

    for (auto n: a2){
        if (set.erase(n) > 0){ // if n exists in set, then 1 is returned and n is erased; otherwise, 0.
            std::pair<uint32_t, uint32_t> arr;
            uint32_t ar1, ar2;
            for(auto e: adj_in[v1]){
              if(n == e.vid) ar1 = e.eid;
            }
            for(auto e: adj_out[v2]){
              if(n == e.vid) ar2 = e.eid;
            }
            arr = std::make_pair(ar1, ar2);
            triangles.push_back(arr);
        } 
    }
    return triangles;
  }

  // Counting flows: two overloads, with or without k

  std::vector<std::pair<uint32_t, uint32_t>>
  GetFlows(const uint32_t eid, const int32_t k) const {
    ASSERT_MSG(UINT32_MAX != edge_info_[eid].first, "invalid edge ID");
    const uint32_t v1 = edge_info_[eid].first;
    const uint32_t v2 = edge_info_[eid].second;
    ASSERT(v1 < n_ && v2 < n_);
    std::vector<std::pair<uint32_t, uint32_t>> triangles;
    
    // guarantee that counting # of vertices instead of triangles
    std::vector<bool> v_in_flow;
    for(int i =0; i< n_;i++){
      v_in_flow.push_back(false);
    }
    // find common neighbors for flow triangles
    size_t p1 = 0, p2 = 0;

    // case 1: in-neighbors of v1 and in-neighbors of v2
    while(p1 < adj_in[v1].size() && p2 < adj_in[v2].size()){
      if(adj_in[v1][p1].vid == adj_in[v2][p2].vid){
        if(!v_in_flow[adj_in[v2][p2].vid] 
          && k_[adj_in[v1][p1].eid] >= k 
          && k_[adj_in[v2][p2].eid] >= k) {
          triangles.push_back({adj_in[p1][v1].eid, adj_in[v2][p2].eid});
          v_in_flow[adj_in[v2][p2].vid] = true;
        }
        ++p1; ++p2;
      } else if (adj_in[v1][p1].vid < adj_in[v2][p2].vid){
        ++p1;
      } else {
        ++p2;
      }
    }

    // case 2: out and out
    p1 = 0, p2 =0;
    while(p1 < adj_out[v1].size() && p2 < adj_out[v2].size()){
      if(adj_out[v1][p1].vid == adj_out[v2][p2].vid){
        if(!v_in_flow[adj_out[v2][p2].vid]
          && k_[adj_out[v1][p1].eid] >= k
          && k_[adj_out[v2][p2].eid] >= k){
          triangles.push_back({adj_out[p1][v1].eid, adj_out[v2][p2].eid});
          v_in_flow[adj_out[v2][p2].vid] = true;
        }
        ++p1; ++p2;
      } else if (adj_out[v1][p1].vid < adj_out[v2][p2].vid){
        ++p1;
      } else {
        ++p2;
      }
    }

    // case 3: out and in
    p1 = 0, p2 =0;
    while(p1 < adj_out[v1].size() && p2 < adj_in[v2].size()){
      if(adj_out[v1][p1].vid == adj_in[v2][p2].vid){
          if(!v_in_flow[adj_in[v2][p2].vid]
            && k_[adj_out[v1][p1].eid] >= k 
            && k_[adj_in[v2][p2].eid] >= k) {
            triangles.push_back({adj_out[p1][v1].eid, adj_in[v2][p2].eid});
            v_in_flow[adj_in[v2][p2].vid] = true;
          }
          ++p1; ++p2;
      } else if (adj_out[v1][p1].vid < adj_in[v2][p2].vid){
        ++p1;
      } else {
        ++p2;
      }
    }

    return triangles;

  }


  // insert an edge (v1, v2) and return its edge ID;
  // the adjacency arrays affected should be sorted later by calling Rectify

  uint32_t DiLazyInsert(const uint32_t v1, const uint32_t v2) {
    ASSERT_MSG(v1 < n_ && v2 < n_, "invalid insertion");
    ASSERT_MSG(m_ + 1 <= l_, "# of edges exceeded");
    // the ID of (v1, v2)
    const uint32_t eid = free_edges_.back();
    free_edges_.pop_back();
    free_[eid] = false;
    ++m_;
    edge_info_[eid] = {v1, v2};
    // insert the edge to the adjacency arrays
    adj_out[v1].push_back({v2, eid});
    adj_in[v2].push_back({v1, eid});
    return eid;
  }


  // sort the adjacency arrays

  void DiRectify() {
    for (uint32_t v = 0; v < n_; ++v) {
      std::sort(adj_in[v].begin(), adj_in[v].end(),
                [](const ArrayEntry& ae1, const ArrayEntry& ae2) {
                  return ae1.vid < ae2.vid;
                });
      std::sort(adj_out[v].begin(), adj_out[v].end(),
                [](const ArrayEntry& ae1, const ArrayEntry& ae2) {
                  return ae1.vid < ae2.vid;
                });
      // no duplicates
      for (size_t i = 1; i < adj_in[v].size(); ++i) {
        ASSERT_MSG(adj_in[v][i].vid > adj_in[v][i - 1].vid,
                   "duplicate edges found for adj_in");
      }
      for (size_t i = 1; i < adj_out[v].size(); ++i) {
        ASSERT_MSG(adj_out[v][i].vid > adj_out[v][i - 1].vid,
                   "duplicate edges found for adj_out");
      }
    }
  }


  // insert an edge (v1, v2) and return its edge ID
  
  uint32_t DiInsert(const uint32_t v1, const uint32_t v2) {
    // ASSERT_MSG(v1 < n_ && v2 < n_, "invalid insertion");
    ASSERT_MSG(m_ + 1 <= l_, "# of edges exceeded");
    // the ID of (v1, v2)
    const uint32_t eid = free_edges_.back();
    // insert the edge to the adjacency arrays
    size_t p1 = 0;
    while (p1 < adj_out[v1].size() && adj_out[v1][p1].vid < v2) ++p1;
    ASSERT_MSG(p1 != adj_out[v1].size() ? adj_out[v1][p1].vid > v2 : true,
               "duplicate insertion for adj_out");
    adj_out[v1].insert(adj_out[v1].begin() + p1, {v2, eid});
    size_t p2 = 0;
    while (p2 < adj_in[v2].size() && adj_in[v2][p2].vid < v1) ++p2;
    ASSERT_MSG(p2 != adj_in[v2].size() ? adj_in[v2][p2].vid > v1 : true,
               "duplicate insertion for adj_in");
    adj_in[v2].insert(adj_in[v2].begin() + p2, {v1, eid});
    // update other information
    edge_info_[eid] = {v1, v2};
    free_edges_.pop_back();
    free_[eid] = false;
    ++m_;
    return eid;
  }
  
  
  // remove the edge with ID eid

  
  void DiRemove(const uint32_t eid) {
    const uint32_t v1 = edge_info_.at(eid).first;
    const uint32_t v2 = edge_info_.at(eid).second;
    ASSERT_MSG(UINT32_MAX != v1 && UINT32_MAX != v2, "invalid deletion");
    // update information
    free_[eid] = true;
    free_edges_.push_back(eid);
    edge_info_[eid] = {UINT32_MAX, UINT32_MAX};
    // remove the edge from the adjacency arrays
    size_t p1 = 0;
    while (adj_out[v1][p1].vid != v2) ++p1;
    adj_out[v1].erase(adj_out[v1].begin() + p1);
    size_t p2 = 0;
    while (adj_in[v2][p2].vid != v1) ++p2;
    adj_in[v2].erase(adj_in[v2].begin() + p2);
    // decrease the # of edges
    --m_;
  }
  
  
  // whether the ID eid is valid
  bool Contain(const uint32_t eid) const {
    return UINT32_MAX != edge_info_.at(eid).first;
  }
  // accessors
  uint32_t n() const { return n_; }
  uint32_t m() const { return m_; }
  uint32_t l() const { return l_; }

 private:
  
  // the maximum # of edges this strcture can hold;
  // the current implementation requires l_ < 2^29
  const uint32_t l_;
  // the # of vertices
  const uint32_t n_;
  // the # of edges at the moment; m_ <= l_
  uint32_t m_;
  // the truss numbers of the edges
  const std::vector<int32_t>& k_;
  // free_[i] = true if edge ID i can be allocated; free_.size() == l_
  std::vector<bool> free_;
  // the set of available edge IDs, i.e., the set of IDs i with free_[i] = true
  std::vector<uint32_t> free_edges_;
  // adjacency arrays
  std::vector<std::vector<ArrayEntry>> adj_;
  std::vector<std::vector<ArrayEntry>> adj_in;
  std::vector<std::vector<ArrayEntry>> adj_out;
  // edge_info_[i] records the endpoints of the edge with ID i
  std::vector<EdgT> edge_info_;



};


}  // namespace truss_maint

#endif
