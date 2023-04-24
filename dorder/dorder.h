#ifndef TRUSS_MAINT_ORDER_H_
#define TRUSS_MAINT_ORDER_H_

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "dgraph.h"

namespace truss_maint {
using std::int32_t;
using std::uint32_t;
class Order final {
 public:
  // ctors and dtors
  // param n: the # of vertices
  // param l: for the ease of implementation, a graph can have at most l edges;
  //          that is, no more edges can be inserted if it has l edges already;
  //          the space complexity is \Theta(l + n)
  //          TODO: remove this constraint
  // param fn: the file name
  Order(const uint32_t n, const uint32_t l, const std::string& fn, uint32_t curf);
  Order(const Order&) = delete;
  Order& operator=(const Order&) = delete;
  ~Order() {}
  // near bounded insertion
  // unit edge insertion

  void DiInsert(const std::vector<EdgT>& nedges);
  // batch edge insertion if |\Delta G| >= 0.01|G|

  void DiBatchInsert(const std::vector<EdgT>& nedges);
  // bounded removal

  void DiRemove(const uint32_t v1, const uint32_t v2);

  void DiBatchRemove(const std::vector<EdgT>& redges);
  void LoadIndex(const std::string& fn, uint32_t m, uint32_t n, uint32_t curf);
  // debug
  void Debug() const;
  void Check(const std::string& fn) const;
  // accessors
  uint32_t l() const { return l_; }
  uint32_t n() const { return n_; }
  Graph g() const { return g_; }
  std::vector<int32_t> k() const { return k_; }

 private:
  struct ListNode final {
    uint32_t rem;
    uint32_t ext;
    uint32_t prev;
    uint32_t next;
  };
  struct OMNode final {
    uint64_t tag;
    uint32_t prev;
    uint32_t next;
  };
  //
  
  // list maintenance
  void ListInsert(const uint32_t e1, const uint32_t e2);
  void ListRemove(const uint32_t e, uint32_t& head, uint32_t& tail);
  // order maintenace
  void OMLoad();
  void OMInsert(const uint32_t e1, const uint32_t e2);
  void OMRemove(const uint32_t e);
  bool OMPred(const uint32_t e1, const uint32_t e2) const;
  // heap maintenance
  void HPInit();
  void HPUp(const uint32_t h, const uint32_t e);
  void HPDown(const uint32_t h, const uint32_t e);
  void HPInsert(const uint32_t e);
  void HPDelete(const uint32_t e);
  uint32_t HPTop() const;
  // members
  const uint32_t l_;
  const uint32_t n_;
  // graph and trussnesses
  Graph g_;
  std::vector<int32_t> k_;
  // basic structures
  std::vector<bool> chg_;
  std::vector<bool> new_;
  std::vector<uint32_t> head_;
  std::vector<uint32_t> tail_;
  std::vector<uint32_t> ts_;
  std::vector<uint32_t> s_;
  std::vector<ListNode> node_;
  // data members for order maintenance
  const uint32_t om_grp_ub_ = 30;
  uint32_t om_avail_;
  std::vector<OMNode> om_nodes_;
  std::vector<uint64_t> om_tag_;
  std::vector<uint32_t> om_grp_;
  std::vector<uint32_t> om_cnt_;
  // data members for heap maintenance
  std::vector<uint32_t> hp_tbl_;
  std::vector<uint32_t> hp_pos_;
  // auxiliary array
  std::vector<uint32_t>& rank_ = hp_pos_;
};

void mainDOrder(uint32_t n, uint32_t m, 
                std::vector<std::vector<ArrayEntry>>& adj_out,
                std::vector<uint32_t>& fs_,
                uint32_t& maxf,
                const std::string old_index_file, 
                const std::string update_file, 
                const std::string ground_truth_file, 
                const std::string final_file,
                const std::string f_index_file,
                const std::string op);

void mainFlow(const std::string final_file,
              std::vector<std::vector<ArrayEntry>>& adj_out,
              std::vector<uint32_t>& fs_,
              uint32_t& maxf);

}  // namespace truss_maint

#endif
