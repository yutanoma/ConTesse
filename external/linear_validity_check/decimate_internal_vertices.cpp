#include "decimate_internal_vertices.h"

#include <Eigen/Core>
#include <tuple>
#include <vector>
#include <queue>
#include <unordered_map>
#include <cstdint>
#include <cassert>

namespace utils {
namespace details {
  using Row2 = Eigen::RowVector2d;
  
  struct HE { int v=-1, twin=-1, next=-1, prev=-1, f=-1; }; // head vertex id
  struct Face { int he=-1; };
  struct Vert { int he=-1, val=0; bool alive=true, boundary=false; };

  struct Mesh {
      const Eigen::MatrixXd &P; // geometry (Nx2)
      std::vector<HE>   he;
      std::vector<Face> face;
      std::vector<Vert> vert;

      explicit Mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) : P(V)
      {
          const int nV = (int)V.rows();
          const int nF = (int)F.rows();
          vert.assign(nV, Vert{});
          face.assign(nF, Face{});
          he.resize(3*nF);

          auto dkey = [](int a, int b)->uint64_t { return (uint64_t(uint32_t(a))<<32) | uint32_t(b); };

          for(int fi=0; fi<nF; ++fi){
              int a = F(fi,0), b = F(fi,1), c = F(fi,2);
              int e0 = 3*fi+0, e1 = 3*fi+1, e2 = 3*fi+2;
              he[e0].v = b; he[e1].v = c; he[e2].v = a;
              he[e0].f = he[e1].f = he[e2].f = fi;
              he[e0].next = e1; he[e1].next = e2; he[e2].next = e0;
              he[e0].prev = e2; he[e1].prev = e0; he[e2].prev = e1;
              face[fi].he = e0;
              vert[a].val++; vert[b].val++; vert[c].val++;
              if(vert[a].he<0) vert[a].he = e0; // tail(e2)=a
              if(vert[b].he<0) vert[b].he = e1; // tail(e0)=b
              if(vert[c].he<0) vert[c].he = e2; // tail(e1)=c
          }
          std::unordered_map<uint64_t,int> map_dir;
          map_dir.reserve(3*nF*2);
          auto tail_of = [&](int e)->int { return he[he[e].prev].v; };

          for(int e=0; e<(int)he.size(); ++e){
              int t = tail_of(e), h = he[e].v;
              map_dir.emplace(dkey(t,h), e);
          }
          for(int e=0; e<(int)he.size(); ++e){
              int t = tail_of(e), h = he[e].v;
              auto it = map_dir.find(dkey(h,t));
              if(it!=map_dir.end()) he[e].twin = it->second;
          }
          // boundary flags via a rotation around each vertex
          for(int u=0; u<nV; ++u){
              if(vert[u].he<0){ vert[u].boundary=true; continue; }
              bool isB=false;
              int start = vert[u].he, e = start;
              for(;;){
                  if(he[e].twin<0){ isB=true; break; }
                  int e_next_around_u = he[he[e].twin].next;
                  if(e_next_around_u==start) break;
                  e = e_next_around_u;
              }
              vert[u].boundary = isB;
          }
      }

      static inline double orient(const Row2& a, const Row2& b, const Row2& c){
          return (b.x()-a.x())*(c.y()-a.y()) - (b.y()-a.y())*(c.x()-a.x());
      }

      bool convex_quad_flip_ok(int a,int u,int b,int v, double eps=0.0) const {
          const Row2 A=P.row(a), U=P.row(u), B=P.row(b), Vv=P.row(v);
          if(orient(U,Row2(P.row(v)),A) <= eps) return false;
          if(orient(Row2(P.row(v)),U,B) <= eps) return false;
          double o1 = orient(U,Row2(P.row(v)),A);
          double o2 = orient(U,Row2(P.row(v)),B);
          double o3 = orient(A,B,U);
          double o4 = orient(A,B,Row2(P.row(v)));
          if(!(o1*o2 < 0 && o3*o4 < 0)) return false;
          if(orient(U,A,B) <= eps) return false;
          if(orient(Row2(P.row(v)),B,A) <= eps) return false;
          return true;
      }

      inline int tail_of(int e) const { return he[he[e].prev].v; }

      bool flip_edge(int e){
          int et = he[e].twin; if(et<0) return false;
          int fL = he[e].f, fR = he[et].f; if(fL<0||fR<0) return false;

          int e_uv_a = he[e].next;
          int e_a_u  = he[e].prev;
          int et_v_b = he[et].next;
          int et_b_v = he[et].prev;

          int u = tail_of(e);
          int v = he[e].v;
          int a = he[e_uv_a].v;
          int b = he[et_v_b].v;

          if(!convex_quad_flip_ok(a,u,b,v)) return false;

          he[e].v  = b;   // now (u->b)
          he[et].v = a;   // now (v->a)

          he[e].next = et_b_v;    he[et_b_v].prev = e;
          he[et_b_v].next = e_a_u; he[e_a_u].prev = et_b_v;
          he[e_a_u].next = e;     he[e].prev     = e_a_u;
          he[e].f = he[et_b_v].f = he[e_a_u].f = fL;
          face[fL].he = e;

          he[et].next = e_uv_a;    he[e_uv_a].prev = et;
          he[e_uv_a].next = et_v_b; he[et_v_b].prev = e_uv_a;
          he[et_v_b].next = et;    he[et].prev     = et_v_b;
          he[et].f = he[e_uv_a].f = he[et_v_b].f = fR;
          face[fR].he = et;

          auto fix_rep = [&](int x, int candidate){
              if(!vert[x].alive) return;
              if(vert[x].he<0) { vert[x].he=candidate; return; }
              if(tail_of(vert[x].he)!=x) vert[x].he = candidate;
          };
          fix_rep(u, e);
          fix_rep(v, et);
          fix_rep(a, e_uv_a);
          fix_rep(b, et_v_b);
          return true;
      }

      void neighbors_ccw(int u, std::vector<int>& nbrs, std::vector<int>& out_he) const {
          nbrs.clear(); out_he.clear();
          if(vert[u].he<0) return;
          int start = vert[u].he;
          if(tail_of(start)!=u){
              if(he[start].twin>=0) start = he[he[start].twin].next;
          }
          int e = start;
          for(;;){
              nbrs.push_back(he[e].v);
              out_he.push_back(e);
              if(he[e].twin<0) break; // boundary reached
              int enext = he[he[e].twin].next;
              if(enext==start) break;
              e = enext;
          }
      }

      // allow_boundary==true permits removing a vertex marked as boundary
      bool remove_valence3(int u, bool allow_boundary){
          if(!vert[u].alive || vert[u].val!=3) return false;
          if(vert[u].boundary && !allow_boundary) return false;

          std::vector<int> nbrs, outhe;
          neighbors_ccw(u, nbrs, outhe);
          if((int)nbrs.size()!=3) return false;
          int a=nbrs[0], b=nbrs[1], c=nbrs[2];
          int e01=outhe[0], e12=outhe[1], e20=outhe[2];

          int f1 = he[e01].f; // (u,a,b)
          int f2 = he[e12].f; // (u,b,c)
          int f3 = he[e20].f; // (u,c,a)
          if(f1<0||f2<0||f3<0) return false;

          auto next_of = [&](int e)-> int { return he[e].next; };
          int h_ab = next_of(e01);      // (a->b)
          int h_bc = next_of(e12);      // (b->c)
          int h_ca = next_of(e20);      // (c->a)

          he[h_ab].f = he[h_bc].f = he[h_ca].f = f1;
          he[h_ab].prev = h_ca; he[h_ab].next = h_bc;
          he[h_bc].prev = h_ab; he[h_bc].next = h_ca;
          he[h_ca].prev = h_bc; he[h_ca].next = h_ab;
          face[f1].he = h_ab;

          auto kill_face = [&](int f){ face[f].he=-1; };
          auto kill_he   = [&](int h){ he[h].f=-1; he[h].next=-1; he[h].prev=-1; };

          kill_face(f2); kill_face(f3);
          kill_he(e01); kill_he(e12); kill_he(e20);
          kill_he(he[e01].prev); kill_he(he[e12].prev); kill_he(he[e20].prev);

          vert[u].alive=false; vert[u].he=-1;
          vert[a].val--; vert[b].val--; vert[c].val--;

          auto ensure_rep = [&](int x, int keep_he){
              if(!vert[x].alive) return;
              if(vert[x].he<0 || he[vert[x].he].f<0 || tail_of(vert[x].he)!=x)
                  vert[x].he = keep_he;
          };
          ensure_rep(a, h_ab);
          ensure_rep(b, h_bc);
          ensure_rep(c, h_ca);
          return true;
      }
  };
}
  
  std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>
  decimate_internal_vertices(const Eigen::MatrixXd &V,
                             const Eigen::MatrixXi &F,
                             const Eigen::VectorXi &is_decimating)
  {
      // --- build mesh
      assert(is_decimating.size() == V.rows() &&
             "is_decimating must be |V|-long");
      using namespace utils::details;
      Mesh M(V,F);
  
      auto eligible = [&](int u)->bool {
          // Eligible if not boundary OR explicitly requested to decimate.
          return M.vert[u].alive && (!M.vert[u].boundary || (is_decimating[u]!=0));
      };
  
      // Work queue of eligible vertices
      std::queue<int> Q;
      const int nV = (int)V.rows();
      std::vector<char> inQ(nV, 0);
      for(int u=0; u<nV; ++u){
          if(eligible(u)) { Q.push(u); inQ[u]=1; }
      }
  
      // Driver
      while(!Q.empty()){
          int u = Q.front(); Q.pop(); inQ[u]=0;
          if(!eligible(u)) continue;
  
          // Reduce valence via interior flips only
          bool local_progress = true;
          while(M.vert[u].alive && M.vert[u].val>3 && local_progress){
              local_progress = false;
              std::vector<int> nbrs, outhe;
              M.neighbors_ccw(u, nbrs, outhe);
              for(int e : outhe){
                  if(!M.vert[u].alive || M.vert[u].val<=3) break;
                  int et = M.he[e].twin;
                  if(et<0) continue; // cannot flip a boundary edge
                  if(M.flip_edge(e)){
                      local_progress = true;
                      // push a small affected neighborhood
                      int v = M.he[e].v;
                      int a = M.he[M.he[e].next].v;
                      int b = M.he[M.he[et].next].v;
                      for(int x : {u,v,a,b}){
                          if(x>=0 && x<nV && eligible(x) && !inQ[x]){ Q.push(x); inQ[x]=1; }
                      }
                      break;
                  }
              }
          }
  
          // If valence==3, remove (allow even if boundary when requested)
          if(M.vert[u].alive && M.vert[u].val==3){
              bool allow_boundary = (is_decimating[u]!=0);
              // capture neighbors before deletion to re-enqueue
              std::vector<int> nbrs,outhe;
              M.neighbors_ccw(u, nbrs, outhe);
              if((int)nbrs.size()==3 && M.remove_valence3(u, allow_boundary)){
                  for(int x : nbrs){
                      if(x>=0 && x<nV && eligible(x) && !inQ[x]){ Q.push(x); inQ[x]=1; }
                  }
              }
          }
      }
  
      // --- pack outputs (O(N))
      std::vector<int> old2new(nV, -1);
      std::vector<int> keep;
      keep.reserve(nV);
      for(int i=0;i<nV;++i){
          if(M.vert[i].alive){
              old2new[i] = (int)keep.size();
              keep.push_back(i);
          }
      }
  
      std::vector<Eigen::Vector3i> Fout;
      Fout.reserve(M.face.size());
      auto tail_of = [&](int e)->int { return M.he[M.he[e].prev].v; };
      for(size_t fi=0; fi<M.face.size(); ++fi){
          int h = M.face[fi].he;
          if(h<0) continue;
          int e0=h, e1=M.he[e0].next, e2=M.he[e1].next;
          int a = tail_of(e0);
          int b = tail_of(e1);
          int c = tail_of(e2);
          if(a<0||b<0||c<0) continue;
          if(!M.vert[a].alive || !M.vert[b].alive || !M.vert[c].alive) continue;
          if(old2new[a]<0 || old2new[b]<0 || old2new[c]<0) continue;
          Row2 A2=V.row(a), B2=V.row(b), C2=V.row(c);
          if(M.orient(A2,B2,C2)<=0) std::swap(b,c);
          Fout.emplace_back(Eigen::Vector3i(old2new[a], old2new[b], old2new[c]));
      }
  
      Eigen::MatrixXd newV((int)keep.size(), V.cols());
      Eigen::VectorXi new2old((int)keep.size());
      for(int i=0;i<(int)keep.size(); ++i){
          newV.row(i) = V.row(keep[i]);
          new2old[i] = keep[i];
      }
      Eigen::MatrixXi newF((int)Fout.size(), 3);
      for(int i=0;i<(int)Fout.size(); ++i) newF.row(i)=Fout[i];
  
      return {newV, newF, new2old};
  }
  
} // namespace utils
