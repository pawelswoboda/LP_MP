#ifndef LP_MP_CUT_PACKING_HXX
#define LP_MP_CUT_PACKING_HXX

#include <vector>
#include <array>
#include <algorithm>
#include "two_dimensional_variable_array.hxx"
#include "union_find.hxx"
#include "help_functions.hxx"

namespace LP_MP {

	// possibly templatize for support for sister pointers and masking edges out
	template<typename EDGE_INFORMATION, bool SUPPORT_SISTER=true, bool SUPPORT_MASKING=false>
	class graph {
		using type = graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>;
		public:
		class edge_type {
			public:
			const std::size_t tail() const { return sister_->head(); }
			const std::size_t head() const { return head_; }
			const edge_type& sister() const { return *sister_; }

			EDGE_INFORMATION& edge() { return edge_; }
			const EDGE_INFORMATION& edge() const { return edge_; }

			bool operator<(const edge_type& o) const { return head() < o.head(); }
			edge_type& operator=(const EDGE_INFORMATION& o) { edge() = o; }

			std::size_t head_;
			edge_type* sister_;
			EDGE_INFORMATION edge_;
		};

		// compute sorted adjacency list representation of graph given by edges.
		// TO DO: initialize EDGE_INFORMATION from EDGE_ITERATOR (by additional optional lambda?)
		template<typename EDGE_ITERATOR>
		graph(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end)
		{
			std::vector<std::size_t> adjacency_list_count;
			// first determine size for adjacency_list
			for(auto edge_it=edge_begin; edge_it!=edge_end; ++edge_it) {
				const auto i = (*edge_it)[0];
				const auto j = (*edge_it)[1];
				adjacency_list_count.resize(std::max({i+1,j+1,adjacency_list_count.size()}));
				adjacency_list_count[i]++;
				adjacency_list_count[j]++; 
			}

			edges_.resize(adjacency_list_count.begin(), adjacency_list_count.end());
			std::fill(adjacency_list_count.begin(), adjacency_list_count.end(), 0);

			for(auto edge_it=edge_begin; edge_it!=edge_end; ++edge_it) {
				const auto i = (*edge_it)[0];
				const auto j = (*edge_it)[1];
				edges_[i][adjacency_list_count[i]].head_ = j;
				adjacency_list_count[i]++;
				edges_[j][adjacency_list_count[j]].head_ = i;
				adjacency_list_count[j]++;
			}

#pragma omp parallel for schedule(guided)
			for(std::size_t i=0; i<adjacency_list_count.size(); ++i) {
				std::sort(edges_[i].begin(), edges_[i].end());
			}

			// set sister pointers
			std::fill(adjacency_list_count.begin(), adjacency_list_count.end(), 0);
			for(std::size_t i=0; i<edges_.size(); ++i) {
				for(auto edge_it=edges_[i].begin(); edge_it!=edges_[i].end(); ++edge_it) {
					if(edge_it->head() > i) {
						const auto head = edge_it->head();
						auto* sister = &edges_(head, adjacency_list_count[head]++);
						edge_it->sister_ = sister;
						sister->sister_ = &(*edge_it);
					} 
				} 
			}

			// check that graph is simple
			for(std::size_t i=0; i<edges_.size(); ++i) {
				assert(std::unique(edges_[i].begin(), edges_[i].end(), [](const edge_type& e1, const edge_type& e2){ return e1.head() == e2.head(); }) == edges_[i].end());
			}

			// check that all sister pointers have been set correctly 
			for(std::size_t i=0; i<edges_.size(); ++i) {
				for(std::size_t e=0; e<edges_[i].size(); ++e) {
					assert(edges_(i,e).tail() == i);
					assert(edges_(i,e).sister_->sister_ == &edges_(i,e));
				}
			}
		}

		template<typename LAMBDA>
		void for_each_edge(LAMBDA f) const
		{
			for(std::size_t i=0; i<edges_.size(); ++i) {
				for(auto edge_it=edges_[i].begin(); edge_it!=edges_[i].end(); ++edge_it) {
					const std::size_t j = edge_it->head();
					if(i<j) {
						f(i, edge_it->head(), edge_it->edge());
					}
				}
			}

		}

		bool edge_present(const std::size_t i, const std::size_t j) const
		{
			assert(i != j && std::max(i,j) < no_nodes());
			struct comparator {
				bool operator()(const std::size_t i, const edge_type& e) const { return i < e.head(); }
				bool operator()(const edge_type& e, const std::size_t i) const { return e.head() < i; }
			};
			auto edge_it = std::lower_bound(begin(i), end(i), j, comparator{});
			return edge_it != end(i);

		}
		const EDGE_INFORMATION& edge(const std::size_t i, const std::size_t j) const
		{
			assert(edge_present(i,j));
			struct comparator {
				bool operator()(const std::size_t i, const edge_type& e) const { return i < e.head(); }
				bool operator()(const edge_type& e, const std::size_t i) const { return e.head() < i; }
			};
			auto edge_it = std::lower_bound(begin(i), end(i), j, comparator{});
			assert(edge_it->head() == j);
			return edge_it->edge(); 
		}

		std::size_t no_nodes() const { return edges_.size(); }
		std::size_t no_edges(const std::size_t i) const { return edges_[i].size(); }
		std::size_t edge_no(edge_type* e) const
		{
			const std::size_t i = e->tail();
			assert(std::distance(edges_[i].begin(), e) < no_edges(i));
			return std::distance(edges_[i].begin(), e); 
		}
		std::size_t edge_no(const edge_type& e) const
		{
			const std::size_t i = e.tail();
			assert(std::distance(const_cast<const edge_type*>(edges_[i].begin()), &e) < no_edges(i));
			return std::distance(const_cast<const edge_type*>(edges_[i].begin()), &e); 
		}
		
		auto begin(const std::size_t i) const { return edges_[i].begin(); }
		auto end(const std::size_t i) const { return edges_[i].end(); }

		// enumerate all triangles and call f on each. f expects the three node indices (i<j<k) of triangles (sorted) and references to edges in lexicographical order (ij,ik,jk)
		template<typename LAMBDA>
		void for_each_triangle(LAMBDA f) const
		{
			struct triangle_intersection_type {
				triangle_intersection_type(std::size_t _k, edge_type* _ik, edge_type* _jk) : k(_k), ik(_ik), jk(_jk) {}
				std::size_t k;
				edge_type* ik;
				edge_type* jk;
			};
			auto edge_intersection_merge = [](edge_type& e1, edge_type& e2) -> triangle_intersection_type { 
				assert(e1.head() == e2.head());
				const auto k = e1.head();
				return triangle_intersection_type(k, &e1, &e2);
			};
			auto edge_intersection_sort = [](const edge_type& e1, const edge_type& e2) { return e1.head() < e2.head(); };

			std::vector<triangle_intersection_type> common_nodes;
			for(std::size_t i=0; i<no_nodes(); ++i) {
				for(auto edge_it=begin(i); edge_it!=end(i); ++edge_it) {
					const auto j = edge_it->head();
					if(i<j) {
						// Now find all neighbors of both i and j to see where the triangles are
						common_nodes.clear();
						auto intersects_iter_end = set_intersection_merge(
								edges_[i].begin(), edges_[i].end(), edges_[j].begin(), edges_[j].end(),
							       	std::back_inserter(common_nodes), edge_intersection_sort, edge_intersection_merge);

						for(const triangle_intersection_type t : common_nodes) {
							const auto k = t.k;

							// Since a triplet shows up three times as an edge plus a node, we only consider it for the case when i<j<k 
							if(!(j<k)) { continue; }

							assert(t.ik->tail() == i && t.jk->tail() == j);
							f(i,j,k, *edge_it, *(t.ik), *(t.jk));
						}
					}
				}
			}
		}

		template<typename ENDPOINT_ITERATOR, typename LAMBDA>
		void enumerate_2_hops(ENDPOINT_ITERATOR endpoint_begin, ENDPOINT_ITERATOR endpoint_end, LAMBDA f)
		{
			// do roughly the same as for enumerating triangles.
		}

		// we enumerate quadrangles with the method C4 from "Arboricity and subgraph listin algorithms" by Norishige Chiba and Takao Nishizeki
		template<typename LAMBDA>
		void for_each_quadrangle(LAMBDA f)
		{
			struct node_degree { std::size_t i, degree; };
			std::vector<node_degree> node_degrees(no_nodes());
			for(std::size_t i=0; i<no_nodes(); ++i) {
				node_degrees[i].i = i;
				node_degrees[i].degree = no_edges(i);
			}
			std::sort(node_degrees.begin(), node_degrees.end(), [](const node_degree& n1, const node_degree& n2) { return n1.degree > n2.degree; });

			// structure for recording masked edges
			two_dim_variable_array<unsigned char> edge_mask(edges_);
			edge_mask.set(true);

			std::vector<std::vector<std::size_t>> U(no_nodes());
			std::vector<std::size_t> relevant_U_entries;
			std::vector<std::size_t> nonempty_U_entries;

			for(const auto& n : node_degrees) {
				const std::size_t v = n.i;
				for(std::size_t c_v=0; c_v<no_edges(v); ++c_v) {
					if(edge_mask(v,c_v) == false) { continue; }
					const auto u = edges_(v,c_v).head();
					assert(u != v);
					for(std::size_t c_u=0; c_u<no_edges(u); ++c_u) {
						if(edge_mask(u,c_u) == false) { continue; }
						const std::size_t w = edges_(u,c_u).head();
						if(w != v) {
							U[w].push_back(u); 
							if(U[w].size() == 1) { nonempty_U_entries.push_back(w); }
							if(U[w].size() == 2) { relevant_U_entries.push_back(w); }
						} 
					}
				}
				for(const std::size_t w : relevant_U_entries) {
					assert(U[w].size() >= 2);
					for(std::size_t i=0; i<U[w].size(); ++i) {
						for(std::size_t j=i+1; j<U[w].size(); ++j) {
							f(v, w, U[w][i], U[w][j]);
						}
					}
				}
				relevant_U_entries.clear();
				for(const auto w : nonempty_U_entries) {
					U[w].clear();
				}
				nonempty_U_entries.clear();
				for(const auto& l : U) { assert(l.empty()); } 

				// mask all edges with v as its endpoint
				for(std::size_t c=0; c<no_edges(v); ++c) {
					const auto u = edges_(v,c).head();
					edge_mask(v,c) = false;
					edge_mask(u, edge_no( edges_(v,c).sister() )) = false;
				}
			}
		}


		// return contracted graph and mapping from original nodes to contracted nodes
		template<typename EDGE_ITERATOR, typename MERGE_FUNC>
		std::tuple<graph, std::vector<std::size_t>> contract(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, MERGE_FUNC merge_func) const
		{
			union_find uf(no_nodes());
			for(auto edge_it=edge_begin; edge_it!=edge_end; ++edge_it) {
				const auto i = (*edge_it)[0];
				const auto j = (*edge_it)[1];
				uf.merge(i,j);
			}
			auto contracted_ids = uf.get_contiguous_ids();
			struct contracted_edge_type : public std::array<std::size_t,2> { 
				contracted_edge_type(const std::size_t i, const std::size_t j, EDGE_INFORMATION e) : std::array<std::size_t,2>({i,j}), edge(e) {}
				EDGE_INFORMATION edge; 
			};
			std::vector<contracted_edge_type> contracted_edges;

			for_each_edge( [&](const std::size_t i, const std::size_t j, const auto& edge) {
					const auto contracted_i = contracted_ids[ uf.find(i) ];
					const auto contracted_j = contracted_ids[ uf.find(j) ];
					if(contracted_i != contracted_j) {
						contracted_edge_type contracted_edge(contracted_i, contracted_j, edge);
						contracted_edges.push_back(contracted_edge); 
					}
					});

			// TO DO: faster sorting by first sorting for tail ids via bucket sort, and then separately sorting for head ids.
			// merge contracted edges
			auto edge_sort = [](const auto e1, const auto e2) { return std::lexicographical_compare(e1.begin(), e1.end(), e2.begin(), e2.end()); };
			std::sort(contracted_edges.begin(), contracted_edges.end(), edge_sort);
			std::vector<contracted_edge_type> unique_contracted_edges;
			for(std::size_t i=0; i<contracted_edges.size()-1; ++i) {
				if(contracted_edges[i][0] == contracted_edges[i+1][0] && contracted_edges[i][0] == contracted_edges[i+1][0]) {
					contracted_edges[i+1].edge = merge_func(contracted_edges[i].edge, contracted_edges[i+1].edge);
				} else {
					unique_contracted_edges.push_back(contracted_edges[i]);
				} 
			}
			unique_contracted_edges.push_back(contracted_edges.back()); 

			type contracted_graph(unique_contracted_edges.begin(), unique_contracted_edges.end());
			
			std::vector<std::size_t> contraction_mapping(no_nodes());
			for(std::size_t i=0; i<no_nodes(); ++i) {
				contraction_mapping[i] = contracted_ids[ uf.find(i) ];
			}

			return std::make_tuple(contracted_graph, contraction_mapping); 
		}

		template<typename EDGE_ITERATOR>
		auto contract(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end) const
		{
			auto standard_merge_func = [](EDGE_INFORMATION& e1, EDGE_INFORMATION& e2) { return e1; };
			return contract(edge_begin, edge_end, standard_merge_func);
		}


		// bfs search
		/*
		struct bfs_data {
			struct item { std::size_t parent; std::size_t flag; };
			bfs_data(const graph& g) 
			{
				d.resize(g.size());
				for(INDEX i=0; i<d.size(); ++i) {
					d[i].flag = 0;
				}
				flag1 = 0;
				flag2 = 1;
			}
			void reset() 
			{
				visit.clear();
				flag1 += 2;
				flag2 += 2; 
			}
			item& operator[](const INDEX i) { return d[i]; }

			void label1(const INDEX i) { d[i].flag = flag1; }
			void label2(const INDEX i) { d[i].flag = flag2; }
			bool labelled(const INDEX i) const { return Labelled1(i) || Labelled2(i); }
			bool labelled1(const INDEX i) const { return d[i].flag == flag1; }
			bool labelled2(const INDEX i) const { return d[i].flag == flag2; }

			std::size_t& parent(const INDEX i) { return d[i].parent; }
			std::size_t parent(const INDEX i) const { return d[i].parent; }

			std::tuple<REAL,std::vector<INDEX>> TracePath(const INDEX i1, const INDEX i2, const REAL cost_i1i2) const 
			{
				REAL c = cost_i1i2;
				std::vector<INDEX> path({i1});
				INDEX j=i1;
				while(Parent(j) != j) {
					c = std::min(c, Cost(j) );
					j = Parent(j);
					path.push_back(j);
				}
				std::reverse(path.begin(),path.end());
				path.push_back(i2);
				j=i2;
				while(Parent(j) != j) {
					c = std::min(c, Cost(j) );
					j = Parent(j);
					path.push_back(j);
				}

				return std::move(std::make_tuple(c,std::move(path)));
			}

			//static auto no_mask_op = [](const INDEX i, const INDEX j, const REAL weight) { return true; };
			static bool no_mask_op(const INDEX, const INDEX, const REAL) { return true;}

			// do bfs with thresholded costs and iteratively lower threshold until enough cycles are found
			// only consider edges that have cost equal or larger than th
			std::tuple<REAL, std::vector<INDEX>> FindPath(const INDEX startNode, const INDEX endNode, const Graph& g, const REAL th = 0.0)
			{
				return FindPath(startNode, endNode, g, th, no_mask_op);
			}

			template<typename MASK_OP>
				std::tuple<REAL,std::vector<INDEX>> 
				FindPath(
						const INDEX startNode, const INDEX endNode, const Graph& g, const REAL th,
						MASK_OP mask_op
					)
				{
					assert(startNode != endNode);
					assert(startNode < g.size() && endNode < g.size());
					Reset();
					visit.push_back({startNode, 0});
					Label1(startNode);
					Parent(startNode) = startNode;
					Cost(startNode) = std::numeric_limits<REAL>::infinity();
					visit.push_back({endNode, 0});
					Label2(endNode);
					Parent(endNode) = endNode;
					Cost(endNode) = std::numeric_limits<REAL>::infinity();

					while(!visit.empty()) {
						const INDEX i = visit.front()[0];
						const INDEX distance = visit.front()[1];
						visit.pop_front();

						assert(g[i].begin() < g[i].end());

						if(Labelled1(i)) {
							for(auto* a=g[i].begin(); a->cost>=th && a!=g[i].end(); ++a) { 
								auto* head = a->head;
								const INDEX j = g[head];

								if(mask_op(i,j,a->cost)) {

									if(!Labelled(j)) {
										visit.push_back({j, distance+1});
										Parent(j) = i;
										Cost(j) = a->cost;
										Label1(j);
									} else if(Labelled2(j)) { // shortest path found
										// trace back path from j to endNode and from i to startNode
										return std::move(TracePath(i,j, a->cost));
									}

								}
							}
						} else {
							assert(Labelled2(i));
							for(auto* a=g[i].begin(); a->cost>=th && a!=g[i].end(); ++a) { 
								auto* head = a->head;
								const INDEX j = g[head];

								if(mask_op(i,j,a->cost)) {

									if(!Labelled(j)) {
										visit.push_back({j, distance+1});
										Parent(j) = i;
										Cost(j) = a->cost;
										Label2(j);
									} else if(Labelled1(j)) { // shortest path found
										// trace back path from j to endNode and from i to startNode
										return std::move(TracePath(i,j, a->cost));
									}

								}
							}
						}


					}
					return std::make_tuple(-std::numeric_limits<REAL>::infinity(),std::vector<INDEX>(0));
				}

			private:
			std::vector<item> d;
			std::deque<std::array<std::size_t,2>> visit; // node number, distance from start or end 
			std::size_t flag1, flag2;
		};
	*/

		private:
		two_dim_variable_array<edge_type> edges_;
	};

} // namespace LP_MP

#endif // LP_MP_CUT_PACKING_HXX
