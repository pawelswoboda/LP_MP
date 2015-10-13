#include "LP_MP.h"
#include "graph_matching.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <exception>
#include <stdio.h>
#include <algorithm>

using namespace std;
using namespace LP_MP;

std::vector<int> ComputeAssignment(string filename, bool load_neib = true)
{
   std::vector<std::vector<INDEX> >  leftGraph, rightGraph;
   std::vector<std::tuple<INDEX,INDEX,REAL> > leftAssignment, rightAssignment;
   GraphMatching gm;

   fstream fs{filename};
   if(!fs.is_open()) { 
      cout << "Cannot open " << filename << "\n";
      exit(1);
   }

   string LINE;
   INDEX N[2], A, E, _A = 0, _E = 0;
   REAL cost;
   while(getline(fs,LINE)) {
      stringstream s{LINE};
      char b;
      s >> b;
      if (b == 'p') {
         s >> N[0] >> N[1] >> A >> E;
         if(s.fail()) {
            cout << filename << ": wrong format1: " << N[0] << ", " << N[1] << ", " << A << ", " << E << "\n";
            exit(1); 
         }
         leftGraph.resize(N[0]);
         rightGraph.resize(N[1]);

         for(INDEX i=0; i<leftGraph.size(); i++) gm.AddLeftNode(i,0,1);
         for(INDEX i=0; i<rightGraph.size(); i++) gm.AddRightNode(i,0,1);

         leftAssignment.resize(A);
         rightAssignment.resize(A);

      } else if (b == 'a') {
         INDEX a, i0, i1;
         s >> a >> i0 >> i1 >> cost;
         if(s.fail() || a!=_A++ || _A>A || i0>=N[0] || i1>=N[1]) { 
            assert(false); throw runtime_error(filename + " wrong format2!\n"); 
         }
         gm.AddAssignmentCost(i0,i1,cost);

         leftAssignment[a] = make_tuple(i0,i1,0.5*cost);
         rightAssignment[a] = make_tuple(i1,i0,0.5*cost);
      } else if (LINE[0] == 'e') {
         if(LINE[1] ==  ' ') {
            INDEX a, b;
            s >> a >> b >> cost;
            if (s.fail() || a>=A || b>=A || a==b) {
               assert(false); 
               cout << filename << ": wrong format3!\n"; 
               exit(1); 
            }
            assert(std::abs(cost) < 1000000);

            const INDEX l1 = std::get<0>(leftAssignment[a]);
            const INDEX r1 = std::get<1>(leftAssignment[a]);
            const INDEX l2 = std::get<0>(leftAssignment[b]);
            const INDEX r2 = std::get<1>(leftAssignment[b]);
            assert(l1!=l2);
            assert(r1!=r2);
            gm.AddPairwiseCost(l1,r1,l2,r2,cost);
         } else if(LINE[1] == '2') {
            INDEX l1, l2, r1, r2;
            if(LINE[2] == 'l') {
               char quad_pot;
               s >> quad_pot; // == '2'
               char side;
               s >> side;
               s >> l1 >> r1 >> l2 >> r2 >> cost;
               gm.AddLeftPairwiseCost(l1,r1,l2,r2,cost);
            } else if(LINE[2] == 'r') {
               char quad_pot;
               s >> quad_pot; // == '2'
               char side;
               s >> side;
               s >> l1 >> r1 >> l2 >> r2 >> cost;
               gm.AddRightPairwiseCost(r1,l1,r2,l2,cost);
            } else if(LINE[2] == ' ') {
               char quad_pot;
               s >> quad_pot;
               s >> l1 >> r1 >> l2 >> r2 >> cost;
               gm.AddPairwiseCost(l1,r1,l2,r2,cost);
            }
            if(s.fail() || l1>=N[0] || l2>=N[0] || r1>=N[1] || r2>=N[1] ) { 
               assert(false); 
               std::cout << l1 << ", " << r1 << ", " << l2 << ", " << r2 << ", cost = " << cost << "\n";
               cout << filename << ": wrong format3!\n";
               exit(1); 
            }
         } else if(LINE[1] == '3') {
            assert(false); // not yet ready
            INDEX a, b, c;
            s >> a >> b >> c >> cost;
            if (s.fail() || a>=A || b>=A || c>=A || a==b) { 
               assert(false); 
               cout << filename << ": wrong format3!\n";
               exit(1); 
            }
            //assert(cost != 0.0);

            const INDEX leftNode1 = std::get<0>(leftAssignment[a]);
            const INDEX rightNode1 = std::get<1>(leftAssignment[a]);
            const INDEX leftNode2 = std::get<0>(leftAssignment[b]);
            const INDEX rightNode2 = std::get<1>(leftAssignment[b]);
            const INDEX leftNode3 = std::get<0>(leftAssignment[c]);
            const INDEX rightNode3 = std::get<1>(leftAssignment[c]);
            gm.AddTernaryCost(leftNode1,leftNode2,leftNode3,rightNode1,rightNode2,rightNode3,cost);
         } else {
            throw std::runtime_error("Wrong format6");
         }
      } else if( LINE[0] == 'i') {
         INDEX a,b;
         char side;
         s >> side;
         if(side == '1') {
            s >> a >> b >> cost;
            if(s.fail() || a>=N[0] || b>=N[0] || a==b) { 
               assert(false); 
               throw std::runtime_error("Wrong format5l"); 
            }
            gm.AddLeftPottsTerm(a,b,cost);
            std::cout << "Potts cost = " << cost << std::endl;
         } else if(side == '2') {
            s >> a >> b >> cost;
            if(s.fail() || a>=N[0] || b>=N[0] || a==b) { 
               assert(false); 
               throw std::runtime_error("Wrong format5r"); 
            }
            gm.AddRightPottsTerm(a,b,cost);
            std::cout << "Potts cost = " << cost << std::endl;
         } else {
            assert(false);
            throw std::runtime_error("Wrong format for Potts");
         }
      } else if (LINE[0] == 'n') {
         throw std::runtime_error("neighbors not supported yet");
         if (load_neib)
         {
            INDEX r, i, j;
            s >> r >> i >> j;
            if (s.fail() || r>1 ||  i>=N[r] || j>=N[r] || i==j) { 
               cout << filename << ": wrong format4!\n";
               exit(1); 
            }
         }
      }



   }

   fs.close();
      
   cout << "problem from " << filename << " loaded (N0=" << N[0] << ", N1=" << N[1] << ", A=" << A << ", E=" << E << ")\n";

   gm.Solve(750);
   std::cout << "Primal bound = " << gm.primalBound() << std::endl;

   return gm.GetLeftAssignment();
}

int main(int argc, char * argv[])
{
   if(argc != 2) {
      std::cout << "1 argument must be supplied" << std::endl;
      exit(1);
   }
   std::string fileName(argv[1]);
   std::cout.precision(9);

   std::vector<int> assignment = ComputeAssignment(fileName, true);

   std::string file_out(argv[1]);
   //std::replace(file_out.begin(), file_out.end(), ".", "_output.");
   file_out.replace(file_out.find_last_of("."), 1, "_assignment.");
   std::cout << "Writing assignment to file " << file_out << std::endl;
   std::ofstream file_out_stream(file_out);
   for(INDEX i=0; i<assignment.size(); i++) {
      file_out_stream << assignment[i] << " ";
   }
   file_out_stream << std::endl;
   file_out_stream.close();
}
