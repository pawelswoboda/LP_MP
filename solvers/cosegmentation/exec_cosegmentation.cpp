#include "LP_MP.h"
#include "cosegmentation.h"

#include <fstream>
#include <sstream>
#include <string>

using namespace std;
using namespace LP_MP;

std::vector<int> ComputeAssignment(string filename, bool load_neib = true)
{
   Cosegmentation cs;

   fstream fs{filename};
   if(!fs.is_open()) { 
      cout << "Cannot open " << filename << "\n";
      exit(1);
   }

   string LINE;
   REAL cost;
   INDEX N[2] = {0,0};

   // first run through file and record number of left and right nodes
   while(getline(fs,LINE)) {
      stringstream s{LINE};
      char b;
      s >> b;
      if(LINE[0] == 'a') {
         INDEX a, i0, i1;
         s >> a >> i0 >> i1 >> cost;
         N[0] = std::max(N[0], i0);
         N[1] = std::max(N[1], i1);
      } else if( LINE[0] == 'i') {
         INDEX a,b;
         char side;
         s >> side;
         if(side == '1') {
            s >> a >> b >> cost;
            N[0] = std::max(N[0],a);
            N[0] = std::max(N[0],b);
         } else if(side == '2') {
            s >> a >> b >> cost;
            N[1] = std::max(N[1],a);
            N[1] = std::max(N[1],b);
         }
      }
   }

   ++N[0];
   ++N[1];

   cs.SetNumberLeftNodes(N[0]);
   cs.SetNumberRightNodes(N[1]);

   fs.close();
   fs.open(filename);

   while(getline(fs,LINE)) {
      stringstream s{LINE};
      char b;
      s >> b;
      if (b == 'a') {
         INDEX a, i0, i1;
         s >> a >> i0 >> i1 >> cost;
         if(s.fail() || i0>N[0] || i1>N[1]) { 
            assert(false); throw runtime_error(filename + " wrong format2!\n"); 
         }
         cs.AddAssignmentCost(i0,i1,cost);

      } else if( LINE[0] == 'i') {
         INDEX a,b;
         char side;
         s >> side;
         if(side == '1') {
            s >> a >> b >> cost;
            if(s.fail() || a>N[0] || b>N[0] || a==b) { 
               assert(false); 
               throw std::runtime_error("Wrong format5l"); 
            }
            cs.AddLeftPottsTerm(a,b,cost);
         } else if(side == '2') {
            s >> a >> b >> cost;
            if(s.fail() || a>N[0] || b>N[0] || a==b) { 
               assert(false); 
               throw std::runtime_error("Wrong format5r"); 
            }
            cs.AddRightPottsTerm(a,b,cost);
         } else {
            assert(false);
            throw std::runtime_error("Wrong format for Potts");
         }
      }
   }

   fs.close();
      
   cout << "problem from " << filename << " loaded\n";

   cs.Solve(750);
   std::cout << "Primal bound = " << cs.primalBound() << std::endl;

   auto leftSeg = cs.GetLeftSegmentation();
   for(INDEX i=0; i<leftSeg.size(); ++i) {
      std::cout << leftSeg[i] << ",";
   }
   std::cout << "\n";

   auto rightSeg = cs.GetRightSegmentation();
   for(INDEX i=0; i<rightSeg.size(); ++i) {
      std::cout << rightSeg[i] << ",";
   }
   std::cout << "\n";

   return std::vector<int>(0);
}

int main(int argc, char * argv[])
{
   Cosegmentation cs;

   // test problem
   /*
   cs.AddLeftNode(0).AddLeftNode(1).AddRightNode(0).AddRightNode(1);
   cs.AddLeftPottsTerm(0,1,1.0);
   cs.AddRightPottsTerm(0,1,1.0);

   cs.AddAssignmentCost(0,0,0.0);
   cs.AddAssignmentCost(0,1,-2.0);
   cs.AddAssignmentCost(1,0,-3.0);
   cs.AddAssignmentCost(1,1,-4.0);
   cs.Solve(1000);
   */

   if(argc != 2) {
      std::cout << "1 argument must be supplied" << std::endl;
      exit(1);
   }
   std::string fileName(argv[1]);
   std::cout.precision(9);

   std::vector<int> assignment = ComputeAssignment(fileName, true);
   return 0;
}
