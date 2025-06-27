/*
 *  siteNode.cpp
 *
 *  Created by Mark Kirkpatrick on 7/28/08.
 *  Modified to include CSV output, additional functions for snptree compatibility,
 *  and branch-length calculation functions.
 *
 *  See "sitenode.h" for documentation.
 */

 #include <iostream>
 using std::cout;
 using std::endl;
 #include "sitenode.h"
 #include "argnode.h"
 #include "ran_mk.h"
 #include <fstream>
 #include <sstream>
 #include <algorithm>
 #include <limits>
 #include <cstdlib>
 
 //---------------------------------------------------------------------
 // Helper function: returns true if sitePos lies within any segment in segList.
 static bool siteQuery(double sitePos, vector<Segment> segList)
 {
     for (int i = 0; i < segList.size(); ++i)
     {
         if ((segList.at(i).L < sitePos && sitePos < segList.at(i).R) ||
             (segList.at(i).L == sitePos))
         {
             return true;
         }
     }
     return false;
 }
  
 //---------------------------------------------------------------------
 // Descend the ARG until a node is reached that has more than one descendant carrying the site.
 static shared_ptr<ARGNode> getNextSiteNode(double sitePos, shared_ptr<ARGNode> argRoot)
 {
     int nSiteDesc = 0;
     int uniqueDesc = 0;
     unsigned long nARGDesc = argRoot->getNDescendants();
     for (unsigned long i = 0; i < nARGDesc; i++) {
         if (siteQuery(sitePos, argRoot->getDescendantSegmentVector(i))) {
             nSiteDesc++;
             uniqueDesc = i;
         }
     }
     if (nSiteDesc > 1) {
         return argRoot;
     } else if (nSiteDesc == 1 && argRoot->getTime() != 0) {
         return getNextSiteNode(sitePos, argRoot->getDescendantNode(uniqueDesc));
     }
     return argRoot;
 }
  
 //---------------------------------------------------------------------
 // Null constructor.
 SiteNode::SiteNode()
 {
 }
  
 //---------------------------------------------------------------------
 // Constructor for a gene tree at a given site position, starting at an ARG node.
 SiteNode::SiteNode(double sitePos, shared_ptr<ARGNode> argRoot)
 {
     // Descend the ARG to find the appropriate node for this site.
     shared_ptr<ARGNode> argSiteRoot = getNextSiteNode(sitePos, argRoot);
      
     // Set basic attributes.
     sitePosition = sitePos;
     nodeNumber = argSiteRoot->getNodeNumber();
     time = argSiteRoot->getTime();
     context = argSiteRoot->getContext();
     // Save the original context.
     original_context = context;
      
     // Build descendant SiteNodes for every branch carrying this site.
     unsigned long argSiteRootNDesc = argSiteRoot->getNDescendants();
     for (unsigned long i = 0; i < argSiteRootNDesc; i++) {
         if (siteQuery(sitePos, argSiteRoot->getDescendantSegmentVector(i))) {
             shared_ptr<SiteNode> d(new SiteNode(sitePos, argSiteRoot->getDescendantNode(i)));
             descendant.push_back(d);
         }
     }
 }
  
 //---------------------------------------------------------------------
 // Destructor.
 SiteNode::~SiteNode()
 {
     // (Optional cleanup or debug output)
 }
  
 //---------------------------------------------------------------------
 // Returns a vector of base states along the branch using the Jukes-Cantor model.
 vector<Base> SiteNode::mutate_jc(Base base0, double mu, double t)
 {
     vector<Base> baseVec(1, base0);
     double timeLeft = t;
     do {
         double waitTime = randexp(mu);
         if (timeLeft < waitTime)
             return baseVec;
         else {
             int newBase = (baseVec.back() + randint(1, 3)) % 4;
             baseVec.push_back(Base(newBase));
         }
         timeLeft -= waitTime;
     } while (timeLeft > 0);
     return baseVec;
 }
  
 //---------------------------------------------------------------------
 // Put mutations on descendant branches.
 void SiteNode::mutateDescendants(Base base0, double mu)
 {
     static unsigned long nCalls;
     if (nCalls == 0) {
         baseState.push_back(base0);
         nCalls++;
     }
     for (unsigned long iDesc = 0; iDesc < descendant.size(); iDesc++) {
         vector<Base> baseVec = mutate_jc(base0, mu, time - descendant[iDesc]->getTime());
         descendant[iDesc]->putBaseStates(baseVec);
         descendant[iDesc]->mutateDescendants(descendant[iDesc]->getBaseState(), mu);
     }
 }
  
 //---------------------------------------------------------------------
 // Set the base states.
 void SiteNode::putBaseStates(vector<Base> baseVec)
 {
     baseState = baseVec;
 }
  
 //---------------------------------------------------------------------
 // Print all data of this SiteNode (structured for debugging).
 void SiteNode::printAllData(std::ostream &out, int depth)
 {
     std::string indent(depth * 4, ' ');
     out << indent << "# SiteNode:\n";
     out << indent << "    ├─ sitePosition: " << sitePosition << "\n";
     out << indent << "    ├─ nodeNumber: " << nodeNumber << "\n";
     out << indent << "    ├─ time: " << time << "\n";
     out << indent << "    ├─ branchL: " << branchL << "\n";
     out << indent << "    ├─ branchL_Informative: " << branchL_Informative << "\n";
     out << indent << "    ├─ context: { pop=" << context.pop 
         << ", inv=" << context.inversion << " }\n";
     out << indent << "    ├─ original_context: { pop=" << original_context.pop
         << ", inv=" << original_context.inversion << " }\n";
     out << indent << "    ├─ baseState(s): [";
     for (size_t i = 0; i < baseState.size(); ++i) {
         out << baseState[i];
         if (i < baseState.size() - 1)
             out << ", ";
     }
     out << "]\n";
     out << indent << "    └─ Number of descendants: " << descendant.size() << "\n";
     for (size_t i = 0; i < descendant.size(); ++i) {
         if (descendant[i])
             descendant[i]->printAllData(out, depth + 1);
         else
             out << indent << "    (null descendant)\n";
     }
 }
  
 //---------------------------------------------------------------------
 // Print the tree structure of the gene tree.
 void SiteNode::printTreeStructure(std::ostream &out, const std::string &prefix, bool isLast)
 {
     out << prefix << (isLast ? "└─" : "├─");
     out << "Node #" << nodeNumber << " (time=" << time << ", pop=" 
         << context.pop << ", inv=" << context.inversion << ")\n";
     std::string childPrefix = prefix + (isLast ? "  " : "│ ");
     for (size_t i = 0; i < descendant.size(); ++i) {
         bool lastChild = (i == descendant.size() - 1);
         if (descendant[i])
             descendant[i]->printTreeStructure(out, childPrefix, lastChild);
     }
 }
  
 //---------------------------------------------------------------------
 // Recursively gather pointers to all SiteNode objects in this gene tree.
 void SiteNode::gatherAllNodes(std::vector<SiteNode*>& allNodes)
 {
     allNodes.push_back(this);
     for (size_t i = 0; i < descendant.size(); i++) {
         if (descendant[i])
             descendant[i]->gatherAllNodes(allNodes);
     }
 }
  
 //---------------------------------------------------------------------
 // Write all nodes in the gene tree to a CSV file.
 // CSV columns:
 // nodeNumber,time,sitePosition,branchL,branchL_Informative,context_pop,context_inv,
 // original_context_pop,original_context_inv,baseStates,descendantNodeNumbers
 void SiteNode::printAllNodesSorted(const std::vector<SiteNode*>& allNodes, const std::string &filename)
 {
     std::vector<SiteNode*> sortedNodes = allNodes;
     std::sort(sortedNodes.begin(), sortedNodes.end(), [](SiteNode* a, SiteNode* b) {
          return a->getNodeNumber() < b->getNodeNumber();
     });
  
     std::ofstream outFile(filename.c_str());
     if (!outFile.is_open()) {
          std::cerr << "Error: Could not open file " << filename << " for writing.\n";
          return;
     }
     outFile << "nodeNumber,time,sitePosition,branchL,branchL_Informative,context_pop,context_inv,original_context_pop,original_context_inv,baseStates,descendantNodeNumbers\n";
     for (auto *node : sortedNodes) {
          // Skip dummy nodes
          if (node->getNodeNumber() == std::numeric_limits<unsigned long>::max() || node->getTime() < 0)
              continue;
              
          std::stringstream line;
          line << node->getNodeNumber() << "," 
               << node->getTime() << ","
               << node->sitePosition << ","
               << node->branchL << ","
               << node->branchL_Informative << ","
               << node->getContext().pop << ","
               << node->getContext().inversion << ","
               << node->original_context.pop << ","
               << node->original_context.inversion << ",";
  
          // Concatenate baseStates separated by semicolons.
          {
              std::stringstream baseSS;
              for (size_t i = 0; i < node->baseState.size(); ++i) {
                  baseSS << node->baseState[i];
                  if (i < node->baseState.size() - 1)
                      baseSS << ";";
              }
              line << baseSS.str() << ",";
          }
          // Concatenate descendant node numbers separated by semicolons.
          {
              std::stringstream descSS;
              bool firstDesc = true;
              for (size_t i = 0; i < node->descendant.size(); ++i) {
                  if (node->descendant[i]) {
                      unsigned long descNum = node->descendant[i]->getNodeNumber();
                      if (descNum != std::numeric_limits<unsigned long>::max()) {
                          if (!firstDesc)
                              descSS << ";";
                          descSS << descNum;
                          firstDesc = false;
                      }
                  }
              }
              line << descSS.str();
          }
          outFile << line.str() << "\n";
     }
     outFile.close();
 }
  
 //---------------------------------------------------------------------
 // New wrapper function: writeCSV.
 // This function gathers all nodes from this gene tree and writes them to a CSV file.
 void SiteNode::writeCSV(const std::string &filename)
 {
     std::vector<SiteNode*> nodes;
     gatherAllNodes(nodes);
     printAllNodesSorted(nodes, filename);
 }
  
 //---------------------------------------------------------------------
 // Additional functions required by snptree.h:
  
 // Returns the total length of the gene tree by summing branchL recursively.
 double SiteNode::getTotalLength(double runtot)
 {
     double runTotal = runtot;
     for (size_t i = 0; i < descendant.size(); i++) {
         runTotal += descendant[i]->branchL;
         runTotal = descendant[i]->getTotalLength(runTotal);
     }
     return runTotal;
 }
  
 // Returns the total informative length (using branchL_Informative) recursively.
 double SiteNode::getTotalLength_Informative(double runtot)
 {
     double runTotal = runtot;
     for (size_t i = 0; i < descendant.size(); i++) {
         runTotal += descendant[i]->branchL_Informative;
         runTotal = descendant[i]->getTotalLength_Informative(runTotal);
     }
     return runTotal;
 }
  
 // Recursively finds a node where the running total branch length exceeds the target.
 snpHit SiteNode::getSNPhit(double target, snpHit x)
 {
     snpHit c = x;
     for (size_t i = 0; i < descendant.size(); i++) {
         if (!c.found) {
             c.runningTot += descendant[i]->branchL;
             if (c.runningTot > target) {
                 c.found = true;
                 c.chosen = descendant[i];
                 return c;
             } else {
                 c = descendant[i]->getSNPhit(target, c);
             }
         }
     }
     return c;
 }
  
 // Recursively outputs SNP information. For terminal nodes (time==0), marks index nodeNumber as 1.
 vector<unsigned> SiteNode::outputSNPs(vector<unsigned> s)
 {
     vector<unsigned> sample = s;
     if (time == 0) {
         if (nodeNumber < sample.size())
             sample[nodeNumber] = 1;
     }
     for (size_t i = 0; i < descendant.size(); i++) {
         sample = descendant[i]->outputSNPs(sample);
     }
     return sample;
 }
  
 // A simple getter for branchL.
 double SiteNode::getBranchLength()
 {
     return branchL;
 }
  
 //---------------------------------------------------------------------
 // Private helper: setBranch.
 // Sets branchL as the difference between the provided ancestral time and this node's time.
 void SiteNode::setBranch(double ancTime)
 {
     if (time != -1)
         branchL = ancTime - time;
     else
         branchL = 0;
     if (branchL < 0)
         cout << "Error in setBranch: negative branch length\n";
 }
  
 // Private helper: setBranch_informative.
 // Sets branchL_Informative similarly (using a condition on time).
 void SiteNode::setBranch_informative(double ancTime)
 {
     if (time > 1)
         branchL_Informative = ancTime - time;
     else
         branchL_Informative = 0;
     if (branchL_Informative < 0)
         cout << "Error in setBranch_informative: negative branch length\n";
 }
  
 //---------------------------------------------------------------------
 // Recursively calculates branch lengths for the gene tree.
 // If nCall is 0, sets this node's branchL to 0; then updates descendants.
 int SiteNode::calcBranchLengths(int nCall)
 {
     int nCalls = nCall;
     if(nCalls == 0) {
         branchL = 0;
         nCalls++;
     }
     for (size_t i = 0; i < descendant.size(); i++) {
         descendant[i]->setBranch(time);
         descendant[i]->calcBranchLengths(nCalls);
     }
     return nCalls;
 }
  
 // Recursively calculates informative branch lengths for the gene tree.
 int SiteNode::calcBranchLengths_informative(int nCall)
 {
     int nCalls = nCall;
     if(nCalls == 0) {
         branchL_Informative = 0;
         nCalls++;
     }
     for (size_t i = 0; i < descendant.size(); i++) {
         descendant[i]->setBranch_informative(time);
         descendant[i]->calcBranchLengths_informative(nCalls);
     }
     return nCalls;
 }
  
 // End of sitenode.cpp
 