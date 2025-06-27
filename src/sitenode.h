/*
 *  sitenode.h
 *
 *  Created by Mark Kirkpatrick on 7/28/08, modified to include CSV output,
 *  additional functions for snptree compatibility, and branch-length calculations.
 *
 *  The SiteNode class represents a node in the gene tree for a site.
 *  It contains the node’s attributes (such as nodeNumber, time, branch lengths,
 *  context, and base states) and pointers to descendant nodes.
 *
 *  New functionality:
 *    - gatherAllNodes(): recursively collects pointers to all nodes in the gene tree.
 *    - printAllNodesSorted(): sorts the collected nodes and writes them to a CSV file.
 *    - writeCSV(): a wrapper that gathers nodes and writes CSV output.
 *
 *  Additional functions (required by snptree.h):
 *    - getTotalLength(), getTotalLength_Informative(), getSNPhit(), outputSNPs(), and getBranchLength().
 *
 *  Also added:
 *    - calcBranchLengths(int) and calcBranchLengths_informative(int) to compute branch lengths.
 *    - Private helper functions: setBranch(double) and setBranch_informative(double)
 *
 *  CSV columns:
 *    nodeNumber,time,sitePosition,branchL,branchL_Informative,
 *    context_pop,context_inv,original_context_pop,original_context_inv,
 *    baseStates,descendantNodeNumbers
 *
 */

 #ifndef SITENODE_
 #define SITENODE_
 
 #include "typedefs.h"
 #include "argnode.h"
 #include <vector>
 #include <string>
 #include <memory>
 
 using std::vector;
 using std::string;
 using std::shared_ptr;
 
 struct snpHit;
 
 class SiteNode {
 public:
	 // Member variables for node attributes:
	 double sitePosition;         // Position of the site (typically in Morgans)
	 unsigned long nodeNumber;    // Unique node identifier
	 double time;                 // Time (age) of this node
	 double branchL;              // Branch length from its ancestor to this node
	 double branchL_Informative;  // Informative branch length (e.g., excluding terminal segments)
	 Context context;             // Current context (e.g., population, inversion status)
	 Context original_context;    // Original context (saved at node creation)
	 vector< shared_ptr<SiteNode> > descendant; // Pointers to descendant SiteNodes
	 vector<Base> baseState;      // Base states along the ancestral branch
 
	 // Constructors and destructor:
	 SiteNode();
	 SiteNode(double sitePos, shared_ptr<ARGNode> argRoot);
	 ~SiteNode();
 
	 // Mutation and branch-length functions:
	 vector<Base> mutate_jc(Base base0, double mu, double t);
	 void mutateDescendants(Base base0, double mu);
	 void putBaseStates(vector<Base> baseVec);
	 // These functions are referenced by snptree.h:
	 int calcBranchLengths(int nCall);
	 int calcBranchLengths_informative(int nCall);
 
	 // Accessors:
	 unsigned long getNodeNumber() { return nodeNumber; }
	 double getTime() { return time; }
	 Context getContext() { return context; }
	 Base getBaseState() { return baseState.back(); }
 
	 // New functions for CSV export:
	 void gatherAllNodes(vector<SiteNode*>& allNodes);
	 void printAllNodesSorted(const vector<SiteNode*>& allNodes, const string &filename);
	 void writeCSV(const string &filename);
 
	 // Additional functions required by snptree.h:
	 double getTotalLength(double runtot);
	 double getTotalLength_Informative(double runtot);
	 snpHit getSNPhit(double target, snpHit x);
	 vector<unsigned> outputSNPs(vector<unsigned> s);
	 double getBranchLength();
 
	 // Optional debugging output functions:
	 void printAllData(std::ostream &out, int depth);
	 void printTreeStructure(std::ostream &out, const string &prefix, bool isLast);
 
 private:
	 // Private helper functions for branch length calculations.
	 void setBranch(double ancTime);
	 void setBranch_informative(double ancTime);
 };
 
 struct snpHit {
	 double runningTot;
	 bool found;
	 shared_ptr<SiteNode> chosen;
 };
 
 #endif /* SITENODE_ */
 