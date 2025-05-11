/*
 * runningmax.cc
 *
 *  created: Sep 2017
 *   author: Daniel A. 2023
 */

/*
* An adjustment on abstracting both the subsystems and the interconnection based on a data-driven approach.
*/

#include <iostream>
#include <array>
#include <cmath>
#include <algorithm>

#include "gurobi_c++.h"

/* SCOTS header */
#include "scots.hh"
// #include "cuddObj.hh"
// #include "SymbolicSet.hh"
#include "SymbolicSets.hh"
// #include "SymbolicModelGrowthBound.hh"

/* time profiling */
#include "TicToc.hh"
/* memory profiling */
#include <sys/time.h>
#include <sys/resource.h>
struct rusage usage;
using namespace std;

double const samples = 150;
double h_eta = 1.0;
double h_eta2 = h_eta/samples;
double exo_eta = 1; 
double rho_y = 1;
double L_x = 2*1.0;
double L_w = 1.0;
double const bias = L_x*h_eta2 + L_w*exo_eta; // 2*L_x*h_eta2 + L_w*exo_eta

//* state space dim */
const int state_dim=1;
/* input space dim */
const int control_dim=1;
/* exog space dim*/
const int exog_dim = 1;
/* input space of system is a cartesian product*/
const int input_dim = state_dim + control_dim + exog_dim; 
/* Create N identical systems */
const int N = 32;

const int inter_dim = N - 2; 

/*
 * data types for the state space elements and input space
 * elements used ll uniform grid
 */
using state_type = std::array<double,state_dim>;
using control_type = std::array<double,control_dim>;
using exog_type  = std::array<double,exog_dim>;
using input_type  = std::array<double,input_dim>;
using inter_type = std::array<double, inter_dim>;

/* Data types for the interconnection relation */
using prod_state_type = std::array<double, state_dim * N>;
using prod_exog_type = std::array<double, exog_dim * N>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

/* Generate max functions; the interconnection model */ 
  void max2(std::array<double,2> ll, std::array<double, 2> ur, std::array<double, 1> &o_ll, std::array<double,1> &o_ur){
    o_ll[0] = std::max(ll[0], ll[1]); 
    o_ur[0] = std::max(ur[0], ur[1]);
  }

/* Lasso reg approximation of the interconnection model */ 
  void max3(std::array<double,2> ll, std::array<double, 2> ur, std::array<double, 1> &o_ll, std::array<double,1> &o_ur){
    double b1 = 0.4999882697947213;
    double b2 = 0.4999882697947213;
    double L_M = 1.0; // interconnection map lipschitz constant
    double eps_hat1 = L_M*rho_y + fabs(max(ll[0], ll[1]) - b1*(ll[0]+rho_y) - b2*(ll[1]+rho_y));
    double eps_hat2 = L_M*rho_y + fabs(max(ur[0], ur[1]) - b1*(ur[0]+rho_y) - b2*(ur[1]+rho_y));
    o_ll[0] = b1*ll[0] + b2*ll[1] + eps_hat1;
    o_ur[0] = b1*ur[0] + b2*ur[1] + eps_hat2;
  }
///////////////////////////////////////

/* Saturation function */
inline double saturate(double x, double lb, double ub){
  if (x < lb){
    return lb;
  } else if(x > ub){
    return ub;
  }
  else{
    return x;
  }
}


void print_support(const Cudd& mgr, const BDD& x){
  std::vector< unsigned int >  indices = mgr.SupportIndices({x});
  for (size_t i = 0; i < indices.size(); i++){
    std::cout << indices[i] << " ";
  }
  std::cout << std::endl;
}

////////////////////////////////////////////////////////////////////////
// Traverse the first BDD and write DOT statements to dotFile1
void traverseBDD1(DdManager* manager, DdNode* node, std::ofstream& dotFile1) {
    if (Cudd_IsConstant(node))
        return;

    // Write the current node to the DOT file
    // dotFile1 << "node_" << node << " [label=\"" << node << "\"];\n";

    DdNode* low = Cudd_E(node);
    DdNode* high = Cudd_T(node);

    if (!Cudd_IsConstant(low)) {
        // Write the edge from the current node to the low child
        dotFile1 << node << ",->," << low << ",0\n";
        // Recursively traverse the low child
        traverseBDD1(manager, low, dotFile1);
    }

    if (!Cudd_IsConstant(high)) {
        // Write the edge from the current node to the high child
        dotFile1 << node << ",->," << high << ",1\n";
        // Recursively traverse the high child
        traverseBDD1(manager, high, dotFile1);
    }
};

//////////////////////////////////////////////////////////////////////

int main() {
  /* to measure time */
  TicToc tt;
  /* cudd manager */
  Cudd mgr;
  mgr.AutodynEnable(CUDD_REORDER_SIFT_CONVERGE);
  
  double K = .75;
  /* Dynamics for individual subsystem */ 
  auto dynamics = [K](state_type& x, const control_type u, const exog_type w) -> void {
    x[0] = std::min(K*(x[0] + u[0]), w[0] + 1);
    x[0] = saturate(x[0], 0, 32.0);
  };

  /* Takes an input box and computes an overapproximating box. Both are represented with their lower left
  and upper right corners.
  */
  auto sys_overapprox = [dynamics](const input_type i_ll, const input_type i_ur, state_type& o_ll, state_type& o_ur){
    // data-driven approach
    GRBEnv env = GRBEnv(true);
    //remove all logs
    env.set(GRB_IntParam_OutputFlag, 0);
    env.set(GRB_IntParam_LogToConsole,0);
    //start environment
    env.start();
    
    // Create an empty model
    GRBModel model = GRBModel(env);
    // Create variables
    /* 
    theta_1 = [m11];       
    */
    // diagonal values (can be anything) theta_1
    GRBVar m11 = model.addVar(-200, 200, 0, GRB_CONTINUOUS, "m11");
    // boundary values for diagonals (for optimisation)
    GRBVar u11 = model.addVar(0, 200, 0, GRB_CONTINUOUS, "u11");

    // theta_2
    GRBVar o1 = model.addVar(0, 200, 0, GRB_CONTINUOUS, "o1");
    
    // Set objective: minimize c^T * theta(concatenated theta_1 and theta_2 by columns)
    model.setObjective(u11 + o1, GRB_MINIMIZE);
    
    // Transitions modelled as from pre to post
    // rather than from first to next
          
    //setup other values
    state_type init_pre;
    state_type init_post;
    state_type other_pre;
    state_type other_post;
    state_type lhs;
    state_type rhs;
    double r1, r2;
    // Add constraints by for loop
    for (int i = 0; i < samples; i++)
    {        
      r1 = ((double) rand()) / (double) RAND_MAX;
      r2 = ((double) rand()) / (double) RAND_MAX;
      r1 = (2*r1+1)* h_eta2;
      //random value +/- h of r for each dimension
      other_pre[0] = i_ur[0] + r1;

      r2 = (2*r2+1) * h_eta2;
      //random value +/- h of r for each dimension
      init_pre[0] = i_ll[0] + r2;

      init_post = init_pre;
      other_post = other_pre;
      //updates value _post with actual post value
      dynamics(other_post, {i_ur[1]}, {i_ur[2]});
      dynamics(init_post, {i_ll[1]}, {i_ll[2]});

      for (int j = 0; j < state_dim; j++)
        {
            rhs[j] = abs(other_post[j] - init_post[j]);
            lhs[j] = abs(other_pre[j] - init_pre[j]);
        }
                
      model.addConstr(-(lhs[0]*m11 + o1) + bias <= -rhs[0]);
        
    }
    //absolute values for optimising diagonals of matrix
    model.addConstr( m11 - u11 <= 0);
    model.addConstr(-m11 - u11 <= 0);
                
    // Optimize model
    model.optimize();
    //return r values based on optimised growth bound
    double h_eta1 = h_eta/2;
    state_type post_ll;
    state_type post_ur;
    double kap1 = m11.get(GRB_DoubleAttr_X)*h_eta1 + o1.get(GRB_DoubleAttr_X);
    post_ll[0] = std::max(0.0, init_post[0] + kap1);
    post_ur[0] = std::max(0.0, other_post[0] + kap1);
    o_ll = post_ll;
    o_ur = post_ur;
  };

  /* State spaces */
  std::vector<scots::SymbolicSets> ss_pre; ss_pre.resize(N);
  scots::SymbolicSets pre_product = scots::SymbolicSets();
  std::vector<scots::SymbolicSets> ss_post; ss_post.resize(N);
  scots::SymbolicSets post_product = scots::SymbolicSets();
  state_type s_lb={{0}};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={{32}};
  /* grid node distance diameter */
  state_type s_eta={{1.0}};
  for (int i = 0; i < N; i++){
    ss_pre[i] = scots::SymbolicSets(mgr, state_dim,s_lb,s_ub,s_eta);
    ss_post[i] = scots::SymbolicSets(mgr, state_dim,s_lb,s_ub,s_eta);
    // cartesian product of X_i's
    pre_product = scots::SymbolicSets(pre_product, ss_pre[i]);
    post_product = scots::SymbolicSets(post_product, ss_post[i]);
  }
  std::cout << "Pre State Product Information" << std::endl;
  pre_product.print_info(1);
  
  /*Input spaces*/
  std::vector<scots::SymbolicSets> ss_control; ss_control.resize(N);
  scots::SymbolicSets control_product = scots::SymbolicSets();
  /* lower bounds of the hyper rectangle */
  control_type i_lb={{0}};
  /* upper bounds of the hyper rectangle */
  control_type i_ub={{7}};
  /* grid node distance diameter */
  control_type i_eta={{1.0}};
  for(int i = 0; i < N; i++){
    ss_control[i] = scots::SymbolicSets(mgr, control_dim,i_lb,i_ub,i_eta);
    // cartesian product of U_i's
    control_product = scots::SymbolicSets(control_product, ss_control[i]);
  }
  std::cout << "Control Product Information" << std::endl;
  control_product.print_info(1);

  /* Exogenous spaces */
  scots::SymbolicSets ss_exog;
  exog_type e_lb = {{0}};
  exog_type e_ub = {{32}};
  exog_type e_eta = {{exo_eta}};
  ss_exog = scots::SymbolicSets(mgr, exog_dim,e_lb,e_ub,e_eta);
  
  tt.tic();
  /* 
  Intermediate Variables representing maxima of small sets of variables
  */
  std::cout << "\nIntermediate Variables" << std::endl;
  std::vector<scots::SymbolicSets> ss_intermed; // intermediate layers ll tree
  ss_intermed.resize(N-2);

  std::vector<scots::FunctionDependency>inter_deps;
  inter_deps.resize(N-1); 

  std::cout << "\n\nInterconnection Abstraction" << std::endl;
  BDD interconnection = mgr.bddOne();
  tt.tic();
  for (int i = 0; i < N-1; i++){
    std::cout << "Layer " << i << std::endl;
    std::array<double, 1> inter_lb, inter_ub, inter_eta;
    inter_lb[0] = 0;
    inter_ub[0] = 32;
    inter_eta[0] = 1.0;
    if (i < N-2){
      ss_intermed[i] = scots::SymbolicSets(mgr, 1, inter_lb, inter_ub, inter_eta);
    }

    if (i == 0){ // takes two states, outputs to intermed
      inter_deps[i] = scots::FunctionDependency({ss_pre[i], ss_pre[i+1]}, {ss_intermed[i]});
      inter_deps[i].set_dependency(ss_intermed[i][0], {ss_pre[i][0], ss_pre[i+1][0]});
    }
    else if (i == N - 2){ // takes state and intermed, output to exog
      inter_deps[i] = scots::FunctionDependency({ss_intermed[i-1], ss_pre[i+1]}, { ss_exog });
      inter_deps[i].set_dependency(ss_exog[0], {ss_intermed[i-1][0], ss_pre[i+1][0]});
    }
    else{ // takes state and intermed, output to intermed
      inter_deps[i] = scots::FunctionDependency({ss_intermed[i-1], ss_pre[i+1]}, {ss_intermed[i]});
      inter_deps[i].set_dependency(ss_intermed[i][0], {ss_intermed[i-1][0], ss_pre[i+1][0]});
    }

    scots::FunctionAbstracter<std::array<double, 2>, std::array<double, 1> > layer(inter_deps[i], max3);
    interconnection &= layer.compute_abstraction(mgr);

  }
  tt.toc();

  /* Declare dependencies for individual systems */
  std::vector<scots::FunctionDependency> sysdeps(N, scots::FunctionDependency());
  for (int i = 0; i < N; i++){
    sysdeps[i] = scots::FunctionDependency({ss_pre[i], ss_control[i], ss_exog},{ss_post[i]});
    sysdeps[i].set_dependency(ss_post[i][0], {ss_pre[i][0], ss_control[i][0], ss_exog[0]});
  }

  /*Compute system abstractions using dependencies*/
  std::vector<scots::FunctionAbstracter<input_type, state_type> > abs_comp(N, scots::FunctionAbstracter<input_type, state_type>());
  std::vector<BDD> abs_systems(N, mgr.bddOne());
  BDD composed_systems = mgr.bddOne();
  std::cout << "\n\nSystem Abstractions" << std::endl;
  tt.tic();
  for (int i = 0; i < N; i++){
    abs_comp[i] = scots::FunctionAbstracter<input_type, state_type>(sysdeps[i], sys_overapprox);
    std::cout << "System " << i << " abstraction \n";
    composed_systems &= abs_comp[i].compute_abstraction(mgr);
  }
  tt.toc();

  std::cout << "\n\nComposing Systems" << std::endl;
  tt.tic();
  BDD monolithic = composed_systems & interconnection;
  tt.toc();

  return 0;
} // close main








