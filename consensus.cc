/*
 * consensus.cc
 * 
 * Abstract 10 tank systems and apply an interconnection  
 *  
 *
 * Synthesize for a consensus objective. 
 * 
 *  created: Oct 24
 *  author: Daniel A.
 * 
 * ### Set membership (FRR) V is intrinsic as the norm ||x-x'|| ###
 * 
 */

#include <iostream>
#include <algorithm>
#include <array>
#include <cmath>
#include "SymbolicSets.hh"
#include "gurobi_c++.h"

/* SCOTS header */
#include "scots.hh"

/* time profiling */
#include "TicToc.hh"
/* memory profiling */
#include <sys/time.h>
#include <sys/resource.h>
struct rusage usage;

// for subsystems
double samples = 150;
double Lx = 0.9083 * 2;
double Lw = 0.9988;

double h_eta = 1.25; // grid parameter for state set
double h_eta2 = h_eta/samples; // sub-grid parameter for each state cell
double exo_eta = 2.5; // grid parameter for internal input
double exo_eta2 = exo_eta/samples; // sub-grid parameter for each internal input cell
double const bias = Lx * h_eta2 + Lw * exo_eta2; //+ Lw * exo_eta2; bias term for SCP for each (x^,u^,w^)

double const etaz = 5;
double const etax = 2.5;
double const etau = 1.5;
double const etaw = 5;

// interconnection map approximation
double param_w1 = 0; // features of M at w level
double param_w2 = 0.079867198176204;
double param_w3a = 0.03993299871286744;
double param_w3b = 0.039932998712867485;
double param_w4 = 0;
double param_w5 = 0.079867198176204;
double param_w6a = 0.019965898981201116;
double param_w6b = 0.019965898981200166;
double param_w6c = 0.019965898981199996;
double param_w6d = 0.01996589898120013;
double param_w7 = 0.079867198176204;
double param_w8 = 0.079867198176204;
double param_w9a = 0.03993299871286744;
double param_w9b = 0.039932998712867485;
double param_w10a = 0.03993299871286744;
double param_w10b = 0.039932998712867485;
const double M[10][10] = {{param_w1,0,0,0,0,0,0,0,0,0}, {param_w2,0,0,0,0,0,0,0,0,0}, {param_w3a,param_w3b,0,0,0,0,0,0,0,0}, {0,0,param_w4,0,0,0,0,0,0,0}, {0,0,0,param_w5,0,0,0,0,0,0}, {0,0,param_w6a,param_w6b,param_w6c,0,param_w6d,0,0,0}, {0,0,0,0,param_w7,0,0,0,0,0}, {0,0,param_w8,0,0,0,0,0,0,0}, {0,0,0,0,0,param_w9a,0,param_w9b,0,0}, {0,0,0,0,0,0,0,param_w10a,param_w10b,0}};

/* state space dim */
const int state_dim=1;
const int state_dim_s=1;
/* input space dim */
const int control_dim=10;
const int control_dim_s=1;
/* exog space dim*/
const int exog_dim = 10;
const int exog_dim_s = 1;
/* input space of system is a cartesian product*/
const int input_dim = state_dim + control_dim + exog_dim; 
const int input_dim_s = state_dim_s + control_dim_s + exog_dim_s; 
/* Create N identical systems */
const int N = 10;

const int inter_dim = 8; 

// for the controller synthesis part
double const Tt = 200;
double const bd = 100;
double const lbd = 0;
double const theta = 20;
double thetabd = 30;

const double tau = 0.5;
/*
 * data types for the state space elements and input space
 * elements used in uniform grid
 */
using state_type_s = std::array<double,state_dim_s>;
using state_type = std::array<double,state_dim>;
using control_type = std::array<double,control_dim>;
using control_type_s = std::array<double,control_dim_s>;
using exog_type  = std::array<double,exog_dim>;
using exog_type_s  = std::array<double,exog_dim_s>;
using input_type  = std::array<double,input_dim>;
using input_type_s  = std::array<double,input_dim_s>;
using inter_type = std::array<double, inter_dim>;

/* Data types for the interconnection relation */
using prod_state_type = std::array<double, state_dim*N>;
using prod_control_type = std::array<double, control_dim>;
using prod_exog_type = std::array<double, exog_dim>;
using prod_intermed_type_lvl2 = std::array<double, inter_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

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

inline double tank(double x, double u, double w){
  double x1 = -tau/2.0 + sqrt(pow(tau,2)/4.0 + x + tau*(u + sqrt(w)));
  double xn = pow(x1, 2);
  return xn; 
}

void print_support(const Cudd& mgr, const BDD& x){
  std::vector< unsigned int >  indices = mgr.SupportIndices({x});
  for (size_t i = 0; i < indices.size(); i++){
    std::cout << indices[i] << " ";
  }
  std::cout << std::endl;
}

int main() {
  /* to measure time */
  TicToc tt;
  /* cudd manager and BDD reordering options */
  Cudd mgr;
  mgr.AutodynEnable(CUDD_REORDER_SIFT_CONVERGE);
  // mgr.AutodynEnable(CUDD_REORDER_RANDOM_PIVOT);
  mgr.EnableReorderingReporting(); 
  
  /* Dynamics for individual subsystem */ 
  auto dynamics = [](state_type& x, const control_type u, const exog_type w) -> void {
    // for (int i = 0; i < 10; i++) {
      x[0] = tank(x[0], u[0], w[0]);
    // }
  };

  /* Takes an input box and computes an overapproximating box. Both are represented with their lower left
  and upper right corners.
  */
  auto sys_overapprox = [dynamics](const input_type_s i_ll, const input_type_s i_ur, state_type_s& o_ll, state_type_s& o_ur){
    // data-driven approach
    try {
    GRBEnv env = GRBEnv(true);
    //remove all logs
    env.set(GRB_IntParam_OutputFlag, 0);
    env.set(GRB_IntParam_LogToConsole, 0);
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
    double h_eta3 = h_eta;//2;
    // Add constraints by for loop
    for (int i = 0; i < samples; i++)
    {        
      r1 = ((double) rand()) / (double) RAND_MAX;
      r2 = ((double) rand()) / (double) RAND_MAX;
      r1 = (2*r1+1)* h_eta3;
      //random value +/- h of r for each dimension
      for (int j = 0; j < state_dim; j++){
        other_pre[j] = i_ur[j] + r1;
      }

      r2 = (2*r2+1) * h_eta3;
      //random value +/- h of r for each dimension
      for (int j = 0; j < state_dim; j++) {
        init_pre[j] = i_ll[j] + r2;
      }

      init_post = init_pre;
      other_post = other_pre;
      
      //updates value _post with actual post value
      // exog_type i_ur_w;
      // exog_type i_ll_w;
      // control_type i_ur_u;
      // control_type i_ll_u;
      
      // for (int j = 0; j < state_dim; j++) {
      //   i_ur_u[j] = i_ur[j+N];
      //   i_ll_u[j] = i_ll[j+N];
      //   i_ur_w[j] = i_ur[j+2*N];
      //   i_ll_w[j] = i_ll[j+2*N];
      // }
      // dynamics(other_post, i_ur_u, i_ur_w);
      // dynamics(init_post, i_ll_u, i_ll_w);  

      dynamics(other_post, {i_ur[1]}, {i_ur[2]});
      dynamics(init_post, {i_ll[1]}, {i_ll[2]});  

      for (int j = 0; j < state_dim; j++)
        {
            rhs[j] = abs(other_post[j] - init_post[j]);
            lhs[j] = abs(other_pre[j] - init_pre[j]);

            // model.addConstr(-(lhs[j]*m11 + o1) + bias <= -rhs[j]);
        }
                
      model.addConstr(-(lhs[0]*m11 + o1) + bias <= -rhs[0]);
    }
    //absolute values for optimising diagonals of matrix
    model.addConstr( m11 - u11 <= 0);
    model.addConstr(-m11 - u11 <= 0);
                
    // Optimize model
    model.optimize();
    //return r values based on optimised growth bound
    double h_eta1 = h_eta;///2;
    state_type post_ll;
    state_type post_ur;
    for (int j = 0; j < state_dim; j++) {
      post_ll[0] = init_post[0] + m11.get(GRB_DoubleAttr_X)*h_eta1 + o1.get(GRB_DoubleAttr_X);
      post_ur[0] = other_post[0] + m11.get(GRB_DoubleAttr_X)*h_eta1 + o1.get(GRB_DoubleAttr_X);
    }
    o_ll = post_ll;
    o_ur = post_ur;
    // for(int i=0; i<3; i++){
    //   cout << i << r[i] << "\n";
    // } 
    
    } catch (GRBException &e) {
        // std::cerr << "Error code = " << e.getErrorCode() << std::endl;
        // std::cerr << e.getMessage() << std::endl;
    } catch(...) {
            // std::cout << "Exception during optimization" << std::endl;
        }
  };

  /* State spaces */
  std::vector<scots::SymbolicSets> ss_pre; ss_pre.resize(N);
  scots::SymbolicSets pre_product = scots::SymbolicSets();
  std::vector<scots::SymbolicSets> ss_post; ss_post.resize(N);
  scots::SymbolicSets post_product = scots::SymbolicSets();
  state_type_s s_lb={{lbd}};
  /* upper bounds of the hyper rectangle */
  state_type_s s_ub={{bd}};
  /* grid node distance diameter */
  state_type s_eta={{etax}};
  for (int i = 0; i < N; i++){
    ss_pre[i] = scots::SymbolicSets(mgr, state_dim_s,s_lb,s_ub,s_eta);
    ss_post[i] = scots::SymbolicSets(mgr, state_dim_s,s_lb,s_ub,s_eta);

    pre_product = scots::SymbolicSets(pre_product, ss_pre[i]);
    post_product = scots::SymbolicSets(post_product, ss_post[i]);
  }
  std::cout << "Pre State Product Information" << std::endl;
  pre_product.print_info(1);
  
  /*Input spaces*/
  std::vector<scots::SymbolicSets> ss_control; ss_control.resize(N);
  scots::SymbolicSets control_product = scots::SymbolicSets();
  /* lower bounds of the hyper rectangle */
  control_type_s i_lb={{0}};
  /* upper bounds of the hyper rectangle */
  control_type_s i_ub={{ 10}};
  /* grid node distance diameter */
  control_type_s i_eta={{etau}};
  for(int i = 0; i < N; i++){
    if (i == 0 || i == 3) {
      ss_control[i] = scots::SymbolicSets(mgr, control_dim_s,i_lb,i_ub,i_eta);
      control_product = scots::SymbolicSets(control_product, ss_control[i]);
    }
    else{
      ss_control[i] = scots::SymbolicSets(mgr, control_dim_s,i_lb,i_lb,i_eta);
      control_product = scots::SymbolicSets(control_product, ss_control[i]);
    }
  }
  std::cout << "Control Product Information" << std::endl;
  control_product.print_info(1);

  // /* Exogenous spaces */
  // scots::SymbolicSets ss_exog;
  // exog_type e_lb = {{0}};
  // exog_type e_ub = {{bd}};
  // exog_type e_eta = {{etaw}};
  // ss_exog = scots::SymbolicSets(mgr, exog_dim,e_lb,e_ub,e_eta);

  /* Exogenous spaces */
  std::vector<scots::SymbolicSets> ss_exog; ss_exog.resize(N);
  scots::SymbolicSets se_product = scots::SymbolicSets();
  exog_type_s e_lb = {{lbd}};
  exog_type_s e_ub = {{bd}};
  exog_type_s e_eta = {{etaw}};
  for (int i = 0; i < N; i++){
    ss_exog[i] = scots::SymbolicSets(mgr, exog_dim_s,e_lb,e_ub,e_eta);
    se_product = scots::SymbolicSets(se_product, ss_exog[i]);
  }

  /* Declare dependencies for individual systems */
  std::vector<scots::FunctionDependency> sysdeps(N, scots::FunctionDependency());
  for (int i = 0; i < N; i++){
    sysdeps[i] = scots::FunctionDependency({ss_pre[i], ss_control[i], ss_exog[i]},{ss_post[i]});
    sysdeps[i].set_dependency(ss_post[i][0], {ss_pre[i][0], ss_control[i][0], ss_exog[i][0]});
  }

  /*Compute system abstractions using dependencies*/
  std::vector<scots::FunctionAbstracter<input_type_s, state_type_s> > abs_comp(N, scots::FunctionAbstracter<input_type_s, state_type_s>());
  std::vector<BDD> abs_systems(N, mgr.bddOne());
  BDD systems = mgr.bddOne(); 
  for (int i = 0; i < N; i++){
    abs_comp[i] = scots::FunctionAbstracter<input_type_s, state_type_s>(sysdeps[i], sys_overapprox);
    tt.tic();
    std::cout << "System " << i << " abstraction ";
    abs_systems[i] = abs_comp[i].compute_abstraction(mgr);
    tt.toc();
  }


  /* 
  Intermediate Variables representing sums of small sets of variables
  */
  std::cout << "\nIntermediate Variables" << std::endl;
  inter_type inter_lb = {{lbd,lbd,lbd,lbd,lbd,lbd,lbd,lbd}};
  inter_type inter_ub = {{bd,bd,bd,bd,bd,bd,bd,bd}};
  inter_type inter_eta = {{etaz, etaz, etaz, etaz, etaz, etaz, etaz, etaz}};
  scots::SymbolicSets ss_inter = scots::SymbolicSets(mgr, inter_dim, inter_lb,inter_ub,inter_eta);
  ss_inter.print_info(1);
  tt.tic();

  /*Interconnection Level 1. */
  std::cout << "Level 1" << std::endl;
  scots::FunctionDependency intermed_dep_1({pre_product}, {ss_inter});
  intermed_dep_1.set_dependency(ss_inter[0], {ss_pre[0][0]});
  intermed_dep_1.set_dependency(ss_inter[1], {ss_pre[0][0], ss_pre[1][0]});
  intermed_dep_1.set_dependency(ss_inter[2], {ss_pre[3][0]});
  intermed_dep_1.set_dependency(ss_inter[3], {ss_pre[2][0], ss_pre[3][0], ss_pre[4][0], ss_pre[6][0]});
  intermed_dep_1.set_dependency(ss_inter[4], {ss_pre[4][0]});
  intermed_dep_1.set_dependency(ss_inter[5], {ss_pre[2][0]});
  intermed_dep_1.set_dependency(ss_inter[6], {ss_pre[5][0], ss_pre[7][0]});
  intermed_dep_1.set_dependency(ss_inter[7], {ss_pre[7][0], ss_pre[8][0]});
  std::cout << intermed_dep_1 << std::endl;
  auto inter_overapprox_1 = [ss_pre](const prod_state_type ll, const prod_state_type ur, inter_type &o_ll, inter_type &o_ur){
    std::vector<double> eta = {etaz}; //ss_pre[0].get_eta();
    prod_state_type center; 
    // Compute center of box from ll and ur corners 
    for (int i = 0; i < N; i++){
      center[i]   = (ll[i] + ur[i])/2.0;
    }
    o_ll[0] = std::max(M[1][0]*center[0] - 2.0*eta[0], lbd);
    o_ur[0] = std::min(M[1][0]*center[0] + 2.0*eta[0], bd);
    o_ll[1] = std::max(M[2][0]*center[0] + M[2][1]*center[1] - 2.0*eta[0], lbd);
    o_ur[1] = std::min(M[2][0]*center[0] + M[2][1]*center[1] + 2.0*eta[0], bd);
    o_ll[2] = std::max(M[4][3]*center[3] - 2.0*eta[0], lbd); 
    o_ur[2] = std::min(M[4][3]*center[3] + 2.0*eta[0], bd);
    o_ll[3] = std::max(M[5][2]*center[2] + M[5][3]*center[3] + M[5][4]*center[4] + M[5][6]*center[6] - 2.0*eta[0], lbd); 
    o_ur[3] = std::min(M[5][2]*center[2] + M[5][3]*center[3] + M[5][4]*center[4] + M[5][6]*center[6] + 2.0*eta[0], bd);
    o_ll[4] = std::max(M[6][4]*center[4] - 2.0*eta[0], lbd); 
    o_ur[4] = std::min(M[6][4]*center[4] + 2.0*eta[0], bd);
    o_ll[5] = std::max(M[7][2]*center[2] - 2.0*eta[0], lbd); 
    o_ur[5] = std::min(M[7][2]*center[2] + 2.0*eta[0], bd);
    o_ll[6] = std::max(M[8][5]*center[5] + M[8][7]*center[7] - 2.0*eta[0], lbd);
    o_ur[6] = std::min(M[8][5]*center[5] + M[8][7]*center[7] + 2.0*eta[0], bd);
    o_ll[7] = std::max(M[9][7]*center[7] + M[9][8]*center[8] - 2.0*eta[0], lbd);
    o_ur[7] = std::min(M[9][5]*center[7] + M[9][7]*center[8] + 2.0*eta[0], bd);
  };
  scots::FunctionAbstracter<prod_state_type, inter_type> inter_abs_1(intermed_dep_1 , inter_overapprox_1);
  BDD abs_inter_1 = inter_abs_1.compute_abstraction(mgr);

  /*Interconnection Level 2.*/
  std::cout << "\nLevel 2" << std::endl;
  scots::FunctionDependency intermed_dep_2({ss_inter},{ss_exog});
  intermed_dep_2.set_dependency(ss_exog[1][0], {ss_inter[0]});
  intermed_dep_2.set_dependency(ss_exog[2][0], {ss_inter[1]});
  intermed_dep_2.set_dependency(ss_exog[4][0], {ss_inter[2]});
  intermed_dep_2.set_dependency(ss_exog[5][0], {ss_inter[3]});
  intermed_dep_2.set_dependency(ss_exog[6][0], {ss_inter[4]});
  intermed_dep_2.set_dependency(ss_exog[7][0], {ss_inter[5]});
  intermed_dep_2.set_dependency(ss_exog[8][0], {ss_inter[6]});
  intermed_dep_2.set_dependency(ss_exog[9][0], {ss_inter[7]});
  std::cout << intermed_dep_2 << std::endl;
  auto inter_overapprox_2 = [pre_product, ss_inter](const prod_intermed_type_lvl2 ll, const prod_intermed_type_lvl2 ur, exog_type &o_ll, exog_type &o_ur){
    std::vector<double> eta = {etaz}; //pre_product.get_eta();
    std::vector<double> mu = ss_inter.get_eta();
    exog_type center;
    for (int i = 0; i < 4; i++){
      center[i]  = (ll[i] + ur[i])/2.0;
    }
    o_ll[0] = 0;
    o_ur[0] = 0;
    o_ll[1] = std::max(center[0] - .5*eta[0], lbd);
    o_ur[1] = std::min(center[0] + .5*eta[0], bd); 
    o_ll[2] = std::max(center[1] - .5*eta[0], lbd);
    o_ur[2] = std::min(center[1] + .5*eta[0], bd); 
    o_ll[3] = 0;
    o_ur[3] = 0;
    o_ll[4] = std::max(center[2] - .5*eta[0], lbd);
    o_ur[4] = std::min(center[2] + .5*eta[0], bd); 
    o_ll[5] = std::max(center[3] - .5*eta[0], lbd);
    o_ur[5] = std::min(center[3] + .5*eta[0], bd);
    o_ll[6] = std::max(center[4] - .5*eta[0], lbd);
    o_ur[6] = std::min(center[4] + .5*eta[0], bd);
    o_ll[7] = std::max(center[5] - .5*eta[0], lbd);
    o_ur[7] = std::min(center[5] + .5*eta[0], bd);
    o_ll[8] = std::max(center[6] - .5*eta[0], lbd);
    o_ur[8] = std::min(center[6] + .5*eta[0], bd); 
    o_ll[9] = std::max(center[7] - .5*eta[0], lbd);
    o_ur[9] = std::min(center[7] + .5*eta[0], bd);
  };
  scots::FunctionAbstracter<prod_intermed_type_lvl2, exog_type> inter_abs_2(intermed_dep_2 , inter_overapprox_2);
  BDD abs_inter_2 = inter_abs_2.compute_abstraction(mgr); 


  /* 
  Declare and abstract interconnection. 
  */

  BDD abs_inter = abs_inter_1 & abs_inter_2;
  print_support(mgr, abs_inter);
  std::cout << (int)(abs_inter == mgr.bddZero()) << std::endl;
  tt.toc();

  /** Construct abstraction of interconnected system.
  First construct system that's a subset of X x W x U x X'
  Then use this system to get a monolithic system X x U x X'
  with nondeterminism from W, which is constrained by values in X
   **/
  std::cout << "Composing smaller systems " << std::endl;
  BDD interconnected_sys = mgr.bddOne();
  tt.tic();
  for (int i = 0; i < N; i++){
    std::cout << i << std::endl;
    interconnected_sys = interconnected_sys & abs_systems[i];
  }
  tt.toc();

  /** Construct Invariant Set on monlithic space **/
  std::cout << "Constructing Target Set" << std::endl;
  BDD target = mgr.bddZero();
  for (int t = thetabd; t < thetabd+1 ; t++){
    /* Predicate functions */
    auto inv_predicate = [t, &ss_pre](const abs_type& idx){
      state_type x;
      ss_pre[0].itox(idx,x);
      /* function returns true if cell associated with x is in invariant set  */
      if ( x[0] <= t + theta && x[0] >= t - theta)
        return true;
      return false;
    };

    /* All systems must have x[0] state close to t */
    BDD allin = mgr.bddOne();
    for (int i = 0; i < N; i++){
      allin &= ss_pre[i].ap_to_bdd(mgr,inv_predicate);
    }
    /* Union over all t's*/
    target |= allin;
  }

  /* Controller synthesis over monolithic system */
  scots::SymbolicSets aux_product = scots::SymbolicSets(se_product, ss_inter); // exogenous and intermediate variables
  scots::InterconnectedEnfPre enf_pre(mgr,interconnected_sys, pre_product, control_product, post_product, aux_product, abs_inter, abs_systems);
  scots::SymbolicSets controller(pre_product,control_product);
  // /* the controller */
  BDD X , XX, C = mgr.bddZero(), newbasin;
  // /* BDD cube for existential abstract inputs */
  const BDD U = control_product.get_cube(mgr);
  const BDD E = aux_product.get_cube(mgr); // all auxiliary variables, exog + intermediate ones
  std::cout << "\nState Space Size: " << pre_product.get_size(mgr,mgr.bddOne()) << std::endl;
  std::cout << "Target size: " << pre_product.get_size(mgr,target) << std::endl << std::endl;
  tt.tic();
  /*Safety objective*/
  std::cout<< "Invariance Controller Synthesis" << std::endl;
  X = mgr.bddZero(); XX = mgr.bddOne();
  for(int i=1; XX != X; i++) { 
    X = XX;
    C = enf_pre(X); // (state, input) controlled pre pairs
    XX = C.ExistAbstract(U*E);
    XX = (X & XX & target) & abs_inter.ExistAbstract(E);
    std::cout << i << "-th winning domain size: " << pre_product.get_size(mgr,XX) << std::endl;
  }
  
  if(write_to_file(mgr, controller, C.ExistAbstract(E),"consensus_inv_controller"))
    std::cout << "Done writing controller to file. \n";

  BDD inv = XX;
  std::cout << "Invariant Set BDD Var ID Support" <<std::endl;
  print_support(mgr,inv.ExistAbstract(U*E));
  std::cout << "Target Set BDD Var ID Support" <<std::endl;
  print_support(mgr,target);
  std::cout << "Controller BDD Var ID Support" << std::endl;
  print_support(mgr,C);
  std::cout << "Inter BDD Var ID Support" << std::endl;
  print_support(mgr, abs_inter); 

  /*Reach objective*/
  std::cout<< "Reachability Controller Synthesis" << std::endl;
  X = mgr.bddOne(); XX = mgr.bddZero();
  for(int i = 1; XX != X; i++){
    std::cout << i << "-th reach basin size: " << pre_product.get_size(mgr,XX) << std::endl;
    X = XX;
    XX = (enf_pre(X) | inv) & abs_inter.ExistAbstract(E);
    std::cout << i << "-th predecessor " << pre_product.get_size(mgr,XX) << std::endl;
    newbasin = pre_product.get_grid_bdd(mgr) & XX & (!(C.ExistAbstract(U*E)));
    std::cout << i << "-th new states " << pre_product.get_size(mgr,newbasin) << std::endl << std::endl;
    XX = XX.ExistAbstract(U*E);
    C = C | newbasin;
  }
  tt.toc();

  /* Print final reach set */
  if (false){
    std::ofstream file;
    file.open("better_consensus_reachable.txt");
    auto a = pre_product.bdd_to_grid_points(mgr, XX);
    for(size_t j = 0; j < a.size(); j++){
      file << a[j] << " ";
      if (j % (state_dim *N) == (state_dim *N)-1)
        file << "\n";
    }
    file.close();
  }

  std::cout << "\nWrite controller to better_consensus_controller.scs \n";
  if(write_to_file(mgr, controller, C.ExistAbstract(E),"better_consensus_controller"))
    std::cout << "Done. \n";

  /* Dynamics for monolithic system */
  auto prod_dynamics = [](prod_state_type &x,  prod_control_type u) {
    x[0] = pow(-tau/2.0 + sqrt(pow(tau,2)/4.0 + x[0] + tau*(u[0])), 2);
    x[1] = pow(-tau/2.0 + sqrt(pow(tau,2)/4.0 + x[1] + tau*(sqrt(x[0]))), 2);
    x[2] = pow(-tau/2.0 + sqrt(pow(tau,2)/4.0 + x[2] + tau*(sqrt(x[0])+sqrt(x[1]))/2.0), 2);
    x[3] = pow(-tau/2.0 + sqrt(pow(tau,2)/4.0 + x[3] + tau*(u[3])), 2); 
    x[4] = pow(-tau/2.0 + sqrt(pow(tau,2)/4.0 + x[4] + tau*(sqrt(x[3]))), 2);
    x[5] = pow(-tau/2.0 + sqrt(pow(tau,2)/4.0 + x[5] + tau*(sqrt(x[2])+sqrt(x[3])+sqrt(x[4])+sqrt(x[7]))/4.0), 2);
    x[6] = pow(-tau/2.0 + sqrt(pow(tau,2)/4.0 + x[6] + tau*(sqrt(x[4]))), 2);
    x[7] = pow(-tau/2.0 + sqrt(pow(tau,2)/4.0 + x[7] + tau*(sqrt(x[2]))), 2);
    x[8] = pow(-tau/2.0 + sqrt(pow(tau,2)/4.0 + x[8] + tau*(sqrt(x[7])+sqrt(x[5]))/2.0), 2);
    x[9] = pow(-tau/2.0 + sqrt(pow(tau,2)/4.0 + x[9] + tau*(sqrt(x[7])+sqrt(x[8]))/2.0), 2);
  };

  /*
  Simulate Dynamics with synthesized controller 

  Two options:
  - active_control == true: Use synthesized controller
  - active_control == false: All inputs u are equal to zero 
  */
  std::cout << "Simulating Active Controller" << std::endl;
  prod_state_type x={5, 6, 65, 7, 7.5, 70, 8, 65.5, 75, 80};
  int u_index;

  std::ofstream file;
  file.open("traj_active1.txt");

  for(int i=0; i<Tt; i++) {
    file << "State: ";
    for(int j = 0; j < N; j++){
      file << x[j] << " ";
    }
    file<< std::endl;

    file << "Getting Control Input" << std::endl;
    auto u = controller.restriction<prod_state_type>(mgr,C,x);
    if (u.size() == 0){
      file << "No valid control. Exiting" << std::endl;
      break;
    }
    u_index = u.size() - (rand() % (u.size()/N))*N;//rand() % (u.size()/control_dim);
    file << "Number of Permitted Actions: " << u.size() << std::endl;
    file << "Input Index: " << u_index << std::endl;
    file << "Input: ";
    for(int j = 0; j < N; j++){
      file << u[u_index+j] << " ";
    } 
    
    prod_dynamics(x,{u[u_index],u[u_index+1],u[u_index+2],u[u_index+3],u[u_index+4],u[u_index+5],u[u_index+6],u[u_index+7],u[u_index+8],u[u_index+9]});
    file<< std::endl << std::endl;

  } // close simulation for loop
  file.close();
  std::cout << "Trajectory written to traj_active1.txt" << std::endl;

  /* Passive Control */
  std::cout << "Simulating Passive Controller" << std::endl;
  x={5, 6, 65, 7, 7.5, 70, 8, 65.5, 75, 80};
  
  file.open("traj_passive1.txt");
  for(int i=0; i<Tt; i++) {
    file << "State: ";
    for(int j = 0; j < N; j++){
      file << x[j] << " ";
    }
    file<< std::endl;
    prod_dynamics(x, {0,0,0,0,0,0,0,0,0,0});
  }
  file.close();
  std::cout << "Trajectory written to traj_passive1.txt" << std::endl;

} // close main

