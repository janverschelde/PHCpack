/** @file */
#include <iostream>
#include <sstream>

#include "path_host.h"
#include "path_gpu.h"
#include "parameter.h"
#include "err_check.h"

#include "init_test.h"
#include "mgs_test.h"
#include "predict_test.h"
#include "eval_test.h"
#include "newton_test.h"
#include "path_test.h"
#include "path_multi_test.h"
#include "pieri_test.h"
#include "witness_set_test.h"
#include "ada_test.h"

#include <ctime>

const std::string DATA = "/home/jan/Problems/GPUdata";

void all_test
  (PolySys& Target_Sys, Workspace& workspace_cpu,
   CPUInstHom& cpu_inst_hom, Parameter path_parameter, 
   CT* sol0, int n_eq, int dim);

void sys_filename(string& start_sys_filename, string& target_sys_filename,\
                  int sys_type, int dim)
{
   std::ostringstream start_filename;
   std::ostringstream target_filename;

   if(sys_type == 0)
   {
      start_filename << DATA << "/cyclic/cyc" << dim << "p1";
      target_filename << DATA << "/cyclic/cyc" << dim << "q1";
   }
   else if(sys_type == 1)
   {
      start_filename << DATA << "/MultiPath/cyclic" << dim << ".start";
      target_filename << DATA << "/MultiPath/cyclic"<< dim << ".target";
   }
   else if(sys_type == 2)
   {
      start_filename << DATA << "/PieriBig1/pieri353start" << dim-32;
      target_filename << DATA << "/PieriBig1/pieri353target" << dim-32;
   }
   else if(sys_type == 3)
   {
      start_filename << DATA << "/MultiPath/game8two.start";
      target_filename << DATA << "/MultiPath/game8two.target";
   }
   else if(sys_type == 4)
   {
      start_filename << DATA << "/MultiPath/pieri44.start";
      target_filename << DATA << "/MultiPath/pieri44.target";
   }
   else if(sys_type == 5)
   {
      start_filename << DATA << "/MultiPath/eq.start";
      target_filename << DATA << "/MultiPath/eq.target";
   }
   else
   {
      start_filename << DATA << "/PieriBig2/pieri364start" << dim-32;
      target_filename << DATA << "/PieriBig2/pieri364target" << dim-32;
   }
   start_sys_filename = start_filename.str();
   target_sys_filename = target_filename.str();

   std::cout << "start system file name : "
             << start_sys_filename << std::endl;
   std::cout << "target system file name : "
             << target_sys_filename << std::endl;
}

int main_test(int test, int dim, int device_option, int n_path)
{
   // Tests
   // 0: Predict Test
   // 1: Evalute Test
   // 2: Modified Gram-Smith Test
   // 3: Modified Gram-Smith Large Test
   // 4: Modified Gram-Smith Large Test for any dim
   // 5: Newton Test
   // 6: Path Tracking Test
   // 7: Path Tracking Reverse Test
   // 8: Witeness Set Test
   // 9: Pieri Test
   // 10: Multiple path Test

   // Initialize parameters
   Parameter path_parameter(N_PREDICTOR, MAX_STEP, MAX_IT, MAX_DELTA_T, \
      MAX_DELTA_T_END, MIN_DELTA_T, ERR_MAX_RES, ERR_MAX_DELTA_X, \
      ERR_MAX_FIRST_DELTA_X, ERR_MIN_ROUND_OFF, MAX_IT_REFINE, \
      ERR_MIN_ROUND_OFF_REFINE, STEP_INCREASE, STEP_DECREASE);

   int sys_type = 1; // Chooose system type

   int n_eq;
   CT* sol0;
   PolySys Target_Sys;
   PolySys Start_Sys;
   PolySolSet start_sol_set;
   CPUInstHom cpu_inst_hom;
   Workspace workspace_cpu;

   if(test == 9) // Pieri test
   {
      int dim_start = 32;
      int dim_end = 36;
      int sys_type = 3;
      Pieri_Test(dim_start, dim_end, path_parameter, sys_type, device_option);
      return 0;
   }
   else if(test == 4) // Initialize workspace only for MGS
   {
      int n_eq = dim;
      cpu_inst_hom.dim = dim;
      cpu_inst_hom.n_eq = n_eq;
      int workspace_size = 5;
      int n_coef = 1;
      int n_constant = 1;
      int n_predictor = 1;
      workspace_cpu.init(workspace_size, n_coef, n_constant,
                         n_eq, dim, n_predictor);
      sol0 = new CT[dim];
   }
   else
   {
      string start_sys_filename;
      string target_sys_filename;
      sys_filename(start_sys_filename, target_sys_filename, sys_type, dim);

		bool ada_test = false;
		if(ada_test==true){
			char start_file[] = "../../Problems/MultiPath/cyclic10.start";
			char target_file[] = "../../Problems/MultiPath/cyclic10.target";
			ada_read_homotopy(start_file, target_file, Start_Sys, Target_Sys, start_sol_set);
			bool init_success = init_homotopy_test_ada(Target_Sys, Start_Sys, \
						   cpu_inst_hom, workspace_cpu, path_parameter.n_predictor);
			if(init_success == false){
				std::cout << "Test Initialization Failed." << std::endl;
				return false;
			}
			sol0 = start_sol_set.get_sol(0);
			std::cout << "----------- ada read Finished -----------" << std::endl;
		}
		else{
			// Initialize homotopy and workspace
			bool init_success = init_homotopy_test(Target_Sys, Start_Sys, start_sol_set, dim, n_eq, sol0, cpu_inst_hom, \
					  workspace_cpu, test, path_parameter.n_predictor, start_sys_filename, target_sys_filename);
			if(init_success == false){
				std::cout << "Test Initialization Failed." << std::endl;
				return false;
			}
		}
	}

	CT t;
	CT* sol_new = NULL;
	switch(test) {
	case 0:
		// Predict Test
		t = CT(0.5,0);
		predict_test(workspace_cpu, cpu_inst_hom, t);
		break;
	case 1:
		// Evaluate Test
		std::cout << "Evaluation Test" << std::endl;
		t = CT(T1(1),T1(0));
		eval_test_classic(workspace_cpu, cpu_inst_hom, sol0, t, Target_Sys, n_path);
		break;
	case 2:
		// Modified Gram-Smith Test
		t = CT(0.1,0);
		workspace_cpu.update_x_t_value(sol0, t);
		mgs_test(workspace_cpu, cpu_inst_hom);
		break;
	case 3:
		// Modified Gram-Smith Test
		t = CT(0.1,0);
		workspace_cpu.update_x_t_value(sol0, t);
		mgs_test_large(workspace_cpu, cpu_inst_hom);
		break;
	case 4:
		// Modified Gram-Smith Large Test for any dim
		t = CT(0.1,0);
		workspace_cpu.update_x_t_value(sol0, t);
		mgs_test_any(workspace_cpu, cpu_inst_hom, device_option, n_path);
		break;
	case 5:
		// Newton Test
		t = CT(0.0001,0);
		workspace_cpu.update_x_t_value(sol0, t);
		newton_test(workspace_cpu, cpu_inst_hom, path_parameter, sol0, t);
		break;
	case 6:
		// Single Path Tracking Test
		t = CT(0.0,0.0);
		workspace_cpu.update_x_t(sol0, t);
		path_test(workspace_cpu, cpu_inst_hom, path_parameter, sol0, sol_new, t, Target_Sys, device_option);
		break;
	case 7:
		// Path Tracking Reverse Test
		t = CT(0.0,0.0);
		workspace_cpu.update_x_t(sol0, t);
		path_test_reverse(workspace_cpu, cpu_inst_hom, path_parameter, sol0, t, device_option);
		break;
	case 8:
		// Witness Set Reverse Test
		t = CT(0.0,0.0);
		workspace_cpu.update_x_t(sol0, t);
		witness_set_test(workspace_cpu, cpu_inst_hom, path_parameter, sol0, t, device_option);
		break;
	case 10:
		// Multiple Path Tracking Test start_sol_set.n_sol
		path_multiple_test(workspace_cpu, cpu_inst_hom, path_parameter,\
				start_sol_set, Target_Sys, device_option, n_path);
		break;
	default:
		all_test(Target_Sys, workspace_cpu, cpu_inst_hom, path_parameter, sol0, n_eq, dim);
	}
	//delete[] sol0;

	return 0;
}

int parse_arguments
 ( int argc, char *argv[], int& test, int& dim, int& device_option, int& n_path)
/* Parses the arguments on the command line.
   Returns 0 if okay, otherwise the function returns 1. */
{
	if(argc < 4){
		cout << argv[0] << " needs 3 parameters: test, dim, device_option" << endl;
		cout << "Test  option: " << std::endl
		     << "  0: Predict Test" << endl \
		     << "  1: Evalute Test" << endl  \
		     << "  2: Modified Gram-Smith Test" << endl  \
		     << "  3: Modified Gram-Smith Large Test" << endl  \
		     << "  4: Modified Gram-Smith Large Test for Any Dim" << endl  \
		     << "  5: Newton Test" << endl  \
		     << "  6: Path Tracking Test" << endl  \
		     << "  7: Path Tracking Reverse Test" << endl  \
		     << "  8: Witeness Set Test" << endl  \
		     << "  9: Pieri Test" << endl  \
		     << "  else: All Test" << endl;
		cout << "Device option: " << std::endl \
			 << "  0. CPU and GPU (only for test 5 Path Test, 8 Pieri Test)" << std::endl \
			 << "  1. CPU         (only for test 5 Path Test, 6: Path Reverse Test, 7 Witeness Set Test, 8 Pieri Test)" << std::endl \
			 << "  2. GPU         (only for test 5 Path Test, 6: Path Reverse Test, 7 Witeness Set Test, 8 Pieri Test)" << std::endl \
			 << "All other test runs both on CPU and GPU." << device_option << std::endl;
		cout << "please try again..." << endl; return 1;
   }
   test = atoi(argv[1]);     // block size
   dim = atoi(argv[2]);    // dimension
   device_option = atoi(argv[3]);     // number of monomials
   if(argc==5){
	   n_path = atoi(argv[4]);
   }

   return 0;
}

int main ( int argc, char *argv[] )
{

   // Initialization of the execution parameters

   int test,dim,device_option;
   int n_path = -1;
   if(parse_arguments(argc,argv,test,dim,device_option, n_path) == 1) return 1;

   main_test(test,dim,device_option,n_path);
}

void all_test(PolySys& Target_Sys, Workspace& workspace_cpu, \
   CPUInstHom& cpu_inst_hom, Parameter path_parameter, CT* sol0, 
   int n_eq, int dim) {

	// Predict Test
	T1* err = new T1[5];

	CT t(0.5,0);
	err[0] = predict_test(workspace_cpu, cpu_inst_hom, t);

	// Evaluate Test
	t = CT(1,0);
	workspace_cpu.update_x_t_value(sol0, t);
	err[1] = eval_test_classic(workspace_cpu, cpu_inst_hom, sol0, t, Target_Sys);

	// Modified Gram-Smith Test
	t = CT(0.1,0);
	workspace_cpu.update_x_t_value(sol0, t);
	err[2] = mgs_test(workspace_cpu, cpu_inst_hom);

	// Newton Test
	t = CT(0.1,0);
	workspace_cpu.update_x_t_value(sol0, t);
	err[3] = newton_test(workspace_cpu, cpu_inst_hom, path_parameter, sol0, t);

	// Path Tracking Test
	t = CT(0.0,0.0);
	workspace_cpu.update_x_t(sol0, t);
	CT* sol_new = NULL;
	bool path_success = path_test(workspace_cpu, cpu_inst_hom, path_parameter, sol0, sol_new, t, Target_Sys);

	std::cout << "--------- Test Error Report ----------" << std::endl;
	std::cout << "Predict : " << err[0] << std::endl;
	std::cout << "Eval    : " << err[1] << std::endl;
	std::cout << "MGS     : " << err[2] << std::endl;
	std::cout << "Newton  : " << err[3] << std::endl;
	std::cout << "Path    : " << path_success << std::endl;

	delete[] err;
}
