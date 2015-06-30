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

void all_test(PolySys& Target_Sys, Workspace& workspace_cpu,
			  CPUInstHom& cpu_inst_hom, Parameter path_parameter, CT* sol0, int n_eq, int dim);


int main_test(int test, int dim, int device_option, int n_path) {
	//int test = 7;
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
	// else: All Test

	Parameter path_parameter(N_PREDICTOR, MAX_STEP, MAX_IT, MAX_DELTA_T, MAX_DELTA_T_END,\
							 MIN_DELTA_T, ERR_MAX_RES, ERR_MAX_DELTA_X, \
			                 ERR_MAX_FIRST_DELTA_X, ERR_MIN_ROUND_OFF, \
			                 MAX_IT_REFINE, ERR_MIN_ROUND_OFF_REFINE, \
			                 STEP_INCREASE, STEP_DECREASE);

	int pr = 2 * sizeof(T1);
	std::cout.precision(pr);

	// Pieri_Test
	if(test == 9){
		int dim_start = 32;
		int dim_end = 36;
		int sys_type = 3;
		Pieri_Test(dim_start, dim_end, path_parameter, sys_type, device_option);
		return 0;
	}

	int sys_type = 1;
	if(test==10){
		//sys_type = 1;
	}
	int n_eq;
	CT* sol0;
	PolySys Target_Sys;
	PolySys Start_Sys;
	PolySolSet start_sol_set;
	CPUInstHom cpu_inst_hom;
	Workspace workspace_cpu;

	bool ada_test = true;
	if(ada_test==true){
		char start_file[] = "../../Problems/MultiPath/cyclic10.start";
		char target_file[] = "../../Problems/MultiPath/cyclic10.target";
                adainit()
		ada_read_homotopy(start_file, target_file, Start_Sys, Target_Sys, start_sol_set);
                adafinal()
		bool init_success = init_test_ada(Target_Sys, Start_Sys, \
				       cpu_inst_hom, workspace_cpu, path_parameter.n_predictor);
		if(init_success == false){
			std::cout << "Test Initialization Failed." << std::endl;
			return false;
		}
		sol0 = start_sol_set.get_sol(0);
		n_eq = 10;
		std::cout << "ada read" << std::endl;
	}
	else if(test != 4){
		bool init_success = init_test(Target_Sys, Start_Sys, start_sol_set, dim, n_eq, sol0, cpu_inst_hom, \
				  workspace_cpu, test, path_parameter.n_predictor, sys_type);
		if(init_success == false){
			std::cout << "Test Initialization Failed." << std::endl;
			return false;
		}
	}
	else{
		int n_eq = dim;
		cpu_inst_hom.dim = dim;
		cpu_inst_hom.n_eq = n_eq;
		int workspace_size = 5;
		int n_coef = 1;
		int n_constant = 1;
		int n_predictor = 1;
		workspace_cpu.init(workspace_size, n_coef, n_constant, n_eq, dim, n_predictor);
		sol0 = new CT[dim];
	}

	if (test == 0) {
		// Predict Test
		CT t(0.5,0);
		predict_test(workspace_cpu, cpu_inst_hom, t);
	}
	else if(test == 1) {
		// Evaluate Test
		std::cout << "Evaluation Test" << std::endl;
		CT t(T1(1),T1(0));
		eval_test_classic(workspace_cpu, cpu_inst_hom, sol0, t, Target_Sys, n_eq, dim, n_path);
	}
	else if(test == 2) {
		// Modified Gram-Smith Test
		CT t = CT(0.1,0);
		workspace_cpu.update_x_t_value(sol0, t);
		mgs_test(workspace_cpu, cpu_inst_hom);
	}
	else if(test == 3) {
		// Modified Gram-Smith Test
		CT t = CT(0.1,0);
		workspace_cpu.update_x_t_value(sol0, t);
		mgs_test_large(workspace_cpu, cpu_inst_hom);
	}
	else if(test == 4) {
		// Modified Gram-Smith Large Test for any dim
		CT t = CT(0.1,0);
		workspace_cpu.update_x_t_value(sol0, t);
		mgs_test_any(workspace_cpu, cpu_inst_hom, device_option, n_path);
	}
	else if(test == 5) {
		// Newton Test
		CT t = CT(0.0001,0);
		workspace_cpu.update_x_t_value(sol0, t);
		newton_test(workspace_cpu, cpu_inst_hom, path_parameter, sol0, t);
	}
	else if(test == 6) {
		// Path Tracking Test
		CT t(0.0,0.0);
		workspace_cpu.update_x_t(sol0, t);
		CT* sol_new = NULL;
		path_test(workspace_cpu, cpu_inst_hom, path_parameter, sol0, sol_new, t, Target_Sys, device_option);
	}
	else if(test == 7) {
		// Path Tracking Reverse Test
		CT t(0.0,0.0);
		workspace_cpu.update_x_t(sol0, t);
		path_test_reverse(workspace_cpu, cpu_inst_hom, path_parameter, sol0, t, device_option);
	}
	else if(test == 8) {
		// Witness Set Reverse Test
		CT t(0.0,0.0);
		workspace_cpu.update_x_t(sol0, t);
		witness_set_test(workspace_cpu, cpu_inst_hom, path_parameter, sol0, t, device_option);
	}
	else if(test == 10) {
		// Multiple Path Tracking Test start_sol_set.n_sol
		path_multiple_test(workspace_cpu, cpu_inst_hom, path_parameter,\
				start_sol_set, Target_Sys, device_option, n_path);
	}
	else if(test == 11){

	}
	else {
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
			  CPUInstHom& cpu_inst_hom, Parameter path_parameter, CT* sol0,  int n_eq, int dim) {

	// Predict Test
	T1* err = new T1[5];

	CT t(0.5,0);
	err[0] = predict_test(workspace_cpu, cpu_inst_hom, t);

	// Evaluate Test
	t = CT(1,0);
	workspace_cpu.update_x_t_value(sol0, t);
	err[1] = eval_test_classic(workspace_cpu, cpu_inst_hom, sol0, t, Target_Sys, n_eq, dim);

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

/*void generate_cyclic_system(PolySys& Target_Sys, PolySys& Start_Sys, PolySolSet& start_sol_set, int dim);

void generate_cyclic_system(PolySys& Target_Sys, PolySys& Start_Sys,
		PolySolSet& start_sol_set, int dim) {
	string* sys_string = string_cyclic(dim);
	for (int i = 0; i < dim; i++) {
		std::cout << sys_string[i] << std::endl;
	}

	VarDict pos_dict_target;
	Target_Sys.read(sys_string, dim, pos_dict_target);

	string x_name = "x";
	string* x_names = x_var(x_name, dim);
	Target_Sys.pos_var = x_names;

	Target_Sys.print();

	VarDict pos_dict;

	std::ostringstream cyclic_filename;
	cyclic_filename << "cyclic" << dim << ".start";

	ifstream myfile(cyclic_filename.str().c_str());

	Start_Sys.read_file(myfile, pos_dict);

	start_sol_set.init(myfile);
	//start_sol_set.print();

	string v_name = "z";
	string* v_names = x_var(v_name, dim);
	Start_Sys.pos_var = v_names;

	Start_Sys.print();
}*/
