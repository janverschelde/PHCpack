/*
 * CPU_instruction_eval.h
 *
 *  Created on: Nov 25, 2014
 *      Author: yxc
 */

#ifndef CPU_INSTRUCTION_EVAL_H_
#define CPU_INSTRUCTION_EVAL_H_

#include "varset.h"
//#include "utilities.h"
#include "poly.h"
#include "workspace_host.h"
#include <sys/time.h>
#include <unistd.h>
#include "path_data.h"

#define warp_size 32

class CPUInstHomCoef{
  public:
	int n_coef; // *2 = coef_orig size
	CT* coef_orig;
	CT alpha;

	CPUInstHomCoef(){
		n_coef = 0;
		coef_orig = NULL;
		alpha = CT(0.0,0);
	}

	CPUInstHomCoef(MonSet* hom_monset, int total_n_mon, int n_monset, int n_constant, CT alpha){
		init(hom_monset, total_n_mon, n_monset, n_constant, alpha);
	}

	~CPUInstHomCoef(){
		std::cout << "Delete CPUInstHomCoef" << std::endl;
		delete[] coef_orig;
	}

	void init(MonSet* hom_monset, int total_n_mon, int n_monset, int n_constant, CT alpha);

	void print();

	void eval(const CT t, CT* coef, int reverse=0);

	void update_alpha(CT alpha=CT(0.0,0.0));
};

class CPUInstHomMon{
  public:
	int level; // size of n_mon_level
	int* n_mon_level;
	int n_mon; // size of pos_start, sum of n_mon_level
	int* mon_pos_start;
	int mon_pos_size;
	unsigned short* mon_pos;
	int flops_multiple;

	CPUInstHomMon(){
		level = 0;
		n_mon_level = NULL;
		n_mon = 0;
		mon_pos_start = NULL;
		mon_pos_size = 0;
		mon_pos = NULL;
		flops_multiple = 0;
	}

	CPUInstHomMon(MonSet* hom_monset, int n_monset, int total_n_mon, int n_constant){
		CPUInstHomMon();
		init(hom_monset, n_monset, total_n_mon, n_constant);
	}

	~CPUInstHomMon(){
		std::cout << "Delete CPUInstHomMon" << std::endl;
		delete[] n_mon_level;
		delete[] mon_pos_start;
		delete[] mon_pos;
	}

	void init(MonSet* hom_monset, int n_monset, int total_n_mon, int n_constant);

	void eval(const CT* sol, CT* mon, CT* coef);

	void print();
};

class CPUInstHomMonBlock{
public:
	int n_mon;
	int BS;
	int NB;
	int mon_pos_block_size;
	int* mon_pos_start_block;
	unsigned short* mon_pos_block;
	unsigned short* max_var_block;
	int n_mon_single;
	unsigned short* mon_single_pos_block;

	CPUInstHomMonBlock(){
		n_mon = 0;
		n_mon_single = 0;
		BS = 0;
		NB = 0;
		mon_pos_block_size = 0;
		mon_pos_start_block = NULL;
		mon_pos_block = NULL;
		max_var_block = NULL;
		mon_single_pos_block = NULL;
	}

	CPUInstHomMonBlock(CPUInstHomMon& orig, int BS){
		init(orig, BS);
	}

	~CPUInstHomMonBlock(){
		delete[] max_var_block;
		delete[] mon_pos_start_block;
		delete[] mon_pos_block;
	}

	void init(CPUInstHomMon& orig, int BS);

	void print();
};

class CPUInstHomSumBlock{
  public:
	int n_sum; // size of sum_start
	int n_sum_levels;
	int* n_sum_level;
	int* n_sum_level_rest;
	int* sum_pos_start;
	int sum_pos_size;
	int* sum_pos;

	CPUInstHomSumBlock(){
		n_sum = 0;
		n_sum_levels = 0;
		n_sum_level = NULL;
		n_sum_level_rest = NULL;
		sum_pos_start = NULL;
		sum_pos_size = 0;
		sum_pos = NULL;
	}

	CPUInstHomSumBlock(MonSet* hom_monset, int n_monset, const int* mon_pos_start,
			             int dim, int n_eq, int n_constant, int n_mon0, int* mon_pos_start_block){
		init(hom_monset, n_monset, mon_pos_start, dim, n_eq, n_constant, n_mon0, mon_pos_start_block);
	}

	~CPUInstHomSumBlock(){
		std::cout << "Delete CPUInstHomSum" << std::endl;
		delete[] sum_pos_start;
		delete[] sum_pos;
	}

	void init(MonSet* hom_monset, int n_monset, const int* mon_pos_start,
			             int dim, int n_eq, int n_constant, int n_mon0, int* mon_pos_start_block);

	void eval(CT* sum, CT* matrix);

	void print();
};



class CPUInstHomSum{
  public:
	int n_sum; // size of sum_start
	int n_sum_levels;
	int* n_sum_level;
	int* n_sum_level_rest;
	int* sum_pos_start;
	int sum_pos_size;
	int* sum_pos;
	int n_sum_zero;
	int* sum_zeros;

	CPUInstHomSum(){
		n_sum = 0;
		n_sum_levels = 0;
		n_sum_level = NULL;
		n_sum_level_rest = NULL;
		sum_pos_start = NULL;
		sum_pos_size = 0;
		sum_pos = NULL;
		n_sum_zero = 0;
		sum_zeros = NULL;
	}

	CPUInstHomSum(MonSet* hom_monset, int n_monset, const int* mon_pos_start, \
			             int dim, int n_eq, int n_constant){
		init(hom_monset, n_monset, mon_pos_start, dim, n_eq, n_constant);
	}

	~CPUInstHomSum(){
		std::cout << "Delete CPUInstHomSum" << std::endl;
		delete[] sum_pos_start;
		delete[] sum_pos;
	}

	void init(MonSet* hom_monset, int n_monset, const int* mon_pos_start, \
			             int dim, int n_eq, int n_constant);

	void eval(CT* sum, CT* matrix);

	void print();
};

class CPUInstHomEq{
public:
	int n_eq;
	int* n_mon_eq;
	int* eq_pos_start;
	int n_mon_total;
	int* mon_pos_start_eq;
	int n_pos_total;
	unsigned short * mon_pos_eq;
	CT* coef;

	CPUInstHomEq(){
		n_eq = 0;
		n_mon_eq = NULL;
		eq_pos_start = NULL;
		n_mon_total = 0;
		mon_pos_start_eq = NULL;
		n_pos_total = 0;
		mon_pos_eq = NULL;
		coef = NULL;
	}

	CPUInstHomEq(MonSet* hom_monset, int n_monset, \
			         int n_eq, int n_constant){
		init(hom_monset, n_monset, n_eq, n_constant);
	}

	void init(MonSet* hom_monset, int n_monset, \
	         int n_eq, int n_constant);

	void print();
};

class CPUInstHom{
public:
    CPUInstHomCoef CPU_inst_hom_coef;
    CPUInstHomMon CPU_inst_hom_mon;
    CPUInstHomSum CPU_inst_hom_sum;
    CPUInstHomMonBlock CPU_inst_hom_block;
    CPUInstHomSumBlock CPU_inst_hom_sum_block;

    CPUInstHomEq CPU_inst_hom_eq;

    int n_constant;
    int dim;
    int n_eq;
    int n_coef;
    int n_predictor;

    // Record timing for both CPU and GPU
    Path path_data;
	double timeSec_Path_CPU;
	double timeSec_Path_GPU;

	bool success_CPU;
	bool success_GPU;
	CT t_CPU;
	CT t_GPU;
	T1 max_residual;
	T1 max_delta_x;

	int n_step_CPU;
	int n_step_GPU;
	int n_point_CPU;
	int n_point_GPU;
	int n_eval_CPU;
	int n_eval_GPU;
	int n_mgs_CPU;
	int n_mgs_GPU;

    void init(MonSet* hom_monset, int n_monset, \
    		  int n_constant, int total_n_mon, int dim, int n_eq, int n_predictor, CT alpha);

	void init(PolySys& Target_Sys, PolySys& Start_Sys, int dim, int n_eq, int n_predictor, CT alpha);

    CPUInstHom(){
    	n_constant = 0;
    	dim = 0;
    	n_eq = 0;
    	n_coef = 0;
    	n_predictor = 0;
    	// For record only
    	timeSec_Path_CPU = 0;
    	timeSec_Path_GPU = 0;
    	success_CPU = 0;
    	success_GPU = 0;
    	t_CPU = CT(0.0,0.0);
    	t_GPU = CT(0.0,0.0);
    	max_residual = 0.0;
    	max_delta_x = 0.0;
    	n_step_CPU = 0;
    	n_step_GPU = 0;
    	n_point_CPU = 0;
    	n_point_GPU = 0;
    	n_eval_CPU = 0;
    	n_eval_GPU = 0;
    	n_mgs_CPU = 0;
    	n_mgs_GPU = 0;
    }

    CPUInstHom(MonSet* hom_monset, int n_monset, int n_constant, int total_n_mon, int dim, int n_eq, int n_predictor, CT alpha)
	{
    	// For record only
    	timeSec_Path_CPU = 0;
    	timeSec_Path_GPU = 0;
    	success_CPU = 0;
    	success_GPU = 0;
    	t_CPU = CT(0.0,0.0);
    	t_GPU = CT(0.0,0.0);
    	max_residual = 0.0;
    	max_delta_x = 0.0;
    	n_step_CPU = 0;
    	n_step_GPU = 0;
    	n_point_CPU = 0;
    	n_point_GPU = 0;
    	n_eval_CPU = 0;
    	n_eval_GPU = 0;
    	n_mgs_CPU = 0;
    	n_mgs_GPU = 0;
    	init(hom_monset, n_monset, n_constant, total_n_mon, dim, n_eq, n_predictor, alpha);
	}

    ~CPUInstHom(){
    }

    void print();

    void init_workspace(Workspace& workspace_cpu);

    void eval(Workspace& workspace_cpu, const CT* sol, const CT t, int reverse=0);

    void update_alpha(CT alpha=CT(0.0,0.0));

};


#endif /* CPU_INSTRUCTION_EVAL_H_ */
