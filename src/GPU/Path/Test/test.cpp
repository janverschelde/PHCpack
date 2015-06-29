#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string.h>

#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"
#include "jump_track.h"

#include "poly.h"

extern "C" void adainit(void);
extern "C" void adafinal(void);

void var_name(char* x_string, int x_string_len, string*& x_names){
	int n_var = 1;
	for(int pos=0; pos<x_string_len; pos++){
		if(x_string[pos]==' '){
			n_var++;
		}
	}
	std::cout << "n_var = " << n_var << std::endl;
	x_names = new string[n_var];

	int begin=0;
	string mystring = x_string;
	int var_idx = 0;
	for(int pos=0; pos<x_string_len; pos++){
		if(x_string[pos]==' '){
			x_names[var_idx] = mystring.substr(begin, pos-begin);
			begin = pos+1;
			var_idx++;
		}
	}
	x_names[var_idx] = mystring.substr(begin, x_string_len-begin);
}

void ada_read_sys(PolySys& sys)
{
	int fail;
	std::cout << "testing reading and writing a system" << std::endl;
	//fail = syscon_read_system();
	std::cout << "the system is .." << std::endl;
	fail = syscon_write_system();

	// Get variable names
	int s_dim = 80;
	char *s = (char*) calloc(80,sizeof(char));
	fail = syscon_string_of_symbols(&s_dim, s);
	string* x_names;
	var_name(s, s_dim, x_names);

	int dim = 4;
	int i = 1;
	double c[2];
	int d[dim];

	int n_eq = 0;
	fail = syscon_number_of_polynomials(&n_eq);

	sys.n_eq = n_eq;
	sys.dim  = dim;
	sys.eq_space = new PolyEq[n_eq];
	sys.pos_var = x_names;

	PolyEq* tmp_eq = sys.eq_space;

	for(int i=1; i<n_eq+1; i++){
		int nt;
		fail = syscon_number_of_terms(i,&nt);
		//std::cout << "  #terms in polynomial " << i << " : " << nt << std::endl;
		tmp_eq->n_mon = nt;
		tmp_eq->dim = dim;
		for(int j=1; j<=nt; j++)
		{
			fail = syscon_retrieve_term(i,j,dim,d,c);
			//std::cout << c[0] << " " << c[1] << std::endl;
			//for (int k=0; k<n; k++) std::cout << " " << d[k];
			//std::cout << std::endl;
			bool constant_term = true;
			for (int k=0; k<dim; k++){
				if(d[k]!=0){
					constant_term = false;
				}
			}

			if(constant_term==true){
				tmp_eq->n_mon--;
				tmp_eq->constant += CT(c[0],c[1]);
				//std::cout << "constant " << c[0] \
				          << " " << c[1] << std::endl;
			}
			else{
				PolyMon* a = new PolyMon(dim,d,c);
				tmp_eq->mon.push_back(a);
			}
		}

		tmp_eq->print(x_names);

		sys.eq.push_back(tmp_eq);
		tmp_eq++;
	}

	sys.print();
}

void ada_read_sols(PolySys& start_sys, PolySolSet& sols){
	int fail, len;
	//fail = copy_start_solutions_to_container();
	fail = solcon_number_of_solutions(&len);
	printf("Number of start solutions : %d\n",len);
	int dim=start_sys.dim;
	sols.dim = dim;
	double sol[2*dim+5];

	for(int sol_idx=1; sol_idx<len+1; sol_idx++){
		int mm;
		solcon_retrieve_solution(dim,sol_idx,&mm,sol);
		std::cout << sol[0] << " " << sol[1] << std::endl;

		for(int var_idx=0; var_idx<4; var_idx++){
			std::cout << sol[2+2*var_idx] << " " << sol[2+2*var_idx] << std::endl;
		}
		std::cout << sol[2+2*dim] \
				  << " " << sol[3+2*dim] \
				  << " " << sol[4+2*dim] << std::endl;

		PolySol* tmp_sol = new PolySol(dim, sol[0], sol[1], sol+2, \
				sol[2+2*dim], sol[3+2*dim], sol[4+2*dim]);
		//tmp_sol->print_info(start_sys.pos_var);
		sols.add_sol(tmp_sol);
	}
	sols.print_info(start_sys.pos_var);
}

void ada_read_homotopy(char* start_file, char* target_file, \
		PolySys& start_sys, PolySys& target_sys, PolySolSet& sols){
	int fail;

	std::cout << "Target " << std::endl;
	fail = read_standard_target_system_from_file(strlen(target_file), target_file);
	ada_read_sys(target_sys);

	std::cout << "Start " << std::endl;
	fail = read_standard_start_system_from_file(strlen(start_file),start_file);
	ada_read_sys(start_sys);
	ada_read_sols(start_sys, sols);
}

int main ( int argc, char *argv[] )
{
	adainit();

	char start_file[] = "../../Problems/cyclic/cyc4p1";
	char target_file[] = "../../Problems/cyclic/cyc4q1";

	PolySys target_sys;
	PolySys start_sys;
	PolySolSet sols;

	ada_read_homotopy(start_file, target_file, \
			start_sys, target_sys, sols);

	adafinal();
}
