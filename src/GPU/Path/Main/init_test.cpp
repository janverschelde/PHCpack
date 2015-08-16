/* init_test.cpp, created on Feb 8, 2015 by yxc with edits by jv */

#include "init_test.h"

const std::string DATA = "/home/jan/Problems/GPUdata";

bool read_homotopy_file
 ( PolySys& Target_Sys, PolySys& Start_Sys, int dim, int& n_eq,
   string start_sys_filename, string target_sys_filename, PolySolSet* sol_set )
{
   ifstream myfile(start_sys_filename.c_str());
   ifstream myfile_target(target_sys_filename.c_str());

   if(myfile.is_open() == false)
   {
      std::cout << "Start System File is unable to open." << std::endl;
      return false;
   }
   if(myfile_target.is_open() == false)
   {
      std::cout << "Target System File is unable to open." << std::endl;
      return false;
   }
   VarDict pos_dict;

   Start_Sys.read_file(myfile,pos_dict);

   string x_name = "x";
   string* x_names = x_var(x_name, dim);
   Start_Sys.pos_var = x_names;

   if(dim < 8) 
   {
      std::cout << "Start Sys" << std::endl;
      Start_Sys.print();
   }
   if(sol_set != NULL)
   {
      sol_set -> init(myfile, pos_dict);
   }
   Target_Sys.read_file(myfile_target,pos_dict);

   string x_name_target = "x";
   string* x_names_target = x_var(x_name_target,dim);
   Target_Sys.pos_var = x_names_target;

   if(dim < 8)
   {
      std::cout << "Target Sys" << std::endl;
      Target_Sys.print();
   }
   n_eq = Start_Sys.n_eq;

   return true;
}

CT read_gamma ( int dim )
{
   std::ostringstream cyclic_filename_gamma;
   cyclic_filename_gamma << DATA << "/MultiPath/cyclic10.gamma";

   string filename_gamma = cyclic_filename_gamma.str();
   ifstream myfile_gamma(filename_gamma.c_str());

   if(myfile_gamma.is_open() == false)
   {
      std::cout << "Can't open gamma file." << std::endl;
      return CT(0.0,0.0);
   }
   double gamma_real;
   myfile_gamma >> gamma_real;
   double gamma_imag;
   myfile_gamma >> gamma_imag;

   return CT(gamma_real,gamma_imag);
}

void init_cpu_inst_workspace
 ( PolySys& Target_Sys, PolySys& Start_Sys,
   int dim, int n_eq, int n_predictor, 
   CPUInstHom& cpu_inst_hom, Workspace& workspace_cpu, int test )
{
   CT alpha;
   bool read_gamma_file = true;

   if(read_gamma_file)
   {
      alpha = read_gamma(10);
   }
   else
   {
      // srand(time(NULL));
      srand(1);
      int r = rand();
      T1 tmp = T1(r);
      alpha = CT(sin(tmp),cos(tmp));
      // alpha = CT(1,0);
      // alpha = CT(-2.65532737234004E-02,-9.99647399663787E-01);
   }
   // Fix gamma for evaluation test
   if(test == 1)
   {
      alpha = CT(1,0);
      //alpha = CT(1,0);
   }

   cpu_inst_hom.init(Target_Sys,Start_Sys,dim,n_eq,n_predictor,alpha);
   //cpu_inst_hom.print();
   cpu_inst_hom.init_workspace(workspace_cpu);

}

bool init_homotopy_test
 ( PolySys& Target_Sys, PolySys& Start_Sys, PolySolSet& sol_set, int dim,
   int& n_eq, CT*& sol0, CPUInstHom& cpu_inst_hom, Workspace& workspace_cpu,
   int test, int n_predictor,
   string start_sys_filename, string target_sys_filename )
{
   bool read_success = read_homotopy_file
      (Target_Sys, Start_Sys, dim, n_eq,
       start_sys_filename, target_sys_filename, &sol_set);

   if(read_success == false)
   {
      return false;
   }
   // sol0 = rand_val_complex_n(dim);

   if(sol_set.n_sol == 0) 
   {
      std::cout << "Error: No solution in start system file." << std::endl;
      return false;
   }
   sol0 = sol_set.get_sol(0);

   /*std::cout << "---------- Start Solution Top 10 ----------" << std::endl;
     for(int i = 0; i < min(dim, 10); i++)
     {
        std::cout << i << " " << sol0[i];
     }*/

   init_cpu_inst_workspace
      (Target_Sys,Start_Sys,dim,n_eq,n_predictor,
       cpu_inst_hom,workspace_cpu,test);

   return true;
}

bool init_homotopy_test_ada
 ( PolySys& Target_Sys, PolySys& Start_Sys, 
   CPUInstHom& cpu_inst_hom, Workspace& workspace_cpu, int n_predictor )
{
   int dim = Target_Sys.dim;
   int n_eq = Target_Sys.n_eq;
   int test = 2;
   init_cpu_inst_workspace
      (Target_Sys,Start_Sys,dim,n_eq,n_predictor,
       cpu_inst_hom,workspace_cpu,test);
   return true;
}
