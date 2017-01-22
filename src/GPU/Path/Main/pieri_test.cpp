/* This file "pieri_test.cpp" defines the function with prototype in
 * the file "pieri_test.h".
 * The file was written on 8 February 2015 by Xiangcheng Yu. */

#include "pieri_test.h"

void write_complex_array ( string file_name, CT* array, int dim )
{
   ofstream tfile(file_name.c_str());

   int pr = 2 * sizeof(T1);
   tfile.precision(pr+10);
   tfile << fixed;

   for(int i=0; i<dim; i++)
   {
      tfile << array[i].real << " "
            << array[i].imag << std::endl;
   }
   tfile.close();
};

CT* read_complex_array ( string file_name, int dim )
{
   CT* array = new CT[dim];
   ifstream tfile(file_name.c_str());
   for(int i=0; i<dim; i++)
   {
      array[i] = get_complex_number(tfile);
   }
   tfile.close();
   return array;
}

string sol_filename ( int dim, int sys_type = 2 )
{
   std::ostringstream sol_filename;

   if(sys_type == 0)
   {
      sol_filename << "../Problems/cyclic/cyc" << dim << "q1" << ".sol";
   }
   else if(sys_type == 2)
   {
      sol_filename << "../Problems/PieriBig1/pieri353target"
                   << dim-32 << ".sol";
   }
   else if(sys_type == 3)
   {
      sol_filename << "../Problems/PieriBig2/pieri364target"
                   << dim-32 << ".sol";
   }
   return sol_filename.str();
}

void pieri_sys_filename
 ( string& Start_Sys_filename, string& Target_Sys_file_name,
   int sys_type, int dim)
{
   std::ostringstream cyclic_filename;
   std::ostringstream cyclic_filename_target;

   if(sys_type == 0)
   {
      cyclic_filename << "../Problems/cyclic/cyc" << dim << "p1";
      cyclic_filename_target << "../Problems/cyclic/cyc" << dim << "q1";
   }
   else if(sys_type == 1)
   {
      cyclic_filename << "../Problems/MultiPath/cyclic7.start";
      cyclic_filename_target << "../Problems/MultiPath/cyclic7.target";
   }
   else if(sys_type == 2)
   {
      cyclic_filename << "../Problems/PieriBig1/pieri353start" << dim-32;
      cyclic_filename_target << "../Problems/PieriBig1/pieri353target"
                             << dim-32;
   }
   else
   {
      cyclic_filename << "../Problems/PieriBig2/pieri364start" << dim-32;
      cyclic_filename_target << "../Problems/PieriBig2/pieri364target"
                             << dim-32;
   }
   Start_Sys_filename = cyclic_filename.str();
   Target_Sys_file_name = cyclic_filename_target.str();
}

void init_test_pieri
 ( PolySys& Target_Sys, PolySys& Start_Sys, PolySolSet& sol_set,
   int dim, int& n_eq, CT*& sol0, CPUInstHom& cpu_inst_hom,
   Workspace& workspace_cpu, int test, int n_predictor, int sys_type = 2)
{
   string Start_Sys_filename;
   string Target_Sys_file_name;

   pieri_sys_filename(Start_Sys_filename,Target_Sys_file_name,sys_type,dim);

   if(dim == 32)
   {
      read_homotopy_file
         (Target_Sys,Start_Sys,dim,n_eq,Start_Sys_filename,
          Target_Sys_file_name,&sol_set);
      sol0 = sol_set.get_sol(0);
   }
   else
   {
      read_homotopy_file
         (Target_Sys,Start_Sys,dim,n_eq,Start_Sys_filename,
          Target_Sys_file_name);
      if(sol0 == NULL)
      {
         string pieri_sol_filename = sol_filename(dim-1, sys_type);
         sol0 = read_complex_array(pieri_sol_filename, dim-1);
      }
      CT* sol_new = new CT[dim];
      for(int i=0; i<dim-1; i++) sol_new[i] = sol0[i];
      sol_new[dim-1] = CT(0.0,0);
      delete[] sol0;
      sol0 = sol_new;
   }
   init_cpu_inst_workspace
      (Target_Sys,Start_Sys,dim,n_eq,n_predictor,
       cpu_inst_hom,workspace_cpu,test);
}

void Pieri_Test
 ( int dim_start, int dim_end, Parameter path_parameter,
   int sys_type, int device_option )
{
   CT* sol0 = NULL;
   int test = 7;

   std::ostringstream file_name;
   file_name << "../Problems/pieri_timing_" << sys_type-1;
   ofstream tfile(file_name.str().c_str());

   int pr = 2 * sizeof(T1);
   tfile.precision(pr+10);
   tfile << fixed;

   int dim = dim_start;

   for(dim=dim_start; dim<dim_end; dim++)
   {
      int n_eq;
      PolySys Target_Sys;
      PolySys Start_Sys;
      PolySolSet sol_set;
      CPUInstHom cpu_inst_hom;
      Workspace workspace_cpu;

      init_test_pieri
         (Target_Sys,Start_Sys,sol_set,dim,n_eq,sol0,cpu_inst_hom,
          workspace_cpu,test,path_parameter.n_predictor,sys_type);

      // Path Tracking Test
      CT t(0.0,0.0);
      workspace_cpu.init_x_t_idx();
      workspace_cpu.update_x_t(sol0, t);
      CT* sol_new = NULL;
      bool success = path_test
         (workspace_cpu,cpu_inst_hom,path_parameter,sol0,sol_new,t,
          Target_Sys,device_option);
      if(success == 0)
      {
         std::cout << "Pieri Test Fail!" << std::endl;
         std::cout << "dim = " << dim << std::endl;
         break;
      }
      for(int i=0; i<dim; i++) sol0[i] = sol_new[i];

      // free gpu solution memory
      if(device_option==0 or device_option==2) free(sol_new);
      /* CT* f_val = Target_Sys.eval(sol0);
         for(int i=0; i<n_eq; i++)
         {
             std::cout << i << " " << f_val[i] << std::endl;
         }*/

      string pieri_sol_filename = sol_filename(dim, sys_type);
      write_complex_array(pieri_sol_filename, sol0, dim);

      /* delete[] sol0;
         sol0 = read_complex_array(pieri_sol_filename.str(), dim);
         for(int i=0; i<dim; i++)
         {
            std::cout << sol0[i];
         }*/
      tfile << dim << " " << cpu_inst_hom.success_CPU
            << " " << cpu_inst_hom.success_GPU
            << " " << cpu_inst_hom.n_step_CPU
            << " " << cpu_inst_hom.n_step_GPU
            << " " << cpu_inst_hom.timeSec_Path_CPU
            << " " << cpu_inst_hom.timeSec_Path_GPU << std::endl;
   }
   if(dim == dim_end)
   {
      std::cout << "Pieri Test Success!" << std::endl;
      std::cout << "dim = " << dim << std::endl;
   }
   tfile.close();
}
