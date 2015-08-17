/* mgs_test.cpp created on Feb 8, 2015 by yxc with edits by jv */

#include "mgs_test.h"

void rand_matrix ( CT* matrix, int dim, int n_eq, int n_path )
{
   srand(time(NULL));
   // srand(1);
   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      for(int i=0; i<dim+1; i++)
      {
         for(int j=0; j<n_eq; j++)
         {
            /*if(i>=j)
              {
                 workspace_cpu.matrix[i*n_eq+j] = CT(1+i-j,1);
              }
              if(i<=j)
              {
                 workspace_cpu.matrix[i*n_eq+j] = CT(1-i+j,0);
              }
              if(i==j)
              {
                 workspace_cpu.matrix[i*n_eq+j] = CT((i+j)%5+1,0);
              }
              else
              {
                 workspace_cpu.matrix[i*n_eq+j] = CT(0.0,0);
                 workspace_cpu.matrix[i*n_eq+j] = CT(1,0);
              }
            */
            double temp = rand()/((double) RAND_MAX)*2*M_PI;
            matrix[path_idx*n_eq*(dim+1)+i*n_eq+j] = CT(cos(temp),sin(temp));
            // workspace_cpu.matrix[i*n_eq+j] = CT((i+j)%5+1,0);
            // std::cout << i << " " << j << " "
            //           << workspace_cpu.matrix[i*n_eq+j];
            // workspace_cpu.matrix[i*n_eq+j].imag = 0;
         }
      }
   }
}

T1 right_hand_check ( CT* tmp_matrix, CT* sol, int n_eq, int dim )
{
   CT* tmp_right = new CT[n_eq];

   for(int i = 0; i < n_eq; i++)
   {
      tmp_right[i] = CT(0.0,0);
      for(int j=0; j<dim; j++)
      {
         tmp_right[i] += sol[j]*tmp_matrix[j*n_eq+i];
      }
   }
   T1 err1 = err_check_workspace(tmp_right,tmp_matrix + dim * n_eq,n_eq);
   delete[] tmp_right;
   return err1;
}

T1 mgs_test_any
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom,
   int device_option, int n_path )
{
   std::cout << "--------------- Modified Gram-Smith Test ----------"
             << std::endl;

   if(n_path < 0)
   {
      std::cout << "Default n_path = 1000" << std::endl;
      n_path = 1000;
   }
   else
   {
      std::cout << "n_path = " << n_path << std::endl;
   }
   bool cpu_test = 0;
   bool gpu_test = 0;
   if(device_option == 0)
   {
      cpu_test = 1;
      gpu_test = 1;
      std::cout << "CPU + GPU Testing..." << std::endl;
   }
   else if(device_option == 1)
   {
      cpu_test = 1;
      gpu_test = 0;
      std::cout << "CPU Testing..." << std::endl;
   }
   else if(device_option == 2)
   {
      cpu_test = 0;
      gpu_test = 1;
      std::cout << "GPU Testing..." << std::endl;
   }
   CT* x = workspace_cpu.x;
   CT t = *(workspace_cpu.t);
   CT** V = (workspace_cpu.V);
   CT** R = (workspace_cpu.R);
   CT* sol_cpu = (workspace_cpu.sol);

   int dim = cpu_inst_hom.dim;
   int n_eq = cpu_inst_hom.n_eq; // to be removed

   std::cout << "cpu_inst_hom.dim  = " << cpu_inst_hom.dim << std::endl;
   std::cout << "cpu_inst_hom.n_eq = " << cpu_inst_hom.n_eq << std::endl;

   // cpu_inst_hom.eval(workspace_cpu, x, t);

   CT** sol_gpu;
   CT** matrix_gpu_q;
   CT** matrix_gpu_r;

   // Save tmp_matrix for right hand side check
   CT* tmp_matrix = new CT[n_path*n_eq*(dim+1)];
   rand_matrix(tmp_matrix, dim, n_eq, n_path);

   for(int i = 0; i < n_eq * (dim + 1); i++) 
   {
      workspace_cpu.matrix[i] = tmp_matrix[i];
   }
   // workspace_cpu.print_result(dim,n_eq);

   CT*** V_tmp = new CT**[n_path];
   CT**  V_col_tmp = new CT*[n_path*(dim+1)];
   CT*   V_data = new CT[n_path*(dim+1)*n_eq];

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      V_tmp[path_idx] = V_col_tmp + path_idx*(dim+1);
      for(int col_idx=0; col_idx<dim+1; col_idx++)
      {
         V_col_tmp[col_idx + path_idx*(dim+1)] = V_data
            + path_idx*(dim+1)*n_eq + col_idx*n_eq;
      }
   }
   for(int i = 0; i < n_path*n_eq*(dim + 1); i++)
   {
      V_data[i] = tmp_matrix[i];
   }
   CT*** R_tmp = new CT**[n_path];
   CT**  R_col_tmp = new CT*[n_path*(dim+1)];
   CT*   R_data = new CT[n_path*(dim+1)*(dim+1)];

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      R_tmp[path_idx] = R_col_tmp + path_idx*(dim+1);
      for(int col_idx=0; col_idx<dim+1; col_idx++)
      {
         R_col_tmp[col_idx + path_idx*(dim+1)] = R_data
            + path_idx*(dim+1)*n_eq + col_idx*(dim+1);
      }
   }
   CT** sol_cpu_tmp = new CT*[n_path];
   CT* sol_cpu_data = new CT[n_path*dim];
   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      sol_cpu_tmp[path_idx] = sol_cpu_data + path_idx*dim;
   }
   // GPU MGS
   // Run GPU test first, cpu test writes directly on V_data
   if(gpu_test == 1)
   {  // GPU MGS
      if(n_path == 1)
      {
         sol_gpu = new CT*[1];
         matrix_gpu_q = new CT*[1];
         matrix_gpu_r = new CT*[1];
         GPU_MGS(cpu_inst_hom,sol_gpu[0],matrix_gpu_q[0],matrix_gpu_r[0],
            V_data);
      }
      else
      {
         GPU_MGS_Mult(cpu_inst_hom,sol_gpu,matrix_gpu_q,matrix_gpu_r,
            V_data,n_path);
      }
   }
   // CPU MGS
   if(cpu_test == 1)
   {
      struct timeval start, end;
      long seconds, useconds;
      gettimeofday(&start, NULL);
      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         CPU_mgs2qrls(V_tmp[path_idx],R_tmp[path_idx],sol_cpu_tmp[path_idx],
            n_eq,dim + 1);
      }
      gettimeofday(&end, NULL);
      seconds  = end.tv_sec  - start.tv_sec;
      useconds = end.tv_usec - start.tv_usec;
      double timeMS_MGS_CPU = seconds*1000 + useconds/1000.0;
      double timeSec_MGS_CPU = timeMS_MGS_CPU/1000;
      std::cout << "CPU Time MS  = " << timeMS_MGS_CPU << std::endl;
      std::cout << "CPU Time Sec = " << timeSec_MGS_CPU << std::endl;
   }
   T1 err1 = 0;
   T1 err2 = 0;
   T1 err3 = 0;
   T1 err4 = 0;

   if(cpu_test == 1 && gpu_test == 1)
   {
      std::cout << "----------- Solution Check CPU vs GPU-------------"
                << std::endl;
      for(int path_idx = 0; path_idx < n_path; path_idx++)
      {
         T1 tmp_err = err_check_workspace(sol_cpu_tmp[path_idx],
                         sol_gpu[path_idx], cpu_inst_hom.dim);
         if(err2 < tmp_err)
         {
            err2 = tmp_err;
         }
      }
      std::cout << "err = " << err2 << std::endl;

      /*std::cout << "----------- Q Check CPU vs GPU-------------"
                  << (cpu_inst_hom.dim + 1) * cpu_inst_hom.n_eq << std::endl;
        for(int path_idx = 0; path_idx < n_path; path_idx++)
        {
           T1 tmp_err = err_check_workspace_matrix
              (V_data+path_idx*n_eq*(dim + 1), matrix_gpu_q[path_idx],
               cpu_inst_hom.n_eq,(cpu_inst_hom.dim+1));
           if(err3 < tmp_err)
           {
              err3 = tmp_err;
           }
        }
        std::cout << "err = " << err3 << std::endl;
        std::cout << "--------- R Check CPU vs GPU -----------" << std::endl;
        for(int path_idx = 0; path_idx < n_path; path_idx++)
        {
           T1 tmp_err = err_check_r(R_tmp[path_idx],matrix_gpu_r[path_idx],
                           cpu_inst_hom.dim, 0);
           if(err4 < tmp_err)
           {
              err4 = tmp_err;
           }
        }
        std::cout << "err = " << err4 << std::endl;
      */
   }
   std::cout << "----------- Right hand side check -------------" << std::endl;
   if(cpu_test == 1)
   {
      std::cout << "CPU ";
      err1 = 0.0;
      for(int path_idx = 0; path_idx < n_path; path_idx++)
      {
         T1 tmp_err = right_hand_check(tmp_matrix+path_idx*n_eq*(dim+1),
                         sol_cpu_tmp[path_idx], n_eq, dim);
         if(err1 < tmp_err)
         {
            err1 = tmp_err;
         }
      }
      std::cout << err1 << std::endl;
   }
   if(gpu_test == 1)
   {
      std::cout << "GPU ";
      err1 = 0.0;
      for(int path_idx = 0; path_idx < n_path; path_idx++)
      {
         T1 tmp_err = right_hand_check(tmp_matrix+path_idx*n_eq*(dim+1),
                         sol_gpu[path_idx], n_eq, dim);
         if(err1 < tmp_err)
         {
            err1 = tmp_err;
         }
      }
      std::cout << err1 << std::endl;
   }
   delete[] V_tmp;
   delete[] V_col_tmp;
   delete[] V_data;
   delete[] R_tmp;
   delete[] R_col_tmp;
   delete[] R_data;
   delete[] tmp_matrix;

   if(gpu_test == 1)
   {
      free(sol_gpu[0]);
      delete[] sol_gpu;
      // free(matrix_gpu_q[path_idx]);
      // free(matrix_gpu_r[path_idx]);
      // delete[] matrix_gpu_q;
      // delete[] matrix_gpu_r;
   }
   return max(err1, err2);
}

T1 mgs_test_large
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom )
{
   std::cout << "--------------- Modified Gram-Smith Test ----------"
             << std::endl;
   CT* x = workspace_cpu.x;
   CT t = *(workspace_cpu.t);
   CT** V = (workspace_cpu.V);
   CT** R = (workspace_cpu.R);
   CT* sol_cpu = (workspace_cpu.sol);
   int dim = cpu_inst_hom.dim;
   int n_eq = cpu_inst_hom.n_eq; // to be removed

   cpu_inst_hom.eval(workspace_cpu, x, t);

   CT* sol_gpu = NULL;
   CT* matrix_gpu_q = NULL;
   CT* matrix_gpu_r = NULL;

   /*for(int i=0; i<dim+1; i++)
     {
        for(int j=0; j<n_eq; j++)
        {
           if(i>=j)
           {
              workspace_cpu.matrix[i*n_eq+j] = CT(1+i-j,1);
           }
           if(i<=j)
           {
              workspace_cpu.matrix[i*n_eq+j] = CT(1-i+j,0);
           }
           if(i==j)
           {
              workspace_cpu.matrix[i*n_eq+j] = CT((i+j)%5+1,0);
           }
           else
           {
              workspace_cpu.matrix[i*n_eq+j] = CT(0.0,0);
              workspace_cpu.matrix[i*n_eq+j] = CT(1,0);
           }
           workspace_cpu.matrix[i*n_eq+j] = CT((i+j)%5+1,0);
           std::cout << i << " " << j << " " << workspace_cpu.matrix[i*n_eq+j];
           workspace_cpu.matrix[i*n_eq+j].imag = 0;
        }
     }
   */
   // Save matrix for right hand side check
   CT* tmp_matrix = new CT[n_eq*(dim+1)];
   for(int i = 0; i < n_eq * (dim + 1); i++)
   {
      tmp_matrix[i] = workspace_cpu.matrix[i];
   }
   // workspace_cpu.print_result(dim,n_eq);
   // GPU MGS
   GPU_MGS(cpu_inst_hom,sol_gpu,matrix_gpu_q,matrix_gpu_r,
      workspace_cpu.matrix);
   // CPU MGS
   CPU_mgs2qrls(V, R, sol_cpu, n_eq, dim + 1);

   std::cout << "----------- Right hand side check -------------" << std::endl;
   CT* tmp_right = new CT[n_eq];
   for(int i = 0; i < n_eq; i++)
   {
      tmp_right[i] = CT(0.0,0);
      for(int j=0; j<dim; j++)
      {
         tmp_right[i] += sol_cpu[j]*tmp_matrix[j*n_eq+i];
      }
   }
   T1 err1 = err_check_workspace(tmp_right,tmp_matrix + dim * n_eq,
                cpu_inst_hom.n_eq);
   delete[] tmp_matrix;
   delete[] tmp_right;
   T1 err2 = 0;
   T1 err3 = 0;
   std::cout << "----------- Solution Check CPU vs GPU-------------"
             << std::endl;
   err2 = err_check_workspace(sol_cpu, sol_gpu, cpu_inst_hom.dim);
   cout << "         sol[0] = " << sol_cpu[0];
   cout << "     sol_gpu[0] = " << sol_gpu[0];
   cout << "       v_cpu[0] = " << V[0][0];
   cout << "       v_gpu[0] = " << matrix_gpu_q[0];
   cout << "       v_cpu[1] = " << V[0][1];
   cout << "       v_gpu[1] = " << matrix_gpu_q[1];
   std::cout << "----------- Q Check CPU vs GPU-------------"
             << (cpu_inst_hom.dim + 1) * cpu_inst_hom.n_eq << std::endl;
   err3 = err_check_workspace_matrix(workspace_cpu.matrix,matrix_gpu_q,
             cpu_inst_hom.n_eq,(cpu_inst_hom.dim+1));
   std::cout << "----------- R Check CPU vs GPU -------------" << std::endl;
   err_check_r(R, matrix_gpu_r, cpu_inst_hom.dim, 0);

   free(sol_gpu);
   free(matrix_gpu_q);
   free(matrix_gpu_r);
   return max(err1, err2);
}

T1 mgs_test ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom)
{
   std::cout << "--------------- Modified Gram-Smith Test ----------"
             << std::endl;
   CT* x = workspace_cpu.x;
   CT t = *(workspace_cpu.t);
   CT** V = (workspace_cpu.V);
   CT** R = (workspace_cpu.R);
   CT* sol_cpu = (workspace_cpu.sol);
   int dim = cpu_inst_hom.dim;
   int n_eq = cpu_inst_hom.n_eq; // to be removed

   cpu_inst_hom.eval(workspace_cpu, x, t);

   CT* sol_gpu = NULL;
   CT* matrix_gpu_q = NULL;
   CT* matrix_gpu_r = NULL;

   // Save matrix for right hand side check
   CT* tmp_matrix = new CT[n_eq*(dim+1)];
   for(int i = 0; i < n_eq * (dim + 1); i++)
   {
      tmp_matrix[i] = workspace_cpu.matrix[i];
   }
   // workspace_cpu.print_result(dim,n_eq);
   // GPU MGS
   GPU_MGS(cpu_inst_hom,sol_gpu,matrix_gpu_q,matrix_gpu_r,
      workspace_cpu.matrix);
   // CPU MGS
   CPU_mgs2qrls(V, R, sol_cpu, n_eq, dim + 1);

   std::cout << "----------- Right hand side check -------------" << std::endl;
   CT* tmp_right = new CT[n_eq];
   for(int i = 0; i < n_eq; i++)
   {
      tmp_right[i] = CT(0.0,0);
      for(int j=0; j<dim; j++)
      {
         tmp_right[i] += sol_cpu[j]*tmp_matrix[j*n_eq+i];
      }
   }
   T1 err1 = err_check_workspace(tmp_right,tmp_matrix + dim * n_eq,
                                 cpu_inst_hom.n_eq);
   delete[] tmp_matrix;
   delete[] tmp_right;
   std::cout << "----------- Solution Check CPU vs GPU-------------"
             << std::endl;
   T1 err2 = err_check_workspace(sol_cpu, sol_gpu, cpu_inst_hom.dim);

   std::cout << "----------- Q Check CPU vs GPU -------------"
             << (cpu_inst_hom.dim + 1) * cpu_inst_hom.n_eq << std::endl;
   T1 err3 = err_check_workspace_matrix(workspace_cpu.matrix, matrix_gpu_q,
                cpu_inst_hom.n_eq,(cpu_inst_hom.dim + 1));

   std::cout << "----------- R Check CPU vs GPU -------------" << std::endl;
   err_check_r(R, matrix_gpu_r, dim);

   cout << "         sol[0] = " << sol_cpu[0];
   cout << "     sol_gpu[0] = " << sol_gpu[0];
   cout << "       v_cpu[0] = " << V[0][0];
   cout << "       v_gpu[0] = " << matrix_gpu_q[0];
   cout << "       v_cpu[1] = " << V[0][1];
   cout << "       v_gpu[1] = " << matrix_gpu_q[1];

   free(sol_gpu);
   free(matrix_gpu_q);
   free(matrix_gpu_r);

   return max(err1, err2);
}
