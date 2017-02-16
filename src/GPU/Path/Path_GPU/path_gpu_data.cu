#ifndef __GPU_DATA_CU_
#define __GPU_DATA_CU_

#include "path_gpu_data.h"

int get_NB ( int n_job, int BS, int n_thread_per_job )
{
   return (n_job*n_thread_per_job - 1)/BS + 1;
}

dim3 get_grid ( int NB, int n_path )
{
   int NBy = 1;
   int NBx = NB;
   while(NBx > 65535)
   {
      NBy++;
      NBx = (NB-1)/NBy + 1;
   }
   return dim3(NBx,NBy,n_path);
}

dim3 get_grid ( int n_job, int BS, int n_path, int n_thread_per_job )
{
   int NB = get_NB(n_job, BS, n_thread_per_job);
   return get_grid(NB, n_path);
}

void GPUWorkspace::init_workspace_eq ( int n_pos_total_eq, int n_path )
{
   size_t size_workspace_eq = n_pos_total_eq*n_path*sizeof(GT);
   cudaMalloc((void **)&workspace_eq, size_workspace_eq);
}

void GPUWorkspace::init_matrix ( int dim, int n_eq )
{
   V = matrix;
   cudaMalloc((void **)&R, n_matrix_R*sizeof(GT));
   std::cout << "dim = " << dim << " " << dim*sizeof(GT) << std::endl;
   cudaMalloc((void **)&sol, dim*sizeof(GT));

   int rows = n_eq;
   int cols = dim+1;
   int row_block = (rows-1)/matrix_block_row+1;
   cudaMalloc((void**)&P, row_block*matrix_block_pivot_col*cols*sizeof(GT));
}

void GPUWorkspace::init_V_value ( CT* V_cpu )
{
   GT* V_host = (GT*)malloc(n_path*n_matrix*sizeof(GT));
   // std::cout << "----------------n_eq = " << n_eq << " dim = "
   //           << dim << std::endl;
   for(int i=0; i<n_path*n_matrix; i++)
   {
      comp1_qd2gqd(&V_cpu[i],&V_host[i]);
   }
   cudaMemcpy(V, V_host, n_path*n_matrix*sizeof(GT), cudaMemcpyHostToDevice);
}

void GPUWorkspace::init_x_t ( int dim, int n_predictor )
{
   x_t_idx = 0;
   this->n_predictor = n_predictor;
   this->dim = dim;
   n_array = n_predictor + 1;
   cudaMalloc((void **)&x_array, n_array*dim*sizeof(GT));
   cudaMalloc((void **)&t_array, n_array*sizeof(GT));

   x = x_array;
   t = t_array;

   cudaMalloc((void **)&one_minor_t, sizeof(GT));
}

void GPUWorkspace::update_x_t_idx()
{
   x_last = x;
   t_last = t;
   x_t_idx = (x_t_idx+1)%n_array;
   cudaMemcpy(x_t_idx_mult, &x_t_idx, sizeof(int), cudaMemcpyHostToDevice);
   x = x_array + dim*x_t_idx;
   t = t_array + x_t_idx;
}

void GPUWorkspace::update_x_t_idx_all ( int* x_t_idx_host )
{
   size_t x_t_idx_mult_size = n_path*sizeof(int);
   cudaMemcpy(x_t_idx_mult, x_t_idx_host,
              x_t_idx_mult_size, cudaMemcpyHostToDevice);
   // update x and t pointer
   if(n_path==1)
   {
      x = x_array + dim*x_t_idx_host[0];
      t = t_array + x_t_idx_host[0];
   }
}

int* GPUWorkspace::get_x_t_idx_all()
{
   size_t x_t_idx_mult_size = n_path*sizeof(int);
   int* x_t_idx_host0 = (int*)malloc(x_t_idx_mult_size);
   cudaMemcpy(x_t_idx_host0, x_t_idx_mult,
              x_t_idx_mult_size, cudaMemcpyDeviceToHost);
   return x_t_idx_host0;
}

void GPUWorkspace::update_t_value ( CT cpu_t )
{
   GT* host_t = (GT*)malloc(sizeof(GT));
   GT* host_one_minor_t = (GT*)malloc(sizeof(GT));
   GT* tmp_t = t;
   GT* tmp_one_minor_t = one_minor_t;
   for(int i=0; i<n_path; i++)
   {
      CT cpu_one_minor_t(1.0-cpu_t.real, -cpu_t.imag);
      cpu_one_minor_t *= alpha;
      comp1_qd2gqd(&cpu_t, host_t);
      cudaMemcpy(tmp_t, host_t, sizeof(GT), cudaMemcpyHostToDevice);

      comp1_qd2gqd(&cpu_one_minor_t, host_one_minor_t);
      cudaMemcpy(tmp_one_minor_t, host_one_minor_t,
                 sizeof(GT), cudaMemcpyHostToDevice);
      tmp_t += workspace_size;
      tmp_one_minor_t += workspace_size;
   }
   free(host_t);
   free(host_one_minor_t);
}

void GPUWorkspace::update_t_value_array ( CT* cpu_t, int* x_t_idx_host )
{
   GT* tmp_t = t_array;
   GT* tmp_one_minor_t = one_minor_t;

   GT* host_t = (GT*)malloc(n_path*(n_predictor+1)*sizeof(GT));
   GT* host_one_minor_t = (GT*)malloc(n_path*sizeof(GT));
   GT* host_t_tmp = host_t;
   CT* cpu_t_tmp = cpu_t;

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      int x_t_idx_tmp = x_t_idx_host[path_idx];
      for(int pred_idx=0; pred_idx<(n_predictor+1); pred_idx++)
      {
         comp1_qd2gqd(cpu_t_tmp+pred_idx, host_t_tmp+pred_idx);
      }
      CT cpu_one_minor_t(1.0-cpu_t_tmp[x_t_idx_tmp].real,
                            -cpu_t_tmp[x_t_idx_tmp].imag);
      cpu_one_minor_t *= alpha;
      comp1_qd2gqd(&cpu_one_minor_t, host_one_minor_t+path_idx);
      cpu_t_tmp += (n_predictor+1);
      host_t_tmp += (n_predictor+1);
   }
   cudaMemcpy(tmp_t, host_t,
              n_path*(n_predictor+1)*sizeof(GT), cudaMemcpyHostToDevice);
   cudaMemcpy(tmp_one_minor_t, host_one_minor_t,
              n_path*sizeof(GT), cudaMemcpyHostToDevice);
   free(host_t);
   free(host_one_minor_t);
}

void GPUWorkspace::update_t_value_mult(CT* cpu_t)
{
   GT* host_t = (GT*)malloc(sizeof(GT));
   GT* host_one_minor_t = (GT*)malloc(sizeof(GT));
   GT* tmp_t = t_array;
   GT* tmp_one_minor_t = one_minor_t;
   CT* cpu_t_tmp = cpu_t;

   size_t x_t_idx_mult_size = n_path*sizeof(int);
   int* x_t_idx_mult_host = (int*)malloc(x_t_idx_mult_size);
   cudaMemcpy(x_t_idx_mult_host, x_t_idx_mult,
              x_t_idx_mult_size, cudaMemcpyDeviceToHost);

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      // std::cout << "x_t_idx_mult_host = "
      //           << x_t_idx_mult_host[path_idx] << std::endl;
      // std::cout << cpu_t[path_idx];
      comp1_qd2gqd(cpu_t_tmp, host_t);
      cudaMemcpy(tmp_t+x_t_idx_mult_host[path_idx], host_t,
                 sizeof(GT), cudaMemcpyHostToDevice);
      CT cpu_one_minor_t(1.0-cpu_t_tmp->real, - cpu_t_tmp->imag);
      cpu_one_minor_t *= alpha;
      comp1_qd2gqd(&cpu_one_minor_t, host_one_minor_t);
      cudaMemcpy(tmp_one_minor_t, host_one_minor_t,
                 sizeof(GT), cudaMemcpyHostToDevice);
      tmp_t += workspace_size;
      tmp_one_minor_t += workspace_size;
      cpu_t_tmp++;
   }
   free(host_t);
   free(host_one_minor_t);
}

void GPUWorkspace::update_t_value_inverse ( CT cpu_one_minor_t )
{
   CT cpu_t(1.0-cpu_one_minor_t.real, -cpu_one_minor_t.imag);
   cpu_one_minor_t *= alpha;

   GT* host_t = (GT*)malloc(sizeof(GT));
   comp1_qd2gqd(&cpu_t, host_t);
   cudaMemcpy(t, host_t, sizeof(GT), cudaMemcpyHostToDevice);
   free(host_t);

   GT* host_one_minor_t = (GT*)malloc(sizeof(GT));
   comp1_qd2gqd(&cpu_one_minor_t, host_one_minor_t);
   cudaMemcpy(one_minor_t, host_one_minor_t,
              sizeof(GT), cudaMemcpyHostToDevice);
   free(host_one_minor_t);
}

void GPUWorkspace::update_x_value ( CT* cpu_sol0 )
{
   GT* host_sol0 = (GT*)malloc(n_path*dim*sizeof(GT));
   for(int i=0; i<n_path*dim; i++)
   {
      comp1_qd2gqd(&cpu_sol0[i],&host_sol0[i]);
   }
   size_t sol0_size = dim*sizeof(GT);
   GT* tmp_x_array = x_array;
   GT* tmp_host_sol0 = host_sol0;
   for(int i=0; i<n_path; i++)
   {
      cudaMemcpy(tmp_x_array, tmp_host_sol0,
                 sol0_size, cudaMemcpyHostToDevice);
      tmp_x_array += workspace_size;
      tmp_host_sol0 += dim;
   }
   free(host_sol0);
}

void GPUWorkspace::update_x_value_array ( CT* cpu_sol0 )
{
   size_t sol0_size = n_path*(n_predictor+1)*dim*sizeof(GT);
   GT* host_sol0 = (GT*)malloc(sol0_size);
   for(int i=0; i<n_path*(n_predictor+1)*dim; i++)
   {
      comp1_qd2gqd(&cpu_sol0[i],&host_sol0[i]);
   }
   cudaMemcpy(x_array, host_sol0, sol0_size, cudaMemcpyHostToDevice);
   free(host_sol0);
}

void GPUWorkspace::update_x_value_mult ( CT* cpu_sol0 )
{
   GT* host_sol0 = (GT*)malloc(n_path*dim*sizeof(GT));
   for(int i=0; i<n_path*dim; i++)
   {
      comp1_qd2gqd(&cpu_sol0[i],&host_sol0[i]);
   }
   size_t sol0_size = dim*sizeof(GT);
   GT* tmp_x_array = x_array;
   GT* tmp_host_sol0 = host_sol0;
   for(int i=0; i<n_path; i++)
   {
      cudaMemcpy(tmp_x_array, tmp_host_sol0,
                 sol0_size, cudaMemcpyHostToDevice);
      tmp_x_array += workspace_size;
      tmp_host_sol0 += dim;
   }
   free(host_sol0);
}

void GPUWorkspace::update_x_t_value ( CT* cpu_sol0, CT cpu_t )
{
   update_x_value(cpu_sol0);
   update_t_value(cpu_t);
}

void GPUWorkspace::update_x_t_value_array
 ( CT* cpu_sol0, CT* cpu_t, int* x_t_idx_host )
{
   update_x_value_array(cpu_sol0);
   update_t_value_array(cpu_t, x_t_idx_host);
}

void GPUWorkspace::update_x_t_value_mult ( CT* cpu_sol0, CT* cpu_t )
{
   update_x_value_mult(cpu_sol0);
   update_t_value_mult(cpu_t);
}

void GPUWorkspace::update_x_mult_vertical ( CT* cpu_sol0, int* x_t_idx_host )
{
   CT* start_sol_gpu = (CT*)malloc(n_path*dim*sizeof(CT));
   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      CT* tmp_sol
         = cpu_sol0+path_idx*dim*(n_predictor+1)+dim*x_t_idx_host[path_idx];
      for(int var_idx=0; var_idx<dim; var_idx++)
      {
         start_sol_gpu[path_idx*dim+var_idx] = tmp_sol[var_idx];
      }
   }
   // Vertical multiple x
   GT* host_sol0 = (GT*)malloc(n_path*dim*sizeof(GT));
   for(int i=0; i<n_path; i++)
   {
      for(int j=0; j<dim; j++)
      {
         comp1_qd2gqd(&start_sol_gpu[i*dim+j],&host_sol0[i+j*n_path]);
      }
   }
   cudaMemcpy(x_mult, host_sol0,
              n_path*dim*sizeof(GT), cudaMemcpyHostToDevice);

   // horizontal multiple x
   /*
    for(int i=0; i<n_path; i++)
    {
       for(int j=0; j<dim; j++)
       {
          comp1_qd2gqd(&start_sol_gpu[i*dim+j],&host_sol0[i*dim+j]);
       }
    }
    cudaMemcpy(x_mult, host_sol0,
               n_path*dim*sizeof(GT), cudaMemcpyHostToDevice);
    */

   free(host_sol0);
   free(start_sol_gpu);
}

void GPUWorkspace::update_x_mult_horizontal ( CT* cpu_sol0 )
{
   // Vertical multiple x
   GT* host_sol0 = (GT*)malloc(n_path*dim*sizeof(GT));
   // horizontal multiple x
   for(int i=0; i<n_path; i++)
   {
      for(int j=0; j<dim; j++)
      {
         comp1_qd2gqd(&cpu_sol0[i*dim+j],&host_sol0[i*dim+j]);
      }
   }
   cudaMemcpy(x_mult, host_sol0,
              n_path*dim*sizeof(GT), cudaMemcpyHostToDevice);
   free(host_sol0);
}

void GPUWorkspace::update_t_value_mult2 ( CT* cpu_t, int* x_t_idx_host )
{
   size_t t_mult_size = n_path*sizeof(CT);
   CT* t_gpu = (CT*)malloc(t_mult_size);
   CT* one_minor_t_gpu = (CT*)malloc(t_mult_size);
   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      CT tmp_t = cpu_t[path_idx*(n_predictor+1)+x_t_idx_host[path_idx]];
      // std::cout << "tmp_t = " << tmp_t << std::endl;
      CT cpu_one_minor_t(1.0-tmp_t.real, -tmp_t.imag);
      cpu_one_minor_t *= alpha;
      t_gpu[path_idx] = tmp_t;
      one_minor_t_gpu[path_idx] = cpu_one_minor_t;
   }
   GT* host_t = (GT*)malloc(t_mult_size);
   GT* host_one_minor_t = (GT*)malloc(t_mult_size);
   for(int i=0; i<n_path; i++)
   {
      comp1_qd2gqd(&t_gpu[i],&host_t[i]);
      comp1_qd2gqd(&one_minor_t_gpu[i],&host_one_minor_t[i]);
   }
   cudaMemcpy(t_mult, host_t, t_mult_size, cudaMemcpyHostToDevice);
   cudaMemcpy(newton_t_mult, host_t, t_mult_size, cudaMemcpyHostToDevice);
   cudaMemcpy(one_minor_t, host_one_minor_t,
              t_mult_size, cudaMemcpyHostToDevice);

   free(t_gpu);
   free(one_minor_t_gpu);
   free(host_t);
   free(host_one_minor_t);
}

void GPUWorkspace::update_x_t ( CT* cpu_sol0, CT cpu_t )
{
   update_x_t_value(cpu_sol0, cpu_t);
   update_x_t_idx();
}

CT* GPUWorkspace::get_coef_mult()
{
   size_t coef_mult_size = n_path*n_coef*sizeof(GT);
   GT* coef_mult_host = (GT*)malloc(coef_mult_size);
   cudaMemcpy(coef_mult_host, coef, coef_mult_size, cudaMemcpyDeviceToHost);

   CT* coef_mult_gpu = (CT*)malloc(coef_mult_size);
   for(int i=0; i<n_coef*n_path; i++)
   {
      comp1_gqd2qd(&coef_mult_host[i], &coef_mult_gpu[i]);
   }
   free(coef_mult_host);

   return coef_mult_gpu;
}

CT* GPUWorkspace::get_mon_mult()
{
   size_t mon_mult_size = n_path*mon_pos_size*sizeof(GT);
   GT* mon_mult_host = (GT*)malloc(mon_mult_size);
   cudaMemcpy(mon_mult_host, mon, mon_mult_size, cudaMemcpyDeviceToHost);

   CT* mon_mult_gpu = (CT*)malloc(mon_mult_size);
   for(int i=0; i<mon_pos_size*n_path; i++)
   {
      comp1_gqd2qd(&mon_mult_host[i], &mon_mult_gpu[i]);
   }
   free(mon_mult_host);

   return mon_mult_gpu;
}

CT* GPUWorkspace::get_matrix ( int path_idx )
{
   GT* host_matrix = (GT*)malloc(n_matrix*sizeof(GT));
   cudaMemcpy(host_matrix, matrix+path_idx*n_matrix,
              n_matrix*sizeof(GT), cudaMemcpyDeviceToHost);

   CT* gpu_matrix = (CT*)malloc(n_matrix*sizeof(GT));
   for(int i=0; i<n_matrix; i++) comp1_gqd2qd(&host_matrix[i], &gpu_matrix[i]);

   free(host_matrix);

   return gpu_matrix;
}

CT* GPUWorkspace::get_matrix()
{
   size_t matrix_mult_size = n_path*n_matrix*sizeof(GT);
   GT* matrix_mult_host = (GT*)malloc(matrix_mult_size);
   cudaMemcpy(matrix_mult_host, matrix,
              matrix_mult_size, cudaMemcpyDeviceToHost);

   CT* matrix_mult_gpu = (CT*)malloc(matrix_mult_size);
   for(int i=0; i<n_matrix*n_path; i++)
   {
      comp1_gqd2qd(&matrix_mult_host[i], &matrix_mult_gpu[i]);
   }
   free(matrix_mult_host);
   return matrix_mult_gpu;
}

CT** GPUWorkspace::get_matrix_mult()
{
   CT* gpu_matrix_mult = get_matrix();
   CT**gpu_matrix_all = new CT*[n_path];
   gpu_matrix_all[0] = new CT[n_path*n_matrix];
   for(int path_idx=1; path_idx<n_path; path_idx++)
   {
      gpu_matrix_all[path_idx] = gpu_matrix_all[path_idx-1]+n_matrix;
   }
   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      for(int i=0; i<n_matrix; i++)
      {
         gpu_matrix_all[path_idx][i] = gpu_matrix_mult[n_path*i+path_idx];
      }
   }
   return gpu_matrix_all;
}

CT* GPUWorkspace::get_workspace ( int path_idx )
{
   // XXXXX it is wrong, not path_idx
   int n_arrays =  n_coef+mon_pos_size;
   GT* host_workspace = (GT*)malloc(n_arrays*sizeof(GT));
   cudaMemcpy(host_workspace, coef,
              (n_arrays)*sizeof(GT), cudaMemcpyDeviceToHost);
   CT* gpu_workspace = (CT*)malloc((n_arrays)*sizeof(GT));
   for(int i=0; i<n_arrays; i++)
      comp1_gqd2qd(&host_workspace[i], &gpu_workspace[i]);

   free(host_workspace);

   return gpu_workspace;
}

CT* GPUWorkspace::get_workspace()
{
   GT* host_workspace = (GT*)malloc(size_GT_arrays);
   cudaMemcpy(host_workspace, GT_arrays,
              size_GT_arrays, cudaMemcpyDeviceToHost);
   CT* gpu_workspace = (CT*)malloc(size_GT_arrays);

   for(int i=0; i<n_GT_arrays; i++)
   {
      comp1_gqd2qd(&host_workspace[i], &gpu_workspace[i]);
   }
   free(host_workspace);

   return gpu_workspace;
}

CT** GPUWorkspace::get_workspace_mult()
{
   CT* gpu_workspace_mult = get_workspace();
   CT** gpu_workspace_all = new CT*[n_path];
   int n_workspace= n_coef+mon_pos_size;
   gpu_workspace_all[0] = new CT[n_path*n_workspace];

   for(int path_idx=1; path_idx<n_path; path_idx++)
   {
      gpu_workspace_all[path_idx] = gpu_workspace_all[path_idx-1]+n_workspace;
   }
   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      for(int i=0; i<n_workspace; i++)
      {
         gpu_workspace_all[path_idx][i]
            = gpu_workspace_mult[n_path*i+path_idx];
      }
   }
   return gpu_workspace_all;
}

CT* GPUWorkspace::get_matrix_r ( int path_idx )
{
   GT* host_matrix_r = (GT*)malloc(n_matrix_R*sizeof(GT));
   cudaMemcpy(host_matrix_r, R+path_idx*n_matrix_R,
              n_matrix_R*sizeof(GT), cudaMemcpyDeviceToHost);

   CT* gpu_matrix_r = (CT*)malloc(n_matrix_R*sizeof(GT));
   for(int i=0; i<n_matrix_R; i++)
      comp1_gqd2qd(&host_matrix_r[i], &gpu_matrix_r[i]);

   free(host_matrix_r);

   return gpu_matrix_r;
}

void GPUWorkspace::print_matrix_r()
{
   CT* R = get_matrix_r();
   int tmp_idx_r = 0;
   for(int i=0; i<dim+1; i++)
   {
      for(int j=0; j<dim+1-i; j++)
      {
         if(i!=0 or j!= dim)
         {
            std::cout << dim-i << " " << j << " " << R[tmp_idx_r];
         }
         tmp_idx_r++;
      }
   }
}

CT* GPUWorkspace::get_x()
{
   size_t x_size = dim*sizeof(GT);
   GT* x_host = (GT*)malloc(x_size);
   cudaMemcpy(x_host, x, x_size, cudaMemcpyDeviceToHost);

   CT* x_cpu = (CT*)malloc(x_size);
   for(int i=0; i<dim; i++) comp1_gqd2qd(&x_host[i], &x_cpu[i]);

   free(x_host);

   return x_cpu;
}

CT** GPUWorkspace::get_mult_x_horizontal()
{
   CT** x_all_cpu = new CT*[n_path];
   size_t x_size = dim*sizeof(GT);

   GT* x_host = (GT*)malloc(n_path*x_size);

   cudaMemcpy(x_host, x_mult, n_path*x_size, cudaMemcpyDeviceToHost);

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      CT* x_cpu = (CT*)malloc(x_size);
      for(int var_idx=0; var_idx<dim; var_idx++)
      {
         comp1_gqd2qd(&x_host[path_idx*dim+var_idx], &x_cpu[var_idx]);
      }
      x_all_cpu[path_idx] = x_cpu;
   }
   free(x_host);

   return x_all_cpu;
}

CT* GPUWorkspace::get_x_array()
{
   GT* x_array_host = (GT*)malloc(n_array*dim*sizeof(GT));
   cudaMemcpy(x_array_host, x_array,
              n_array*dim*sizeof(GT), cudaMemcpyDeviceToHost);

   CT* x_array_ct = (CT*)malloc(n_array*dim*sizeof(GT));
   for(int i=0; i<n_array*dim; i++)
      comp1_gqd2qd(&x_array_host[i], &x_array_ct[i]);

   free(x_array_host);

   return x_array_ct;
}

CT** GPUWorkspace::get_x_all()
{
   CT** x_all_cpu = new CT*[n_path];

   size_t x_t_idx_mult_size = n_path*sizeof(int);
   int* x_t_idx_mult_host = (int*)malloc(x_t_idx_mult_size);
   cudaMemcpy(x_t_idx_mult_host, x_t_idx_mult,
              x_t_idx_mult_size, cudaMemcpyDeviceToHost);

   GT* x_host = (GT*)malloc(dim*sizeof(GT));

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      std::cout << "x_t_idx_mult_host[path_idx] = "
                << x_t_idx_mult_host[path_idx] << std::endl;
      cudaMemcpy(x_host, x_array+path_idx*dim*n_array
                 +dim*x_t_idx_mult_host[path_idx],
                 dim*sizeof(GT), cudaMemcpyDeviceToHost);
      CT* x_cpu = (CT*)malloc(dim*sizeof(GT));
      for(int i=0; i<dim; i++)
      {
         comp1_gqd2qd(&x_host[i], &x_cpu[i]);
         // std::cout << i << " " << sol_cpu[i] << std::endl;
      }
      x_all_cpu[path_idx] = x_cpu;
   }
   free(x_host);

   return x_all_cpu;
}

CT** GPUWorkspace::get_x_last_all()
{
   CT** x_all_cpu = new CT*[n_path];

   size_t x_t_idx_mult_size = n_path*sizeof(int);
   int* x_t_idx_mult_host = (int*)malloc(x_t_idx_mult_size);
   cudaMemcpy(x_t_idx_mult_host, x_t_idx_mult,
              x_t_idx_mult_size, cudaMemcpyDeviceToHost);

   GT* x_host = (GT*)malloc(dim*sizeof(GT));

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      // std::cout << "x_t_idx_mult_host[path_idx] = "
      //           << x_t_idx_mult_host[path_idx] << std::endl;
      cudaMemcpy(x_host, x_array+path_idx*dim*n_array
         +dim*((x_t_idx_mult_host[path_idx]+n_predictor)%(n_predictor+1)),
         dim*sizeof(GT), cudaMemcpyDeviceToHost);
      CT* x_cpu = (CT*)malloc(dim*sizeof(GT));
      for(int i=0; i<dim; i++)
      {
         comp1_gqd2qd(&x_host[i], &x_cpu[i]);
         // std::cout << i << " " << sol_cpu[i] << std::endl;
      }
      x_all_cpu[path_idx] = x_cpu;
   }
   free(x_host);

   return x_all_cpu;
}

void GPUWorkspace::print_x()
{
   CT* x_ct = get_x();

   for(int i=0; i<dim; i++)
   {
      std::cout << i << " "  << x_ct[i];
   }
   free(x_ct);
}


void GPUWorkspace::print_x_mult ( int path_idx_one )
{
   CT** x_all_cpu = get_x_all();

   if(path_idx_one==-1)
   {
      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         std::cout << "path_idx = " << path_idx << std::endl;
         for(int i=0; i<dim; i++)
         {
            std::cout << i << " "  << x_all_cpu[path_idx][i];
         }
      }
   }
   else
   {
      for(int i=0; i<dim; i++)
      {
         std::cout << i << " "  << x_all_cpu[path_idx_one][i];
      }
   }
   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      free(x_all_cpu[path_idx]);
   }
   free(x_all_cpu);
}

void GPUWorkspace::print_t_mult ( int path_idx_one )
{
   CT* t_all_cpu = new CT[n_path];

   size_t x_t_idx_mult_size = n_path*sizeof(int);
   int* x_t_idx_mult_host = (int*)malloc(x_t_idx_mult_size);
   cudaMemcpy(x_t_idx_mult_host, x_t_idx_mult,
              x_t_idx_mult_size, cudaMemcpyDeviceToHost);

   GT* t_host = (GT*)malloc(n_path*sizeof(GT));

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      // std::cout << "x_t_idx_mult_host[path_idx] = "
      //           << x_t_idx_mult_host[path_idx] << std::endl;
      cudaMemcpy(t_host+path_idx,
         t_array+path_idx*(n_predictor+1)+x_t_idx_mult_host[path_idx],
         sizeof(GT), cudaMemcpyDeviceToHost);
      // cudaMemcpy(t_host+path_idx,
      //    t_array+path_idx*workspace_size+x_t_idx_mult_host[path_idx],
      //    sizeof(GT), cudaMemcpyDeviceToHost);
      comp1_gqd2qd(&t_host[path_idx], &t_all_cpu[path_idx]);
   }
   if(path_idx_one == -1)
   {
      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         std::cout << path_idx << " " << t_all_cpu[path_idx];
      }
   }
   else
   {
      std::cout << path_idx_one << " " << t_all_cpu[path_idx_one];
   }
   free(t_all_cpu);
   free(x_t_idx_mult_host);
   free(t_host);
}

void GPUWorkspace::print_delta_t_mult ( int path_idx_one )
{
   CT* delta_t_mult_cpu = new CT[n_path];

   GT* delta_t_mult_host = (GT*)malloc(n_path*sizeof(GT));

   cudaMemcpy(delta_t_mult_host, delta_t_mult,
              n_path*sizeof(GT), cudaMemcpyDeviceToHost);

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      // std::cout << "x_t_idx_mult_host[path_idx] = "
      //           << x_t_idx_mult_host[path_idx] << std::endl;
      comp1_gqd2qd(&delta_t_mult_host[path_idx], &delta_t_mult_cpu[path_idx]);
   }
   if(path_idx_one == -1)
   {
      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         std::cout << path_idx << " " << delta_t_mult_cpu[path_idx];
      }
   }
   else
   {
      std::cout << path_idx_one << " " << delta_t_mult_cpu[path_idx_one];
   }
   free(delta_t_mult_cpu);
   free(delta_t_mult_host);
}

void GPUWorkspace::print_x_last_mult()
{
   CT** x_all_cpu = get_x_last_all();

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      for(int i=0; i<dim; i++)
      {
         std::cout << i << " "  << x_all_cpu[path_idx][i];
      }
      free(x_all_cpu[path_idx]);
   }
   free(x_all_cpu);
}

void GPUWorkspace::print_x_array()
{
   CT* x_array_ct = get_x_array();
   for(int i=0; i<dim; i++)
   {
      for(int j=0; j<n_array; j++)
      {
         std::cout << i << " " << j << " " << x_array_ct[j*dim+i];
      }
      std::cout << std::endl;
   }
   free(x_array_ct);
}

CT* GPUWorkspace::get_t_array()
{
   GT* t_array_host = (GT*)malloc(n_array*sizeof(GT));
   cudaMemcpy(t_array_host, t_array,
              n_array*sizeof(GT), cudaMemcpyDeviceToHost);

   CT* t_array_ct = (CT*)malloc(n_array*sizeof(GT));
   for(int i=0; i<n_array; i++) comp1_gqd2qd(&t_array_host[i], &t_array_ct[i]);

   free(t_array_host);

   return t_array_ct;
}

void GPUWorkspace::print_t_array()
{
   CT* t_array_ct = get_t_array();
   for(int j=0; j<n_array; j++)
   {
      std::cout << j << " " << t_array_ct[j];
   }
   std::cout << std::endl;
   free(t_array_ct);
}


CT* GPUWorkspace::get_f_val()
{
   size_t f_val_size = n_eq*sizeof(GT);
   GT* f_val_host = (GT*)malloc(f_val_size);
   cudaMemcpy(f_val_host, f_val, f_val_size, cudaMemcpyDeviceToHost);

   CT* f_val_cpu = (CT*)malloc(f_val_size);
   for(int i=0; i<dim; i++) comp1_gqd2qd(&f_val_host[i], &f_val_cpu[i]);

   free(f_val_host);

   return f_val_cpu;
}

void GPUWorkspace::print_f_val()
{
   CT* f_val_cpu = get_f_val();
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
   {
      std::cout << eq_idx << " " << f_val_cpu[eq_idx];
   }
   free(f_val_cpu);
}

CT* GPUWorkspace::get_x_last()
{
   size_t x_size = dim*sizeof(GT);
   GT* x_host = (GT*)malloc(x_size);
   cudaMemcpy(x_host, x_last, x_size, cudaMemcpyDeviceToHost);

   CT* x_cpu = (CT*)malloc(x_size);
   for(int i=0; i<dim; i++) comp1_gqd2qd(&x_host[i], &x_cpu[i]);

   free(x_host);
   return x_cpu;
}

CT* GPUWorkspace::get_sol ( int path_idx )
{
   GT* sol_host = (GT*)malloc(dim*sizeof(GT));
   cudaMemcpy(sol_host, sol+path_idx*dim,
              dim*sizeof(GT), cudaMemcpyDeviceToHost);
   CT* sol_cpu = (CT*)malloc(dim*sizeof(GT));
   for(int i=0; i<dim; i++)
   {
      comp1_gqd2qd(&sol_host[i], &sol_cpu[i]);
   }
   free(sol_host);

   return sol_cpu;
}

CT* GPUWorkspace::get_sol_array()
{
   size_t size_sol_mult = n_path*dim*sizeof(GT);
   GT* host_sol_mult = (GT*)malloc(size_sol_mult);
   cudaMemcpy(host_sol_mult, sol, size_sol_mult, cudaMemcpyDeviceToHost);
   CT* gpu_sol_mult = (CT*)malloc(size_sol_mult);
   for(int i=0; i<dim*n_path; i++)
   {
      comp1_gqd2qd(&host_sol_mult[i], &gpu_sol_mult[i]);
   }
   free(host_sol_mult);

   return gpu_sol_mult;
}

CT** GPUWorkspace::get_sol_mult()
{
   CT* sol_cpu_all = get_sol_array();
   CT** sol_cpu = new CT*[n_path];
   sol_cpu[0] = sol_cpu_all;

   for(int path_idx=1; path_idx<n_path; path_idx++)
   {
      sol_cpu[path_idx] = sol_cpu[path_idx-1] + dim;
   }
   return sol_cpu;
}

T1 GPUWorkspace::sol_norm()
{
   GT* sol_host = (GT*)malloc(dim*sizeof(GT));
   cudaMemcpy(sol_host, sol, dim*sizeof(GT), cudaMemcpyDeviceToHost);

   CT tmp_sol;
   comp1_gqd2qd(&sol_host[0], &tmp_sol);
   T1 max_delta = tmp_sol.real*tmp_sol.real+tmp_sol.imag*tmp_sol.imag;

   for(int k=1; k<dim; k++)
   {
      comp1_gqd2qd(&sol_host[0], &tmp_sol);
      T1 tmp_delta = tmp_sol.real*tmp_sol.real+tmp_sol.imag*tmp_sol.imag;
      
      if(tmp_delta>max_delta)
      {
         max_delta = tmp_delta;
      }
   }
   return max_delta;
}

void GPUWorkspace::init_x_t_predict_test()
{
   std::cout
      << "--------- Initialize x and t value for Testing only --------- "
      <<std::endl;
   CT* x_array_cpu = (CT*)malloc(n_predictor*dim*sizeof(CT));
   CT* t_array_cpu = (CT*)malloc(n_predictor*sizeof(CT));

   for(int i=0; i<n_predictor; i++)
   {
      for(int j=0; j<dim; j++)
      {
         x_array_cpu[i*dim + j] = CT(((i+1)*(i+1)*(i+1)+5),0);
         std::cout << i << " " << j << " " << x_array_cpu[i*dim + j];
      }
      t_array_cpu[i] = CT(i+1,0);
      std::cout << i << " " << t_array_cpu[i];
   }
   GT* x_array_host = (GT*)malloc(n_predictor*dim*sizeof(GT));

   for(int i=0; i<dim*n_predictor; i++)
   {
      comp1_qd2gqd(&x_array_cpu[i],&x_array_host[i]);
   }
   GT* t_array_host = (GT*)malloc(n_predictor*sizeof(GT));

   for(int i=0; i<n_predictor; i++)
   {
      comp1_qd2gqd(&t_array_cpu[i],&t_array_host[i]);
   }
   cudaMemcpy(x_array+dim, x_array_host,
              n_predictor*dim*sizeof(GT), cudaMemcpyHostToDevice);
   cudaMemcpy(t_array+1, t_array_host,
              n_predictor*sizeof(GT), cudaMemcpyHostToDevice);

   free(x_array_cpu);
   free(t_array_cpu);
   free(x_array_host);
   free(t_array_host);
}

void GPUInst::init_predict()
{
   predict_BS = 32;
   predict_grid = get_grid(dim,predict_BS,n_path);
}

void GPUInst::init_coef ( const CPUInstHomCoef& cpu_inst_coef )
{
   n_coef = cpu_inst_coef.n_coef;
   size_t coef_size = n_coef*2*sizeof(GT);

   GT* host_coef = (GT *)malloc(coef_size);
   CT* tmp_cpu_coef = cpu_inst_coef.coef_orig;
   for(int i=0; i<n_coef*2; i++)
   {
      comp1_qd2gqd(&tmp_cpu_coef[i],&host_coef[i]);
   }
   cudaMalloc((void **)&coef, coef_size);
   cudaMemcpy(coef, host_coef, coef_size, cudaMemcpyHostToDevice);
   free(host_coef);

   coef_BS = 64;
   coef_grid = get_grid(n_coef, coef_BS, n_path);

   alpha = cpu_inst_coef.alpha;
}

void GPUInst::init_mon ( const CPUInstHomMonBlock& cpu_inst_mon_block )
{
   n_mon_block = cpu_inst_mon_block.n_mon;
   mon_block_grid = get_grid(n_mon_block, BS_Mon_Align, n_path);
   BS_mon_block = cpu_inst_mon_block.BS;
   NB_mon_block = cpu_inst_mon_block.NB;
   // mon_pos_block_size = cpu_inst_mon_block.mon_pos_block_size;
   mon_pos_size = cpu_inst_mon_block.mon_pos_block_size;

   n_mon_single = cpu_inst_mon_block.n_mon_single;
   size_t mon_pos_single_size = 2*n_mon_single*sizeof(unsigned short);
   cudaMalloc((void **)&mon_single_pos_block, mon_pos_single_size);
   cudaMemcpy(mon_single_pos_block, cpu_inst_mon_block.mon_single_pos_block,
              mon_pos_single_size, cudaMemcpyHostToDevice);

   mon_level0_BS = 32;
   mon_global_BS = 64;
   mon_level_BS = shmemsize/2;

   size_t mon_pos_start_block_size = NB_mon_block*sizeof(int);
   cudaMalloc((void **)&mon_pos_start_block, mon_pos_start_block_size);
   cudaMemcpy(mon_pos_start_block, cpu_inst_mon_block.mon_pos_start_block,
              mon_pos_start_block_size, cudaMemcpyHostToDevice);

   size_t mon_pos_block_size
      = cpu_inst_mon_block.mon_pos_block_size*sizeof(unsigned short);
   cudaMalloc((void **)&mon_pos_block, mon_pos_block_size);
   cudaMemcpy(mon_pos_block, cpu_inst_mon_block.mon_pos_block,
              mon_pos_block_size, cudaMemcpyHostToDevice);
}

void GPUInst::init_mon ( const CPUInstHomMon& cpu_inst_mon )
{
   level = cpu_inst_mon.level;
   n_mon_level = cpu_inst_mon.n_mon_level;
   n_mon = cpu_inst_mon.n_mon;

   // std::cout << "n_mon = " << n_mon << " level = " << level << std::endl;
   // for(int i=0; i<level; i++)
   // {
   //    std::cout << "level " << i << " " << n_mon_level[i] << std::endl;
   // }

   size_t mon_pos_start_size = n_mon*sizeof(int);
   cudaMalloc((void **)&mon_pos_start, mon_pos_start_size);
   cudaMemcpy(mon_pos_start, cpu_inst_mon.mon_pos_start,
              mon_pos_start_size, cudaMemcpyHostToDevice);

   mon_pos_size = cpu_inst_mon.mon_pos_size;
   cudaMalloc((void **)&mon_pos, mon_pos_size*sizeof(unsigned short));
   cudaMemcpy(mon_pos, cpu_inst_mon.mon_pos,
              mon_pos_size*sizeof(unsigned short), cudaMemcpyHostToDevice);

   mon_level_grid = new dim3[level];
   // rest monomial after this level use global method
   n_mon_level_rest = new int[level];
   mon_level_grid_rest = new dim3[level];

   mon_level0_BS = 32;
   mon_global_BS = 64;
   mon_level_grid[0] = get_grid(n_mon_level[0],mon_level0_BS,n_path);
   n_mon_level_rest[0] = n_mon -n_mon_level[0];
   mon_level_grid_rest[0] = get_grid(n_mon_level_rest[0],mon_global_BS,n_path);

   mon_level_grid[1] = get_grid(n_mon_level[1],mon_global_BS,n_path);
   n_mon_level_rest[1] = n_mon_level_rest[0] -n_mon_level[1];
   mon_level_grid_rest[1] = get_grid(n_mon_level_rest[1],mon_global_BS,n_path);
   // for level 1
   // for dd cyclic 96
   // 64 is the best
   // 32 128 has similar result
   // 16 is bad

   mon_level_BS = shmemsize/2;
   int n_thread_per_job = 2;
   for(int i=2; i<level; i++)
   {
      mon_level_grid[i]
         = get_grid(n_mon_level[i],mon_level_BS,n_path, n_thread_per_job);
      n_mon_level_rest[i] = n_mon_level_rest[i-1] -n_mon_level[i];
      mon_level_grid_rest[i]
         = get_grid(n_mon_level_rest[i],mon_global_BS,n_path, n_thread_per_job);
      n_thread_per_job *= 2;
   }

   n_mon_global = n_mon - n_mon_level[0];
   mon_global_grid = get_grid(n_mon_global,mon_global_BS,n_path);
}

void GPUInst::init_base ( const CPUInstHomMon& cpu_inst_mon )
{
   base_table_size = 0;
   int* base_table_start_host = new int[dim];

   std::cout << "max_deg_base " << std::endl;
   for(int var_idx=0; var_idx<dim; var_idx++)
   {
      std::cout << var_idx << " " << cpu_inst_mon.max_deg_base[var_idx]
                << std::endl;
      base_table_start_host[var_idx] = base_table_size;
      base_table_size += cpu_inst_mon.max_deg_base[var_idx];
   }
   std::cout << "base_table_size = " << base_table_size << std::endl;
   for(int var_idx=0; var_idx<dim; var_idx++)
   {
      std::cout << var_idx << " " << base_table_start_host[var_idx]
                << std::endl;
   }
   cudaMalloc((void **)&max_deg_base, dim*sizeof(int));
   cudaMemcpy(max_deg_base, cpu_inst_mon.max_deg_base,
              dim*sizeof(int), cudaMemcpyHostToDevice);

   cudaMalloc((void **)&base_table_start, dim*sizeof(int));
   cudaMemcpy(base_table_start, base_table_start_host,
              dim*sizeof(int), cudaMemcpyHostToDevice);

   cudaMalloc((void **)&mon_exp, mon_pos_size*sizeof(unsigned short));
   cudaMemcpy(mon_exp, cpu_inst_mon.mon_exp,
              mon_pos_size*sizeof(unsigned short), cudaMemcpyHostToDevice);

   n_mon_base_start = cpu_inst_mon.n_mon_base_start;
   n_mon_base = n_mon - n_mon_base_start;
   std::cout << "n_mon = " << n_mon << std::endl;
   std::cout << "n_mon_base = " << n_mon_base << std::endl;
   std::cout << "n_mon_base_start = " << n_mon_base_start << std::endl;
}

void GPUInst::init_sum
 ( const CPUInstHomSumBlock& cpu_inst_sum,
   const CPUInstHomSum& cpu_inst_sum_orig )
{
   n_sum = cpu_inst_sum.n_sum;
   n_sum_levels = cpu_inst_sum.n_sum_levels;
   n_sum_level= cpu_inst_sum.n_sum_level;
   n_sum_level_rest = cpu_inst_sum.n_sum_level_rest;

   sum_BS = 64;
   sum_grid = get_grid(n_sum,sum_BS,n_path);

   sum_level_grid = new dim3[n_sum_levels];
   sum_level_grid_rest = new dim3[n_sum_levels];

   int n_thread_per_job = 1;

   sum_level_grid[0]
      = get_grid(n_sum_level[0], sum_BS, n_path, n_thread_per_job);
   sum_level_grid_rest[0]
      = get_grid(n_sum, sum_BS, n_path, n_thread_per_job);
   n_thread_per_job *= 2;

   for(int i=1; i<n_sum_levels; i++)
   {
      sum_level_grid[i]
         = get_grid(n_sum_level[i], sum_BS, n_path, n_thread_per_job);
      sum_level_grid_rest[i]
         = get_grid(n_sum_level_rest[i-1], sum_BS, n_path, n_thread_per_job);
      n_thread_per_job *= 2;
   }

   size_t sum_pos_start_size = n_sum*sizeof(int);
   cudaMalloc((void **)&sum_pos_start, sum_pos_start_size);
   cudaMemcpy(sum_pos_start, cpu_inst_sum.sum_pos_start,
              sum_pos_start_size, cudaMemcpyHostToDevice);

   size_t sum_pos_size = cpu_inst_sum.sum_pos_size*sizeof(int);
   cudaMalloc((void **)&sum_pos, sum_pos_size);
   cudaMemcpy(sum_pos, cpu_inst_sum.sum_pos,
              sum_pos_size, cudaMemcpyHostToDevice);

   n_sum_zero = 0;
   sum_zeros = NULL;
   if(cpu_inst_sum_orig.n_sum_zero>0)
   {
      n_sum_zero = cpu_inst_sum_orig.n_sum_zero;
      cudaMalloc((void **)&sum_zeros, n_sum_zero*sizeof(int));
      cudaMemcpy(sum_zeros, cpu_inst_sum_orig.sum_zeros,
                 n_sum_zero*sizeof(int), cudaMemcpyHostToDevice);
   }
}

void GPUInst::init_sum ( const CPUInstHomSum& cpu_inst_sum )
{
   n_sum = cpu_inst_sum.n_sum;
   n_sum_levels = cpu_inst_sum.n_sum_levels;
   n_sum_level= cpu_inst_sum.n_sum_level;
   n_sum_level_rest = cpu_inst_sum.n_sum_level_rest;
   sum_BS = 32;
   sum_grid = get_grid(n_sum,sum_BS,n_path);

   sum_level_grid = new dim3[n_sum_levels];
   sum_level_grid_rest = new dim3[n_sum_levels];

   int n_thread_per_job = 1;

   sum_level_grid[0]
      = get_grid(n_sum_level[0], sum_BS, n_path, n_thread_per_job);
   sum_level_grid_rest[0]
      = get_grid(n_sum, sum_BS, n_path, n_thread_per_job);
   n_thread_per_job *= 2;

   for(int i=1; i<n_sum_levels; i++)
   {
      sum_level_grid[i]
         = get_grid(n_sum_level[i], sum_BS, n_path, n_thread_per_job);
      sum_level_grid_rest[i]
         = get_grid(n_sum_level_rest[i-1], sum_BS, n_path, n_thread_per_job);
      n_thread_per_job *= 2;
   }

   size_t sum_pos_start_size = n_sum*sizeof(int);
   cudaMalloc((void **)&sum_pos_start, sum_pos_start_size);
   cudaMemcpy(sum_pos_start, cpu_inst_sum.sum_pos_start,
              sum_pos_start_size, cudaMemcpyHostToDevice);

   size_t sum_pos_size = cpu_inst_sum.sum_pos_size*sizeof(int);
   cudaMalloc((void **)&sum_pos, sum_pos_size);
   cudaMemcpy(sum_pos, cpu_inst_sum.sum_pos,
              sum_pos_size, cudaMemcpyHostToDevice);

   n_sum_zero = 0;
   sum_zeros = NULL;
   // std::cout << "cpu_inst_sum.n_sum_zero::::"
   //           << cpu_inst_sum.n_sum_zero << std::endl;
   if(cpu_inst_sum.n_sum_zero>0)
   {
      n_sum_zero = cpu_inst_sum.n_sum_zero;
      cudaMalloc((void **)&sum_zeros, n_sum_zero*sizeof(int));
      cudaMemcpy(sum_zeros, cpu_inst_sum.sum_zeros,
                 n_sum_zero*sizeof(int), cudaMemcpyHostToDevice);
   }
}

void GPUInst::init_workspace( const CPUInstHom& cpu_inst )
{
   // std::cout << "n_workspace1 = " << n_workspace1 << std::endl;
   // std::cout << "cpu_inst.CPU_inst_hom_block.mon_pos_block_size = "
   //           << cpu_inst.CPU_inst_hom_block.mon_pos_block_size << std::endl;
   // n_workspace = n_coef 
   //    + cpu_inst.CPU_inst_hom_block.mon_pos_block_size + n_mon_level[0]*2;
   // n_workspace = n_coef + cpu_inst.CPU_inst_hom_mon.mon_pos_size;
   n_constant = cpu_inst.n_constant;
};

void GPUInst::init_eq ( const CPUInstHomEq& cpu_inst_eq )
{
   cudaMalloc((void **)&eq_pos_start, n_eq*sizeof(int));
   cudaMemcpy(eq_pos_start, cpu_inst_eq.eq_pos_start,
              n_eq*sizeof(int), cudaMemcpyHostToDevice);

   n_mon_total_eq = cpu_inst_eq.n_mon_total;
   cudaMalloc((void **)&mon_pos_start_eq, n_mon_total_eq*sizeof(int));
   cudaMemcpy(mon_pos_start_eq, cpu_inst_eq.mon_pos_start_eq,
              n_mon_total_eq*sizeof(int), cudaMemcpyHostToDevice);

   GT* coef_eq_host = (GT*)malloc(2*n_mon_total_eq*sizeof(GT));
   for(int coef_idx=0; coef_idx<2*n_mon_total_eq; coef_idx++)
   {
      // std::cout << coef_idx << " " << cpu_inst_eq.coef[coef_idx];
      comp1_qd2gqd(&cpu_inst_eq.coef[coef_idx],&coef_eq_host[coef_idx]);
   }
   cudaMalloc((void **)&coef_eq, 2*n_mon_total_eq*sizeof(GT));
   cudaMemcpy(coef_eq, coef_eq_host, 2*n_mon_total_eq*sizeof(GT),
              cudaMemcpyHostToDevice);
   free(coef_eq_host);

   n_pos_total_eq = cpu_inst_eq.n_pos_total;
   cudaMalloc((void **)&mon_pos_eq, n_pos_total_eq*sizeof(unsigned short));
   cudaMemcpy(mon_pos_eq, cpu_inst_eq.mon_pos_eq,
              n_pos_total_eq*sizeof(unsigned short), cudaMemcpyHostToDevice);
};

#endif /*__GPU_DATA_CU__*/
