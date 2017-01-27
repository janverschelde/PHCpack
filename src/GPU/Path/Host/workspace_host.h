// The workspace stores the data during the path tracking.

#ifndef WORKSPACE_HOST_H_
#define WORKSPACE_HOST_H_

template <class ComplexType>
class Workspace
{
   public:

      ComplexType* all;  // combination of coef and mon
      ComplexType* coef; // coef evaluation, start from workspace
      ComplexType* mon;  // monomial evaluation
      ComplexType* sum;  // part of coef and mon

      ComplexType* matrix;

      ComplexType* sol;  // for modified Gram-Schmidt
      ComplexType** R;
      ComplexType** V;

      int dim;           // array for predictor
      int n_eq;
      int n_array;
      int n_predictor;
      int x_t_idx_last;
      int x_t_idx;

      ComplexType** x_array;
      ComplexType* x;
      ComplexType* x_last;

      ComplexType* t_array;
      ComplexType* t;
      ComplexType* t_last;

      ComplexType** deg_table;

      int path_idx;

      Workspace()
      {
         all = NULL;
         coef = NULL;
         mon = NULL;
         sum = NULL;
         matrix = NULL;
         sol = NULL;
         R = NULL;
         V = NULL;
         x_array = NULL;
         x_last = NULL;
         x = NULL;
         t_array = NULL;
         t_last = NULL;
         t = NULL;
         dim = 0;
         n_eq = 0;
         n_array = 0;
         n_predictor = 0;
         x_t_idx = 0;
         x_t_idx_last = 0;
         path_idx = 0;
         deg_table = NULL;
      }

      void init ( int workspace_size, int n_coef, int n_constant,
                  int n_eq, int dim, int n_predictor, int* max_deg_base=NULL );

      Workspace ( int workspace_size, int coef_size,
                  int workspace_constant_size, int n_eq, int dim,
                  int n_predictor, int* max_deg_base=NULL )
      {
         init(workspace_size, coef_size, workspace_constant_size, n_eq,
              dim, n_predictor, max_deg_base);
      }

      ~Workspace()
      {
         if(all != NULL) delete[] all;
         if(matrix != NULL) delete[] matrix;
         if(V != NULL) delete[] V;
         if(R != NULL) delete[] R;
         if(sol != NULL) delete[] sol;
         if(deg_table != NULL)
         {
            delete[] deg_table[0];
            delete[] deg_table;
         }
         // more to add ... careful with memory leaks ...
         // delete[] tmp_x;
      }

      void init_x_t ( int dim, int n_predictor );

      void init_deg_table ( int* max_deg_base );

      void init_x_t_idx();

      void init_x_t_predict_test();

      void update_x_t ( ComplexType* cpu_sol0, ComplexType cpu_t );

      void update_x_t_value ( ComplexType* cpu_sol0, ComplexType cpu_t );

      void update_x_value ( ComplexType* cpu_sol0 );

      void update_t_value ( ComplexType cpu_t );

      void update_x_t_idx();

      void print_t();

      void print_x();

      void print_x_last();

      void print_x_array();

      void print_coef();

      void print_result();

      void print_matrix_r()
      {
         for(int i=0; i<dim+1; i++)
         {
            for(int j=0; j<dim+1; j++)
            {
               std::cout << i << " " << j << " " << R[i][j];
            }
         }
      }

      void copy_x_last ( ComplexType* result )
      {
         for(int i=0; i<dim; i++)
         {
            result[i] = x_last[i];
         }
      }
};

#include "workspace_host.tpp"

#endif /* WORKSPACE_CPU_H_ */
