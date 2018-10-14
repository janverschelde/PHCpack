// The file workspace.h defines the class Workspace,
// to store the data during the path tracking.

#ifndef __WORKSPACE_HOST_H__
#define __WORKSPACE_HOST_H__

template <class ComplexType>
class Workspace
{
   public:

      ComplexType* all;  // combination of coef and mon
      ComplexType* coef; // coef evaluation, start from workspace
      ComplexType* mon;  // monomial evaluation
      ComplexType* sum;  // part of coef and mon

      ComplexType* matrix;

      ComplexType* sol;  // solution for modified Gram-Schmidt (mgs)
      ComplexType* rhs;  // right hand size for mgs
      ComplexType** R;   // upper triangular form computed by mgs
      ComplexType** V;   // basis vectors in mgs

      int dim;           // ambient dimension equals number of variables
      int n_eq;          // number of equations
      int n_array;       // number of past solutions stored
      int n_predictor;   // number of points in predictor

      int x_t_idx_last;
      int x_t_idx;

      ComplexType** x_array;
      ComplexType* x;
      ComplexType* x_last;

      ComplexType* t_array;
      ComplexType* t;
      ComplexType* t_last;

      ComplexType** deg_table; // powers of coordinates

      ComplexType* div_diff4pred; // divided differences, size n_predictor
      ComplexType* t_array4pred;  // auxiliary array for predictor
      ComplexType* t_diff4pred;   // auxiliary array for predictor

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
         div_diff4pred = NULL;
         t_array4pred = NULL;
         t_diff4pred = NULL;
      }

      void init ( int workspace_size, int n_coef, int n_constant,
                  int n_eq, int dim, int n_predictor, int* max_deg_base=NULL,
                  int verbose=0 );
      /*
       * Initializes the workspace.
       *
       * ON ENTRY :
       *   workspace_size  is the size of the workspace;
       *   n_coef          number of coefficients;
       *   n_constant      number of constants;
       *   n_eq            number of equations in the system;
       *   dim             number of variables in the system;
       *   n_predictor     number of predictor points;
       *   max_deg_base    largest degrees of the common factor, equals
       *                   NULL if all monomials are products of variables;
       *   verbose         verbose flag, 0 for no output.
       */

      Workspace ( int workspace_size, int coef_size,
                  int workspace_constant_size, int n_eq, int dim,
                  int n_predictor, int* max_deg_base=NULL, int verbose=0 )
      /*
       * Wrapper to the initialization method of the workspace.
       */
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
         if(rhs != NULL) delete[] rhs;
         if(deg_table != NULL)
         {
            delete[] deg_table[0];
            delete[] deg_table;
         }
         if(div_diff4pred != NULL) delete[] div_diff4pred;
         if(t_array4pred != NULL) delete[] t_array4pred;
         if(t_diff4pred != NULL) delete[] t_diff4pred;
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
            for(int j=0; j<dim+1; j++)
               std::cout << i << " " << j << " " << R[i][j];
      }

      void copy_x_last ( ComplexType* result )
      {
         for(int i=0; i<dim; i++) result[i] = x_last[i];
      }
};

#include "workspace_host.tpp"

#endif /* __WORKSPACE_CPU_H__ */
