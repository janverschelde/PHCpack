/* The path_data.h contains the definition of the class Path,
 * to store all the solutions on a path. */

#ifndef __PATH_DATA_H__
#define __PATH_DATA_H__

#include <vector>
#include "utilities.h"
#include "path_step.h"

using namespace std;

template <class ComplexType, class RealType>
class Path
{
   public:

      int dim;
      int n_step;
      bool success;
      ComplexType* start_pt;
      ComplexType* end_pt;
      vector<PathStep<ComplexType, RealType>*> steps;

      Path()
      {
         dim = 0;
         n_step = 0;
         success = false;
         start_pt = NULL;
         end_pt = NULL;
      }

      ~Path()
      {
         for(typename vector<PathStep<ComplexType, RealType>*>::iterator
             it = steps.begin(); 
             it != steps.end(); ++it) 
         {
            delete (*it);
         }
         if(start_pt != NULL)
         {
            delete[] start_pt;
            start_pt = NULL;
         }
         if(end_pt != NULL)
         {
            delete[] end_pt;
            end_pt = NULL;
         }
      }

      Path(int dim, ifstream& path_file)
      {
         this->dim = dim;
         n_step = 0;
         success = false;
         start_pt = NULL;
         end_pt = NULL;
         read_phc_file(path_file);
      }

      void read_phc_file(ifstream& path_file);

      void add_step(ifstream& path_file);

      void add_step_empty();

      void update_step_predict_pt(ComplexType* predict_pt);

      void update_step_correct_pt(ComplexType* correct_pt);

      void update_step_t(ComplexType delta_t, ComplexType t);

      void add_iteration(double max_delta_x, double r_max_delta_x);

      void update_iteration_res(double residual_a, double residual_r);

      void print_t();

      void print();

      void print_phc();

      void compare(Path& that);

      void add_start_pt(ComplexType* start_pt);

      void add_end_pt(ComplexType* end_pt);

      void update_success(bool success);

      void clear();
};

#include "path_data.tpp"

#endif /* __PATH_DATA_H__ */
