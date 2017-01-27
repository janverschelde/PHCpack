/* path_data.h created on Feb 20, 2015 by yxc with edits by jv */

#ifndef PATH_DATA_H_
#define PATH_DATA_H_

#include <vector>
#include "utilities.h"

using namespace std;

struct correct_iteration
{
   double correct_a;
   double correct_r;
   double residual_a;
   double residual_r;
};

template <class ComplexType, class RealType>
class PathStep
{
   public:

      int dim;
      bool success;
      ComplexType t;
      RealType delta_t;
      ComplexType* predict_pt;
      ComplexType* correct_pt;
      vector<correct_iteration> correct_it;
      int n_it;

      PathStep()
      {
         dim = 0;
         success = false;
         t = ComplexType(0.0,0.0);
         delta_t = 0.0;
         predict_pt = NULL;
         correct_pt = NULL;
         n_it =0;
      }

      PathStep ( int dim )
      {
         this->dim = dim;
         success = false;
         t = ComplexType(0.0,0.0);
         delta_t = 0.0;
         predict_pt = NULL;
         correct_pt = NULL;
         n_it =0;
      }

      ~PathStep()
      {
         if(predict_pt != NULL) delete[] predict_pt;
         if(correct_pt != NULL) delete[] correct_pt;
      }

      PathStep ( int dim, ifstream& path_file )
      {
         this->dim = dim;
         t = ComplexType(0.0,0.0);
         delta_t = 0.0;
         n_it =0;
         success = false;
         read_phc_file(path_file);
      }

      void read_phc_file(ifstream& path_file);

      void read_predict_pt(ifstream& path_file);

      void read_correct_pt(ifstream& path_file);

      void read_correct_it(ifstream& path_file);

      void update_predict_pt(ComplexType* predict_pt);

      void update_correct_pt(ComplexType* correct_pt);

      void update_t(ComplexType delta_t, ComplexType t);

      void add_iteration(double max_delta_x, double r_max_delta_x);

      void update_iteration_res(double residual_a, double residual_r);

      void print();

      void print_it();
};

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

#endif /* PATH_DATA_H_ */
