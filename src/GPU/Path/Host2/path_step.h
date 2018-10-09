/* The file path_step.h contains the definition of the class PathStep,
 * to store the data of one step along a solution path. */

#ifndef __PATH_STEP_H__
#define __PATH_STEP_H__

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

#include "path_step.tpp"

#endif /* __PATH_STEP_H__ */
