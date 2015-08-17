/* path_data.h created on Feb 20, 2015 by yxc with edits by jv */

#ifndef PATH_DATA_H_
#define PATH_DATA_H_

#include <vector>
#include "utilities.h"
#include "DefineType_Host.h"

using namespace std;

struct correct_iteration
{
   double correct_a;
   double correct_r;
   double residual_a;
   double residual_r;
};

class PathStep
{
   public:

      int dim;
      bool success;
      CT t;
      T1 delta_t;
      CT* predict_pt;
      CT* correct_pt;
      vector<correct_iteration> correct_it;
      int n_it;

      PathStep()
      {
         dim = 0;
         success = false;
         t = CT(0.0,0.0);
         delta_t = 0.0;
         predict_pt = NULL;
         correct_pt = NULL;
         n_it =0;
      }

      PathStep ( int dim )
      {
         this->dim = dim;
         success = false;
         t = CT(0.0,0.0);
         delta_t = 0.0;
         predict_pt = NULL;
         correct_pt = NULL;
         n_it =0;
      }

      ~PathStep()
      {
         delete[] predict_pt;
         delete[] correct_pt;
      }

      PathStep ( int dim, ifstream& path_file )
      {
         this->dim = dim;
         t = CT(0.0,0.0);
         delta_t = 0.0;
         n_it =0;
         success = false;
         read_phc_file(path_file);
      }

      void read_phc_file(ifstream& path_file);

      void read_predict_pt(ifstream& path_file);

      void read_correct_pt(ifstream& path_file);

      void read_correct_it(ifstream& path_file);

      void update_predict_pt(CT* predict_pt);

      void update_correct_pt(CT* correct_pt);

      void update_t(CT delta_t, CT t);

      void add_iteration(double max_delta_x, double r_max_delta_x);

      void update_iteration_res(double residual_a, double residual_r);

      void print();

      void print_it();
};

class Path
{
   public:

      int dim;
      int n_step;
      bool success;
      CT* start_pt;
      CT* end_pt;
      vector<PathStep*> steps;

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
         for(vector<PathStep*>::iterator it = steps.begin(); 
             it!=steps.end(); ++it) 
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

      void update_step_predict_pt(CT* predict_pt);

      void update_step_correct_pt(CT* correct_pt);

      void update_step_t(CT delta_t, CT t);

      void add_iteration(double max_delta_x, double r_max_delta_x);

      void update_iteration_res(double residual_a, double residual_r);

      void print_t();

      void print();

      void print_phc();

      void compare(Path& that);

      void add_start_pt(CT* start_pt);

      void add_end_pt(CT* end_pt);

      void update_success(bool success);

      void clear();
};

#endif /* PATH_DATA_H_ */
