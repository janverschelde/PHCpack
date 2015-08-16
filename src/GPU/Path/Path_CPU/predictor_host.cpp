/* predictor.cpp, created on Dec 6, 2014 by yxc with edits by jv */

#include "predictor_host.h"

void predictor_newton
 ( CT** x_array, CT* t_array_orig, int x_t_idx, int np_predictor, int dim )
/*
    Newton predictor
    Reference: pre_newton1.m
*/
{
   CT* div_diff = new CT[np_predictor];
   CT* t_array = new CT[np_predictor];
   CT* t_diff = new CT[np_predictor];
   int k = 0;
   for(int np_idx = x_t_idx+1; np_idx < np_predictor+1; np_idx++)
   {
      t_array[k++] = t_array_orig[np_idx];
   }
   for(int np_idx = 0; np_idx < x_t_idx; np_idx++)
   {
      t_array[k++] = t_array_orig[np_idx];
   }
    /*for(int i=0; i<np_predictor; i++){
    	std::cout << "t" << i << " " << t_array[i];
    }*/
   // std::cout << "Predict Result" << std::endl;
   CT t_new = t_array_orig[x_t_idx];
   for(int i=0; i<np_predictor; i++)
   {
      t_diff[i] = t_new - t_array[i];
   }
   for(int dim_idx=0; dim_idx < dim; dim_idx++)
   {
      // Copy initial X value to divide difference
      k=0;
      for(int np_idx = x_t_idx+1; np_idx < np_predictor+1; np_idx++)
      {
         div_diff[k] = x_array[np_idx][dim_idx];
         k++;
      }
      for(int np_idx = 0; np_idx < x_t_idx; np_idx++)
      {
         div_diff[k] = x_array[np_idx][dim_idx];
         k++;
      }

      /*for(int i=0; i<np_predictor; i++)
        {
           std::cout << i << " " << div_diff[i];
        }*/
      // compute divided differences
      for(int i = 1; i < np_predictor; i++)
      {
         for(int j = np_predictor-1; j >= i; j--)
         {
            div_diff[j] = (div_diff[j] - div_diff[j-1])
                         /(t_array[j]-t_array[j-i]); // need to be computed
         }
      }

      /*for(int i = 0; i < np_predictor; i++)
        {
           std::cout << div_diff[i] << " ";
        }
        std::cout << std::endl;*/

      // Compute predict point
      CT x_tmp = 0.0;
      for(int i=np_predictor-1; i > 0; i--)
      {
         x_tmp = (x_tmp + div_diff[i]) * t_diff[i-1];
      }
      // Put X back
      x_array[x_t_idx][dim_idx] = x_tmp + div_diff[0];
      // std::cout << dim_idx << " " << x_array[x_t_idx][dim_idx];
   }
}
