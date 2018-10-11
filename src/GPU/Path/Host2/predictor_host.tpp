// The file predictor_host.tpp contains the definitions of the functions
// with prototypes in the file predictor_host.h.

template <class ComplexType>
void predictor_divdif
 ( ComplexType** x_array, ComplexType* t_array_orig,
   int x_t_idx, int np_predictor, int dim,
   ComplexType* div_diff, ComplexType* t_array, ComplexType* t_diff )
{
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
   ComplexType t_new = t_array_orig[x_t_idx];
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
      ComplexType x_tmp = 0.0;
      for(int i=np_predictor-1; i > 0; i--)
      {
         x_tmp = (x_tmp + div_diff[i]) * t_diff[i-1];
      }
      // Put X back
      x_array[x_t_idx][dim_idx] = x_tmp + div_diff[0];
      // std::cout << dim_idx << " " << x_array[x_t_idx][dim_idx];
   }
}

template <class ComplexType>
void predictor_newton
 ( ComplexType** x_array, ComplexType* t_array_orig,
   int x_t_idx, int np_predictor, int dim )
{
   ComplexType* div_diff = new ComplexType[np_predictor];
   ComplexType* t_array = new ComplexType[np_predictor];
   ComplexType* t_diff = new ComplexType[np_predictor];

   predictor_divdif
     (x_array,t_array_orig,x_t_idx,np_predictor,dim,
      div_diff,t_array,t_diff);

   delete[] div_diff;
   delete[] t_array;
   delete[] t_diff;
}
