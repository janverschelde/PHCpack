// The file cpuinsthomsumblock.tpp provides the definitions of the methods in
// the class CPUInstHomSumBlock, with prototypes in cpuinsthomsumblock.h.

template <class ComplexType>
void CPUInstHomSumBlock<ComplexType>::init
 ( MonSet<ComplexType>* hom_monset, int n_monset,
   const int* mon_pos_start, int dim, int n_eq, int n_constant,
   int n_mon_single, int* mon_pos_start_block, int verbose )
{
   if(verbose > 0)
   {
      std::cout << "dim = " << dim << " n_eq = " << n_eq << std::endl;
   }
   // Step 1: count number of terms to sum in Jacobian matrix
   int* n_sums_loc = new int[n_eq*(dim+1)];
   for(int i=0; i<n_eq*(dim+1); i++)
   {
      n_sums_loc[i] = 0;
   }
   int** n_sums = new int*[n_eq];
   int* n_sums_tmp = n_sums_loc;
   for(int i=0; i<n_eq; i++)
   {
      n_sums[i] = n_sums_tmp;
      n_sums_tmp += dim+1;
   }
   MonSet<ComplexType>* tmp_hom_monset = hom_monset;
   for(int set_idx=0; set_idx<n_monset; set_idx++)
   {
      for(int mon_idx=0; mon_idx<tmp_hom_monset->get_n_mon(); mon_idx++)
      {
         int tmp_eq_idx = tmp_hom_monset->get_eq_idx(mon_idx);
         n_sums[tmp_eq_idx][dim] += 1;
         for(int k=0; k<tmp_hom_monset->get_n(); k++)
         {
            n_sums[tmp_eq_idx][tmp_hom_monset->get_pos(k)] += 1;
         }
      }
      tmp_hom_monset++;
   }
   // Step 2: Count number of sums for certain number of terms
   //         total number of terms to sum
   //         max number of terms to sum
   int max_n_sums = 0;
   for(int i=0; i<n_eq; i++)
   {
      for(int j=0; j<dim+1; j++)
      {
         if(n_sums[i][j] > max_n_sums)
         {
            max_n_sums = n_sums[i][j];
         }
      }
   }
   int* n_sums_count = new int[max_n_sums+1];
   for(int i=0; i<max_n_sums+1; i++)
   {
      n_sums_count[i] = 0;
   }
   for(int i=0; i<n_eq; i++)
   {
      for(int j=0; j<dim+1; j++)
      {
         n_sums_count[n_sums[i][j]]++;
      }
   }
   // Step 3: Sum level
   n_sum = (dim+1)*n_eq - n_sums_count[0];
   n_sum_levels = log2ceil(max_n_sums);
   n_sum_level = new int[n_sum_levels];
   n_sum_level_rest = new int[n_sum_levels];

   for(int i=0; i<n_sum_levels; i++)
   {
      n_sum_level[i] = 0;
      n_sum_level_rest[i] = 0;
   }
   n_sum_level[0] = n_sums_count[1];
   int tmp_level_size = 4;
   int tmp_level = 1;

   for(int i=2; i<max_n_sums+1; i++)
   {
      if(tmp_level_size < i)
      {
         tmp_level_size *= 2;
         tmp_level++;
      }
      n_sum_level[tmp_level] += n_sums_count[i];
   }
   if(verbose > 0)
   {
      std::cout << "n_sum = " << n_sum << std::endl;
   }
   n_sum_level_rest[0] = n_sum - n_sum_level[0];
   if(verbose > 0)
   {
      std::cout << 0 << " " << n_sum_level[0]
                << " " << n_sum_level_rest[0] << std::endl;
   }
   for(int i=1; i<n_sum_levels; i++)
   {
      n_sum_level_rest[i] = n_sum_level_rest[i-1] - n_sum_level[i];
      if(verbose > 0)
      {
         std::cout << i << " " << n_sum_level[i]
                   << " " << n_sum_level_rest[i] << std::endl;
      }
   }
   // Step 4: sum start
   sum_pos_start = new int[n_sum];
   int tmp_idx = 0;
   int last_length = 0;
   for(int i=1; i<max_n_sums+1; i++)
   {
      for(int j=0; j<n_sums_count[i]; j++)
      {
         if(tmp_idx == 0)
         {
            sum_pos_start[0] = 0;
         }
         else
         {
            sum_pos_start[tmp_idx] = sum_pos_start[tmp_idx-1] + last_length;
         }
         tmp_idx++;
         last_length = i+2;
      }
   }
   // Step 5: Start pos of sums
   int* n_sums_start = new int[max_n_sums+1];
   n_sums_start[0] = 0;
   n_sums_start[1] = 0;
   for(int i=2; i<max_n_sums+1; i++)
   {
      n_sums_start[i] = n_sums_start[i-1] + n_sums_count[i-1]*(1+i);
   }
   sum_pos_size
    = n_sums_start[max_n_sums] + n_sums_count[max_n_sums]*(2+max_n_sums);

   if(verbose > 0)
   {
      std::cout << "sum_pos_size = " << sum_pos_size << std::endl;
   }
   int* sum_pos_start_loc = new int[n_eq*(dim+1)];
   for(int i=0; i<n_eq*(dim+1); i++)
   {
      sum_pos_start_loc[i] = 0;
   }
   int** sum_pos_start_matrix = new int*[n_eq];
   int* sum_pos_start_matrix_tmp = sum_pos_start_loc;
   for(int i=0; i<n_eq; i++)
   {
      sum_pos_start_matrix[i] = sum_pos_start_matrix_tmp;
      sum_pos_start_matrix_tmp += dim+1;
   }
   sum_pos = new int[sum_pos_size];
   for(int i=0; i<sum_pos_size; i++)
   {
      sum_pos[i] = 0;
   }
   for(int i=0; i<n_eq; i++)
   {
      for(int j=0; j<dim+1; j++)
      {
         int tmp_n = n_sums[i][j];
         if(tmp_n > 0)
         {
            int tmp_start = n_sums_start[tmp_n];
            sum_pos[tmp_start] = tmp_n;
            sum_pos_start_matrix[i][j] = tmp_start+1;
            sum_pos[tmp_start+tmp_n+1] = j*n_eq + i;
            n_sums_start[tmp_n] += tmp_n+2;
         }
      }
   }
   tmp_hom_monset = hom_monset;
   for(int i=0; i<tmp_hom_monset->get_n_mon(); i++)
   {
      int tmp_eq_idx = tmp_hom_monset->get_eq_idx(i);
      sum_pos[sum_pos_start_matrix[tmp_eq_idx][dim]] = i;
      sum_pos_start_matrix[tmp_eq_idx][dim]++;
   }
   tmp_hom_monset = hom_monset+1;
   int mon_idx = 0;
   for(int i=1; i<n_monset; i++)
   {
      for(int j=0; j<tmp_hom_monset->get_n_mon(); j++)
      {
         int n_var = tmp_hom_monset->get_n();
         int tmp_pos;
         int bidx;
         int tidx;
         int tmp_start_block;
         int tmp_mon_idx;
         if(n_var < 2)
         {
            tmp_pos = mon_pos_start[mon_idx]+n_constant;
         }
         else
         {
            tmp_mon_idx = mon_idx - n_mon_single;
            bidx = tmp_mon_idx/warp_size;
            tidx = tmp_mon_idx - bidx*warp_size;
            tmp_start_block = mon_pos_start_block[bidx];
            tmp_start_block = tmp_start_block + n_constant;
            tmp_pos = tmp_start_block + tidx;
         }
         int tmp_eq_idx = tmp_hom_monset->get_eq_idx(j);
         // Value
         sum_pos[sum_pos_start_matrix[tmp_eq_idx][dim]] = tmp_pos;
         tmp_pos++;
         sum_pos_start_matrix[tmp_eq_idx][dim]++;
         n_sums[tmp_eq_idx][dim] += 1;
         // Derivative
         for(int k=0; k<n_var; k++)
         {
            if(n_var < 2)
            {
               sum_pos[sum_pos_start_matrix[tmp_eq_idx]
                      [tmp_hom_monset->get_pos(k)]] = tmp_pos;
               tmp_pos++;
            }
            else
            {
               sum_pos[sum_pos_start_matrix[tmp_eq_idx]
                  [tmp_hom_monset->get_pos(k)]] = \
                  tmp_start_block+(k+1)*warp_size+tidx;
            }
            sum_pos_start_matrix[tmp_eq_idx][tmp_hom_monset->get_pos(k)]++;
         }
         mon_idx++;
      }
      tmp_hom_monset++;
   }
   delete[] n_sums;
   delete[] n_sums_loc;
   delete[] n_sums_count;
   delete[] n_sums_start;
   delete[] sum_pos_start_loc;
   delete[] sum_pos_start_matrix;
}

template <class ComplexType>
void CPUInstHomSumBlock<ComplexType>::eval
 ( ComplexType* sum, ComplexType* matrix )
{
   for(int i=0; i<n_sum; i++)
   {
      int tmp_start = sum_pos_start[i];
      int* tmp_pos = sum_pos+tmp_start;
      int tmp_n = *(tmp_pos++);
      ComplexType tmp = sum[*tmp_pos++];

      for(int j=1; j<tmp_n; j++) tmp += sum[*tmp_pos++];

      matrix[*tmp_pos] = tmp;
   }
}

template <class ComplexType>
void CPUInstHomSumBlock<ComplexType>::print()
{
   std::cout << "n_sum = " << n_sum << std::endl;
   std::cout << "sum_pos_size = " << sum_pos_size << std::endl;
   for(int i=0; i<n_sum; i++)
   {
      int tmp_start = sum_pos_start[i];
      int* tmp_pos = sum_pos+tmp_start;
      int tmp_n = *(tmp_pos++);
      std::cout << "i = " << i << " n = " << tmp_n << ", ";

      for(int j=0; j<tmp_n; j++) std::cout << *tmp_pos++ << " ";

      std::cout << "   sum_pos_start = " << tmp_start << " output = "
                << *tmp_pos++ << std::endl;
   }
}
