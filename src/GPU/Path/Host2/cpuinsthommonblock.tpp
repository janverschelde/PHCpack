// The file cpuinsthommonblock.tpp provides the definitions of the methods
// in the class CPUInstHomMonBlock, with prototypes in cpuinsthommonblock.h.

template <class ComplexType>
void CPUInstHomMonBlock<ComplexType>::init
 ( CPUInstHomMon<ComplexType>& orig, int BS, int verbose )
{
   n_mon = orig.n_mon;
   this->BS = BS;
   n_mon_single = 0;
   for(int i=0; i<n_mon; i++)
      if(orig.mon_pos[orig.mon_pos_start[i]]==1) n_mon_single++;
   
   mon_single_pos_block = new unsigned short[2*n_mon_single];

   unsigned short* tmp_mon_single_pos_block = mon_single_pos_block;
   for(int i=0; i<n_mon; i++)
   {
      if(orig.mon_pos[orig.mon_pos_start[i]]==1)
      {
         *tmp_mon_single_pos_block++ = 1;
	 *tmp_mon_single_pos_block++ = orig.mon_pos[orig.mon_pos_start[i]+1];
      }
   }

   n_mon -= n_mon_single;
   NB = (n_mon-1)/BS + 1;
   max_var_block = new unsigned short[NB];
   for(int i=0; i<n_mon; i++)
   {
      int bidx = i/BS;
      int tidx = i - bidx*BS;
      unsigned short n_var = orig.mon_pos[orig.mon_pos_start[i+n_mon_single]];
      if(tidx == 0)
      {
         max_var_block[bidx] = n_var;
      }
      else
      {
         if(n_var > max_var_block[bidx]) max_var_block[bidx] = n_var;
      }
   }
   mon_pos_start_block = new int[NB];
   mon_pos_block_size = 2*n_mon_single;
   for(int i=0; i<NB; i++)
   {
      mon_pos_start_block[i] = mon_pos_block_size;
      mon_pos_block_size += BS*(max_var_block[i]+1);
   }
   if(verbose > 0)
   {
      std::cout << "mon_pos_block_size = " << mon_pos_block_size << std::endl;
   }
   mon_pos_block = new unsigned short[mon_pos_block_size];
   unsigned short* tmp_mon_pos_block = mon_pos_block;
   for(int i=0; i<n_mon; i++)
   {
      if(orig.mon_pos[orig.mon_pos_start[i]]==1)
      {
         *tmp_mon_pos_block++ = 1;
         *tmp_mon_pos_block++ = orig.mon_pos[orig.mon_pos_start[i]+1];
      }
   }
   unsigned short* tmp_mon_pos = orig.mon_pos+orig.mon_pos_start[n_mon_single];
   for(int i=0; i<NB; i++)
   {
      unsigned short* tmp_mon_pos_block = mon_pos_block 
         + mon_pos_start_block[i];
      for(int j=0; j<BS; j++)
      {
         int mon_idx = i*BS+j;
         if(mon_idx<n_mon)
         {
            unsigned short n_var = *tmp_mon_pos++;
            tmp_mon_pos_block[j] = n_var;

            for(int k=0; k<n_var; k++)
               tmp_mon_pos_block[(k+1)*BS+j] = *tmp_mon_pos++;
         }
      }
   }
}

template <class ComplexType>
void CPUInstHomMonBlock<ComplexType>::print()
{
   std::cout << "BS = " << BS << std::endl;
   std::cout << "NB = " << NB << std::endl;
   for(int i=0; i<NB; i++)
   {
      std::cout << "BS " << i << " n_var = " << max_var_block[i] \
                << " start = " << mon_pos_start_block[i] << std::endl;
      unsigned short* tmp_mon_pos_block = mon_pos_block
         + mon_pos_start_block[i];
      for(int j=0; j<BS; j++)
      {
         int mon_idx = i*BS+j;
         if(mon_idx<n_mon)
         {
            unsigned short n_var = tmp_mon_pos_block[j];
            std::cout << mon_idx << " n_var=" << n_var;
            for(int k=0; k<n_var; k++)
            {
               std::cout << " " << tmp_mon_pos_block[(k+1)*BS+j];
            }
            std::cout << std::endl;
         }
      }
   }
}
