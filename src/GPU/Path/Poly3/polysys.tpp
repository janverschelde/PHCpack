template <class ComplexType, class RealType>
ComplexType* PolySys<ComplexType,RealType>::eval ( const ComplexType* x_val )
{
   ComplexType* val = new ComplexType[n_eq];

   for(int i=0; i<n_eq; i++) val[i] = eq[i]->eval(x_val);

   return val;
}

template <class ComplexType, class RealType>
ComplexType* PolySys<ComplexType,RealType>::eval
 ( const ComplexType* x_val, ComplexType** deri_val )
{
   ComplexType* val = new ComplexType[n_eq];

   for(int i=0; i<n_eq; i++) val[i] = eq[i]->eval(x_val, deri_val[i]);

   return val;
}

template <class ComplexType, class RealType>
void PolySys<ComplexType,RealType>::balance_eq ( const ComplexType* x_val )
{
   ComplexType* f_val = eval(x_val);

   for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
      eq[eq_idx]->constant -= f_val[eq_idx];
}

template <class ComplexType, class RealType>
ComplexType** PolySys<ComplexType,RealType>::allocate_deg_table ( void )
{
   int n_total_deg = 0;
   for(int var_idx=0; var_idx<dim; var_idx++)
       n_total_deg += max_deg_base[var_idx];

   ComplexType** deg_table = new ComplexType*[dim];

   deg_table[0] = new ComplexType[n_total_deg];
   for(int var_idx=1; var_idx<dim; var_idx++)
      deg_table[var_idx]=deg_table[var_idx-1]+max_deg_base[var_idx-1];

   return deg_table;
}

template <class ComplexType, class RealType>
void PolySys<ComplexType,RealType>::compute_deg_table
 ( const ComplexType* x_val, ComplexType** deg_table )
{
   for(int var_idx=0; var_idx<dim; var_idx++)
   {
      if(max_deg_base[var_idx]>0)
      {
         ComplexType* tmp_deg_table = deg_table[var_idx];
         ComplexType tmp_var = x_val[var_idx];
         tmp_deg_table[0] = tmp_var;
         for(int deg_idx=1; deg_idx<max_deg_base[var_idx]; deg_idx++)
            tmp_deg_table[deg_idx] = tmp_deg_table[deg_idx-1]*tmp_var;
      }
   }
}

template <class ComplexType, class RealType>
ComplexType** PolySys<ComplexType,RealType>::eval_deg
 ( const ComplexType* x_val )
{
   ComplexType** deg_table = allocate_deg_table();

   compute_deg_table(x_val,deg_table);

   // print degree table
   // for(int var_idx=0; var_idx<dim; var_idx++)
   // {
   //    for(int deg_idx=0; deg_idx<max_deg_base[var_idx]; deg_idx++)
   //    {
   //       std::cout << var_idx << " " << deg_idx
   //                 << " " << deg_table[var_idx][deg_idx];
   //    }
   // }

   return deg_table;
}

template <class ComplexType, class RealType>
void PolySys<ComplexType,RealType>::eval
 ( const ComplexType* x_val, ComplexType* f_val, ComplexType** deri_val )
{
   if(eval_base)
   {
      ComplexType** deg_table = eval_deg(x_val);
      for(int i=0; i<n_eq; i++)
         f_val[i] = eq[i]->eval(x_val, deri_val[i], deg_table);

      delete[] deg_table[0];
      delete[] deg_table;
   }
   else
      for(int i=0; i<n_eq; i++)
         f_val[i] = eq[i]->eval(x_val, deri_val[i]);
}

template <class ComplexType, class RealType>
void PolySys<ComplexType,RealType>::eval
 ( const ComplexType* x_val, ComplexType* f_val, ComplexType** deri_val,
   ComplexType** deg_table )
{
   if(eval_base)
      for(int i=0; i<n_eq; i++)
         f_val[i] = eq[i]->eval(x_val, deri_val[i], deg_table);
   else
      for(int i=0; i<n_eq; i++)
         f_val[i] = eq[i]->eval(x_val, deri_val[i]);
}

template <class ComplexType, class RealType>
void PolySys<ComplexType,RealType>::print()
{
   cout << "dim = " << dim << ", n_eq = " << n_eq << endl;
   for(int i=0; i<n_eq; i++)
   {
       cout << "f"<< i << "=" << endl;
       eq[i]->print(pos_var);
   }
}

template <class ComplexType, class RealType>
void PolySys<ComplexType, RealType>::read
 ( const string& sys_string, VarDict& pos_dict, int verbose )
{
   if(verbose > 0) cout << "entering PolySys.read() ..." << endl;

   int l = sys_string.length();
   n_eq = 1;

   int eq_first = 0;

   LinkList<int>* eq_list = new LinkList<int>();

   eq_list->append(eq_first);

   int ii = eq_first;
   while(ii<l)
   {
      if(sys_string[ii] == ';')
      {
         ii++;
         while(ii<l)
         {
            if(sys_string[ii] != ' ')
            {
               n_eq++;
               eq_list->append(ii);
               break;
            }
            ii++;
         }
      }
      ii++;
   }
   eq_list->append(l);

   // eq = new PolyEq[n_eq];

   LinkNode<int>* tmp = eq_list->header();

   int eq_start = tmp->data;
   for(int i=0; i<n_eq; i++)
   {
      PolyEq<ComplexType,RealType>* new_eq = new PolyEq<ComplexType,RealType>();
      tmp = tmp->next;
      int eq_end = tmp->data;
      if(verbose > 0)
      {
         cout << n_eq << " "<< eq_start <<" " << eq_end<< endl;
         string new_eq = sys_string.substr(eq_start, eq_end-eq_start);
         cout << new_eq << endl;
      }
      new_eq->read(sys_string, pos_dict, eq_start, eq_end, verbose);
      eq.push_back(new_eq);
      eq_start = eq_end;
   }
   eq_list->destroy();
   delete eq_list;

   dim = pos_dict.n_job;
   pos_var = pos_dict.reverse();
   // dim = pos_dict.n_job;

   if(verbose > 0) cout << "... leaving polysys.read()" << endl;
}

template <class ComplexType, class RealType>
void PolySys<ComplexType,RealType>::read
 ( const string* sys_string, int n_eq, VarDict& pos_dict, int verbose )
{
   eq.reserve(n_eq);
   this->n_eq = n_eq;
   dim = n_eq;

   // Here is confusion, I should use either vector or array, 
   // but not to mix them
   eq_space = new PolyEq<ComplexType,RealType>[n_eq];
   PolyEq<ComplexType,RealType>* tmp_eq_space = eq_space;
   for(int i=0; i< n_eq; i++)
   {
      tmp_eq_space->read(sys_string[i], pos_dict, verbose);
      tmp_eq_space->dim = dim;
      eq.push_back(tmp_eq_space);
      tmp_eq_space++;
   }
   pos_var = pos_dict.reverse();
   update_max_deg_base();
}

template <class ComplexType, class RealType>
void PolySys<ComplexType,RealType>::read_file
 ( const string& file_name, int verbose )
{
   if(verbose > 0) cout << "entering PolySys.read_file() ..." << endl;

   VarDict pos_dict;
   ifstream myfile (file_name.c_str());
   if(myfile.is_open())
   {
      read_file(myfile, pos_dict, verbose);
      myfile.close();
   }
   else
   {
      cout << "There is no such file." << endl;
   }
   // dim = pos_dict.n_job;
   pos_var = pos_dict.reverse();
   /*
   std::cout << "dim = " << dim << std::endl;
   for(int var_idx=0; var_idx<dim; var_idx++)
   {
      std::cout << var_idx << " " << pos_var[var_idx] << std::endl;
   }*/

   if(verbose > 0) cout << "... leaving PolySys.read_file()" << endl;
}

template <class ComplexType, class RealType>
void PolySys<ComplexType,RealType>::read_file
 ( ifstream& myfile, VarDict& pos_dict, int verbose )
{
   if(verbose > 0) cout << "entering the second read_file ..." << endl;

   string line;
   getline(myfile,line);

   while(true)
   {
      if(not_empty_line(line)!=-1) break;

      getline(myfile,line);
   }

   std::vector<std::string> line_parts = split(line, ' ');
   if(line_parts.size() == 2)
   {
      dim = atoi(line_parts[1].c_str());
      n_eq = atoi(line_parts[0].c_str());
   }
   else
   {
      dim = atoi(line_parts[0].c_str());
      n_eq = dim;
   }
   eq.reserve(n_eq);
   eq_space = new PolyEq<ComplexType,RealType>[n_eq];
   PolyEq<ComplexType,RealType>* tmp_eq_space = eq_space;
   /*for(int i=0; i< n_eq; i++)
     {
        eq.push_back(tmp_eq_space);
        tmp_eq_space->dim = dim;
        tmp_eq_space++;
     }*/

   string tmp_line;
   line = "";
   int n_eq_file = 0;
   while( getline(myfile,tmp_line) )
   {
      int not_empty = not_empty_line(tmp_line);

      if(not_empty != -1)
      {
         int finish = 0;
         for(unsigned int j=not_empty; j<tmp_line.length(); j++)
         {
            char c =tmp_line[j];
            if(c ==';') finish = 1;
         }
         line += tmp_line;
         if(finish)
         {
            // std::cout << line << std::endl << std::endl;
            tmp_eq_space->read(line, pos_dict, verbose);
            tmp_eq_space->dim = dim;
            eq.push_back(tmp_eq_space);
            tmp_eq_space++;
            n_eq_file++;
            // std::cout << "n_eq_file = " << n_eq_file << std::endl;
            /* the following is not an error
            if(n_eq_file >= n_eq)
            {
               cout << "Error: too many lines" << endl;
               break;
            }*/
            line = "";
         }
      }
   }
   pos_var = pos_dict.reverse();

   update_max_deg_base();

   if(verbose > 0) cout << "... leaving the second read_file" << endl;
}

template <class ComplexType, class RealType>
void PolySys<ComplexType,RealType>::update_max_deg_base()
{
   max_deg_base = new int[dim];
   for(int i=0; i<dim; i++)
   {
      max_deg_base[i] = 0;
   }
   for(int i=0; i<n_eq; i++)
   {
      // cout << "eq " << i << endl;
      eq[i]->update_max_deg(max_deg_base);
   }

   for(int var_idx=0; var_idx<dim; var_idx++)
   {
      max_deg_base[var_idx] -= 1;
      if(max_deg_base[var_idx] > 0)
      {
         eval_base = true;
      }
   }
}

template <class ComplexType, class RealType>
void PolySys<ComplexType,RealType>::gpu_mon
 ( int& dim, int& level, int& workspace_size, int*& workspace_level,
   int*& n_mon_level, int& pos_size, unsigned short*& pos, int*& pos_level,
   int& sum_level, int*& n_sum_level,
   int& total_n_sum, int& sum_array_size, int*& sum_start, int*& sum_array )
{
   dim = this->dim;                  
    
   int max_n_var = 0;
   int total_n_mon = 0;

   // Count the largest n_var and the total n_mon
   for(int i=0; i<n_eq; i++)
   {
      PolyEq<ComplexType,RealType>* eq_tmp = eq[i];
      for(int j=0; j<eq_tmp->n_mon; j++)
      {
         int new_n_var = eq_tmp->mon[j]->n_var;
         if(new_n_var > max_n_var) max_n_var = new_n_var;
         total_n_mon++;
      }
   }
   // cout << "max_n_var = " << max_n_var << endl;
    
   // Put monomials index into vector ordered by n_var
   vector<int_idx>* n_var_list = new vector<int_idx>[max_n_var+1];
   for(int i=0; i<n_eq; i++)
   {
      PolyEq<ComplexType,RealType>* eq_tmp = eq[i];
      for(int j=0; j<eq_tmp->n_mon; j++)
      {
         int new_n_var = eq_tmp->mon[j]->n_var;
         n_var_list[new_n_var].push_back(int_idx(i,j));
         // cout << new_n_var << " " << i << " " << j << endl;
      }
   }
    
   /*for(int i=0; i< max_n_var+1; i++)
     {
        cout << i << " size = " << n_var_list[i].size() << endl;
        for(vector<int_idx>::iterator it = n_var_list[i].begin();
            it!=n_var_list[i].end(); ++it)
        {
           cout << (*it).eq_idx << " " << (*it).mon_idx << endl;
        }
     }*/
    
   // compute level
   level = 0;
   for(int i=1; i<max_n_var; i<<=1) level++;

   cout << " level = " << level << endl;
    
   // Count monomials on each level and count n_vars
   n_mon_level = new int[level + 1];
    
   for(int i=0; i<=level; i++)
   {
      n_mon_level[i] = 0;
      int m_size = 1<<i;
      for(int j= m_size/2+1; j<= min(max_n_var, m_size); j++)
      {
         int size_tmp = n_var_list[j].size();
         n_mon_level[i] += size_tmp;
      }
   }
    
   // Count pos on each level and total pos
   pos_level = new int[level + 1];
   pos_size = 0;
   workspace_level = new int[level + 1];
   workspace_size = 0;
   pos_level[0] = 0;
   workspace_level[0] = 0;
   int pos_level_size = 0;
   for(int i=0; i<=level; i++)
   {
      pos_level[i] = pos_level_size;
      workspace_level[i] = pos_level_size;
      pos_level_size = n_mon_level[i] * ((1<<i)+1);
      // cout << i << " pos_level_size = " << pos_level_size 
      //           << " level_start = " << pos_level[i] << endl;
      pos_size += pos_level_size;
      workspace_size += pos_level_size;
   }
   // cout << "pos_size = " << pos_size << " workspace_size = "
   //      << workspace_size << endl;
    
   // Put data into pos
   pos = new unsigned short[pos_size];    
   unsigned short* pos_tmp = pos;
   for(int i=0; i<=level; i++)
   {
      int m_size = 1<<i;
      // cout << "m_size = " << m_size << endl;
      for(int j= m_size/2+1; j<= min(max_n_var, m_size); j++)
      {
         // cout << "size = " << j << endl;
         for(vector<int_idx>::iterator it = n_var_list[j].begin();
             it!=n_var_list[j].end(); ++it)
         {
            int eq_idx = (*it).eq_idx;
            int mon_idx = (*it).mon_idx;
            int* mon_pos = eq[eq_idx] -> mon[mon_idx] -> pos;
            *pos_tmp++ = j;
            for(int k=0; k<j; k++)
            {
               pos_tmp[k] = mon_pos[k];
               // cout << pos_tmp[k] << " ";
            }
            // cout << endl;
            pos_tmp += m_size;
            // cout << m_size << endl;
         }
      }
   }
    
   // print pos
   pos_tmp = pos;
   for(int i=1; i<=level; i++)
   {
      int m_size = 1<<i;
      // cout << "level = " << i << " m_size = " << m_size
      //      << " level_start = " << pos_level[i] << endl;
      pos_tmp += pos_level[i];
      for(int j=0; j<n_mon_level[i]; j++)
      {
         //cout << m_size << endl;
         unsigned short* pos_start = pos_tmp + (m_size+1)*j;
         int n_var_tmp = *pos_start++;
         /*cout << "size = " << n_var_tmp << "   ";
           for(int k=0; k<n_var_tmp; k++) cout << pos_start[k] << " ";
           cout << endl;*/
      }
   }
    
   // Get location of each monomials in pos
   // cout << "n_mon = " << n_mon << endl;
   int* mon_pos_data = new int[total_n_mon];
   int** mon_pos = new int*[n_eq];
   mon_pos[0] = mon_pos_data;
   for(int i=1; i<n_eq; i++) mon_pos[i] = mon_pos[i-1] + eq[i-1] -> n_mon;
    
   int mon_loc = 0;
   for(int i=0; i<=level; i++)
   {
      int m_size = 1<<i;
      // cout << "m_size = " << m_size << endl;
      for(int j= m_size/2+1; j<= min(max_n_var, m_size); j++)
      {
         // cout << "size = " << j << endl;
         for(vector<int_idx>::iterator it = n_var_list[j].begin();
             it!=n_var_list[j].end(); ++it)
         {
            int eq_idx = (*it).eq_idx;
            int mon_idx = (*it).mon_idx;
            mon_pos[eq_idx][mon_idx] = mon_loc;
            mon_loc += m_size + 1;
         }
      }
   }
    
   /*for(int i=0; i<n_eq; i++)
     {
        cout << "eq " << i << endl;
        for(int j=0; j<eq[i] -> n_mon; j++){
            cout << mon_pos[i][j] << " "; 
        }
        cout << endl;
     }*/
    
   int* sum_size_data = new int[n_eq*(dim+1)];
   for(int i=0; i<n_eq*(dim+1); i++) sum_size_data[i] = 0;
   int** sum_size = new int*[n_eq];
   sum_size[0] = sum_size_data;
   for(int i=1; i<n_eq; i++) sum_size[i] = sum_size[i-1] + dim + 1;
    
   //cout << "Equation size:";
   /*for(int i=0; i<n_eq; i++)
     {
        //cout << sum_size[n_eq][i] << " ";
     }
     //cout << endl;*/
    
   // Equation Sum Size and Derivative Sum Size
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
   {
      sum_size[eq_idx][dim] = eq[eq_idx] -> n_mon; // +1 for constant
      // cout << "eq " << eq_idx << " eq_size " << eq[eq_idx] -> n_mon << endl;
      for(int mon_idx=0; mon_idx<eq[eq_idx] -> n_mon; mon_idx++)
      {
         int mon_start = mon_pos[eq_idx][mon_idx];
         int mon_size = pos[mon_start];
         // cout << "eq " << eq_idx << " mon_idx " << mon_idx
         //      << " mon_size " << mon_size << endl;
         for(int i=1; i<mon_size+1; i++)
            (sum_size[eq_idx][pos[mon_start+i]])++;
      }
   }
    
   // Largest number in Sum
   int max_sum_size = 0;
   total_n_sum = 0;
   int total_n_sum_pos = 0;
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
   {
      // cout << "eq " << eq_idx << "   ";
      for(int x_idx=0; x_idx < dim+1; x_idx++)
      {
         int tmp_sum_size = sum_size[eq_idx][x_idx];
         if(tmp_sum_size>max_sum_size)
         {
            max_sum_size = tmp_sum_size;
         }
         if(tmp_sum_size > 0)
         {
            total_n_sum++;
            total_n_sum_pos += tmp_sum_size;
         }
         // cout << sum_size[eq_idx][x_idx] << " ";
      }
      // cout << endl;
   }
   cout << "max_sum_size    = " << max_sum_size << endl;
   cout << "total_n_sum     = " << total_n_sum  << endl;
   cout << "total_n_sum_pos = " << total_n_sum_pos  << endl;
    
   // Put index into Sum Vector
   vector<int_idx>* sum_size_list = new vector<int_idx>[max_sum_size+1];
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
   {
      for(int x_idx=0; x_idx<dim+1; x_idx++)
      {
         sum_size_list[sum_size[eq_idx][x_idx]].push_back
            (int_idx(eq_idx, x_idx));
      }
   }

   // Generate n_sum_level
   sum_level = 0;
   for(int i=1; i<max_sum_size; i<<=1)
   {
      sum_level ++;
      // cout << "i = " << i << " level = " << level << endl;
   }
   cout << "sum_level = " << sum_level << endl;
    
   // Count n_sum on each level
   n_sum_level = new int[sum_level + 1];
    
   for(int i=0; i<=sum_level; i++)
   {
      n_sum_level[i] = 0;
      int m_size = 1<<i;
      for(int j= m_size/2+1; j<= min(max_sum_size, m_size); j++)
      {
         int size_tmp = sum_size_list[j].size();
         n_sum_level[i] += size_tmp;
      }
   }

   /*for(int i=0; i<=sum_level; i++)
     {
        cout << i << " " << n_sum_level[i] << endl;
     }*/

   // Generate Sum Start
   int* sum_start_loc_data = new int[n_eq*(dim+1)];
   for(int i=0; i<n_eq*(dim+1); i++) sum_start_loc_data[i] = 0;
    
   int** sum_start_loc = new int*[n_eq+1];
   sum_start_loc[0] = sum_start_loc_data;
   for(int i=1; i<n_eq; i++)
       sum_start_loc[i] = sum_start_loc[i-1] + dim + 1;
    
   sum_array_size = total_n_sum*2+total_n_sum_pos;
   sum_array = new int[sum_array_size]; // need to return
   for(int i=0; i<total_n_sum*2+total_n_sum_pos; i++) sum_array[i] = 0;

   sum_start = new int[total_n_sum+1]; // need to be return
   int* sum_start_tmp = sum_start;
   *sum_start_tmp = 0;
    
   for(int i=0; i<=max_sum_size; i++)
   {
      // cout << "size = " << i+2  << " n = "
      //      << sum_size_list[i].size() << endl;
      for(vector<int_idx>::iterator it = sum_size_list[i].begin();
          it!=sum_size_list[i].end(); ++it)
      {
         // cout << "   " << (*it).eq_idx << " " << (*it).mon_idx << endl;
         int eq_idx = (*it).eq_idx;
         int x_idx = (*it).mon_idx;
         sum_start_loc[eq_idx][x_idx] = *sum_start_tmp;
         sum_array[*sum_start_tmp+1] = eq_idx + x_idx*dim;
         *(sum_start_tmp+1) = (*sum_start_tmp) + i+2;
         sum_start_tmp++;
      }
   }
    
   /*for(int i=0; i<total_n_sum; i++)
     {
        cout << i << " " << sum_start[i] << endl;
     }*/
    
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
   {
      //cout << "eq " << eq_idx << endl;
      for(int mon_idx=0; mon_idx < eq[eq_idx]->n_mon; mon_idx++)
      {
         int mon_start = mon_pos[eq_idx][mon_idx];
         int mon_size = pos[mon_start];
         // cout << "eq " << eq_idx << " mon_idx " << mon_idx
         //      << " mon_size " << mon_size << endl;
         for(int i=1; i<mon_size+1; i++)
         {
            int x_idx = pos[mon_start+i];
            int deri_start = sum_start_loc[eq_idx][x_idx];
            int current_deri_size = sum_array[deri_start];
            sum_array[deri_start+current_deri_size+2] = mon_start+i;
            sum_array[deri_start]++;
         }
      }
   }
    
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
   {
      // cout << "eq[eq_idx]->n_mon = " << eq[eq_idx]->n_mon << endl;
      int deri_start = sum_start_loc[eq_idx][dim];
      // cout << "deri_start = " << deri_start << endl;
      for(int mon_idx=0; mon_idx < eq[eq_idx]->n_mon; mon_idx++)
      {
         int mon_start = mon_pos[eq_idx][mon_idx];
         int current_deri_size = sum_array[deri_start];
         // cout << "current_deri_size = " << current_deri_size << endl;
         sum_array[deri_start+current_deri_size+2] = mon_start;
         sum_array[deri_start]++;
      }
   }
    
   /*cout << "sum array" << endl;
     for(int i=0; i<total_n_sum; i++)
     {
        int sum_start_tmp = sum_start[i];
        for(int j=0; j<sum_array[sum_start_tmp]+2; j++)
        {
            cout << sum_array[sum_start_tmp+j] << " ";
        }
        cout << endl;
     }*/
            
   delete[] n_var_list;
   delete[] sum_size;
   delete[] sum_size_data;
   delete[] mon_pos_data;
   delete[] mon_pos;
   delete[] sum_size_list;
   delete[] sum_start_loc;
   delete[] sum_start_loc_data;
}
