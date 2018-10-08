// The templated definitions are included in polysolset.h.

template <class ComplexType, class RealType>
void PolySolSet<ComplexType,RealType>::init ( ifstream& myfile )
{
   n_sol = 0;
   dim = 0;

   string prefix_two[2] = {"START SOLUTIONS : ", "THE SOLUTIONS :"};
   read_until_line(myfile, prefix_two, 2);

   myfile >> n_sol;
   myfile >> dim;

   string prefix = "======";
   read_until_line(myfile, prefix);

   sols.reserve(n_sol);

   //std::cout << "n_sol = " << n_sol << std::endl;
   for(int i=0; i<n_sol; i++)
   {
      // std::cout << i << std::endl;
      PolySol<ComplexType,RealType>* tmp_sol
         = new PolySol<ComplexType,RealType>(myfile, dim);
      sols.push_back(tmp_sol);
   }
}

template <class ComplexType, class RealType>
void PolySolSet<ComplexType,RealType>::init
 ( ifstream& myfile, VarDict& pos_dict )
{
   n_sol = 0;
   dim = 0;

   string prefix_two[2] = {"START SOLUTIONS : ", "THE SOLUTIONS :"};
   read_until_line(myfile, prefix_two, 2);

   myfile >> n_sol;
   myfile >> dim;

   string prefix = "======";
   read_until_line(myfile, prefix);

   sols.reserve(n_sol);

   // std::cout << "n_sol = " << n_sol << std::endl;
   for(int i=0; i<n_sol; i++)
   {
      // std::cout << i << std::endl;
      PolySol<ComplexType,RealType>* tmp_sol
         = new PolySol<ComplexType,RealType>(myfile, dim, pos_dict);
      sols.push_back(tmp_sol);
   }
}

template <class ComplexType, class RealType>
bool PolySolSet<ComplexType,RealType>::find_same_sol
 ( PolySol<ComplexType,RealType>* tmp_sol )
{
   for(int i=0; i<n_sol; i++)
   {
      if(*tmp_sol == *sols[i])
      {
         return true;
      }
   }
   return false;
}

template <class ComplexType, class RealType>
int PolySolSet<ComplexType,RealType>::count_same_sol
 ( PolySol<ComplexType,RealType>* tmp_sol )
{
   int n_same_sol = 0;
   for(int i=0; i<n_sol; i++)
   {
      if(*tmp_sol == *sols[i])
      {
         n_same_sol++;
      }
   }
   return n_same_sol;
}

template <class ComplexType, class RealType>
bool PolySolSet<ComplexType,RealType>::add_diff_sol ( ComplexType* new_sol )
{
   PolySol<ComplexType,RealType>* tmp_sol
      = new PolySol<ComplexType,RealType>(new_sol, dim);
   if(find_same_sol(tmp_sol)==true)
   {
      return false;
   }
   std::cout << "Add New Solution" << std::endl;
   sols.push_back(tmp_sol);
   n_sol++;
   return true;
}

template <class ComplexType, class RealType>
void PolySolSet<ComplexType,RealType>::add_sol
 ( PolySol<ComplexType,RealType>* tmp_sol )
{
   sols.push_back(tmp_sol);
   n_sol++;
}

template <class ComplexType, class RealType>
void PolySolSet<ComplexType,RealType>::change_sol
 ( int idx, ComplexType* coords )
{
   PolySol<ComplexType,RealType>* idxsol = sols[idx];
   for(int k=0; k<dim; k++) idxsol->sol[k] = coords[k];
}

template <class ComplexType, class RealType>
void PolySolSet<ComplexType,RealType>::change_solt
 ( int idx, ComplexType* coords, ComplexType* tval )
{
   PolySol<ComplexType,RealType>* idxsol = sols[idx];
   for(int k=0; k<dim; k++) idxsol->sol[k] = coords[k];
   idxsol->t = *tval;
}

template <class ComplexType, class RealType>
void PolySolSet<ComplexType,RealType>::add_sol
 ( ComplexType* new_sol, RealType max_residual, RealType max_delta_x,
   int path_idx, string path_info )
{
   PolySol<ComplexType,RealType>* tmp_sol
      = new PolySol<ComplexType,RealType>
           (new_sol,dim,max_residual,max_delta_x,path_idx,path_info);
   add_sol(tmp_sol);
}

template <class ComplexType, class RealType>
void PolySolSet<ComplexType,RealType>::print()
{
   std::cout << "dim   = " << dim << std::endl
             << "n_sol = " << n_sol << std::endl;

   for(int i=0; i<n_sol; i++) sols[i]->print();
}

template <class ComplexType, class RealType>
void PolySolSet<ComplexType,RealType>::print_info ( string* pos_var )
{
   std::cout << "dim   = " << dim << std::endl
             << "n_sol = " << n_sol << std::endl;

   for(int i=0; i<n_sol; i++) sols[i]->print_info(pos_var);
}

template <class ComplexType, class RealType>
void PolySolSet<ComplexType,RealType>::print_short()
{
   std::cout << "n_sol = " << n_sol << std::endl;
   for(int i=0; i<n_sol; i++)
      std::cout << i << " " << sols[i]->info << " x[0] = " << sols[i]->sol[0];
}

template <class ComplexType, class RealType>
ComplexType* PolySolSet<ComplexType,RealType>::get_sol ( int idx )
{
   return sols[idx]->get_sol();
}

template <class ComplexType, class RealType>
void PolySolSet<ComplexType,RealType>::get_solt
 ( int idx, ComplexType* sol, ComplexType *t )
{
   PolySol<ComplexType,RealType>* solution = sols[idx];

   for(int i=0; i<solution->dim; i++) sol[i] = solution->sol[i];

   *t = solution->t;
}

template <class ComplexType, class RealType>
void PolySolSet<ComplexType,RealType>::sort_set()
{
   sort(sols.begin(), sols.end(), compare_sol);
}

template <class ComplexType, class RealType>
void PolySolSet<ComplexType,RealType>::compare
 ( PolySolSet<ComplexType,RealType>& that )
{
   vector<PolySol<ComplexType,RealType>*> that_sols(that.sols);
   vector<PolySol<ComplexType,RealType>*> this_sols(this->sols);

   sort(this_sols.begin(), this_sols.end(), compare_sol);
   sort(that_sols.begin(), that_sols.end(), compare_sol);

   typename vector<PolySol<ComplexType,RealType>*>::iterator this_pointer
      = this_sols.begin();
   typename vector<PolySol<ComplexType,RealType>*>::iterator that_pointer
      = that_sols.begin();

   vector<PolySol<ComplexType,RealType>*> this_sols_only;
   vector<PolySol<ComplexType,RealType>*> that_sols_only;

   int same_sol = 0;
   int this_sol = 0;
   int that_sol = 0;

   while(true)
   {
      if((**this_pointer)==(**that_pointer))
      {
         this_pointer++;
         that_pointer++;
         same_sol++;
      }
      else
      {
         if((**this_pointer)<(**that_pointer))
         {
            this_sols_only.push_back(*this_pointer);
            this_pointer++;
            this_sol++;
         }
         else
         {
            /* std::cout << this_pointer - this_sols.begin() << std::endl;
               std::cout << that_pointer - that_sols.begin() << std::endl;
               (*this_pointer)->print();
               (*that_pointer)->print();
             */
            that_sols_only.push_back(*that_pointer);
            that_pointer++;
            that_sol++;
         }
      }
      if(this_pointer==this_sols.end() || that_pointer==that_sols.end())
      {
         break;
      }
   }

   if(this_pointer != this_sols.end())
   {
      this_sol += this_sols.end() - this_pointer;
      while(this_pointer < this_sols.end())
      {
         this_sols_only.push_back(*this_pointer);
         this_pointer++;
      }
   }

   if(that_pointer != that_sols.end())
   {
      that_sol += that_sols.end() - that_pointer;
      while(that_pointer < that_sols.end())
      {
         that_sols_only.push_back(*that_pointer);
         that_pointer++;
      }
   }
   std::cout << "same_sol = " << same_sol << std::endl;

   if(this_sol > 0) std::cout << "this_sol = " << this_sol << std::endl;
   if(that_sol > 0) std::cout << "that_sol = " << that_sol << std::endl;

   that_pointer = that_sols_only.begin();
   int i=0;
   while(that_pointer < that_sols_only.end())
   {
      bool find_same_this = find_same_sol(*that_pointer);
      int that_n_same_sol = that.count_same_sol(*that_pointer);
      std::cout << std::setw(4) << i++ << " "
                << std::setw(8) << (*that_pointer)->path_idx << " "
                << std::setw(4) << find_same_this << " "
                << std::setw(4) << that_n_same_sol;
      (*that_pointer)->print_short();
      that_pointer++;
   }
}
