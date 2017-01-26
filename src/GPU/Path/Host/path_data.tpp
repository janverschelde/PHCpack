// definitions of the templated code with prototypes in path_data.h

template <class ComplexType, class RealType>
ComplexType* read_pt ( ifstream& path_file, int dim )
{
   ComplexType* pt = new ComplexType[dim];
   string tmp_line;
   for(int i=0; i<dim; i++)
   {
      getline(path_file,tmp_line, ':');
      pt[i] = get_complex_number<ComplexType>(path_file);
      getline(path_file,tmp_line);
   }
   return pt;
}

template <class ComplexType, class RealType>
void print_pt ( ComplexType* pt, int dim )
{
   int pr = 2 * sizeof(RealType);
   std::cout << setprecision(pr);
   std::cout << scientific;

   for(int i=0; i<dim; i++)
   {
      std::cout << i << ": ";
      if(pt[i].real > 0)
      {
         std::cout << "+";
      }
      std::cout << pt[i].real;
      std::cout<< "  ";
      if(pt[i].imag > 0)
      {
         std::cout << "+";
      }
      std::cout << pt[i].imag << std::endl;
   }
}

template <class ComplexType, class RealType>
void PathStep<ComplexType,RealType>::read_phc_file ( ifstream& path_file )
{
   string tmp_line;
   path_file >> delta_t;
   getline(path_file,tmp_line, ':');
   path_file >> t.real;
   path_file >> t.imag;
   getline(path_file,tmp_line);
   getline(path_file,tmp_line);
   predict_pt = read_pt<ComplexType,RealType>(path_file, dim);
   read_correct_it(path_file);
   read_until_line(path_file, "the solution");
   correct_pt = read_pt<ComplexType,RealType>(path_file, dim);
}

template <class ComplexType, class RealType>
void PathStep<ComplexType,RealType>::read_correct_it ( ifstream& path_file )
{
   string tmp_line;
   getline(path_file,tmp_line, ':');
   while(tmp_line == "correction (a&r)")
   {
      correct_iteration tmp_it;
      path_file >> tmp_it.correct_a;
      path_file >> tmp_it.correct_r;
      getline(path_file,tmp_line, ':');
      path_file >> tmp_it.residual_a;
      path_file >> tmp_it.residual_r;
      correct_it.push_back(tmp_it);
      n_it++;
      getline(path_file,tmp_line);
      getline(path_file,tmp_line, ':');
   }
}

template <class ComplexType, class RealType>
void PathStep<ComplexType,RealType>::update_predict_pt ( ComplexType* predict_pt )
{
   this-> predict_pt = new ComplexType[dim];

   for(int i=0; i<dim; i++) this->predict_pt[i] = predict_pt[i];
}

template <class ComplexType, class RealType>
void PathStep<ComplexType,RealType>::update_correct_pt ( ComplexType* correct_pt )
{
   this->correct_pt = new ComplexType[dim];

   for(int i=0; i<dim; i++) this->correct_pt[i] = correct_pt[i];
}

template <class ComplexType, class RealType>
void PathStep<ComplexType,RealType>::update_t ( ComplexType delta_t, ComplexType t )
{
   this->delta_t = delta_t.real;
   this->t = t;
}

template <class ComplexType, class RealType>
void PathStep<ComplexType, RealType>::add_iteration
 ( double max_delta_x, double r_max_delta_x )
{
   correct_iteration tmp_it;
   tmp_it.correct_a = max_delta_x;
   tmp_it.correct_r = r_max_delta_x;
   tmp_it.residual_a = 0;
   tmp_it.residual_r = 0;
   correct_it.push_back(tmp_it);
   n_it++;
}

template <class ComplexType, class RealType>
void PathStep<ComplexType,RealType>::update_iteration_res
 ( double residual_a, double residual_r )
{
   correct_it[n_it-1].residual_a = residual_a;
   correct_it[n_it-1].residual_r = residual_r;
}

template <class ComplexType, class RealType>
void PathStep<ComplexType,RealType>::print()
{
   std::cout << scientific;
   std::cout << setprecision(5);
   std::cout << delta_t  << "  t :  " << t.real << "   " << t.imag
             << std::endl;
   print_pt<ComplexType,RealType>(predict_pt, dim);

   print_it();
   std::cout << "the solution for t :" << std::endl;
   if(correct_pt != NULL)
   {
      print_pt<ComplexType,RealType>(correct_pt, dim);
   }
   else
   {
      print_pt<ComplexType,RealType>(predict_pt, dim);
   }
   std::cout << "== err :  " << std::endl;
}

template <class ComplexType, class RealType>
void PathStep<ComplexType,RealType>::print_it()
{
   std::cout << scientific << setprecision(3);

   for(int i=0; i<n_it; i++)
   {
      std::cout << "correct (a&r):  " << correct_it[i].correct_a
                << "  " << correct_it[i].correct_r
                << " residue (a&r):  " << correct_it[i].residual_a
                << "  " << correct_it[i].residual_r
                << std::endl;
   }
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::read_phc_file ( ifstream& path_file )
{
   string output_symbol = "OUTPUT INFORMATION";
   read_until_line(path_file, output_symbol);
   string sol_symbol = "the solution";
   read_until_line(path_file, sol_symbol);
   start_pt = read_pt<ComplexType,RealType>(path_file, dim);
   string step_symbol = "== err : ";
   string corrector_symbol = "correction (a&r)";

   string tmp_line;
   while(true)
   {
      read_until_line(path_file, step_symbol);
      getline(path_file,tmp_line, ':');
      if(tmp_line != "step ")
      {
         if(tmp_line == "correction (a&r)")
         {
            read_until_line(path_file, step_symbol);
            read_until_line(path_file, step_symbol);
            getline(path_file,tmp_line, ':');
            if(tmp_line != "step ")
            {
               break;
            }
         }
         else
         {
            break;
         }
      }
      add_step(path_file);
   }
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::add_step ( ifstream& path_file )
{
   PathStep<ComplexType,RealType>* new_step
      = new PathStep<ComplexType,RealType>(dim,path_file);

   steps.push_back(new_step);
   n_step++;
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::add_step_empty()
{
   PathStep<ComplexType,RealType>* new_step
      = new PathStep<ComplexType,RealType>(dim);

   steps.push_back(new_step);
   n_step++;
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::update_step_predict_pt
 ( ComplexType* predict_pt )
{
   steps[n_step-1]->update_predict_pt(predict_pt);
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::update_step_correct_pt
 ( ComplexType* correct_pt )
{
   steps[n_step-1]->update_correct_pt(correct_pt);
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::update_step_t
 ( ComplexType delta_t, ComplexType t )
{
   steps[n_step-1]->update_t(delta_t,t);
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::add_iteration
 ( double max_delta_x, double r_max_delta_x )
{
   steps[n_step-1]->add_iteration(max_delta_x, r_max_delta_x);
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::update_iteration_res
 ( double residual_a, double residual_r )
{
   steps[n_step-1]->update_iteration_res(residual_a, residual_r);
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::print_t()
{
   std::cout << "n_step = " << n_step << std::endl;
   for(int i=0; i<n_step; i++)
   {
      std::cout << i << " " << steps[i]->t << " " << steps[i]->delta_t
                << " " << 0.0 << std::endl;
   }
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::print()
{
   std::cout << "n_step = " << n_step << std::endl;
   for(int i=0; i<n_step; i++)
   {
      std::cout << "Step " << i << std::endl;
      steps[i]->print();
   }
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::print_phc()
{
   std::cout << "n_step = " << n_step
             << " success = " << success << std::endl;
   std::cout << "OUTPUT INFORMATION DURING CONTINUATION :" << std::endl;
   if(start_pt != NULL)
   {
      std::cout << "the solution for t :" << std::endl;
      print_pt<ComplexType,RealType>(start_pt,dim);
      std::cout << "== err :  " << std::endl;
   }
   for(int i=0; i<n_step; i++)
   {
      std::cout << "step :  ";
      steps[i]->print();
   }
   if(end_pt != NULL)
   {
      std::cout << "the solution for t :" << std::endl;
      print_pt<ComplexType,RealType>(end_pt,dim);
      std::cout << "== err :  " << std::endl;
   }
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::compare ( Path& that )
{
   std::cout << "n_step = " << n_step << std::endl;
   for(int i=0; i<n_step; i++)
   {
      if((abs(steps[i]->t.real - that.steps[i]->t.real) > 1E-10)
       ||(abs(steps[i]->t.imag - that.steps[i]->t.imag) > 1E-10))
      {
         std::cout << i << " " << steps[i]->t.real << " "
                   << that.steps[i]->t.real << std::endl;
         break;
      }
   }
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::add_start_pt ( ComplexType* start_pt )
{
   this->start_pt = new ComplexType[dim];

   for(int i=0; i<dim; i++) this->start_pt[i] = start_pt[i];
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::add_end_pt ( ComplexType* end_pt )
{
   this->end_pt = new ComplexType[dim];

   for(int i=0; i<dim; i++) this->end_pt[i] = end_pt[i];
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::update_success ( bool success )
{
   this->success = success;
}

template <class ComplexType, class RealType>
void Path<ComplexType,RealType>::clear()
{
   n_step = 0;
   success = false;
   start_pt = NULL;
   end_pt = NULL;
   for(typename vector<PathStep<ComplexType,RealType>*>::iterator
       it = steps.begin(); it!=steps.end(); ++it) 
   {
      delete (*it);
   }
   steps.clear();
}
