/* path_data.cpp created on Feb 21, 2015 by yxc with edits by jv */

#include "path_data.h"

CT* read_pt ( ifstream& path_file, int dim )
{
   CT* pt= new CT[dim];
   string tmp_line;
   for(int i=0; i<dim; i++)
   {
      getline(path_file,tmp_line, ':');
      pt[i] = get_complex_number(path_file);
      getline(path_file,tmp_line);
   }
   return pt;
}

void print_pt ( CT* pt, int dim )
{
   int pr = 2 * sizeof(T1);
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

void PathStep::read_phc_file ( ifstream& path_file )
{
   string tmp_line;
   path_file >> delta_t;
   getline(path_file,tmp_line, ':');
   path_file >> t.real;
   path_file >> t.imag;
   getline(path_file,tmp_line);
   getline(path_file,tmp_line);
   predict_pt = read_pt(path_file, dim);
   read_correct_it(path_file);
   read_until_line(path_file, "the solution");
   correct_pt = read_pt(path_file, dim);
}

void PathStep::read_correct_it ( ifstream& path_file )
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

void PathStep::update_predict_pt ( CT* predict_pt )
{
   this-> predict_pt = new CT[dim];
   for(int i=0; i<dim; i++)
   {
      this->predict_pt[i] = predict_pt[i];
   }
}

void PathStep::update_correct_pt ( CT* correct_pt )
{
   this->correct_pt = new CT[dim];
   for(int i=0; i<dim; i++)
   {
      this->correct_pt[i] = correct_pt[i];
   }
}

void PathStep::update_t ( CT delta_t, CT t )
{
   this->delta_t = delta_t.real;
   this->t = t;
}

void PathStep::add_iteration ( double max_delta_x, double r_max_delta_x )
{
   correct_iteration tmp_it;
   tmp_it.correct_a = max_delta_x;
   tmp_it.correct_r = r_max_delta_x;
   tmp_it.residual_a = 0;
   tmp_it.residual_r = 0;
   correct_it.push_back(tmp_it);
   n_it++;
}

void PathStep::update_iteration_res ( double residual_a, double residual_r )
{
   correct_it[n_it-1].residual_a = residual_a;
   correct_it[n_it-1].residual_r = residual_r;
}


void PathStep::print()
{
   std::cout << scientific;
   std::cout << setprecision(5);
   std::cout << delta_t  << "  t :  " << t.real << "   " << t.imag
             << std::endl;
   print_pt(predict_pt, dim);

   print_it();
   std::cout << "the solution for t :" << std::endl;
   if(correct_pt != NULL)
   {
      print_pt(correct_pt, dim);
   }
   else
   {
      print_pt(predict_pt, dim);
   }
   std::cout << "== err :  " << std::endl;
}

void PathStep::print_it()
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

void Path::read_phc_file ( ifstream& path_file )
{
   string output_symbol = "OUTPUT INFORMATION";
   read_until_line(path_file, output_symbol);
   string sol_symbol = "the solution";
   read_until_line(path_file, sol_symbol);
   start_pt = read_pt(path_file, dim);
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

void Path::add_step ( ifstream& path_file )
{
   PathStep* new_step = new PathStep(dim,path_file);
   steps.push_back(new_step);
   n_step++;
}

void Path::add_step_empty()
{
   PathStep* new_step = new PathStep(dim);
   steps.push_back(new_step);
   n_step++;
}

void Path::update_step_predict_pt( CT* predict_pt )
{
   steps[n_step-1]->update_predict_pt(predict_pt);
}

void Path::update_step_correct_pt ( CT* correct_pt )
{
   steps[n_step-1]->update_correct_pt(correct_pt);
}

void Path::update_step_t ( CT delta_t, CT t )
{
   steps[n_step-1]->update_t(delta_t,t);
}

void Path::add_iteration ( double max_delta_x, double r_max_delta_x )
{
   steps[n_step-1]->add_iteration(max_delta_x, r_max_delta_x);
}

void Path::update_iteration_res ( double residual_a, double residual_r )
{
   steps[n_step-1]->update_iteration_res(residual_a, residual_r);
}

void Path::print_t()
{
   std::cout << "n_step = " << n_step << std::endl;
   for(int i=0; i<n_step; i++)
   {
      std::cout << i << " " << steps[i]->t << " " << steps[i]->delta_t
                << " " << 0.0 << std::endl;
   }
}

void Path::print()
{
   std::cout << "n_step = " << n_step << std::endl;
   for(int i=0; i<n_step; i++)
   {
      std::cout << "Step " << i << std::endl;
      steps[i]->print();
   }
}

void Path::print_phc()
{
   std::cout << "n_step = " << n_step
             << " success = " << success << std::endl;
   std::cout << "OUTPUT INFORMATION DURING CONTINUATION :" << std::endl;
   if(start_pt != NULL)
   {
      std::cout << "the solution for t :" << std::endl;
      print_pt(start_pt,dim);
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
      print_pt(end_pt,dim);
      std::cout << "== err :  " << std::endl;
   }
}

void Path::compare ( Path& that )
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

void Path::add_start_pt ( CT* start_pt )
{
   this->start_pt = new CT[dim];
   for(int i=0; i<dim; i++)
   {
      this->start_pt[i] = start_pt[i];
   }
}

void Path::add_end_pt ( CT* end_pt )
{
   this->end_pt = new CT[dim];
   for(int i=0; i<dim; i++)
   {
      this->end_pt[i] = end_pt[i];
   }
}

void Path::update_success ( bool success )
{
   this->success = success;
}

void Path::clear()
{
   n_step = 0;
   success = false;
   start_pt = NULL;
   end_pt = NULL;
   for(vector<PathStep*>::iterator it = steps.begin(); it!=steps.end(); ++it) 
   {
      delete (*it);
   }
   steps.clear();
}
