// contains the definitions for the prototypes in workspace_host.h

template <class ComplexType>
void Workspace<ComplexType>::init
 ( int workspace_size, int n_coef, int n_constant,
   int n_eq, int dim, int n_predictor, int* max_deg_base )
{
   this->dim = dim;
   this->n_eq = n_eq;
   all = new ComplexType[workspace_size];
   coef = all;
   // std::cout << "n_constant = " << n_constant << std::endl;
   // std::cout << "n_coef = " << n_coef << std::endl;
   // std::cout << "workspace_size = " << workspace_size << std::endl;
   mon = coef + n_coef;
   sum = mon - n_constant;
   matrix = new ComplexType[n_eq*(dim + 1)];

   sol = new ComplexType[dim];

   R = new ComplexType*[dim+1];
   R[0] = new ComplexType[(dim+1)*(dim+1)];
   for(int i=0; i<dim; i++)
   {
      R[i+1] = R[i] + dim + 1;
      for(int j=0; j<dim+1; j++) R[i][j].init(0.0,0.0);
   }
   for(int j=0; j<dim+1; j++) R[dim][j].init(0.0,0.0);

   // std::cout << "dim = " << dim << " n_eq = " << n_eq <<std::endl;

   V = new ComplexType*[dim+1];
   V[0] = matrix;
   for(int i=0; i<dim; i++)
   {
      V[i+1] = V[i] + n_eq;
   }

   // tmp_x = new ComplexType[dim];

   this->init_x_t(dim, n_predictor); // added this->
   if(max_deg_base != NULL)
   {
      init_deg_table(max_deg_base);
   }
}

template <class ComplexType>
void Workspace<ComplexType>::init_deg_table ( int* max_deg_base )
{
   int n_total_deg = 0;
   for(int var_idx=0; var_idx<dim; var_idx++)
   {
      n_total_deg += max_deg_base[var_idx];
   }
   deg_table = new ComplexType*[dim];
   // std::cout << "n_total_deg = " << n_total_deg << std::endl;
   deg_table[0] = new ComplexType[n_total_deg];
   for(int var_idx=1; var_idx<dim; var_idx++)
   {
      deg_table[var_idx]=deg_table[var_idx-1]+max_deg_base[var_idx-1];
   }
}

template <class ComplexType>
void Workspace<ComplexType>::init_x_t ( int dim, int n_predictor )
{
   this->dim = dim;
   this->n_predictor = n_predictor;
   this->n_array = n_predictor+1;  // added this->

   this->x_t_idx = 0; // added this->

   this->x_array = new ComplexType*[this->n_array]; // added this->
   this->x_array[0] = new ComplexType[this->n_array*dim]; // added this->
   for(int i=0; i<n_predictor; i++)
   {
      this->x_array[i+1] = this->x_array[i] + dim; // added this->
   }
   this->t_array = new ComplexType[this->n_array]; // added this->

   this->x = this->x_array[0]; // added this->
   this->t = this->t_array;    // added this->
}

template <class ComplexType>
void Workspace<ComplexType>::init_x_t_idx()
{
   std::cout << "inside init_x_t_idx() ..." << std::endl;
   this->x_t_idx = 0; // added this->
   std::cout << "assigning x_array[0] to x ..." << std::endl;
   this->x = this->x_array[0]; // added this->
   std::cout << "after assignement" << std::endl;
   this->t = this->t_array; // added this->
   std::cout << "leaving init_x_t_idx()" << std::endl;
}

template <class ComplexType>
void Workspace<ComplexType>::update_x_t
 ( ComplexType* cpu_sol0, ComplexType cpu_t )
{
   update_x_t_value(cpu_sol0, cpu_t);
   update_x_t_idx();
}

template <class ComplexType>
void Workspace<ComplexType>::update_x_t_value
 ( ComplexType* cpu_sol0, ComplexType cpu_t )
{
   update_x_value(cpu_sol0);
   update_t_value(cpu_t);
}

template <class ComplexType>
void Workspace<ComplexType>::update_x_value ( ComplexType* cpu_sol0 )
{
   for(int i=0; i<dim; i++)
   {
      this->x[i] = cpu_sol0[i]; // added this->
   }
}

template <class ComplexType>
void Workspace<ComplexType>::update_t_value ( ComplexType cpu_t )
{
   *t = cpu_t;
}

template <class ComplexType>
void Workspace<ComplexType>::update_x_t_idx()
{
   x_last = x;
   t_last = t;
   x_t_idx = (x_t_idx+1)%n_array;
   x = x_array[x_t_idx];
   t = t_array + x_t_idx;
}

template <class ComplexType>
void Workspace<ComplexType>::print_t()
{
   for(int i=0; i<n_array; i++)
   {
      std::cout << i << " "<< t_array[i];
   }
}

template <class ComplexType>
void Workspace<ComplexType>::print_coef()
{
   ComplexType* tmp_coef = coef;
   int tmp_idx = 0;
   while(tmp_coef < mon)
   {
      std::cout << tmp_idx++ << " " << *tmp_coef++;
   }
}

template <class ComplexType>
void Workspace<ComplexType>::print_result()
{
   for(int i=0; i<n_eq; i++)
   {
      for(int j=0; j<dim+1; j++)
      {
         std::cout << i << " " << j << " " << matrix[i*(dim+1)+j];
      }
   }
}

template <class ComplexType>
void Workspace<ComplexType>::print_x()
{
   for(int i=0; i<dim; i++)
   {
      std::cout << i << " " << x[i];
   }
}

template <class ComplexType>
void Workspace<ComplexType>::print_x_array()
{
   for(int i=0; i<dim; i++)
   {
      for(int j=0; j<n_array; j++)
      {
         std::cout << i << " " << j << " " << x_array[j][i];
      }
      std::cout << std::endl;
   }
}

template <class ComplexType>
void Workspace<ComplexType>::print_x_last()
{
   for(int i=0; i<dim; i++)
   {
      std::cout << i << " " << x_last[i];
   }
}

template <class ComplexType>
void Workspace<ComplexType>::init_x_t_predict_test()
{
   std::cout << "------- Initializa x and t value for Testing only ---------"
             << std::endl;
   for(int i=1; i<n_array; i++)
   {
      for(int j=0; j<dim; j++)
      {
         x_array[i][j] = ComplexType(i*i*i+5,0);
         // std::cout << i << " " << j << " " << x_array[i][j];
      }
      t_array[i] = ComplexType(i,0);
      // std::cout << i << " " << t_array[i];
   }
}
