// The file cpuinsthom.tpp provides the definitions of the methods in
// the class CPUInstHom, with prototypes in the file cpuinsthom.h.

/* start with definitions of helper functions */

template <class ComplexType, class RealType>
void polysys_mon_set
 ( const PolySys<ComplexType,RealType>& Target_Sys,
   MonIdxSet<ComplexType>* mon_set, bool sys_idx )
{
   PolyEq<ComplexType,RealType>* tmp_eq = Target_Sys.eq_space;

   int mon_idx = 0;
   for(int i=0; i<Target_Sys.n_eq; i++)
   {
      if(tmp_eq->constant.real != 0.0 || tmp_eq->constant.imag != 0.0)
      {
         mon_set[mon_idx] = MonIdxSet<ComplexType>
            (0,NULL,NULL,i,0,sys_idx,tmp_eq->constant);
         mon_idx++;
      }
      for(int j=0; j<tmp_eq->n_mon; j++)
      {
         PolyMon<ComplexType,RealType>* tmp_mon = tmp_eq->mon[j];
         mon_set[mon_idx] = MonIdxSet<ComplexType>
           (tmp_mon->n_var,tmp_mon->pos,tmp_mon->exp,i,j,
            sys_idx,tmp_mon->coef);
         mon_idx++;
      }
      tmp_eq++;
   }
}

template <class ComplexType, class RealType>
MonIdxSet<ComplexType>* polysyshom_monidxset
 ( PolySys<ComplexType,RealType>& Target_Sys,
   PolySys<ComplexType,RealType>& Start_Sys,
   int& total_n_mon, int verbose )
{
   total_n_mon = 0;
   PolyEq<ComplexType,RealType>* tmp_eq = Target_Sys.eq_space;

   for(int i=0; i<Target_Sys.n_eq; i++)
   {
      total_n_mon += tmp_eq->n_mon;
      if(tmp_eq->constant.real != 0.0 || tmp_eq->constant.imag != 0.0)
         total_n_mon++;
      tmp_eq++;
   }
   int total_n_mon_start = total_n_mon;

   tmp_eq= Start_Sys.eq_space;
   for(int i=0; i<Start_Sys.n_eq; i++)
   {
      total_n_mon += tmp_eq->n_mon;
      if(tmp_eq->constant.real != 0.0 || tmp_eq->constant.imag != 0.0)
         total_n_mon++;
      tmp_eq++;
   }
   if(verbose > 0)
      std::cout << "total_n_mon = " << total_n_mon << std::endl;

   MonIdxSet<ComplexType>* mon_set = new MonIdxSet<ComplexType>[total_n_mon];

   polysys_mon_set<ComplexType,RealType>(Target_Sys, mon_set, 0);

   polysys_mon_set<ComplexType,RealType>
      (Start_Sys, mon_set+total_n_mon_start, 1);

   for(int i=0; i<total_n_mon; i++) mon_set[i].sorted();

   // sort all sets
   std::sort(mon_set, mon_set+total_n_mon);

   return mon_set;
}

template <class ComplexType, class RealType>
MonIdxSet<ComplexType>* polysys_monidxset
 ( PolySys<ComplexType,RealType>& Target_Sys, int& total_n_mon, int verbose )
{
   total_n_mon = 0;
   PolyEq<ComplexType,RealType>* tmp_eq = Target_Sys.eq_space;

   for(int i=0; i<Target_Sys.n_eq; i++)
   {
      total_n_mon += tmp_eq->n_mon;
      if(tmp_eq->constant.real != 0.0 || tmp_eq->constant.imag != 0.0)
      {
         total_n_mon++;
      }
      tmp_eq++;
   }
   if(verbose > 0)
      std::cout << "total_n_mon = " << total_n_mon << std::endl;

   MonIdxSet<ComplexType>* mon_set = new MonIdxSet<ComplexType>[total_n_mon];

   polysys_mon_set<ComplexType,RealType>(Target_Sys, mon_set, 0);

   for(int i=0; i<total_n_mon; i++)
   {
      mon_set[i].sorted();
   }
   // sort all sets
   std::sort(mon_set, mon_set+total_n_mon);

   return mon_set;
}

template <class ComplexType>
MonSet<ComplexType>* polysyshom_monset
 ( int total_n_mon, MonIdxSet<ComplexType>* mons, 
   int& n_constant, int& hom_n_mon, int& n_monset, int verbose )
{
   if(verbose > 0)
      std::cout << "total_n_mon = " << total_n_mon << std::endl;

   n_monset = 1;

   // Mark new location
   int* new_type = new int[total_n_mon];
   new_type[0] = 0;

   // Check monomial type
   MonIdxSet<ComplexType> tmp_mon = mons[0];

   for(int i=1; i<total_n_mon; i++)
   {
      if(mons[i] == tmp_mon)
      {
         new_type[i] = 0;
      }
      else
      {
         new_type[i-1] = 1;
         n_monset++;
         tmp_mon = mons[i];
      }
   }
   new_type[total_n_mon-1] = 1;

   if(verbose > 0)
      std::cout << "n_mon_type = " << n_monset << std::endl;

   int* n_mons = new int[n_monset];
   int tmp_n_mon = 0;
   int tmp_eq_idx = -1;
   int monset_idx = 0;
   for(int i=0; i<total_n_mon; i++)
   {
      if(tmp_eq_idx != mons[i].get_eq_idx())
      {
         tmp_n_mon++;
         tmp_eq_idx = mons[i].get_eq_idx();
      }
      if(new_type[i] == 1)
      {
         n_mons[monset_idx] = tmp_n_mon;
         tmp_n_mon = 0;
         tmp_eq_idx = -1;
         monset_idx++;
      }
   }

   MonSet<ComplexType>* hom_monset = new MonSet<ComplexType>[n_monset];
   int mon_idx = 0;
   for(int i=0; i<n_monset; i++)
   {
      hom_monset[i].copy_pos(mons[mon_idx]);
      EqIdxCoef<ComplexType>* tmp_eq_idx_coef
         = new EqIdxCoef<ComplexType>[n_mons[i]];
      for(int j=0; j<n_mons[i]; j++)
      {
         // merge by eq_idx
         if(mons[mon_idx].get_sys_idx() == 0)
         {
            int tmp_eq_idx = mons[mon_idx].get_eq_idx();
            mon_idx++;
            if(mons[mon_idx].get_sys_idx() == 1 
               && tmp_eq_idx == mons[mon_idx].get_eq_idx())
            {
               tmp_eq_idx_coef[j] = EqIdxCoef<ComplexType>
                  (tmp_eq_idx,mons[mon_idx-1].get_coef(),
                   mons[mon_idx].get_coef());
               mon_idx++;
            }
            else
            {
               tmp_eq_idx_coef[j] = EqIdxCoef<ComplexType>
                  (tmp_eq_idx,mons[mon_idx-1].get_coef(),0);
            }
         }
         else
         {
            int tmp_eq_idx = mons[mon_idx].get_eq_idx();
            tmp_eq_idx_coef[j] = EqIdxCoef<ComplexType>
               (tmp_eq_idx,mons[mon_idx].get_coef(),1);
            mon_idx++;
         }
      }
      hom_monset[i].update_eq_idx(n_mons[i],tmp_eq_idx_coef);
   }
   // Get number of constants
   n_constant = 0;
   if(hom_monset[0].get_n() == 0)
      n_constant = hom_monset[0].get_n_mon();

   hom_n_mon = 0;
   for(int i=0; i<n_monset; i++)
      hom_n_mon += hom_monset[i].get_n_mon();

   return hom_monset;
}

template <class ComplexType, class RealType>
MonSet<ComplexType>* hom_monset_generator
 ( PolySys<ComplexType,RealType>& Target_Sys,
   PolySys<ComplexType,RealType>& Start_Sys,
   int& n_monset, int& n_constant,
   int& total_n_mon, int verbose )
{
   int hom_n_mon;

   MonIdxSet<ComplexType>* mons = polysyshom_monidxset<ComplexType>
      (Target_Sys,Start_Sys,hom_n_mon,verbose);
   MonSet<ComplexType>* hom_monset = polysyshom_monset<ComplexType>
      (hom_n_mon,mons,n_constant,total_n_mon,n_monset,verbose);

   return hom_monset;
}

template <class ComplexType, class RealType>
MonSet<ComplexType>* single_monset_generator
 ( PolySys<ComplexType, RealType>& Target_Sys,
   int& n_monset, int& n_constant, int& total_n_mon, int verbose )
{
   int hom_n_mon;
   MonIdxSet<ComplexType>* mons
      = polysys_monidxset<ComplexType>(Target_Sys,hom_n_mon,verbose);
   MonSet<ComplexType>* hom_monset = polysyshom_monset<ComplexType>
      (hom_n_mon, mons,n_constant,total_n_mon,n_monset,verbose);

   return hom_monset;
}

template <class ComplexType, class RealType>
int* get_max_deg_base
 ( PolySys<ComplexType,RealType>& Target_Sys,
   PolySys<ComplexType,RealType>& Start_Sys )
{
   if(Target_Sys.eval_base==false && Start_Sys.eval_base==false) return NULL;

   int dim = Target_Sys.dim;
   int* max_deg_base = new int[dim];

   for(int var_idx=0; var_idx<dim; var_idx++)
      max_deg_base[var_idx] = max(Target_Sys.max_deg_base[var_idx],
                                  Target_Sys.max_deg_base[var_idx]);

   return max_deg_base;
}

template <class ComplexType, class RealType>
int* get_max_deg_base ( PolySys<ComplexType,RealType>& Target_Sys )
{
   if(Target_Sys.eval_base==false) return NULL;

   int dim = Target_Sys.dim;
   int* max_deg_base = new int[dim];

   for(int var_idx=0; var_idx<dim; var_idx++)
      max_deg_base[var_idx] = Target_Sys.max_deg_base[var_idx];

   return max_deg_base;
}

/* end helper functions, begin definition of class members */

template <class ComplexType, class RealType>
void CPUInstHom<ComplexType,RealType>::init
 ( PolySys<ComplexType,RealType>& Target_Sys,
   PolySys<ComplexType,RealType>& Start_Sys,
   int dim, int n_eq, int n_predictor, ComplexType alpha, int verbose )
{
   PED_hom = true;
   int n_constant;
   int total_n_mon;
   int n_monset;
   int n_mon;
   MonSet<ComplexType>* hom_monset = hom_monset_generator
      (Target_Sys,Start_Sys,n_monset,n_constant,total_n_mon,verbose);
   int* max_deg_base 
      = get_max_deg_base<ComplexType,RealType>(Target_Sys, Start_Sys);

   if(verbose > 0)
   {
      std::cout << "n_constant  = " << n_constant << std::endl;
      std::cout << "total_n_mon = " << total_n_mon << std::endl;
      std::cout << "n_monset    = " << n_monset << std::endl;
   }
   init(hom_monset,n_monset,n_constant,total_n_mon,dim,n_eq,
        n_predictor,alpha,max_deg_base,verbose);
}

template <class ComplexType, class RealType>
void CPUInstHom<ComplexType,RealType>::init
 ( PolySys<ComplexType,RealType>& Target_Sys,
   int dim, int n_eq, int n_predictor, ComplexType alpha, int verbose )
{
   PED_hom = false;
   int n_constant;
   int total_n_mon;
   int n_monset;
   int n_mon;
   MonSet<ComplexType>* hom_monset = single_monset_generator
      (Target_Sys,n_monset,n_constant,total_n_mon,verbose);
   int* max_deg_base = get_max_deg_base<ComplexType,RealType>(Target_Sys);

   if(verbose > 0)
   {
      if(Target_Sys.eval_base)
      {
         std::cout << "eval_base is true" << std::endl;
         std::cout << "max_deg_base :";
         for(int idx=0; idx<Target_Sys.dim; idx++)
         std::cout << " " << max_deg_base[idx];
         std::cout << endl;
      }
      else
         std::cout << "eval_base is false" << std::endl;

      std::cout << "n_constant  = " << n_constant << std::endl;
      std::cout << "total_n_mon = " << total_n_mon << std::endl;
      std::cout << "n_monset    = " << n_monset << std::endl;
   }

   init(hom_monset,n_monset,n_constant,total_n_mon,dim,n_eq,n_predictor,
        alpha,max_deg_base,verbose);
}

template <class ComplexType, class RealType>
void CPUInstHom<ComplexType,RealType>::init
 ( MonSet<ComplexType>* hom_monset, int n_monset, int n_constant,
   int total_n_mon, int dim, int n_eq, int n_predictor, ComplexType alpha,
   int* max_deg_base, int verbose )
{
   this->n_constant = n_constant;
   this->dim = dim;
   this->n_eq = n_eq;
   this->n_predictor = n_predictor;
   path_data.dim = dim;

   if(verbose > 0)
   {
      std::cout << "Generating CPU Instruction ..." << std::endl;
      std::cout << "           Coef Instruction ..." << std::endl;
   }
   CPU_inst_hom_coef.init(hom_monset,total_n_mon,n_monset,n_constant,alpha);
   if(verbose > 0)
      std::cout << "           Mon Instruction ..." << std::endl;

   CPU_inst_hom_mon.init
      (hom_monset,n_monset,total_n_mon,n_constant,max_deg_base);
   if(verbose > 0)
      std::cout << "           Sum Instruction ..." << std::endl;

   CPU_inst_hom_sum.init
      (hom_monset,n_monset,CPU_inst_hom_mon.mon_pos_start,dim,n_eq,
       n_constant,verbose);

   this->n_coef = CPU_inst_hom_coef.n_coef;

   if(MON_EVAL_METHOD == 1)
   {
      CPU_inst_hom_block.init(CPU_inst_hom_mon,warp_size,verbose);
      CPU_inst_hom_sum_block.init
         (hom_monset,n_monset,CPU_inst_hom_mon.mon_pos_start,dim,n_eq,
          n_constant,CPU_inst_hom_mon.n_mon_level[0],
          CPU_inst_hom_block.mon_pos_start_block,verbose);
   }
   if(verbose > 0)
      std::cout << "Generating CPU Instruction Finished" << std::endl;
}

template <class ComplexType, class RealType>
void CPUInstHom<ComplexType,RealType>::print()
{
   std::cout << "*************** Coef Instruction ********************"
             << std::endl;
   CPU_inst_hom_coef.print();
   std::cout << "*************** Mon Instruction ********************"
             << std::endl;
   CPU_inst_hom_mon.print();
   std::cout << "*************** Sum Instruction ********************"
             << std::endl;
   CPU_inst_hom_sum.print();
}

template <class ComplexType, class RealType>
void CPUInstHom<ComplexType,RealType>::init_workspace
 ( Workspace<ComplexType>& workspace_cpu, int verbose )
{
   int coef_size = CPU_inst_hom_coef.n_coef;
   int workspace_size = coef_size + CPU_inst_hom_mon.mon_pos_size;
   workspace_cpu.init(workspace_size,coef_size,n_constant,n_eq,dim,
       n_predictor, CPU_inst_hom_mon.max_deg_base,verbose);
}

template <class ComplexType, class RealType>
void CPUInstHom<ComplexType,RealType>::eval
 ( Workspace<ComplexType>& workspace_cpu,
   const ComplexType* sol, const ComplexType t, int reverse )
{
   if(PED_hom == true)
   {
      CPU_inst_hom_coef.eval(t, workspace_cpu.coef, reverse);
   }
   else
   {
      for(int coef_idx=0; coef_idx<CPU_inst_hom_coef.n_coef; coef_idx++)
      {
         workspace_cpu.coef[coef_idx] = CPU_inst_hom_coef.coef_orig[coef_idx];
      }
   }
   CPU_inst_hom_mon.eval(dim,sol,workspace_cpu.mon,workspace_cpu.coef,
                         workspace_cpu.deg_table);
   CPU_inst_hom_sum.eval(workspace_cpu.sum, workspace_cpu.matrix);
}

template <class ComplexType, class RealType>
void CPUInstHom<ComplexType,RealType>::update_alpha ( ComplexType alpha )
{
   CPU_inst_hom_coef.update_alpha(alpha);
}

template <class ComplexType, class RealType>
void CPUInstHom<ComplexType,RealType>::compare_path_cpu_gpu()
{
   std::cout << "n_step_cpu = " << path_data.n_step << std::endl
             << "n_step_gpu = " << path_data_gpu.n_step << std::endl;
   int n_step = min(path_data.n_step, path_data_gpu.n_step);
   for(int i=0; i<n_step; i++)
   {
      if(abs(path_data.steps[i]->t.real-path_data_gpu.steps[i]->t.real)>1E-6)
      {
         std::cout << i << " " << path_data.steps[i]->t.real
                   << " " << path_data.steps[i]->delta_t << std::endl;
         std::cout << i << " " << path_data_gpu.steps[i]->t.real
                   << " " << path_data_gpu.steps[i]->delta_t << std::endl;
         break;
      }
   }
}
