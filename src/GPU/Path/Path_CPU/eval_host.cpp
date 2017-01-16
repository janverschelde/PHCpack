/* eval_host.cpp, created on Dec 6, 2014 by yxc with edits by jv */

#include "eval_host.h"

void polysys_mon_set
 ( const PolySys& Target_Sys, MonIdxSet* mon_set, bool sys_idx )
{
   PolyEq* tmp_eq= Target_Sys.eq_space;

   int mon_idx = 0;
   for(int i=0; i<Target_Sys.n_eq; i++)
   {
      if(tmp_eq->constant.real != 0.0 || tmp_eq->constant.imag != 0.0)
      {
         mon_set[mon_idx] = MonIdxSet(0,NULL,NULL,i,0,sys_idx,
                                      tmp_eq->constant);
         mon_idx++;
      }
      for(int j=0; j<tmp_eq->n_mon; j++)
      {
         PolyMon* tmp_mon = tmp_eq->mon[j];
         mon_set[mon_idx] = MonIdxSet
           (tmp_mon->n_var,tmp_mon->pos,tmp_mon->exp,i,j,
            sys_idx,tmp_mon->coef);
         mon_idx++;
      }
      tmp_eq++;
   }
}

MonIdxSet* polysyshom_monidxset
 ( PolySys& Target_Sys, PolySys& Start_Sys, int& total_n_mon,
   int verbose )
{
   total_n_mon = 0;
   PolyEq* tmp_eq= Target_Sys.eq_space;

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

   MonIdxSet* mon_set = new MonIdxSet[total_n_mon];

   polysys_mon_set(Target_Sys, mon_set, 0);

   polysys_mon_set(Start_Sys, mon_set+total_n_mon_start, 1);

   for(int i=0; i<total_n_mon; i++) mon_set[i].sorted();

   // sort all sets
   std::sort(mon_set, mon_set+total_n_mon);

   return mon_set;
}

MonIdxSet* polysys_monidxset
 ( PolySys& Target_Sys, int& total_n_mon, int verbose )
{
   total_n_mon = 0;
   PolyEq* tmp_eq= Target_Sys.eq_space;

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

   MonIdxSet* mon_set = new MonIdxSet[total_n_mon];

   polysys_mon_set(Target_Sys, mon_set, 0);

   for(int i=0; i<total_n_mon; i++)
   {
      mon_set[i].sorted();
   }
   // sort all sets
   std::sort(mon_set, mon_set+total_n_mon);

   return mon_set;
}

MonSet* polysyshom_monset
 ( int total_n_mon, MonIdxSet* mons, 
   int& n_constant, int& hom_n_mon, int& n_monset, int verbose )
{
   if(verbose > 0)
      std::cout << "total_n_mon = " << total_n_mon << std::endl;

   n_monset = 1;

   // Mark new location
   int* new_type = new int[total_n_mon];
   new_type[0] = 0;

   // Check monomial type
   MonIdxSet tmp_mon = mons[0];

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

   MonSet* hom_monset = new MonSet[n_monset];
   int mon_idx = 0;
   for(int i=0; i<n_monset; i++)
   {
      hom_monset[i].copy_pos(mons[mon_idx]);
      EqIdxCoef* tmp_eq_idx_coef = new EqIdxCoef[n_mons[i]];
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
               tmp_eq_idx_coef[j] = EqIdxCoef
                  (tmp_eq_idx,mons[mon_idx-1].get_coef(),
                   mons[mon_idx].get_coef());
               mon_idx++;
            }
            else
            {
               tmp_eq_idx_coef[j] = EqIdxCoef
                  (tmp_eq_idx,mons[mon_idx-1].get_coef(),0);
            }
         }
         else
         {
            int tmp_eq_idx = mons[mon_idx].get_eq_idx();
            tmp_eq_idx_coef[j] = EqIdxCoef
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

MonSet* hom_monset_generator
 ( PolySys& Target_Sys, PolySys& Start_Sys, int& n_monset, int& n_constant,
   int& total_n_mon, int verbose )
{
   int hom_n_mon;

   MonIdxSet* mons = polysyshom_monidxset
      (Target_Sys,Start_Sys,hom_n_mon,verbose);
   MonSet* hom_monset = polysyshom_monset
      (hom_n_mon,mons,n_constant,total_n_mon,n_monset,verbose);

   return hom_monset;
}

MonSet* single_monset_generator
 ( PolySys& Target_Sys, int& n_monset, int& n_constant, int& total_n_mon,
   int verbose )
{
   int hom_n_mon;
   MonIdxSet* mons = polysys_monidxset(Target_Sys,hom_n_mon,verbose);
   MonSet* hom_monset = polysyshom_monset
      (hom_n_mon, mons,n_constant,total_n_mon,n_monset,verbose);

   return hom_monset;
}

int* get_max_deg_base ( PolySys& Target_Sys, PolySys& Start_Sys )
{
   if(Target_Sys.eval_base==false && Start_Sys.eval_base==false) return NULL;

   int dim = Target_Sys.dim;
   int* max_deg_base = new int[dim];

   for(int var_idx=0; var_idx<dim; var_idx++)
      max_deg_base[var_idx] = max(Target_Sys.max_deg_base[var_idx],
                                  Target_Sys.max_deg_base[var_idx]);

   return max_deg_base;
}

int* get_max_deg_base ( PolySys& Target_Sys )
{
   if(Target_Sys.eval_base==false) return NULL;

   int dim = Target_Sys.dim;
   int* max_deg_base = new int[dim];

   for(int var_idx=0; var_idx<dim; var_idx++)
      max_deg_base[var_idx] = Target_Sys.max_deg_base[var_idx];

   return max_deg_base;
}

void CPUInstHom::init
 ( PolySys& Target_Sys, PolySys& Start_Sys, int dim, int n_eq,
   int n_predictor, CT alpha, int verbose )
{
   PED_hom = true;
   int n_constant;
   int total_n_mon;
   int n_monset;
   int n_mon;
   MonSet* hom_monset = hom_monset_generator
      (Target_Sys,Start_Sys,n_monset,n_constant,total_n_mon,verbose);
   int* max_deg_base = get_max_deg_base(Target_Sys, Start_Sys);

   if(verbose > 0)
   {
      std::cout << "n_constant  = " << n_constant << std::endl;
      std::cout << "total_n_mon = " << total_n_mon << std::endl;
      std::cout << "n_monset    = " << n_monset << std::endl;
   }
   init(hom_monset,n_monset,n_constant,total_n_mon,dim,n_eq,
        n_predictor,alpha,max_deg_base,verbose);
}

void CPUInstHom::init
 ( PolySys& Target_Sys, int dim, int n_eq, int n_predictor, CT alpha,
   int verbose )
{
   PED_hom = false;
   int n_constant;
   int total_n_mon;
   int n_monset;
   int n_mon;
   MonSet* hom_monset = single_monset_generator
      (Target_Sys,n_monset,n_constant,total_n_mon,verbose);
   int* max_deg_base = get_max_deg_base(Target_Sys);

   if(verbose > 0)
   {
      std::cout << "n_constant  = " << n_constant << std::endl;
      std::cout << "total_n_mon = " << total_n_mon << std::endl;
      std::cout << "n_monset    = " << n_monset << std::endl;
   }

   init(hom_monset,n_monset,n_constant,total_n_mon,dim,n_eq,n_predictor,
        alpha,max_deg_base,verbose);
}

void CPUInstHom::init
 ( MonSet* hom_monset, int n_monset, int n_constant, int total_n_mon,
   int dim, int n_eq, int n_predictor, CT alpha, int* max_deg_base,
   int verbose )
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

void CPUInstHomCoef::init
 ( MonSet* hom_monset, int total_n_mon, int n_monset, int n_constant,
   CT alpha, int verbose )
{
   this->alpha = alpha;
   n_coef = total_n_mon;
   coef_orig = new CT[n_coef*2];
   CT* tmp_coef_orig = coef_orig;

   int constant_exist = 0;
   if(n_constant > 0) constant_exist = 1;

   // write start coefficient and target coefficient together
   for(int i=constant_exist; i<n_monset; i++)
      hom_monset[i].write_coef(tmp_coef_orig);

   if(n_constant > 0)
      hom_monset[0].write_coef(tmp_coef_orig);

   // write start coefficient and target coefficient seperately
   tmp_coef_orig = new CT[n_coef*2];
   for(int coef_idx=0; coef_idx<n_coef; coef_idx++)
   {
      tmp_coef_orig[coef_idx] = coef_orig[2*coef_idx];
      tmp_coef_orig[coef_idx+n_coef] = coef_orig[2*coef_idx+1];
   }
   delete[] coef_orig;
   coef_orig = tmp_coef_orig;
}

void CPUInstHomCoef::print()
{
   for(int i=0; i<n_coef; i++)
   {
      std::cout << i << std::endl
                << coef_orig[i]
                << coef_orig[i+n_coef]<< std::endl;
   }
}

void CPUInstHomCoef::eval ( const CT t, CT* coef, int reverse )
{
   CT one_minor_t(1.0- t.real, -t.imag);

   int k = 1;
   CT t_power_k = t;
   CT one_minor_t_power_k = one_minor_t;
   for(int i=1; i<k; i++)
   {
      t_power_k *= t;
      one_minor_t_power_k *= one_minor_t;
   }
   CT t0, t1;
   if(reverse == 0)
   {
      t0 = one_minor_t_power_k*alpha;
      t1 = t_power_k;
   }
   else
   {
      t0 = t_power_k*alpha;
      t1 = one_minor_t_power_k;
   }
   for(int i=0; i<n_coef; i++)
      coef[i] = coef_orig[i+n_coef]*t0 + coef_orig[i]*t1;
}

void CPUInstHomEq::init
 ( MonSet* hom_monset, int n_monset, int n_eq, int n_constant )
{
   this->n_eq = n_eq;
   n_mon_eq = NULL;
   eq_pos_start = NULL;
   n_mon_total = 0;
   mon_pos_start_eq = NULL;
   n_pos_total = 0;
   mon_pos_eq = NULL;
   coef = NULL;

   n_mon_eq = new int[n_eq];
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++) n_mon_eq[eq_idx] = 0;

   for(int set_idx=0; set_idx<n_monset; set_idx++)
   {
      if(hom_monset[set_idx].get_n()!=0)
      {
         int n_mon_set = hom_monset[set_idx].get_n_mon();
         for(int mon_idx=0; mon_idx<n_mon_set; mon_idx++)
         {
            int eq_idx = hom_monset[set_idx].get_eq_idx(mon_idx);
            n_mon_eq[eq_idx]++;
         }
      }
   }
   eq_pos_start = new int[n_eq];
   int tmp_pos=0;
   eq_pos_start[0] = 0;

   for(int eq_idx=0; eq_idx<n_eq-1; eq_idx++)
      eq_pos_start[eq_idx+1] = eq_pos_start[eq_idx] + (n_mon_eq[eq_idx]+1);

   n_mon_total = eq_pos_start[n_eq-1] + n_mon_eq[n_eq-1]+1;

   int* eq_mon_pos_size = new int[n_eq];
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++) eq_mon_pos_size[eq_idx] = 0;

   for(int set_idx=0; set_idx<n_monset; set_idx++)
   {
      if(hom_monset[set_idx].get_n()!=0)
      {
         int n_mon_set = hom_monset[set_idx].get_n_mon();
         for(int mon_idx=0; mon_idx<n_mon_set; mon_idx++)
         {
            int eq_idx = hom_monset[set_idx].get_eq_idx(mon_idx);
            eq_mon_pos_size[eq_idx] += hom_monset[set_idx].get_n() + 1;
         }
      }
   }
   n_pos_total = 0;

   for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
      n_pos_total += eq_mon_pos_size[eq_idx];

   mon_pos_start_eq = new int[n_mon_total];
   coef = new CT[n_mon_total*2];
   mon_pos_start_eq[0] = 0;

   int* eq_mon_idx_tmp = new int[n_eq];
   int* eq_mon_pos_tmp = new int[n_eq];
   eq_mon_pos_tmp[0] = 0;
   eq_mon_idx_tmp[0] = 0;
   for(int eq_idx=0; eq_idx<n_eq-1; eq_idx++)
   {
      eq_mon_pos_tmp[eq_idx+1] = eq_mon_pos_tmp[eq_idx]
                               + eq_mon_pos_size[eq_idx];
   }
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
   {
      eq_mon_idx_tmp[eq_idx] = eq_pos_start[eq_idx];
   }
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
   {
      mon_pos_start_eq[eq_pos_start[eq_idx]] = n_mon_eq[eq_idx];
      eq_mon_idx_tmp[eq_idx]++;
   }
   // init coef
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
   {
      coef[2*eq_pos_start[eq_idx]] = CT(0.0,0.0);
      coef[2*eq_pos_start[eq_idx]+1] = CT(0.0,0.0);
   }
   mon_pos_eq = new unsigned short[n_pos_total];
   for(int set_idx=0; set_idx<n_monset; set_idx++)
   {
      if(hom_monset[set_idx].get_n()!=0)
      {
         int n_mon_set = hom_monset[set_idx].get_n_mon();
         for(int mon_idx=0; mon_idx<n_mon_set; mon_idx++)
         {
            int eq_idx = hom_monset[set_idx].get_eq_idx(mon_idx);
            mon_pos_start_eq[eq_mon_idx_tmp[eq_idx]] = eq_mon_pos_tmp[eq_idx];
            CT* coef_tmp = coef + 2*eq_mon_idx_tmp[eq_idx];
            hom_monset[set_idx].write_coef(coef_tmp, mon_idx);
            int tmp_pos = eq_mon_pos_tmp[eq_idx];
            mon_pos_eq[tmp_pos] = hom_monset[set_idx].get_n();
            unsigned short* tmp_pointer = mon_pos_eq+tmp_pos+1;
            hom_monset[set_idx].write_pos(tmp_pointer);
            eq_mon_pos_tmp[eq_idx] += hom_monset[set_idx].get_n() + 1;
            eq_mon_idx_tmp[eq_idx]++;
         }
      }
      else
      {
         int n_mon_set = hom_monset[set_idx].get_n_mon();
         for(int mon_idx=0; mon_idx<n_mon_set; mon_idx++)
         {
             int eq_idx = hom_monset[set_idx].get_eq_idx(mon_idx);
             CT* coef_tmp = coef + 2*eq_pos_start[eq_idx];
             hom_monset[set_idx].write_coef(coef_tmp, mon_idx);
         }
      }
   }
}

void CPUInstHomEq::print()
{
   std::cout << "CPUInstHomMonEq" << std::endl;
   std::cout << "n_eq = " << n_eq << std::endl;
   std::cout << "n_mon       = " << n_mon_total-n_eq << std::endl;
   std::cout << "n_mon+n_eq  = " << n_mon_total << std::endl;
   std::cout << "n_pos       = " << n_pos_total-(n_mon_total-n_eq) << std::endl;
   std::cout << "n_pos+n_mon = " << n_pos_total << std::endl;
}

void CPUInstHom::print()
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

void CPUInstHom::init_workspace ( Workspace& workspace_cpu )
{
   int coef_size = CPU_inst_hom_coef.n_coef;
   int workspace_size = coef_size + CPU_inst_hom_mon.mon_pos_size;
   workspace_cpu.init(workspace_size,coef_size,n_constant,n_eq,dim,
       n_predictor, CPU_inst_hom_mon.max_deg_base);
}

void CPUInstHom::eval
 ( Workspace& workspace_cpu, const CT* sol, const CT t, int reverse )
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

void CPUInstHom::update_alpha(CT alpha)
{
   CPU_inst_hom_coef.update_alpha(alpha);
}

void CPUInstHom::compare_path_cpu_gpu()
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

void CPUInstHomCoef::update_alpha(CT alpha)
{
   if(alpha.real == 0 && alpha.imag == 0)
   {
      int r = rand();
      T1 tmp = T1(r);
      this->alpha = CT(sin(tmp),cos(tmp));
   }
   else
      this->alpha = alpha;
}

void CPUInstHomMon::init
 ( MonSet* hom_monset, int n_monset, int total_n_mon, int n_constant,
   int* max_deg_base, int verbose )
{
   this->max_deg_base = max_deg_base;
   int max_n_var = hom_monset[n_monset-1].get_n();
   level = 1;
   for(int i=1; i<max_n_var; i<<=1) level++;
   n_mon_level = new int[level];
   for(int i=0; i<level; i++) n_mon_level[i] = 0;
   mon_pos_size = 0;
   n_mon = total_n_mon - n_constant;
   if(verbose > 0)
      std::cout << "n_constant = " << n_constant << std::endl;
   int constant_exist;
   if(n_constant == 0)
      constant_exist = 0;
   else
      constant_exist = 1;
   // Write monomial start position
   int tmp_level = 0;
   int tmp_level_size = 1;
   mon_pos_start = new int[n_mon];

   int mon_idx = 0;
   for(int i=constant_exist; i<n_monset; i++)
   {
      int tmp_n = hom_monset[i].get_n();
      int tmp_n_mon = hom_monset[i].get_n_mon();
      if(tmp_n > 1)
      {
         flops_multiple += tmp_n_mon*3*(tmp_n-1);
      }
      else
      {
         flops_multiple += tmp_n_mon;
      }
      for(int j=0; j<tmp_n_mon; j++)
      {
         mon_pos_start[mon_idx++] = mon_pos_size;
         mon_pos_size += tmp_n+1;
      }
      while(tmp_level_size < tmp_n)
      {
         tmp_level_size *= 2;
         tmp_level++;
      }
      n_mon_level[tmp_level]+=tmp_n_mon;
   }
   // Write position instruction
   mon_pos = new unsigned short[mon_pos_size];
   mon_exp = new unsigned short[mon_pos_size];

   unsigned short* tmp_pos = mon_pos;
   unsigned short* tmp_exp = mon_exp;
   for(int i=constant_exist; i<n_monset; i++)
   {
      // Number of variable each term
      int tmp_n_mon = hom_monset[i].get_n_mon();
      for(int j=0; j<tmp_n_mon; j++)
      {
         int n_var = hom_monset[i].get_n();
         *tmp_pos++ = n_var;
         int base_start = hom_monset[i].get_base_start();
         *tmp_exp++ = base_start;
         // Write position
         hom_monset[i].write_pos(tmp_pos);
         hom_monset[i].write_exp(tmp_exp);
      }
   }
   n_mon_base_start = n_mon;
   for(int mon_idx=0; mon_idx<n_mon; mon_idx++)
   {
      int tmp_mon_pos_start = mon_pos_start[mon_idx];
      int tmp_n_mon = mon_pos[tmp_mon_pos_start];
      if(mon_exp[tmp_mon_pos_start+tmp_n_mon]>1)
      {
         n_mon_base_start = mon_idx;
         break;
      }
   }
   if(verbose > 0)
   {
      std::cout << "*** flops_multiple = " << flops_multiple << std::endl;
      std::cout << "*** flops          = "
                << flops_multiple*4+flops_multiple*2 << std::endl;
   }
}

void CPUInstHomMon::eval
 ( int dim, const CT* x_val, CT* mon, CT* coef, CT** deg_table )
{
   if(deg_table != NULL)
   {
      eval_deg_table(dim, x_val, deg_table);
      eval_base(deg_table, coef);
      for(int j=0; j<n_mon_level[0]; j++)
      {
         int tmp_idx = mon_pos_start[j];
         cpu_speel_with_base0(x_val,mon_pos+tmp_idx,mon_exp+tmp_idx,
            mon+tmp_idx, coef[j]);
      }
      for(int j=n_mon_level[0]; j<n_mon; j++)
      {
         int tmp_idx = mon_pos_start[j];
         cpu_speel_with_base(x_val,mon_pos+tmp_idx,mon_exp+tmp_idx,
            mon+tmp_idx, coef[j]);
      }
   }
   else
   {
      for(int j=0; j<n_mon_level[0]; j++)
      {
         int tmp_idx = mon_pos_start[j];
         cpu_speel0(x_val,mon_pos+tmp_idx,mon+tmp_idx,coef[j]);
      }
      for(int j=n_mon_level[0]; j<n_mon; j++)
      {
         int tmp_idx = mon_pos_start[j];
         cpu_speel(x_val,mon_pos+tmp_idx,mon+tmp_idx,coef[j]);
      }
   }
}


void CPUInstHomMon::eval_deg_table
 ( int dim, const CT* x_val, CT** deg_table )
{
   for(int var_idx=0; var_idx<dim; var_idx++)
   {
      if(max_deg_base[var_idx]>0)
      {
         CT* tmp_deg_table = deg_table[var_idx];
         CT tmp_var = x_val[var_idx];
         tmp_deg_table[0] = tmp_var;
         for(int deg_idx=1; deg_idx<max_deg_base[var_idx]; deg_idx++)
         {
            tmp_deg_table[deg_idx] = tmp_deg_table[deg_idx-1]*tmp_var;
         }
      }
   }
}

void CPUInstHomMon::eval_base ( CT** deg_table, CT* coef )
{
   for(int mon_idx=n_mon_base_start; mon_idx<n_mon; mon_idx++)
   {
      int tmp_idx = mon_pos_start[mon_idx];
      unsigned short* tmp_mon_pos = mon_pos + tmp_idx;
      unsigned short* tmp_mon_exp = mon_exp + tmp_idx;
      int tmp_var_start = *tmp_mon_exp++;
      int tmp_n_var = *tmp_mon_pos++;
      if(tmp_var_start < tmp_n_var)
      {
         CT tmp_val = coef[mon_idx];
         for(int var_idx=tmp_var_start; var_idx<tmp_n_var; var_idx++)
            tmp_val *= deg_table[tmp_mon_pos[var_idx]][tmp_mon_exp[var_idx]-2];
         coef[mon_idx] = tmp_val;
      }
   }
}

void CPUInstHomMon::print()
{
   std::cout << "level = " << level << std::endl;
   for(int i=0; i<level; i++)
   {
      std::cout << i << " " << n_mon_level[i] << std::endl;
   }
   std::cout << "n_mon = " << n_mon << std::endl;
   // Print monomial with position
   for(int i=0; i<n_mon; i++)
   {
      int tmp_n = mon_pos[mon_pos_start[i]];
      int tmp_exp_start = mon_exp[mon_pos_start[i]];
      std::cout << i << " n = " << tmp_n << ": ";
      for(int j=0; j<tmp_n; j++)
         std::cout << mon_pos[mon_pos_start[i]+j+1] << ", ";
      std::cout << " exp_start = " << tmp_exp_start << ": ";
      for(int j=0; j<tmp_n; j++)
         std::cout << mon_exp[mon_pos_start[i]+j+1] << ", ";
      std::cout << " mon_pos_start = " << mon_pos_start[i] << std::endl;
   }
}

void CPUInstHomMonBlock::init ( CPUInstHomMon& orig, int BS, int verbose )
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

void CPUInstHomMonBlock::print()
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

void CPUInstHomSumBlock::init
 ( MonSet* hom_monset, int n_monset, const int* mon_pos_start, int dim,
   int n_eq, int n_constant, int n_mon_single, int* mon_pos_start_block,
   int verbose )
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
   MonSet* tmp_hom_monset = hom_monset;
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

void CPUInstHomSumBlock::eval ( CT* sum, CT* matrix )
{
   for(int i=0; i<n_sum; i++)
   {
      int tmp_start = sum_pos_start[i];
      int* tmp_pos = sum_pos+tmp_start;
      int tmp_n = *(tmp_pos++);
      CT tmp = sum[*tmp_pos++];

      for(int j=1; j<tmp_n; j++) tmp += sum[*tmp_pos++];

      matrix[*tmp_pos] = tmp;
   }
}

void CPUInstHomSumBlock::print()
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

void CPUInstHomSum::init
 ( MonSet* hom_monset, int n_monset, const int* mon_pos_start,
   int dim, int n_eq, int n_constant, int verbose )
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
   MonSet* tmp_hom_monset = hom_monset;

   for(int monset_idx=0; monset_idx<n_monset; monset_idx++)
   {
      for(int j=0; j<tmp_hom_monset->get_n_mon(); j++)
      {
         int tmp_eq_idx = tmp_hom_monset->get_eq_idx(j);
         n_sums[tmp_eq_idx][dim] += 1;
         for(int k=0; k<tmp_hom_monset->get_n(); k++)
         {
            n_sums[tmp_eq_idx][tmp_hom_monset->get_pos(k)] += 1;
         }
      }
      tmp_hom_monset++;
   }
   // Step 2: Count number of sums of certain number of terms
   //         total number of terms to sum
   //         max number of terms to sum
   //         max number of terms
   int max_n_sums = 0;
   for(int i=0; i<n_eq; i++)
   {
      for(int j=0; j<dim+1; j++)
      {
         if(n_sums[i][j] > max_n_sums)
         {
            max_n_sums = n_sums[i][j];
         }
         // std::cout << n_sums[i][j] << " ";
      }
      // std::cout << std::endl;
    }
    if(verbose > 0)
    {
       std::cout << "max_n_sums = " << max_n_sums << std::endl;
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
    // zeros sums
    n_sum_zero = n_sums_count[0];
    sum_zeros = new int[n_sum_zero];
    int sum_zeros_idx = 0;
    for(int var_idx=0; var_idx<dim+1; var_idx++)
    {
       for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
       {
          if(n_sums[eq_idx][var_idx] == 0)
          {
             sum_zeros[sum_zeros_idx++]= eq_idx+var_idx*n_eq;
          }
       }
    }
    if(verbose > 0)
    {
       std::cout << "n_sum_zero = " << n_sum_zero << std::endl;
    }
    // total number of sums
    n_sum = (dim+1)*n_eq - n_sum_zero;
    // Split one sum into multiple threads
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
    }
    // sum start
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
    // Start pos of sums
    int* n_sums_start = new int[max_n_sums+1];
    n_sums_start[0] = 0;
    n_sums_start[1] = 0;
    for(int i=2; i<max_n_sums+1; i++)
    {
       n_sums_start[i] = n_sums_start[i-1] + n_sums_count[i-1]*(1+i);
    }
    sum_pos_size = n_sums_start[max_n_sums]
                 + n_sums_count[max_n_sums]*(2+max_n_sums);

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
             // std::cout << i << " " << j << " " 
             //           << "tmp_start = " << tmp_start << std::endl;
             sum_pos[tmp_start] = tmp_n;
             sum_pos_start_matrix[i][j] = tmp_start+1;
             sum_pos[tmp_start+tmp_n+1] = j*n_eq + i;
             // sum_pos[tmp_start+tmp_n+1] = i*(dim+1) + j;
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
       // std::cout << *tmp_hom_monset;
       for(int j=0; j<tmp_hom_monset->get_n_mon(); j++)
       {
          int tmp_pos = mon_pos_start[mon_idx++]+n_constant;
          int tmp_eq_idx = tmp_hom_monset->get_eq_idx(j);
          // Value
          sum_pos[sum_pos_start_matrix[tmp_eq_idx][dim]] = tmp_pos;
          tmp_pos++;
          sum_pos_start_matrix[tmp_eq_idx][dim]++;
          n_sums[tmp_eq_idx][dim] += 1;
          // Derivative
          for(int k=0; k<tmp_hom_monset->get_n(); k++)
          {
             sum_pos[sum_pos_start_matrix[tmp_eq_idx]
                    [tmp_hom_monset->get_pos(k)]] = tmp_pos;
             tmp_pos++;
             sum_pos_start_matrix[tmp_eq_idx][tmp_hom_monset->get_pos(k)]++;
          }
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

void CPUInstHomSum::eval ( CT* sum, CT* matrix )
{
   // std::cout << "n_sum = " << n_sum << std::endl;

   for(int i=0; i<n_sum_zero; i++)
   {
      matrix[sum_zeros[i]].init(0.0,0.0);
   }
   for(int i=0; i<n_sum; i++)
   {
      int tmp_start = sum_pos_start[i];
      int* tmp_pos = sum_pos+tmp_start;
      int tmp_n = *(tmp_pos++);
      // std::cout << "i = " << i << " n = " << tmp_n << ", ";
      // std::cout << *tmp_pos << " " << sum[*tmp_pos] << " ";
      CT tmp = sum[*tmp_pos++];
      for(int j=1; j<tmp_n; j++)
      {
         // std::cout << *tmp_pos << " " << sum[*tmp_pos] << " ";
         tmp += sum[*tmp_pos++];
      }
      matrix[*tmp_pos] = tmp;
      // std::cout << i << " sum_pos_start = " << tmp_start
      //           << " output = " << *tmp_pos << " " << tmp;
   }
}

void CPUInstHomSum::print()
{
   std::cout << "n_sum = " << n_sum << std::endl;
   std::cout << "sum_pos_size = " << sum_pos_size << std::endl;
   for(int i=0; i<n_sum; i++)
   {
      int tmp_start = sum_pos_start[i];
      int* tmp_pos = sum_pos+tmp_start;
      int tmp_n = *(tmp_pos++);
      std::cout << "sum_idx = " << i << " n = " << tmp_n << ", ";
      for(int j=0; j<tmp_n; j++)
      {
         // std::cout << *tmp_pos++ << " ";
      }
      std::cout << "   sum_pos_start = " << tmp_start 
                << " output = " << sum_pos[tmp_start+tmp_n+1] << std::endl;
   }
}
