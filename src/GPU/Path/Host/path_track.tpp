template <class ComplexType, class RealType>
int manytrack
 ( int verbose, double regamma, double imgamma, int prec,
   PolySys<ComplexType,RealType>& p, PolySys<ComplexType,RealType>& q,
   PolySolSet<ComplexType,RealType>& s )
{
   double tpred,teval,tmgs;
   bool success;
   ComplexType* sol = s.get_sol(0);
   ComplexType alpha,t;
   CPUInstHom<ComplexType,RealType> cpu_inst_hom;
   Workspace<ComplexType> workspace_cpu;
   Parameter pars(prec);

   if(verbose > 0)
   {
      cout << "The first solution on input :" << endl;
      for(int k=0; k<p.dim; k++)
      {
         cout << k << " :";
         cout << setw(prec+8) << scientific << setprecision(prec) << sol[k];
      }
   }
   int fail,n_path;
   n_path = s.n_sol;

   alpha = ComplexType(regamma,imgamma);
   cpu_inst_hom.init(p,q,p.dim,p.n_eq,1,alpha,verbose);
   cpu_inst_hom.init_workspace(workspace_cpu);

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      if(verbose > 0) cout << "tracking path " << path_idx << endl;
      ComplexType* sol0 = s.get_sol(path_idx);
      ComplexType* sol_new = NULL;
      t = ComplexType(0.0,0);
      workspace_cpu.init_x_t_idx();
      workspace_cpu.update_x_t_value(sol0,t);
      workspace_cpu.update_x_t_idx();
      workspace_cpu.path_idx = path_idx;
      tpred = 0.0; teval = 0.0; tmgs = 0.0;
      if(verbose > 0) cout << "calling path_tracker ..." << endl;
      success = path_tracker(workspace_cpu,cpu_inst_hom,pars,
                             tpred,teval,tmgs,0,verbose);
      if(verbose > 0) cout << "done with call to path_tracker." << endl;
      s.change_sol(path_idx,workspace_cpu.x_last);
      if(verbose > 0)
      {
         cout.precision(prec);
         cout << "The solution " << path_idx << endl;
         for(int k=0; k<p.dim; k++)
            cout << k << " :" << setw(prec+8) << workspace_cpu.x_last[k];
      }
   }
   return 0;
}
