template <class ComplexType, class RealType>
int manytrack
 ( int verbose, double regamma, double imgamma, Parameter pars, int prec,
   PolySys<ComplexType,RealType>& p, PolySys<ComplexType,RealType>& q,
   PolySolSet<ComplexType,RealType>& s )
{
   double tpred,teval,tmgs;
   bool success;
   ComplexType* sol = s.get_sol(0);
   ComplexType alpha,t;
   CPUInstHom<ComplexType,RealType> cpu_inst_hom;
   Workspace<ComplexType> workspace_cpu;

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
      s.change_solt(path_idx,workspace_cpu.x_last,workspace_cpu.t_last);
      cout << "t : " << *workspace_cpu.t_last;
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

template <class ComplexType, class RealType>
int crewmanytrack
 ( int crewsize,
   int verbose, double regamma, double imgamma, Parameter pars, int prec,
   PolySys<ComplexType,RealType>& p, PolySys<ComplexType,RealType>& q,
   PolySolSet<ComplexType,RealType>& s )
{
   double tpred,teval,tmgs;
   bool success;
   ComplexType alpha,t;
   alpha = ComplexType(regamma,imgamma);

   CPUInstHom<ComplexType,RealType> *cpu_inst_hom;
   cpu_inst_hom = (CPUInstHom<ComplexType,RealType>*)calloc
      (crewsize,sizeof(CPUInstHom<ComplexType,RealType>));

   Workspace<ComplexType> *workspace_cpu;
   workspace_cpu = (Workspace<ComplexType>*)calloc
      (crewsize,sizeof(Workspace<ComplexType>));

   int fail,n_path;
   n_path = s.n_sol;

   for(int idxworker=0; idxworker<crewsize; idxworker++)
   {
      cpu_inst_hom[idxworker].init(p,q,p.dim,p.n_eq,1,alpha,verbose);
      cpu_inst_hom[idxworker].init_workspace(workspace_cpu[idxworker]);
   }

   ComplexType* points[n_path]; // stores solution coordinates

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      points[path_idx] = s.get_sol(path_idx);
   }
   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      cout << "The coordinates of start solution " << path_idx
           << " :" << endl;
      for(int k=0; k<p.dim; k++) cout << points[path_idx][k];
   }

   typedef struct
   {
      int label;
      CPUInstHom<ComplexType,RealType>* homotopy;
      Workspace<ComplexType>* workdata;
      JobQueue *jobs;
   }
   WorkItem;

   Worker *crew = (Worker*)calloc(crewsize,sizeof(Worker));
   WorkItem *wit = (WorkItem*) calloc(crewsize,sizeof(WorkItem));
   JobQueue jobpaths(n_path);

   for(int idxworker=0; idxworker<crewsize; idxworker++)
   {
      wit[idxworker].label = idxworker;
      wit[idxworker].jobs = &jobpaths;
      wit[idxworker].workdata = &workspace_cpu[idxworker];
      wit[idxworker].homotopy = &cpu_inst_hom[idxworker];
   }
   
   return 0;
}
