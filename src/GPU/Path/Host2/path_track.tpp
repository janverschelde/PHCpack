// The file path_track.tpp defines the functions with prototypes 
// in path_track.h.

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
class WorkItem
{
   public:

      int label;
      int dimsol;
      int prec;
      Parameter* pars;
      CPUInstHom<ComplexType,RealType>* homotopy;
      Workspace<ComplexType>* workdata;
      JobQueue* jobs;        // manages job queue
      ComplexType** points;  // pointer to start solutions
};

template <class ComplexType, class RealType>
void* track_paths ( void* args )
/*
 * A worker in a thread crew tracks path, as defined by the job queue.
 * The args is of type WorkItem. */
{
   WorkItem<ComplexType, RealType> *wit
      = (WorkItem<ComplexType, RealType>*) (args);
   JobQueue *jobs = (JobQueue*) wit->jobs;

   double tpred,teval,tmgs;
   bool success;
   ComplexType t;

   const int workid = wit->label;

   cout << "Worker " << workid << " entered job loop ..." << endl;

   while(1)
   {
      int nextjob = jobs->get_next_job();

      if(nextjob == -1) break;

      cout << "Worker " << workid
           << " has job : " << nextjob
           << " with work " << jobs->work[nextjob] << endl;

      cout << "The start solution at job " << nextjob << " :" << endl;
      for(int k=0; k<wit->dimsol; k++)
      {
         cout << k << " :";
         cout << setw(wit->prec+8) << scientific << setprecision(wit->prec)
              << wit->points[nextjob][k];
      }

      // sleep(jobs->work[nextjob]);
      cout << "calling the path tracker ..." << endl;

      t = ComplexType(0.0,0.0);
      wit->workdata[workid].init_x_t_idx();
      wit->workdata[workid].update_x_t_value(wit->points[nextjob],t);
      wit->workdata[workid].update_x_t_idx();
      wit->workdata[workid].path_idx = nextjob;

      success = path_tracker
         (wit->workdata[workid],wit->homotopy[workid],*(wit->pars),
          tpred,teval,tmgs,0,0);

      cout << "The end solution at job " << nextjob << " :" << endl;
      for(int k=0; k<wit->dimsol; k++)
      {
         wit->points[nextjob][k] = wit->workdata[workid].x_last[k];
         cout << k << " :";
         cout << setw(wit->prec+8) << scientific << setprecision(wit->prec)
              << wit->points[nextjob][k];
      }
      cout << "done with call to path_tracker." << endl;
   }

   return NULL;
}

template <class ComplexType, class RealType>
int track_all_paths ( int crewsize, WorkItem<ComplexType,RealType>* wit )
/*
 * Defines a crew of workers, as many as the value of crewsize,
 * to track all paths, as defined in the array wit.
 * There are as many work items in wit as there are workers. */
{
   Worker *crew = (Worker*)calloc(crewsize,sizeof(Worker));

   cout << "Assigning unique identification to each worker ..." << endl;
   for(int i=0; i<crewsize; i++) crew[i].idn = i;
   cout << "Writing the data at each worker ..." << endl;
   for(int i=0; i<crewsize; i++) crew[i].write();
   cout << "Starting the work crew ..." << endl;
   for(int i=0; i<crewsize; i++)
   {
      cout << "Launching worker " << wit[i].label << endl;

      crew[i].work(&track_paths<ComplexType,RealType>,(void*)&wit[i]);
   }
   cout << "Waiting for the workers to finish ..." << endl;
   for(int i=0; i<crewsize; i++) crew[i].join();

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

   ComplexType** points; // stores solution coordinates
   points = (ComplexType**) calloc(n_path,sizeof(ComplexType*));

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

   Worker *crew = (Worker*)calloc(crewsize,sizeof(Worker));
   WorkItem<ComplexType,RealType> *wit
      = (WorkItem<ComplexType,RealType>*)
            calloc(crewsize,sizeof(WorkItem<ComplexType,RealType>));
   JobQueue jobpaths(n_path);

   for(int idxworker=0; idxworker<crewsize; idxworker++)
   {
      wit[idxworker].label = idxworker;
      wit[idxworker].dimsol = p.dim;
      wit[idxworker].prec = prec;
      wit[idxworker].pars = &pars;
      wit[idxworker].jobs = &jobpaths;
      wit[idxworker].workdata = &workspace_cpu[idxworker];
      wit[idxworker].homotopy = &cpu_inst_hom[idxworker];
      wit[idxworker].points = points;
   }

   fail = track_all_paths<ComplexType,RealType>(crewsize,wit);

   return fail;
}
