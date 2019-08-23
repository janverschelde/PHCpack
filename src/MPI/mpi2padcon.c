/* mpi2padcon.c defines an interactive message passing parallel program
 * to track all solution paths in a static load balancing scheme with
 * Pade approximants as predictor.
 * The user is prompted for the precision, can adjust all parameters,
 * and is asked to provide the file names for the target, start system,
 * and start solutions.  All output is written to file.  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"
#include "jump_track.h"
#include "parallel_phcpack.h"
#include "padcon.h"

extern void adainit ( void );
extern void adafinal ( void );

int prompt_for_precision ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for the working precision,
 *   which is either double, double double, or quad double,
 *   and returns 0, 1, or 2 respectively. */

int prompt_for_verbose ( void );
/*
 * DESCRIPTION :
 *   Ask if extra output is needed during path tracking.
 *   Returns 0 if no, returns 1 otherwise. */

int parameters_broadcast ( int myid, int nbrp, int verbose );
/*
 * DESCRIPTION :
 *   After tuning the parameters for Pade continuation,
 *   the root node broadcasts all values for the parameters
 *   to all nbrp processes.
 *   If verbose, then extra data is written to screen for testing. */

int standard_track_paths
 ( int myid, int nbrp, int nbc, char *outfile, int verbose, int mhom,
   int nvr, int *idz );
/*
 * DESCRIPTION :
 *   For the stored target, start system and start solutions,
 *   worker node with myid tracks all paths in double precision
 *   and writes to the suffixed output file if verbose.
 *   This function is called by standard_run().
 *
 * ON ENTRY :
 *   myid     identification number of the current process,
 *            assumed is that myid > 1
 *   nbrp     number of processes;
 *   nbc      number of characters in the string outfile;
 *   outfile  is the name for the output file, only valid if myid == 0;
 *   verbose  if zero, then the run will be silent,
 *            otherwise extra output will be written to file;
 *   mhom     if zero, then path tracking happens in affine coordinates,
 *            if one, then 1-homogeneous coordinates are applied;
 *            if > 1, then m-homogenization applies, where m = mhom;
 *   nvr      number of (affine) variables;
 *   idz      index representation of the permutation if mhom > 1. */

int dobldobl_track_paths
 ( int myid, int nbrp, int nbc, char *outfile, int verbose, int mhom,
   int nvr, int *idz );
/*
 * DESCRIPTION :
 *   For the stored target, start system and start solutions,
 *   worker node with myid tracks all paths in double double precision
 *   and writes to the suffixed output file if verbose.
 *   This function is called by dobldobl_run().
 *
 * ON ENTRY :
 *   myid     identification number of the current process,
 *            assumed is that myid > 1
 *   nbrp     number of processes;
 *   nbc      number of characters in the string outfile;
 *   outfile  is the name for the output file, only valid if myid == 0;
 *   verbose  if zero, then the run will be silent,
 *            otherwise extra output will be written to file;
 *   mhom     if zero, then path tracking happens in affine coordinates,
 *            if one, then 1-homogeneous coordinates are applied;
 *            if > 1, then m-homogenization applies, where m = mhom;
 *   nvr      number of (affine) variables;
 *   idz      index representation of the permutation if mhom > 1. */

int quaddobl_track_paths
 ( int myid, int nbrp, int nbc, char *outfile, int verbose, int mhom,
   int nvr, int *idz );
/*
 * DESCRIPTION :
 *   For the stored target, start system and start solutions,
 *   worker node with myid tracks all paths in quad double precision
 *   and writes to the suffixed output file if verbose.
 *   This function is called by quaddobl_run().
 *
 * ON ENTRY :
 *   myid     identification number of the current process,
 *            assumed is that myid > 1
 *   nbrp     number of processes;
 *   nbc      number of characters in the string outfile;
 *   outfile  is the name for the output file, only valid if myid == 0;
 *   verbose  if zero, then the run will be silent,
 *            otherwise extra output will be written to file;
 *   mhom     if zero, then path tracking happens in affine coordinates,
 *            if one, then 1-homogeneous coordinates are applied;
 *            if > 1, then m-homogenization applies, where m = mhom;
 *   nvr      number of (affine) variables;
 *   idz      index representation of the permutation if mhom > 1. */

int standard_run
 ( int myid, int nbrp, int nbc, char *outfile, int verbose );
/*
 * DESCRIPTION :
 *   Prompts the user for a target system, start system with solutions,
 *   allows then to tune the parameters interactive, and then runs
 *   the Pade continuation in standard double precision.
 *    
 * ON ENTRY :
 *   myid     identification number of the current process;
 *   nbrp     number of processes;
 *   nbc      number of characters in the string outfile;
 *   outfile  is the name for the output file, only valid if myid == 0;
 *   verbose  if zero, then the run will be silent,
 *            otherwise extra messages will be written to screen. */

int dobldobl_run
 ( int myid, int nbrp, int nbc, char *outfile, int verbose );
/*
 * DESCRIPTION :
 *   Prompts the user for a target system, start system with solutions,
 *   allows then to tune the parameters interactive, and then runs
 *   the Pade continuation in double double precision.
 *    
 * ON ENTRY :
 *   myid     identification number of the current process;
 *   nbrp     number of processes;
 *   nbc      number of characters in the string outfile;
 *   outfile  is the name for the output file, only valid if myid == 0;
 *   verbose  if zero, then the run will be silent,
 *            otherwise extra messages will be written to screen. */

int quaddobl_run
 ( int myid, int nbrp, int nbc, char *outfile, int verbose );
/*
 * DESCRIPTION :
 *   Prompts the user for a target system, start system with solutions,
 *   allows then to tune the parameters interactive, and then runs
 *   the Pade continuation in quad double precision.
 *    
 * ON ENTRY :
 *   myid     identification number of the current process;
 *   nbrp     number of processes;
 *   nbc      number of characters in the string outfile;
 *   outfile  is the name for the output file, only valid if myid == 0;
 *   verbose  if zero, then the run will be silent,
 *            otherwise extra messages will be written to screen. */

int main ( int argc, char *argv[] )
{
   const int verbose = 0; // set to 1 for more info
   int myid,numprocs,precision,fail,ns;
   char *filename;
 
   adainit();
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);

   if(verbose > 0)
   {
      printf("\nHello from process %d!\n", myid);
      MPI_Barrier(MPI_COMM_WORLD);
   }
   if(myid == 0) precision = prompt_for_precision();

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Bcast(&precision,1,MPI_INT,0,MPI_COMM_WORLD);
   if(verbose > 0)
      printf("\nProcess %d has %d as precision.\n", myid, precision);

   MPI_Barrier(MPI_COMM_WORLD);

   filename = (char*)calloc(80,sizeof(char));
   if(myid == 0)
   {
      char nl;

      printf("\nReading the name of the output file ...");
      printf("\nGive a string of characters : "); scanf("%s",filename);
      scanf("%c",&nl); /* skip newline symbol for next reading ...*/
      if(verbose > 0) printf("\nThe output file is \"%s\".\n", filename);
      ns = strlen(filename);
      fail = define_output_file_with_string(ns,filename);
   }
   MPI_Bcast(&ns,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(filename,ns,MPI_CHAR,0,MPI_COMM_WORLD);

   if(verbose > 0)
      printf("\nNode %d has filename \"%s\".\n", myid, filename);

   MPI_Barrier(MPI_COMM_WORLD);

   switch(precision)
   {
      case 0: fail = standard_run(myid,numprocs,ns,filename,verbose);
              break;
      case 1: fail = dobldobl_run(myid,numprocs,ns,filename,verbose);
              break;
      case 2: fail = quaddobl_run(myid,numprocs,ns,filename,verbose);
              break;
      default: printf("Invalid choice of precision.\n");
   }

   MPI_Finalize();
   adafinal();

   return 0;
}

int prompt_for_precision ( void )
{
   int answer;

   printf("\nMENU for the working precision :\n");
   printf("  0. double precision :\n");
   printf("  1. double double precision :\n");
   printf("  2. quad double precision :\n");
   printf("Type 0, 1, or 2 to select the precision : ");
   scanf("%d", &answer);

   return answer;
}

int prompt_for_verbose ( void )
{
   char answer = 'n';

   printf("\nExtra output during path tracking ? (y/n) ");
   do
   {                            // loop to skip remaining newlines
      scanf("%c", &answer);
      if((answer == 'y') || (answer == 'n')) break;
   }
   while(1);

   if(answer == 'y') return 1;

   return 0;
}

int parameters_broadcast ( int myid, int nbrp, int verbose )
{
   int fail;
   double fltval,minstep,maxstep,beta1,beta2,alpha,epsilon,tolzero;
   double gamma[2];
   int numdeg,dendeg,maxiter,maxsteps;

   if(myid == 0)
   {
      fail = padcon_set_default_parameters();
      padcon_tune_homotopy_continuation_parameters();
      fail = padcon_get_homotopy_continuation_parameter(1,gamma);
      fail = padcon_get_homotopy_continuation_parameter(2,&fltval);
      numdeg = (int) fltval;
      fail = padcon_get_homotopy_continuation_parameter(3,&fltval);
      dendeg = (int) fltval;
      fail = padcon_get_homotopy_continuation_parameter(4,&maxstep);
      fail = padcon_get_homotopy_continuation_parameter(5,&minstep);
      fail = padcon_get_homotopy_continuation_parameter(6,&beta1);
      fail = padcon_get_homotopy_continuation_parameter(7,&beta2);
      fail = padcon_get_homotopy_continuation_parameter(8,&alpha);
      fail = padcon_get_homotopy_continuation_parameter(9,&epsilon);
      fail = padcon_get_homotopy_continuation_parameter(10,&tolzero);
      fail = padcon_get_homotopy_continuation_parameter(11,&fltval);
      maxiter = (int) fltval;
      fail = padcon_get_homotopy_continuation_parameter(12,&fltval);
      maxsteps = (int) fltval;
      fail = padcon_write_homotopy_continuation_parameters_to_defined_output();
   }
   MPI_Bcast(gamma,2,MPI_DOUBLE,0,MPI_COMM_WORLD); 
   MPI_Bcast(&numdeg,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&dendeg,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&maxstep,1,MPI_DOUBLE,0,MPI_COMM_WORLD); 
   MPI_Bcast(&minstep,1,MPI_DOUBLE,0,MPI_COMM_WORLD); 
   MPI_Bcast(&beta1,1,MPI_DOUBLE,0,MPI_COMM_WORLD); 
   MPI_Bcast(&beta2,1,MPI_DOUBLE,0,MPI_COMM_WORLD); 
   MPI_Bcast(&alpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD); 
   MPI_Bcast(&epsilon,1,MPI_DOUBLE,0,MPI_COMM_WORLD); 
   MPI_Bcast(&tolzero,1,MPI_DOUBLE,0,MPI_COMM_WORLD); 
   MPI_Bcast(&maxiter,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&maxsteps,1,MPI_INT,0,MPI_COMM_WORLD);
   if(myid != 0)
   {
      fail = padcon_set_homotopy_continuation_parameter(1,gamma);
      fltval = (double) numdeg;
      fail = padcon_set_homotopy_continuation_parameter(2,&fltval);
      fltval = (double) dendeg;
      fail = padcon_set_homotopy_continuation_parameter(3,&fltval);
      fail = padcon_set_homotopy_continuation_parameter(4,&maxstep);
      fail = padcon_set_homotopy_continuation_parameter(5,&minstep);
      fail = padcon_set_homotopy_continuation_parameter(6,&beta1);
      fail = padcon_set_homotopy_continuation_parameter(7,&beta2);
      fail = padcon_set_homotopy_continuation_parameter(8,&alpha);
      fail = padcon_set_homotopy_continuation_parameter(9,&epsilon);
      fail = padcon_set_homotopy_continuation_parameter(10,&tolzero);
      fltval = (double) maxiter;
      fail = padcon_set_homotopy_continuation_parameter(11,&fltval);
      fltval = (double) maxsteps;
      fail = padcon_set_homotopy_continuation_parameter(12,&fltval);
   }
   if(verbose > 0)
      if(myid == 1) padcon_write_homotopy_continuation_parameters();
}

int standard_track_paths
 ( int myid, int nbrp, int nbc, char *outfile, int verbose, int mhom,
   int nvr, int *idz )
{
   char myfile[80];
   char nbr[5];
   int fail,idx,cnt;

   if(nbc == 0)
   {
      fail = padcon_standard_track(0,"",0,verbose,mhom,nvr,idz);
   }
   else
   {
      for(idx=0; idx<nbc; idx++) myfile[idx] = outfile[idx];
      cnt = sprintf(nbr,"%d",myid);
      // printf("the number myid %d = \"%s\", cnt = %d\n",myid,nbr,cnt);
      for(idx=0; idx<nbc; idx++) myfile[nbc+idx] = nbr[idx];
      myfile[nbc+cnt] = '\0';
      if(verbose > 0)
         printf("Node %d will write to file \"%s\".\n",myid,myfile);
      if(myid == 0) // manager uses the defined output file
         fail = padcon_standard_track
                   (nbc+cnt,myfile,0,verbose,mhom,nvr,idz);
      else          // files for workers are local
         fail = padcon_standard_track
                   (nbc+cnt,myfile,1,verbose,mhom,nvr,idz);
   }
   return fail;
}

int dobldobl_track_paths
 ( int myid, int nbrp, int nbc, char *outfile, int verbose, int mhom,
   int nvr, int *idz )
{
   char myfile[80];
   char nbr[5];
   int fail,idx,cnt;

   if(nbc == 0)
   {
      fail = padcon_dobldobl_track(0,"",0,verbose,mhom,nvr,idz);
   }
   else
   {
      for(idx=0; idx<nbc; idx++) myfile[idx] = outfile[idx];
      cnt = sprintf(nbr,"%d",myid);
      // printf("the number myid %d = \"%s\", cnt = %d\n",myid,nbr,cnt);
      for(idx=0; idx<nbc; idx++) myfile[nbc+idx] = nbr[idx];
      myfile[nbc+cnt] = '\0';
      if(verbose > 0)
         printf("Node %d will write to file \"%s\".\n",myid,myfile);
      if(myid == 0) // manager uses the defined output file
         fail = padcon_dobldobl_track
                   (nbc+cnt,myfile,0,verbose,mhom,nvr,idz);
      else          // files for workers are local
         fail = padcon_dobldobl_track
                   (nbc+cnt,myfile,1,verbose,mhom,nvr,idz);
   }
   return fail;
}

int quaddobl_track_paths
 ( int myid, int nbrp, int nbc, char *outfile, int verbose, int mhom,
   int nvr, int *idz )
{
   char myfile[80];
   char nbr[5];
   int fail,idx,cnt;

   if(nbc == 0)
   {
      fail = padcon_quaddobl_track(0,"",0,verbose,mhom,nvr,idz);
   }
   else
   {
      for(idx=0; idx<nbc; idx++) myfile[idx] = outfile[idx];
      cnt = sprintf(nbr,"%d",myid);
      // printf("the number myid %d = \"%s\", cnt = %d\n",myid,nbr,cnt);
      for(idx=0; idx<nbc; idx++) myfile[nbc+idx] = nbr[idx];
      myfile[nbc+cnt] = '\0';
      if(verbose > 0)
         printf("Node %d will write to file \"%s\".\n",myid,myfile);
      if(myid == 0) // manager uses the defined output file
         fail = padcon_quaddobl_track
                   (nbc+cnt,myfile,0,verbose,mhom,nvr,idz);
      else          // files for workers are local
         fail = padcon_quaddobl_track
                   (nbc+cnt,myfile,1,verbose,mhom,nvr,idz);
   }
   return fail;
}

int standard_run
 ( int myid, int nbrp, int nbc, char* outfile, int verbose )
{
   int fail,dim,nbsols,mysolnum,len,*nbpaths,vrbrun,mhom,*idz;
   double startwtime,endwtime,wtime,*time;
   char nlsb;

   startwtime = MPI_Wtime();
   if(myid == 0)
   {
      time = (double*)calloc(nbrp,sizeof(double));
      nbpaths = (int*)calloc(nbrp,sizeof(int));
      fail = read_target_system_without_solutions();
      fail = copy_target_system_to_container();
      fail = syscon_number_of_standard_polynomials(&dim);
      mhom = padcon_prompt_for_multi_homogenization(dim);
      scanf("%c",&nlsb); // swallow new line symbol
      if(mhom == 1)
      {
         fail = syscon_standard_one_homogenization(0);
         fail = copy_container_to_target_system();
         dim = dim + 1;
      }
      else if(mhom > 1)
      {
         idz = (int*)calloc(dim,sizeof(int));
         padcon_define_partition(mhom,dim,idz);
         scanf("%c",&nlsb); // swallow new line symbol
         fail = syscon_standard_multi_homogenization(dim,mhom,idz,0);
         fail = copy_container_to_target_system();
         dim = dim + mhom;
      }
   }
   standard_dimensions_broadcast(myid,&dim,&dim);

   if(verbose > 0) printf("Process %d has dimension %d.\n",myid,dim);

   monomials_broadcast(myid,dim); // broadcast target system

   if(myid != 0) fail = copy_container_to_target_system();

   if(verbose > 0)
      if(myid == 1) fail = write_standard_target_system();

   MPI_Bcast(&mhom,1,MPI_INT,0,MPI_COMM_WORLD);
   if(verbose > 0)
      printf("\nProcess %d has %d for homogeneous.\n", myid, mhom);

   if(myid != 0)
   {
      fail = syscon_initialize_number_of_standard_polynomials(dim);
      if(mhom > 1) idz = (int*)calloc(dim-mhom,sizeof(int));
   }
   else
   {
      fail = read_standard_start_system();
      fail = copy_start_system_to_container();
      if(mhom == 1) 
      {
         fail = syscon_standard_one_homogenization(1);
         fail = copy_container_to_start_system();
      }
      else if(mhom > 1)
      {
         fail = syscon_standard_multi_homogenization(dim-mhom,mhom,idz,1);
         fail = copy_container_to_start_system();
      }
      fail = copy_start_solutions_to_container();
      if(mhom == 1)
      {
         fail = solcon_standard_one_homogenization();
         fail = copy_container_to_start_solutions();
      }
      else if(mhom > 1)
      {
         fail = solcon_standard_multi_homogenization(mhom);
         fail = copy_container_to_start_solutions();
      }
      fail = solcon_number_of_standard_solutions(&nbsols);
      if(verbose > 0) printf("Read %d start solutions.\n",nbsols);
      if(mhom == 1) fail = padcon_add_Z0();
      if(mhom > 1) fail = padcon_add_symbols(mhom);

      fail = write_standard_target_system(); // writes to file
      fail = write_standard_start_system(); // writes to file
      fail = write_start_solutions(); // writes solutions to file
   }
   monomials_broadcast(myid,dim); // broadcast start system

   if(myid != 0) fail = copy_container_to_start_system();

   if(verbose > 0)
      if(myid == 1) fail = write_standard_start_system();

   parameters_broadcast(myid,nbrp,verbose);

   if(myid == 0) vrbrun = prompt_for_verbose();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Bcast(&vrbrun,1,MPI_INT,0,MPI_COMM_WORLD);

   if(mhom > 1)
   {
      MPI_Bcast(idz,dim-mhom,MPI_INT,0,MPI_COMM_WORLD);
      if(verbose > 0)
      {
         printf("\nProcess %d has partition ", myid);
         for(int idx=0; idx<dim-mhom; idx++) printf(" %d",idz[idx]);
         printf("\n");
      }
   }
   MPI_Bcast(&nbsols,1,MPI_INT,0,MPI_COMM_WORLD);
   solutions_distribute(myid,nbsols,dim,nbrp,&mysolnum);

   fail = solcon_number_of_standard_solutions(&len);
   if(verbose > 0) printf("Node %d has %d solutions.\n",myid,len);

   if(myid > 0)
   {
      fail = copy_container_to_start_solutions();
      fail = standard_track_paths
                (myid,nbrp,nbc,outfile,vrbrun,mhom,dim-mhom,idz);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   solutions_collect(myid,nbsols,dim,nbrp,mysolnum);

   if(myid == 0)
   {
      if(mhom == 1)
      {
         fail = solcon_standard_one_affinization();
         fail = copy_target_system_to_container();
         fail = syscon_standard_one_affinization();
      }
      else if(mhom > 1)
      {
         fail = solcon_standard_multi_affinization(dim-mhom,mhom,idz);
         fail = copy_target_system_to_container();
         fail = syscon_standard_multi_affinization(mhom);
      }
      fail = syscon_write_standard_system();
      fail = solcon_write_solution_banner_to_defined_output_file();
      fail = solcon_write_standard_solutions();
   }
   endwtime = MPI_Wtime();
   wtime = endwtime-startwtime;
   MPI_Gather(&wtime,1,MPI_DOUBLE,time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Gather(&mysolnum,1,MPI_INT,nbpaths,1,MPI_INT,0,MPI_COMM_WORLD);

   if(myid == 0)
   {
      nbpaths[0] = nbsols;
      write_time_and_paths_to_defined_output_file(nbrp,time,nbpaths);
      free(time); free(nbpaths);
   }
   return 0;
}

int dobldobl_run
 ( int myid, int nbrp, int nbc, char* outfile, int verbose )
{
   int fail,dim,nbsols,mysolnum,len,*nbpaths,vrbrun,mhom,*idz;
   double startwtime,endwtime,wtime,*time;
   char nlsb;

   startwtime = MPI_Wtime();
   if(myid == 0)
   {
      time = (double*)calloc(nbrp,sizeof(double));
      nbpaths = (int*)calloc(nbrp,sizeof(int));
      fail = read_dobldobl_target_system_without_solutions();
      fail = copy_dobldobl_target_system_to_container();
      fail = syscon_number_of_dobldobl_polynomials(&dim);
      mhom = padcon_prompt_for_multi_homogenization(dim);
      scanf("%c",&nlsb); // swallow new line symbol
      if(mhom == 1)
      {
         fail = syscon_dobldobl_one_homogenization(0);
         fail = copy_dobldobl_container_to_target_system();
         dim = dim + 1;
      }
      else if(mhom > 1)
      {
         idz = (int*)calloc(dim,sizeof(int));
         padcon_define_partition(mhom,dim,idz);
         scanf("%c",&nlsb); // swallow new line symbol
         fail = syscon_dobldobl_multi_homogenization(dim,mhom,idz,0);
         fail = copy_dobldobl_container_to_target_system();
         dim = dim + mhom;
      }
   }
   dobldobl_dimensions_broadcast(myid,&dim,&dim);

   if(verbose > 0) printf("Process %d has dimension %d.\n",myid,dim);

   dobldobl_monomials_broadcast(myid,dim);

   if(myid != 0) fail = copy_dobldobl_container_to_target_system();

   if(verbose > 0)
      if(myid == 1) fail = write_dobldobl_target_system();

   MPI_Bcast(&mhom,1,MPI_INT,0,MPI_COMM_WORLD);
   if(verbose > 0)
      printf("\nProcess %d has %d for homogeneous.\n", myid, mhom);

   if(myid != 0)
   {
      fail = syscon_initialize_number_of_dobldobl_polynomials(dim);
      if(mhom > 1) idz = (int*)calloc(dim-mhom,sizeof(int));
   }
   else
   {
      fail = read_dobldobl_start_system();
      fail = copy_dobldobl_start_system_to_container();
      if(mhom == 1) 
      {
         fail = syscon_dobldobl_one_homogenization(1);
         fail = copy_dobldobl_container_to_start_system();
      }
      else if(mhom > 1)
      {
         fail = syscon_dobldobl_multi_homogenization(dim,mhom-mhom,idz,1);
         fail = copy_dobldobl_container_to_start_system();
      }
      fail = copy_dobldobl_start_solutions_to_container();
      if(mhom == 1)
      {
         fail = solcon_dobldobl_one_homogenization();
         fail = copy_dobldobl_container_to_start_solutions();
      }
      else if(mhom > 1)
      {
         fail = solcon_dobldobl_multi_homogenization(mhom);
         fail = copy_dobldobl_container_to_start_solutions();
      }
      fail = solcon_number_of_dobldobl_solutions(&nbsols);
      if(verbose > 0) printf("Read %d start solutions.\n",nbsols);
      if(mhom == 1) fail = padcon_add_Z0();
      if(mhom > 1) fail = padcon_add_symbols(mhom);

      fail = write_dobldobl_target_system(); // writes to file
      fail = write_dobldobl_start_system(); // writes to file
      fail = write_dobldobl_start_solutions(); // writes solutions to file
   }
   dobldobl_monomials_broadcast(myid,dim); // broadcast start system

   if(myid != 0) fail = copy_dobldobl_container_to_start_system();

   if(verbose > 0)
      if(myid == 1) fail = write_dobldobl_start_system();

   parameters_broadcast(myid,nbrp,verbose);
   if(myid == 0) vrbrun = prompt_for_verbose();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Bcast(&vrbrun,1,MPI_INT,0,MPI_COMM_WORLD);

   if(mhom > 1)
   {
      MPI_Bcast(idz,dim-mhom,MPI_INT,0,MPI_COMM_WORLD);
      if(verbose > 0)
      {
         printf("\nProcess %d has partition ", myid);
         for(int idx=0; idx<dim-mhom; idx++) printf(" %d",idz[idx]);
         printf("\n");
      }
   }
   MPI_Bcast(&nbsols,1,MPI_INT,0,MPI_COMM_WORLD);
   dobldobl_solutions_distribute(myid,nbsols,dim,nbrp,&mysolnum,verbose);

   fail = solcon_number_of_dobldobl_solutions(&len);
   if(verbose > 0) printf("Node %d has %d solutions.\n",myid,len);

   if(myid > 0)
   {
      fail = copy_dobldobl_container_to_start_solutions();
      fail = dobldobl_track_paths
                (myid,nbrp,nbc,outfile,vrbrun,mhom,dim-mhom,idz);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   dobldobl_solutions_collect(myid,nbsols,dim,nbrp,mysolnum);

   if(myid == 0)
   {
      if(mhom == 1)
      {
         fail = solcon_dobldobl_one_affinization();
         fail = copy_dobldobl_target_system_to_container();
         fail = syscon_dobldobl_one_affinization();
      }
      else if(mhom > 1)
      {
         fail = solcon_dobldobl_multi_affinization(dim-mhom,mhom,idz);
         fail = copy_dobldobl_target_system_to_container();
         fail = syscon_dobldobl_multi_affinization(mhom);
      }
      fail = syscon_write_dobldobl_system();
      fail = solcon_write_solution_banner_to_defined_output_file();
      fail = solcon_write_dobldobl_solutions();
   }
   endwtime = MPI_Wtime();
   wtime = endwtime-startwtime;
   MPI_Gather(&wtime,1,MPI_DOUBLE,time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Gather(&mysolnum,1,MPI_INT,nbpaths,1,MPI_INT,0,MPI_COMM_WORLD);

   if(myid == 0)
   {
      nbpaths[0] = nbsols;
      write_time_and_paths_to_defined_output_file(nbrp,time,nbpaths);
      free(time); free(nbpaths);
   }
   return 0;
}

int quaddobl_run
 ( int myid, int nbrp, int nbc, char* outfile, int verbose )
{
   int fail,dim,nbsols,mysolnum,len,*nbpaths,vrbrun,mhom,*idz;
   double startwtime,endwtime,wtime,*time;
   char nlsb;

   startwtime = MPI_Wtime();
   if(myid == 0)
   {
      time = (double*)calloc(nbrp,sizeof(double));
      nbpaths = (int*)calloc(nbrp,sizeof(int));
      fail = read_quaddobl_target_system_without_solutions();
      fail = copy_quaddobl_target_system_to_container();
      fail = syscon_number_of_quaddobl_polynomials(&dim);
      mhom = padcon_prompt_for_multi_homogenization(dim);
      scanf("%c",&nlsb); // swallow new line symbol
      if(mhom == 1)
      {
         fail = syscon_quaddobl_one_homogenization(0);
         fail = copy_quaddobl_container_to_target_system();
         dim = dim + 1;
      }
      else if(mhom > 1)
      {
         idz = (int*)calloc(dim,sizeof(int));
         padcon_define_partition(mhom,dim,idz);
         scanf("%c",&nlsb); // swallow new line symbol
         fail = syscon_dobldobl_multi_homogenization(dim,mhom,idz,0);
         fail = copy_dobldobl_container_to_target_system();
         dim = dim + mhom;
      }
   }
   quaddobl_dimensions_broadcast(myid,&dim,&dim);

   if(verbose > 0) printf("Process %d has dimension %d.\n",myid,dim);

   quaddobl_monomials_broadcast(myid,dim);

   if(myid != 0) fail = copy_quaddobl_container_to_target_system();

   if(verbose > 0)
      if(myid == 1) fail = write_quaddobl_target_system();

   MPI_Bcast(&mhom,1,MPI_INT,0,MPI_COMM_WORLD);
   if(verbose > 0)
      printf("\nProcess %d has %d for homogeneous.\n", myid, mhom);

   if(myid != 0)
   {
      fail = syscon_initialize_number_of_quaddobl_polynomials(dim);
      if(mhom > 1) idz = (int*)calloc(dim-mhom,sizeof(int));
   }
   else
   {
      fail = read_quaddobl_start_system();
      fail = copy_quaddobl_start_system_to_container();
      if(mhom == 1) 
      {
         fail = syscon_quaddobl_one_homogenization(1);
         fail = copy_quaddobl_container_to_start_system();
      }
      else if(mhom > 1)
      {
         fail = syscon_quaddobl_multi_homogenization(dim,mhom,idz,1);
         fail = copy_quaddobl_container_to_start_system();
      }
      fail = copy_quaddobl_start_solutions_to_container();
      if(mhom == 1)
      {
         fail = solcon_quaddobl_one_homogenization();
         fail = copy_quaddobl_container_to_start_solutions();
      }
      else if(mhom > 1)
      {
         fail = solcon_quaddobl_multi_homogenization(mhom);
         fail = copy_quaddobl_container_to_start_solutions();
      }
      fail = solcon_number_of_quaddobl_solutions(&nbsols);
      if(verbose > 0) printf("Read %d start solutions.\n",nbsols);
      if(mhom == 1) fail = padcon_add_Z0();
      if(mhom > 1) fail = padcon_add_symbols(mhom);

      fail = write_quaddobl_target_system(); // writes to file
      fail = write_quaddobl_start_system(); // writes to file
      fail = write_quaddobl_start_solutions(); // writes solutions to file
   }
   quaddobl_monomials_broadcast(myid,dim); // broadcast start system

   if(myid != 0) fail = copy_quaddobl_container_to_start_system();

   if(verbose > 0)
      if(myid == 1) fail = write_quaddobl_start_system();

   parameters_broadcast(myid,nbrp,1);

   if(myid == 0) vrbrun = prompt_for_verbose();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Bcast(&vrbrun,1,MPI_INT,0,MPI_COMM_WORLD);

   if(mhom > 1)
   {
      MPI_Bcast(idz,dim-mhom,MPI_INT,0,MPI_COMM_WORLD);
      if(verbose > 0)
      {
         printf("\nProcess %d has partition ", myid);
         for(int idx=0; idx<dim-mhom; idx++) printf(" %d",idz[idx]);
         printf("\n");
      }
   }
   MPI_Bcast(&nbsols,1,MPI_INT,0,MPI_COMM_WORLD);
   quaddobl_solutions_distribute(myid,nbsols,dim,nbrp,&mysolnum,verbose);

   fail = solcon_number_of_quaddobl_solutions(&len);
   if(verbose > 0) printf("Node %d has %d solutions.\n",myid,len);

   if(myid > 0)
   {
      fail = copy_quaddobl_container_to_start_solutions();
      fail = quaddobl_track_paths
                (myid,nbrp,nbc,outfile,vrbrun,mhom,dim-mhom,idz);
   }
   quaddobl_solutions_collect(myid,nbsols,dim,nbrp,mysolnum);

   if(myid == 0)
   {
      if(mhom == 1) 
      {
         fail = solcon_quaddobl_one_affinization();
         fail = copy_quaddobl_target_system_to_container();
         fail = syscon_quaddobl_one_affinization();
      }
      else if(mhom > 1)
      {
         fail = solcon_dobldobl_multi_affinization(dim-mhom,mhom,idz);
         fail = copy_dobldobl_target_system_to_container();
         fail = syscon_dobldobl_multi_affinization(mhom);
      }
      fail = syscon_write_quaddobl_system();
      fail = solcon_write_solution_banner_to_defined_output_file();
      fail = solcon_write_quaddobl_solutions();
   }
   endwtime = MPI_Wtime();
   wtime = endwtime-startwtime;
   MPI_Gather(&wtime,1,MPI_DOUBLE,time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Gather(&mysolnum,1,MPI_INT,nbpaths,1,MPI_INT,0,MPI_COMM_WORLD);

   if(myid == 0)
   {
      nbpaths[0] = nbsols;
      write_time_and_paths_to_defined_output_file(nbrp,time,nbpaths);
      free(time); free(nbpaths);
   }
   return 0;
}
