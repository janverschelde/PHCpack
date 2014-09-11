#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include "phcpack.h"
#include "unisolvers.h"
#include "giftwrappers.h"
#include "schubert.h"
#include "syscon.h"
#include "solcon.h"
#include "product.h"
#include "lists_and_strings.h"
#include "celcon.h"
#include "witset.h"
#include "mapcon.h"
#include "next_track.h"
#include "structmember.h"

extern void adainit ( void );

int g_ada_initialized = 0;

void initialize(void)
{
   if(!g_ada_initialized)
   {
      adainit();
      g_ada_initialized = 1;
   }
}

/* wrapping functions in phcpack.h starts from here */

static PyObject *py2c_PHCpack_version_string ( PyObject *self, PyObject *args )
{
   int fail;
   int size = 40;
   char s[size];

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = version_string(&size,s);
   s[size] = '\0';
              
   return Py_BuildValue("s",s);
}

static PyObject *py2c_set_seed( PyObject *self, PyObject *args )
{
   int fail,seed;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&seed)) return NULL;

   fail = set_seed(seed);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_standard_target_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = read_standard_target_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_standard_target_system_from_file
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *name;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&name)) return NULL;
   fail = read_standard_target_system_from_file(nc,name);
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_standard_start_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = read_standard_start_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_standard_start_system_from_file
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *name;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&name)) return NULL;   
   fail = read_standard_start_system_from_file(nc,name);
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_dobldobl_target_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = read_dobldobl_target_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_dobldobl_target_system_from_file
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *name;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&name)) return NULL;
   fail = read_dobldobl_target_system_from_file(nc,name);
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_dobldobl_start_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = read_dobldobl_start_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_dobldobl_start_system_from_file
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *name;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&name)) return NULL;   
   fail = read_dobldobl_start_system_from_file(nc,name);
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_quaddobl_target_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = read_quaddobl_target_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_quaddobl_target_system_from_file
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *name;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&name)) return NULL;
   fail = read_quaddobl_target_system_from_file(nc,name);
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_quaddobl_start_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = read_quaddobl_start_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_quaddobl_start_system_from_file
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *name;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&name)) return NULL;   
   fail = read_quaddobl_start_system_from_file(nc,name);
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_define_output_file ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = define_output_file();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_write_standard_target_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = write_standard_target_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_write_standard_start_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = write_standard_start_system();
   
   return Py_BuildValue("i",fail);
}

/* moving systems from and to containers */

static PyObject *py2c_copy_target_system_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_target_system_to_container();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_dobldobl_target_system_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_dobldobl_target_system_to_container();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_quaddobl_target_system_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_quaddobl_target_system_to_container();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_multprec_target_system_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_multprec_target_system_to_container();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_container_to_target_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_container_to_target_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_dobldobl_container_to_target_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_dobldobl_container_to_target_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_quaddobl_container_to_target_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_quaddobl_container_to_target_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_multprec_container_to_target_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_multprec_container_to_target_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_start_system_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_start_system_to_container();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_dobldobl_start_system_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_dobldobl_start_system_to_container();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_quaddobl_start_system_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_dobldobl_start_system_to_container();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_multprec_start_system_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_multprec_start_system_to_container();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_container_to_start_system 
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_container_to_start_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_dobldobl_container_to_start_system 
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_dobldobl_container_to_start_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_quaddobl_container_to_start_system 
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_quaddobl_container_to_start_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_multprec_container_to_start_system 
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_multprec_container_to_start_system();
              
   return Py_BuildValue("i",fail);
}

/* creation of homotopy and tracking all solution paths */

static PyObject *py2c_create_homotopy ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = create_homotopy();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_create_homotopy_with_gamma
 ( PyObject *self, PyObject *args )
{
   int fail;
   double g_re,g_im;

   initialize();
   if(!PyArg_ParseTuple(args,"dd",&g_re,&g_im)) return NULL;   
   fail = create_homotopy_with_given_gamma(g_re,g_im);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_create_dobldobl_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = create_dobldobl_homotopy();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_create_dobldobl_homotopy_with_gamma
 ( PyObject *self, PyObject *args )
{
   int fail;
   double g_re,g_im;

   initialize();
   if(!PyArg_ParseTuple(args,"dd",&g_re,&g_im)) return NULL;   
   fail = create_dobldobl_homotopy_with_given_gamma(g_re,g_im);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_create_quaddobl_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = create_quaddobl_homotopy();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_create_quaddobl_homotopy_with_gamma
 ( PyObject *self, PyObject *args )
{
   int fail;
   double g_re,g_im;

   initialize();
   if(!PyArg_ParseTuple(args,"dd",&g_re,&g_im)) return NULL;   
   fail = create_quaddobl_homotopy_with_given_gamma(g_re,g_im);
              
   return Py_BuildValue("i",fail);
}


static PyObject *py2c_create_multprec_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = create_multprec_homotopy();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_create_multprec_homotopy_with_gamma
 ( PyObject *self, PyObject *args )
{
   int fail;
   double g_re,g_im;

   initialize();
   if(!PyArg_ParseTuple(args,"dd",&g_re,&g_im)) return NULL;   
   fail = create_multprec_homotopy_with_given_gamma(g_re,g_im);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_homotopy ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = clear_homotopy();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_dobldobl_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = clear_dobldobl_homotopy();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_quaddobl_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = clear_quaddobl_homotopy();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_multprec_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = clear_multprec_homotopy();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_write_start_solutions ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = write_start_solutions();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_tune_continuation_parameters
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = tune_continuation_parameters();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_show_continuation_parameters
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = show_continuation_parameters();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_autotune_continuation_parameters
 ( PyObject *self, PyObject *args )
{
   int fail,level,nbdgt;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&level,&nbdgt)) return NULL;   
   fail = autotune_continuation_parameters(level,nbdgt);
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_determine_output_during_continuation
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = determine_output_during_continuation();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solve_by_standard_homotopy_continuation
 ( PyObject *self, PyObject *args )
{
   int fail, nbtasks = 0;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&nbtasks)) return NULL;   
   fail = solve_by_standard_homotopy_continuation(nbtasks);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solve_by_dobldobl_homotopy_continuation
 ( PyObject *self, PyObject *args )
{
   int fail, nbtasks = 0;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&nbtasks)) return NULL;   
   fail = solve_by_dobldobl_homotopy_continuation(nbtasks);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solve_by_quaddobl_homotopy_continuation
 ( PyObject *self, PyObject *args )
{
   int fail, nbtasks = 0;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&nbtasks)) return NULL;   
   fail = solve_by_quaddobl_homotopy_continuation(nbtasks);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solve_by_multprec_homotopy_continuation
 ( PyObject *self, PyObject *args )
{
   int fail,deci;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&deci)) return NULL;   
   fail = solve_by_multprec_homotopy_continuation(deci);

   return Py_BuildValue("i",fail);
}

/* moving solutions from and to containers */

static PyObject *py2c_copy_target_solutions_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_target_solutions_to_container();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_dobldobl_target_solutions_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_dobldobl_target_solutions_to_container();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_quaddobl_target_solutions_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_quaddobl_target_solutions_to_container();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_multprec_target_solutions_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_multprec_target_solutions_to_container();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_container_to_target_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_container_to_target_solutions();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_dobldobl_container_to_target_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_dobldobl_container_to_target_solutions();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_quaddobl_container_to_target_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_quaddobl_container_to_target_solutions();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_multprec_container_to_target_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_multprec_container_to_target_solutions();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_start_solutions_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_start_solutions_to_container();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_dobldobl_start_solutions_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_dobldobl_start_solutions_to_container();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_quaddobl_start_solutions_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_quaddobl_start_solutions_to_container();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_multprec_start_solutions_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_multprec_start_solutions_to_container();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_container_to_start_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_container_to_start_solutions();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_dobldobl_container_to_start_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_dobldobl_container_to_start_solutions();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_quaddobl_container_to_start_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_quaddobl_container_to_start_solutions();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_multprec_container_to_start_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = copy_multprec_container_to_start_solutions();

   return Py_BuildValue("i",fail);
}

/* black box solver, mixed volume calculator, and Newton step */

static PyObject *py2c_solve_system ( PyObject *self, PyObject *args )
{
   int fail,rc;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = solve_system(&rc);
   return Py_BuildValue("i",rc);
}

static PyObject *py2c_solve_Laurent_system ( PyObject *self, PyObject *args )
{
   int silent,fail,rc;

   initialize();
   if (!PyArg_ParseTuple(args,"i",&silent)) return NULL;
   fail = solve_Laurent_system(&rc,silent);
   return Py_BuildValue("i",rc);
}

static PyObject *py2c_mixed_volume ( PyObject *self, PyObject *args )
{
   int stable,fail,mv,smv;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&stable)) return NULL;
   if(stable == 0)
   {
      fail = mixed_volume(&mv);
      return Py_BuildValue("i",mv);
   }
   else
   {
      fail = stable_mixed_volume(&mv,&smv);
      return Py_BuildValue("(i,i)",mv,smv);
   }
}

static PyObject *py2c_standard_deflate ( PyObject *self, PyObject *args )
{
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   {
      int fail = standard_deflate();
      return Py_BuildValue("i",fail);
   }
}

static PyObject *py2c_dobldobl_deflate ( PyObject *self, PyObject *args )
{
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   {
      int fail = dobldobl_deflate();
      return Py_BuildValue("i",fail);
   }
}

static PyObject *py2c_quaddobl_deflate ( PyObject *self, PyObject *args )
{
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   {
      int fail = quaddobl_deflate();
      return Py_BuildValue("i",fail);
   }
}

static PyObject *py2c_standard_Newton_step ( PyObject *self, PyObject *args )
{
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   {
      int fail = standard_Newton_step();
      return Py_BuildValue("i",fail);
   }
}

static PyObject *py2c_dobldobl_Newton_step ( PyObject *self, PyObject *args )
{
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   {
      int fail = dobldobl_Newton_step();
      return Py_BuildValue("i",fail);
   }
}

static PyObject *py2c_quaddobl_Newton_step ( PyObject *self, PyObject *args )
{
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   {
      int fail = quaddobl_Newton_step();
      return Py_BuildValue("i",fail);
   }
}

static PyObject *py2c_multprec_Newton_step ( PyObject *self, PyObject *args )
{
   int fail,decimals;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&decimals)) return NULL;
   fail = multprec_Newton_step(decimals);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_standard_Newton_Laurent_step
 ( PyObject *self, PyObject *args )
{
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   {
      int fail = standard_Newton_Laurent_step();
      return Py_BuildValue("i",fail);
   }
}

static PyObject *py2c_dobldobl_Newton_Laurent_step
 ( PyObject *self, PyObject *args )
{
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   {
      int fail = dobldobl_Newton_Laurent_step();
      return Py_BuildValue("i",fail);
   }
}

static PyObject *py2c_quaddobl_Newton_Laurent_step
 ( PyObject *self, PyObject *args )
{
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   {
      int fail = quaddobl_Newton_Laurent_step();
      return Py_BuildValue("i",fail);
   }
}

static PyObject *py2c_multprec_Newton_Laurent_step
 ( PyObject *self, PyObject *args )
{
   int fail,decimals;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&decimals)) return NULL;
   fail = multprec_Newton_Laurent_step(decimals);

   return Py_BuildValue("i",fail);
}

/* wrapping functions in unisolvers.h starts from here */

static PyObject *py2c_usolve_standard ( PyObject *self, PyObject *args )
{
   int fail,max,nit;
   double eps;

   initialize();
   if(!PyArg_ParseTuple(args,"id",&max,&eps)) return NULL;
   fail = solve_with_standard_doubles(max,eps,&nit);

   return Py_BuildValue("i",nit);
}

static PyObject *py2c_usolve_dobldobl ( PyObject *self, PyObject *args )
{
   int fail,max,nit;
   double eps;

   initialize();
   if(!PyArg_ParseTuple(args,"id",&max,&eps)) return NULL;
   fail = solve_with_double_doubles(max,eps,&nit);

   return Py_BuildValue("i",nit);
}

static PyObject *py2c_usolve_quaddobl ( PyObject *self, PyObject *args )
{
   int fail,max,nit;
   double eps;

   initialize();
   if(!PyArg_ParseTuple(args,"id",&max,&eps)) return NULL;
   fail = solve_with_quad_doubles(max,eps,&nit);

   return Py_BuildValue("i",nit);
}

static PyObject *py2c_usolve_multprec ( PyObject *self, PyObject *args )
{
   int fail,dcp,max,nit;
   double eps;

   initialize();
   if(!PyArg_ParseTuple(args,"iid",&dcp,&max,&eps)) return NULL;
   fail = solve_with_multiprecision(dcp,max,eps,&nit);

   return Py_BuildValue("i",nit);
}

/* wrapping functions in giftwrappers.h starts from here */

static PyObject *py2c_giftwrap_planar ( PyObject *self, PyObject *args )
{
   int fail,nbc_pts;
   char *pts;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nbc_pts,&pts)) return NULL;

   {
      int nbc_hull;
      char hull[10*nbc_pts];

      fail = convex_hull_2d(nbc_pts,pts,&nbc_hull,hull);
      hull[nbc_hull] = '\0';

      return Py_BuildValue("s",hull);
   }
}

static PyObject *py2c_giftwrap_convex_hull ( PyObject *self, PyObject *args )
{
   int fail,nbc_pts;
   char *pts;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nbc_pts,&pts)) return NULL;

   fail = convex_hull(nbc_pts,pts);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_giftwrap_number_of_facets
 ( PyObject *self, PyObject *args )
{
   int fail,dim,number;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;

   fail = number_of_facets(dim,&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_giftwrap_retrieve_facet
 ( PyObject *self, PyObject *args )
{
   int fail,dim,ind,ncp;
   char rep[256];

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&dim,&ind)) return NULL;

   fail = retrieve_facet(dim,ind,&ncp,rep);
   rep[ncp] = '\0';

   return Py_BuildValue("s",rep);
}

static PyObject *py2c_giftwrap_clear_3d_facets
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   

   fail = clear_3d_facets();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_giftwrap_clear_4d_facets
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   

   fail = clear_4d_facets();

   return Py_BuildValue("i",fail);
}

/* wrapping functions in syscon.h starts from here */

static PyObject *py2c_syscon_read_system ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_read_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_read_Laurent_system 
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_read_Laurent_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_read_dobldobl_system 
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_read_dobldobl_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_read_dobldobl_Laurent_system 
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_read_dobldobl_Laurent_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_read_quaddobl_system 
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_read_quaddobl_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_read_quaddobl_Laurent_system 
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_read_quaddobl_Laurent_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_read_multprec_system 
 ( PyObject *self, PyObject *args )
{
   int fail,deci;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&deci)) return NULL;   
   fail = syscon_read_multprec_system(deci);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_read_multprec_Laurent_system 
 ( PyObject *self, PyObject *args )
{
   int fail,deci;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&deci)) return NULL;   
   fail = syscon_read_multprec_Laurent_system(deci);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_random_system ( PyObject *self, PyObject *args )
{
   int fail,n,m,d,c;

   initialize();
   if(!PyArg_ParseTuple(args,"iiii",&n,&m,&d,&c)) return NULL;   
   fail = syscon_random_system(n,m,d,c);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_write_system ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_write_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_write_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_write_Laurent_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_write_dobldobl_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_write_dobldobl_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_write_dobldobl_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_write_dobldobl_Laurent_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_write_quaddobl_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_write_quaddobl_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_write_quaddobl_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_write_quaddobl_Laurent_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_write_multprec_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_write_multprec_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_write_multprec_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_write_multprec_Laurent_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_clear_system ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_clear_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_clear_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_clear_Laurent_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_clear_dobldobl_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_clear_dobldobl_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_clear_dobldobl_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_clear_dobldobl_Laurent_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_clear_quaddobl_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_clear_quaddobl_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_clear_quaddobl_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_clear_quaddobl_Laurent_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_clear_multprec_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_clear_multprec_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_clear_multprec_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_clear_multprec_Laurent_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_number_of_symbols
 ( PyObject *self, PyObject *args )
{
   int nb,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_number_of_symbols(&nb);
              
   return Py_BuildValue("i",nb);
}

static PyObject *py2c_syscon_write_symbols
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_write_symbols();
   printf("\n");
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_string_of_symbols
 ( PyObject *self, PyObject *args )
{
   int nb;
   int fail = syscon_number_of_symbols(&nb);
   int size = 80*nb;
   char s[size];

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_string_of_symbols(&size,s);
              
   return Py_BuildValue("s",s);
}

static PyObject *py2c_syscon_remove_symbol_name
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *s;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&s)) return NULL;
   fail = syscon_remove_symbol_name_from_table(nc,s);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_clear_symbol_table
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_clear_symbol_table();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_number_of_polynomials
 ( PyObject *self, PyObject *args )
{
   int fail,number;

   initialize();
   if (!PyArg_ParseTuple(args,"")) return NULL;
   fail = syscon_number_of_polynomials(&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_syscon_number_of_dobldobl_polynomials
 ( PyObject *self, PyObject *args )
{
   int fail,number;

   initialize();
   if (!PyArg_ParseTuple(args,"")) return NULL;
   fail = syscon_number_of_dobldobl_polynomials(&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_syscon_number_of_quaddobl_polynomials
 ( PyObject *self, PyObject *args )
{
   int fail,number;

   initialize();
   if (!PyArg_ParseTuple(args,"")) return NULL;
   fail = syscon_number_of_quaddobl_polynomials(&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_syscon_number_of_multprec_polynomials
 ( PyObject *self, PyObject *args )
{
   int fail,number;

   initialize();
   if (!PyArg_ParseTuple(args,"")) return NULL;
   fail = syscon_number_of_multprec_polynomials(&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_syscon_number_of_Laurentials
 ( PyObject *self, PyObject *args )
{
   int fail,number;
   static PyObject *a;

   initialize();
   /* if (!PyArg_ParseTuple(args,"i",&number)) return NULL; */
   if (!PyArg_ParseTuple(args,"")) return NULL;
   fail = syscon_number_of_Laurentials(&number);

   /* a = Py_BuildValue("{i:i}",fail,number);	 
   return a; */

   return Py_BuildValue("i",number);
}

static PyObject *py2c_syscon_number_of_dobldobl_Laurentials
 ( PyObject *self, PyObject *args )
{
   int fail,number;

   initialize();
   if (!PyArg_ParseTuple(args,"")) return NULL;
   fail = syscon_number_of_dobldobl_Laurentials(&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_syscon_number_of_quaddobl_Laurentials
 ( PyObject *self, PyObject *args )
{
   int fail,number;

   initialize();
   if (!PyArg_ParseTuple(args,"")) return NULL;
   fail = syscon_number_of_quaddobl_Laurentials(&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_syscon_number_of_multprec_Laurentials
 ( PyObject *self, PyObject *args )
{
   int fail,number;

   initialize();
   if (!PyArg_ParseTuple(args,"")) return NULL;
   fail = syscon_number_of_multprec_Laurentials(&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_syscon_initialize_number
 ( PyObject *self, PyObject *args )
{      
   int fail,dim;

   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;
   initialize();
   fail = syscon_initialize_number(dim);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_initialize_number_of_dobldobl_polynomials
 ( PyObject *self, PyObject *args )
{      
   int fail,dim;

   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;
   initialize();
   fail = syscon_initialize_number_of_dobldobl_polynomials(dim);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_initialize_number_of_quaddobl_polynomials
 ( PyObject *self, PyObject *args )
{      
   int fail,dim;

   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;
   initialize();
   fail = syscon_initialize_number_of_quaddobl_polynomials(dim);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_initialize_number_of_multprec_polynomials
 ( PyObject *self, PyObject *args )
{      
   int fail,dim;

   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;
   initialize();
   fail = syscon_initialize_number_of_multprec_polynomials(dim);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_initialize_number_of_Laurentials
 ( PyObject *self, PyObject *args )
{      
   int fail,dim;

   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;
   initialize();
   fail = syscon_initialize_number_of_Laurentials(dim);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_initialize_number_of_dobldobl_Laurentials
 ( PyObject *self, PyObject *args )
{      
   int fail,dim;

   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;
   initialize();
   fail = syscon_initialize_number_of_dobldobl_Laurentials(dim);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_initialize_number_of_quaddobl_Laurentials
 ( PyObject *self, PyObject *args )
{      
   int fail,dim;

   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;
   initialize();
   fail = syscon_initialize_number_of_quaddobl_Laurentials(dim);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_initialize_number_of_multprec_Laurentials
 ( PyObject *self, PyObject *args )
{      
   int fail,dim;

   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;
   initialize();
   fail = syscon_initialize_number_of_multprec_Laurentials(dim);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_degree_of_polynomial
 ( PyObject *self, PyObject *args )
{
   int fail,equ,deg;

   if(!PyArg_ParseTuple(args,"i",&equ)) return NULL;
   initialize();
   fail = syscon_degree_of_polynomial(equ,&deg);
                 
   return Py_BuildValue("i",deg);
}

static PyObject *py2c_syscon_degree_of_dobldobl_polynomial
 ( PyObject *self, PyObject *args )
{
   int fail,equ,deg;

   if(!PyArg_ParseTuple(args,"i",&equ)) return NULL;
   initialize();
   fail = syscon_degree_of_dobldobl_polynomial(equ,&deg);
                 
   return Py_BuildValue("i",deg);
}

static PyObject *py2c_syscon_degree_of_quaddobl_polynomial
 ( PyObject *self, PyObject *args )
{
   int fail,equ,deg;

   if(!PyArg_ParseTuple(args,"i",&equ)) return NULL;
   initialize();
   fail = syscon_degree_of_quaddobl_polynomial(equ,&deg);
                 
   return Py_BuildValue("i",deg);
}

static PyObject *py2c_syscon_degree_of_multprec_polynomial
 ( PyObject *self, PyObject *args )
{
   int fail,equ,deg;

   if(!PyArg_ParseTuple(args,"i",&equ)) return NULL;
   initialize();
   fail = syscon_degree_of_multprec_polynomial(equ,&deg);
                 
   return Py_BuildValue("i",deg);
}

static PyObject *py2c_syscon_number_of_terms ( PyObject *self, PyObject *args )
{
   int fail,i,number;

   initialize();
   /* if (!PyArg_ParseTuple(args,"ii",&i,&number)) return NULL; */
   if (!PyArg_ParseTuple(args,"i",&i)) return NULL;
   fail = syscon_number_of_terms(i,&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_syscon_number_of_Laurent_terms
 ( PyObject *self, PyObject *args )
{
   int fail,i,number;

   initialize();
   /* if (!PyArg_ParseTuple(args,"ii",&i,&number)) return NULL; */
   if (!PyArg_ParseTuple(args,"i",&i)) return NULL;
   fail = syscon_number_of_Laurent_terms(i,&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_syscon_retrieve_term ( PyObject *self, PyObject *args )
{
   int fail, *exp;
   int i,j,n,k;
   double *c;
   static PyObject *a;

   initialize();
   if(!PyArg_ParseTuple(args, "iiiid", &i,&j,&n,&exp,&c)) return NULL;

   exp = (int *)malloc(n * sizeof(int));
   c   = (double *)malloc(2*sizeof(double));

   fail = syscon_retrieve_term(i,j,n,exp,c);
     
   a = Py_BuildValue("i", fail);	 
   
   free (exp);
   free (c);

   return a;           
}

static PyObject *py2c_syscon_store_polynomial
 ( PyObject *self, PyObject *args )
{      
   int fail,n,nc,k;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iiis",&nc,&n,&k,&p)) return NULL;
   fail = syscon_store_polynomial(nc,n,k,p);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_store_dobldobl_polynomial
 ( PyObject *self, PyObject *args )
{      
   int fail,n,nc,k;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iiis",&nc,&n,&k,&p)) return NULL;
   fail = syscon_store_dobldobl_polynomial(nc,n,k,p);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_store_quaddobl_polynomial
 ( PyObject *self, PyObject *args )
{      
   int fail,n,nc,k;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iiis",&nc,&n,&k,&p)) return NULL;
   fail = syscon_store_quaddobl_polynomial(nc,n,k,p);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_store_multprec_polynomial
 ( PyObject *self, PyObject *args )
{      
   int fail,n,nc,k,dp;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iiiis",&nc,&n,&k,&dp,&p)) return NULL;
   fail = syscon_store_multprec_polynomial(nc,n,k,dp,p);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_load_polynomial
 ( PyObject *self, PyObject *args )
{      
   int fail,nc,k;
   char p[25600];  /* must be computed or retrieved !!!! */
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_load_polynomial(k,&nc,p);
                 
   return Py_BuildValue("s",p);
}

static PyObject *py2c_syscon_load_dobldobl_polynomial
 ( PyObject *self, PyObject *args )
{      
   int fail,nc,k;
   char p[51200];  /* must be computed or retrieved !!!! */
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_load_dobldobl_polynomial(k,&nc,p);
                 
   return Py_BuildValue("s",p);
}

static PyObject *py2c_syscon_load_quaddobl_polynomial
 ( PyObject *self, PyObject *args )
{      
   int fail,nc,k;
   char p[102400];  /* must be computed or retrieved !!!! */
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_load_quaddobl_polynomial(k,&nc,p);
                 
   return Py_BuildValue("s",p);
}

static PyObject *py2c_syscon_load_multprec_polynomial
 ( PyObject *self, PyObject *args )
{      
   int fail,nc,k;
   char p[102400];  /* must be computed or retrieved !!!! */
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_load_multprec_polynomial(k,&nc,p);
                 
   return Py_BuildValue("s",p);
}

static PyObject *py2c_syscon_store_Laurential
 ( PyObject *self, PyObject *args )
{      
   int fail,n,nc,k;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iiis",&nc,&n,&k,&p)) return NULL;
   fail = syscon_store_Laurential(nc,n,k,p);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_store_dobldobl_Laurential
 ( PyObject *self, PyObject *args )
{      
   int fail,n,nc,k;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iiis",&nc,&n,&k,&p)) return NULL;
   fail = syscon_store_dobldobl_Laurential(nc,n,k,p);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_store_quaddobl_Laurential
 ( PyObject *self, PyObject *args )
{      
   int fail,n,nc,k;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iiis",&nc,&n,&k,&p)) return NULL;
   fail = syscon_store_quaddobl_Laurential(nc,n,k,p);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_store_multprec_Laurential
 ( PyObject *self, PyObject *args )
{      
   int fail,n,nc,k,size;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iiiis",&nc,&n,&k,&size,&p)) return NULL;
   fail = syscon_store_multprec_Laurential(nc,n,k,size,p);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_load_standard_Laurential
 ( PyObject *self, PyObject *args )
{      
   int fail,nc,k;
   char p[25600];  /* must be computed or retrieved !!!! */
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_load_standard_Laurential(k,&nc,p);
                 
   return Py_BuildValue("s",p);
}

static PyObject *py2c_syscon_load_dobldobl_Laurential
 ( PyObject *self, PyObject *args )
{      
   int fail,nc,k;
   char p[51200];  /* must be computed or retrieved !!!! */
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_load_dobldobl_Laurential(k,&nc,p);
                 
   return Py_BuildValue("s",p);
}

static PyObject *py2c_syscon_load_quaddobl_Laurential
 ( PyObject *self, PyObject *args )
{      
   int fail,nc,k;
   char p[102400];  /* must be computed or retrieved !!!! */
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_load_quaddobl_Laurential(k,&nc,p);
                 
   return Py_BuildValue("s",p);
}

static PyObject *py2c_syscon_load_multprec_Laurential
 ( PyObject *self, PyObject *args )
{      
   int fail,nc,k;
   char p[102400];  /* must be computed or retrieved !!!! */
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_load_multprec_Laurential(k,&nc,p);
                 
   return Py_BuildValue("s",p);
}

static PyObject *py2c_syscon_total_degree ( PyObject *self, PyObject *args )
{
   int fail,totdeg;
                 
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = syscon_total_degree(&totdeg);
                 
   return Py_BuildValue("i",totdeg);
}

static PyObject *py2c_syscon_standard_drop_variable_by_index
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_standard_drop_variable_by_index(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_standard_drop_variable_by_name
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *s;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&s)) return NULL;
   fail = syscon_standard_drop_variable_by_name(nc,s);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_dobldobl_drop_variable_by_index
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_dobldobl_drop_variable_by_index(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_dobldobl_drop_variable_by_name
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *s;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&s)) return NULL;
   fail = syscon_dobldobl_drop_variable_by_name(nc,s);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_quaddobl_drop_variable_by_index
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_quaddobl_drop_variable_by_index(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_quaddobl_drop_variable_by_name
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *s;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&s)) return NULL;
   fail = syscon_quaddobl_drop_variable_by_name(nc,s);

   return Py_BuildValue("i",fail);
}

/* wrapping functions in solcon.h starts from here */

static PyObject *py2c_solcon_read_solutions ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_read_solutions();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_read_dobldobl_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_read_dobldobl_solutions();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_read_quaddobl_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_read_quaddobl_solutions();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_read_multprec_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_read_multprec_solutions();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_write_solutions ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_write_solutions();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_write_dobldobl_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_write_dobldobl_solutions();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_write_quaddobl_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_write_quaddobl_solutions();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_write_multprec_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_write_multprec_solutions();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_clear_solutions ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_clear_solutions();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_clear_dobldobl_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_clear_dobldobl_solutions();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_clear_quaddobl_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_clear_quaddobl_solutions();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_clear_multprec_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_clear_multprec_solutions();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_open_solution_input_file
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_open_solution_input_file();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_length_solution_string
 ( PyObject *self, PyObject *args )
{
   int fail,i,number;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&i)) return NULL;
   fail = solcon_length_solution_string(i,&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_solcon_length_dobldobl_solution_string
 ( PyObject *self, PyObject *args )
{
   int fail,i,number;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&i)) return NULL;
   fail = solcon_length_dobldobl_solution_string(i,&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_solcon_length_quaddobl_solution_string
 ( PyObject *self, PyObject *args )
{
   int fail,i,number;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&i)) return NULL;
   fail = solcon_length_quaddobl_solution_string(i,&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_solcon_length_multprec_solution_string
 ( PyObject *self, PyObject *args )
{
   int fail,i,number;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&i)) return NULL;
   fail = solcon_length_multprec_solution_string(i,&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_solcon_write_solution_string
 ( PyObject *self, PyObject *args )
{      
   char fail;
   int n,k;
   char *p;
   static PyObject* a;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"ii",&n,&k)) return NULL;
   p = (char*)malloc((k+1)*sizeof(char));
   fail = solcon_write_solution_string(n,k,p);
                 
   a = Py_BuildValue("s",p);

   free(p);
   
   return a;
}

static PyObject *py2c_solcon_write_dobldobl_solution_string
 ( PyObject *self, PyObject *args )
{      
   char fail;
   int n,k;
   char *p;
   static PyObject* a;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"ii",&n,&k)) return NULL;
   p = (char*)malloc((k+1)*sizeof(char));
   fail = solcon_write_dobldobl_solution_string(n,k,p);
                 
   a = Py_BuildValue("s",p);

   free(p);
   
   return a;
}

static PyObject *py2c_solcon_write_quaddobl_solution_string
 ( PyObject *self, PyObject *args )
{      
   char fail;
   int n,k;
   char *p;
   static PyObject* a;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"ii",&n,&k)) return NULL;
   p = (char*)malloc((k+1)*sizeof(char));
   fail = solcon_write_quaddobl_solution_string(n,k,p);
                 
   a = Py_BuildValue("s",p);

   free(p);
   
   return a;
}

static PyObject *py2c_solcon_write_multprec_solution_string
 ( PyObject *self, PyObject *args )
{      
   char fail;
   int n,k;
   char *p;
   static PyObject* a;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"ii",&n,&k)) return NULL;
   p = (char*)malloc((k+1)*sizeof(char));
   fail = solcon_write_multprec_solution_string(n,k,p);
                 
   a = Py_BuildValue("s",p);

   free(p);
   
   return a;
}

static PyObject *py2c_solcon_append_solution_string
 ( PyObject *self, PyObject *args )
{      
   char fail;
   int n,k;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iis",&n,&k,&p)) return NULL;
   fail = solcon_append_solution_string(n,k,p);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_append_dobldobl_solution_string
 ( PyObject *self, PyObject *args )
{      
   char fail;
   int n,k;
   char *p;
   
   initialize();
   if(!PyArg_ParseTuple(args,"iis",&n,&k,&p)) return NULL;
   fail = solcon_append_dobldobl_solution_string(n,k,p);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_append_quaddobl_solution_string
 ( PyObject *self, PyObject *args )
{      
   char fail;
   int n,k;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iis",&n,&k,&p)) return NULL;
   fail = solcon_append_quaddobl_solution_string(n,k,p);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_append_multprec_solution_string
 ( PyObject *self, PyObject *args )
{      
   char fail;
   int n,k;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iis",&n,&k,&p)) return NULL;
   fail = solcon_append_multprec_solution_string(n,k,p);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_number_of_solutions
 ( PyObject *self, PyObject *args )
{      
   int result,n;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   result = solcon_number_of_solutions(&n);
                 
   return Py_BuildValue("i",n);
}

static PyObject *py2c_solcon_number_of_dobldobl_solutions
 ( PyObject *self, PyObject *args )
{      
   int result,n;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   result = solcon_number_of_dobldobl_solutions(&n);
                 
   return Py_BuildValue("i",n);
}

static PyObject *py2c_solcon_number_of_quaddobl_solutions
 ( PyObject *self, PyObject *args )
{      
   int result,n;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   result = solcon_number_of_quaddobl_solutions(&n);
                 
   return Py_BuildValue("i",n);
}

static PyObject *py2c_solcon_number_of_multprec_solutions
 ( PyObject *self, PyObject *args )
{      
   int result,n;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   result = solcon_number_of_multprec_solutions(&n);
                 
   return Py_BuildValue("i",n);
}

static PyObject *py2c_solcon_standard_drop_coordinate_by_index
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = solcon_standard_drop_coordinate_by_index(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_standard_drop_coordinate_by_name
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *s;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&s)) return NULL;
   fail = solcon_standard_drop_coordinate_by_name(nc,s);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_dobldobl_drop_coordinate_by_index
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = solcon_dobldobl_drop_coordinate_by_index(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_dobldobl_drop_coordinate_by_name
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *s;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&s)) return NULL;
   fail = solcon_dobldobl_drop_coordinate_by_name(nc,s);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_quaddobl_drop_coordinate_by_index
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = solcon_quaddobl_drop_coordinate_by_index(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_quaddobl_drop_coordinate_by_name
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *s;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&s)) return NULL;
   fail = solcon_quaddobl_drop_coordinate_by_name(nc,s);

   return Py_BuildValue("i",fail);
}

/* wrapping functions in product.h starts here */

static PyObject *py2c_product_supporting_set_structure 
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = supporting_set_structure();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_product_write_set_structure
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = write_set_structure();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_product_set_structure_string
 ( PyObject *self, PyObject *args )
{
   int fail;
   int size = 1024;
   char s[size];

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = set_structure_string(&size,s);
   s[size] = '\0';
              
   return Py_BuildValue("s",s);
}

static PyObject *py2c_product_parse_set_structure
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *s;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&s)) return NULL;
   fail = parse_set_structure(nc,s);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_product_is_set_structure_supporting
 ( PyObject *self, PyObject *args )
{
   int fail,result;

   initialize();
   if (!PyArg_ParseTuple(args,"i",&result)) return NULL;
   fail = is_set_structure_supporting(&result);

   return Py_BuildValue("i",result);
}

static PyObject *py2c_product_linear_product_root_count
 ( PyObject *self, PyObject *args )
{
   int fail,r;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = linear_product_root_count(&r);

   return Py_BuildValue("i",r);
}

static PyObject *py2c_product_random_linear_product_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = random_linear_product_system();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_product_solve_linear_product_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = solve_linear_product_system();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_product_clear_set_structure
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = clear_set_structure();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_product_m_homogeneous_Bezout_number
 ( PyObject *self, PyObject *args )
{
   int fail,mbz,ncp;
   char partition[256];

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = m_homogeneous_Bezout_number(&mbz,&ncp,partition);

   return Py_BuildValue("(i,s)",mbz,partition);
}

static PyObject *py2c_product_m_partition_Bezout_number
 ( PyObject *self, PyObject *args )
{
   int fail,mbz,ncp;
   char *partition;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&ncp,&partition)) return NULL;
   fail = m_partition_Bezout_number(&mbz,ncp,partition);

   return Py_BuildValue("i",mbz);
}

static PyObject *py2c_product_m_homogeneous_start_system
 ( PyObject *self, PyObject *args )
{
   int fail,mbz,ncp;
   char *partition;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&ncp,&partition)) return NULL;
   fail = m_homogeneous_start_system(ncp,partition);

   return Py_BuildValue("i",mbz);
}

/* wrapping functions in celcon.h starts here */

static PyObject *py2c_celcon_set_type_of_mixture
 ( PyObject *self, PyObject *args )
{
   int fail,r,cnt;
   char *strmix;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&r,&strmix)) return NULL;

   cnt = itemcount(strmix);
   {
      int mix[cnt];

      str2intlist(cnt,strmix,mix);
      fail = celcon_set_type_of_mixture(r,mix);
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_type_of_mixture ( PyObject *self, PyObject *args )
{
   int fail,r,nbc;
   int mix[64];
   char strmix[256];

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;

   fail = celcon_type_of_mixture(&r,mix);
   nbc = intlist2str(r,mix,strmix);

   return Py_BuildValue("s",strmix);
}

static PyObject *py2c_celcon_append_lifted_point
 ( PyObject *self, PyObject *args )
{
   int fail,dim,ind,cnt;
   char *strpoint;

   initialize();
   if(!PyArg_ParseTuple(args,"iis",&dim,&ind,&strpoint)) return NULL;

   cnt = itemcount(strpoint);
   {
      double point[cnt];

      str2dbllist(cnt,strpoint,point);
      fail = celcon_append_lifted_point(dim,ind,point);
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_number_of_cells ( PyObject *self, PyObject *args )
{
   int fail,length;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_number_of_cells(&length);

   return Py_BuildValue("i",length);
}

static PyObject *py2c_celcon_create_random_coefficient_system 
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_create_random_coefficient_system();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_dobldobl_random_coefficient_system 
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_dobldobl_random_coefficient_system();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_quaddobl_random_coefficient_system 
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_quaddobl_random_coefficient_system();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_copy_into_systems_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_copy_into_systems_container();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_copy_into_dobldobl_systems_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_copy_into_dobldobl_systems_container();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_copy_into_quaddobl_systems_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_copy_into_quaddobl_systems_container();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_create_polyhedral_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;
  
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_create_polyhedral_homotopy();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_dobldobl_polyhedral_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;
  
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_dobldobl_polyhedral_homotopy();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_quaddobl_polyhedral_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;
  
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_quaddobl_polyhedral_homotopy();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_solve_start_system
 ( PyObject *self, PyObject *args )
{
   int fail,k,nb;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = celcon_solve_start_system(k,&nb);

   return Py_BuildValue("i",nb);
}

static PyObject *py2c_celcon_solve_dobldobl_start_system
 ( PyObject *self, PyObject *args )
{
   int fail,k,nb;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = celcon_solve_dobldobl_start_system(k,&nb);

   return Py_BuildValue("i",nb);
}

static PyObject *py2c_celcon_solve_quaddobl_start_system
 ( PyObject *self, PyObject *args )
{
   int fail,k,nb;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = celcon_solve_quaddobl_start_system(k,&nb);

   return Py_BuildValue("i",nb);
}

static PyObject *py2c_celcon_track_solution_path
 ( PyObject *self, PyObject *args )
{
   int fail,k,i,otp;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&k,&i,&otp)) return NULL;
   fail = celcon_track_solution_path(k,i,otp);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_track_dobldobl_solution_path
 ( PyObject *self, PyObject *args )
{
   int fail,k,i,otp;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&k,&i,&otp)) return NULL;
   fail = celcon_track_dobldobl_solution_path(k,i,otp);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_track_quaddobl_solution_path
 ( PyObject *self, PyObject *args )
{
   int fail,k,i,otp;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&k,&i,&otp)) return NULL;
   fail = celcon_track_quaddobl_solution_path(k,i,otp);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_copy_target_solution_to_container
 ( PyObject *self, PyObject *args )
{
   int fail,k,i;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&k,&i)) return NULL;
   fail = celcon_copy_target_solution_to_container(k,i);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_copy_target_dobldobl_solution_to_container
 ( PyObject *self, PyObject *args )
{
   int fail,k,i;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&k,&i)) return NULL;
   fail = celcon_copy_target_dobldobl_solution_to_container(k,i);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_copy_target_quaddobl_solution_to_container
 ( PyObject *self, PyObject *args )
{
   int fail,k,i;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&k,&i)) return NULL;
   fail = celcon_copy_target_quaddobl_solution_to_container(k,i);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_permute_system
 ( PyObject *self, PyObject *args )
{
   int fail,k,i;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_permute_system();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_permute_dobldobl_system
 ( PyObject *self, PyObject *args )
{
   int fail,k,i;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_permute_dobldobl_system();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_permute_quaddobl_system
 ( PyObject *self, PyObject *args )
{
   int fail,k,i;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_permute_quaddobl_system();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_celcon_clear_container
 ( PyObject *self, PyObject *args )
{
   int fail,k,i;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_clear_mixed_cell_configuration();

   return Py_BuildValue("i",fail);
}

/* wrapping functions to manipulate algebraic sets */

static PyObject *py2c_embed_system ( PyObject *self, PyObject *args )
{
   int d,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&d)) return NULL;
   fail = embed_system(d);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_standard_cascade_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = create_cascade_homotopy();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_dobldobl_cascade_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = create_dobldobl_cascade_homotopy();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_quaddobl_cascade_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = create_quaddobl_cascade_homotopy();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_set_to_mute ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = set_state_to_silent();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_define_output_file_with_string
 ( PyObject *self, PyObject *args )
{
   int n,fail;
   char *name;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&n,&name)) return NULL;
   fail = define_output_file_with_string(n,name);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_assign_labels ( PyObject *self, PyObject *args )
{
   int n,nbsols,i,j,m,fail;
   double x[2*n+5];

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&n,&nbsols)) return NULL;

   for(i=1; i<=nbsols; i++)
   {
      fail = solcon_retrieve_solution(n,i,&m,x);
      m = i;
      fail = solcon_replace_solution(n,i,m,x);
   }

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_initialize_sampler
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = initialize_sampler(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_initialize_monodromy
 ( PyObject *self, PyObject *args )
{
   int n,d,k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&n,&d,&k)) return NULL;
   fail = initialize_monodromy(n,d,k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_store_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = store_solutions();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_restore_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = restore_solutions();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_track_paths ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = track_paths();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_swap_slices ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = swap_slices();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_new_slices ( PyObject *self, PyObject *args )
{
   int k,n,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&k,&n)) return NULL;
   fail = new_slices(k,n);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_set_trace_slice
 ( PyObject *self, PyObject *args )
{
   int first,fail;
   double r[2];

   initialize();
   if(!PyArg_ParseTuple(args,"i",&first)) return NULL;

   if(first == 1)                  /* determine constant coefficient */
   { 
      r[0] = -1.0; r[1] = 0.0;
   }
   else
   {
      r[0] = +1.0; r[1] = 0.0;
   }
   fail = assign_coefficient_of_slice(0,0,r);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_store_gammas ( PyObject *self, PyObject *args )
{
   int n,i,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&n)) return NULL;

   {
      double re_gamma[n];
      double im_gamma[n];
    
      for(i=0; i<n; i++)
         random_complex(&re_gamma[i],&im_gamma[i]);
      fail = store_gamma(n,re_gamma,im_gamma);
   }

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_permutation_after_loop
 ( PyObject *self, PyObject *args )
{
   int d,fail,nb;
   char *result;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&d)) return NULL;
   {
      int i,permutation[d];
      char s[d*10];

      fail = permutation_after_loop(d,permutation);
   /* printf("the permutation :");
      for(i=0; i<d; i++) printf(" %d",permutation[i]);
      printf("\n"); */

      nb = list2str(d,permutation,s);
      result = (char*)calloc(nb,sizeof(char));
      for(i=0; i<nb; i++) result[i] = s[i];
   }
   /* printf("the resulting string : %s\n",result); */
   return Py_BuildValue("s",result);
}

static PyObject *py2c_factor_update_decomposition
 ( PyObject *self, PyObject *args )
{
   int fail,i,d,nc;
   char *permutation;
   int done = 0;

   initialize();
   if(!PyArg_ParseTuple(args,"iis",&d,&nc,&permutation)) return NULL;
   {
      int nb,perm[d],nf[2];

   /* printf("updating with ");
      for(i=0; i<nc; i++) printf("%c",permutation[i]);
      printf("\n"); */

      nb = str2list(nc,permutation,perm);

   /* printf("after str2list :");
      for(i=0; i<nb; i++) printf(" %d",perm[i]);
      printf("\n"); */

      fail = update_decomposition(d,perm,nf,&done);
   /* printf("number of factors : %d -> %d\n",nf[0],nf[1]); */
   }
   return Py_BuildValue("i",done);
}

static PyObject *py2c_factor_number_of_components
 ( PyObject *self, PyObject *args )
{
   int fail,nf;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;

   fail = number_of_irreducible_factors(&nf);

   return Py_BuildValue("i",nf);
}

static PyObject *py2c_factor_witness_points_of_component
 ( PyObject *self, PyObject *args )
{
   int fail,totdeg,k;
   char *result;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&totdeg,&k)) return NULL;
   {
      int deg,nb,i;
      int w[totdeg];
      char s[10*totdeg];

      fail = witness_points_of_irreducible_factor(k,&deg,w);

      nb = list2str(deg,w,s);
      result = (char*)calloc(nb,sizeof(char));
      for(i=0; i<nb; i++) result[i] = s[i];
   }
   return Py_BuildValue("s",result);
}

static PyObject *py2c_factor_trace_sum_difference
 ( PyObject *self, PyObject *args )
{
   int d,k,nc,fail;
   char *ws;
   double tsd;

   initialize();
   if(!PyArg_ParseTuple(args,"iiis",&d,&k,&nc,&ws)) return NULL;
   {
      int i,nb,witset[d];

      nb = str2list(nc,ws,witset);
   /* printf("the witness points :");
      for(i=0; i<nb; i++) printf(" %d",witset[i]);
      printf("\n"); */

      fail = trace_sum_difference(nb,witset,&tsd);
   /* printf("trace sum difference : %.3e\n",tsd); */
   }
   return Py_BuildValue("d",tsd);
}

static PyObject *py2c_witness_set_of_hypersurface
 ( PyObject *self, PyObject *args )
{
   int fail,nv,nc;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iis",&nv,&nc,&p)) return NULL;
   fail = witness_set_of_hypersurface(nv,nc,p);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_create_diagonal_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,a,b;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&a,&b)) return NULL;
   fail = create_diagonal_homotopy(a,b);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_start_diagonal_cascade_solutions
 ( PyObject *self, PyObject *args )
{
   int fail,a,b;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&a,&b)) return NULL;
   fail = start_diagonal_cascade_solutions(a,b);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_extrinsic_top_diagonal_dimension
 ( PyObject *self, PyObject *args )
{
   int fail,a,b,n1,n2,d;

   initialize();
   if(!PyArg_ParseTuple(args,"iiii",&n1,&n2,&a,&b)) return NULL;
   fail = extrinsic_top_diagonal_dimension(n1,n2,a,b,&d);

   if(fail == 0)
      return Py_BuildValue("i",d);
   else
      return Py_BuildValue("i",fail);
}

static PyObject *py2c_collapse_diagonal ( PyObject *self, PyObject *args )
{
   int fail,k,d;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&k,&d)) return NULL;
   fail = collapse_diagonal(k,d);

   return Py_BuildValue("i",fail);
}

/* wrapping of Pieri and Littlewood-Richardson homotopies */

static PyObject *py2c_schubert_pieri_count ( PyObject *self, PyObject *args )
{
   int fail,m,p,q,r;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&m,&p,&q)) return NULL;
   fail = Pieri_root_count(m,p,q,&r);

   return Py_BuildValue("i",r);
}

static PyObject *py2c_schubert_resolve_conditions
 ( PyObject *self, PyObject *args )
{
   int n,k,nbc,nc,fail,r,vrb;
   char *cond;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiisi",&n,&k,&nbc,&nc,&cond,&vrb)) return NULL;
/*
   printf("the number of characters : %d\n", nc);
   printf("the conditions : %s\n", cond);
   printf("the conditions parsed : ");
*/
   {
      int cds[k*nbc];
      int pos = 0;
      int idx = 0;
      while((idx < k*nbc) && (pos < nc))
      {
         while(cond[pos] == ' ' && pos < nc) pos++;
         if(pos > nc) break;
         cds[idx] = 0;
         while(cond[pos] != ' ')
         {
            if(cond[pos] == '\0') break;
            cds[idx] = cds[idx]*10 + (cond[pos] - '0');
            pos = pos + 1;
            if(pos >= nc) break;
         }
         /* printf(" %d", cds[idx]); */
         idx = idx + 1;
      }
      fail = resolve_Schubert_conditions(n,k,nbc,cds,vrb,&r);
   }
   return Py_BuildValue("i",r);
}

static PyObject *py2c_schubert_littlewood_richardson_homotopies
 ( PyObject *self, PyObject *args )
{
   int i,n,k,nbc,nc,fail,r,vrb,szn;
   char *cond;
   char *name;

   initialize();
   if(!PyArg_ParseTuple
         (args,"iiiisiis",&n,&k,&nbc,&nc,&cond,&vrb,&szn,&name)) return NULL;
/*
   printf("name of the output file : %s\n", name);
   printf("the number of characters : %d\n", nc);
   printf("the conditions : %s\n", cond);
   printf("the conditions parsed : ");
*/
   {
      int cds[k*nbc];
      int pos = 0;
      int idx = 0;
      while((idx < k*nbc) && (pos < nc))
      {
         while(cond[pos] == ' ' && pos < nc) pos++;
         if(pos > nc) break;
         cds[idx] = 0;
         while(cond[pos] != ' ')
         {
            if(cond[pos] == '\0') break;
            cds[idx] = cds[idx]*10 + (cond[pos] - '0');
            pos = pos + 1;
            if(pos >= nc) break;
         }
         /* printf(" %d", cds[idx]); */
         idx = idx + 1;
      }
      const int fgsize = 2*(nbc-2)*n*n;
      double fg[fgsize];
      char stfg[fgsize*24+2];
      fail = Littlewood_Richardson_homotopies(n,k,nbc,cds,vrb,szn,name,&r,fg);
      stfg[0] = '[';
      idx = 1;
      for(i=0; i<fgsize; i++)
      {
         sprintf(&stfg[idx],"%+.16e",fg[i]); idx = idx + 23;
         stfg[idx] = ','; idx = idx + 1;
      }
      stfg[idx-1] = ']';
      stfg[idx] = '\0';
      /* printf("The string with flag coefficients :\n%s\n", stfg); */

      return Py_BuildValue("(i,s)",r,stfg);
   }
}

static PyObject *py2c_schubert_localization_poset
 ( PyObject *self, PyObject *args )
{
   int fail,m,p,q,nc;
   const int buffer_size = 10240; /* must be calculated based on m,p,q */
   char ps[buffer_size];

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&m,&p,&q)) return NULL;
   fail = localization_poset(m,p,q,&nc,ps);

   return Py_BuildValue("s",ps);
}

static PyObject *py2c_schubert_pieri_homotopies
 ( PyObject *self, PyObject *args )
{
   int fail,m,p,q,nc,r;
   char *A;
   char *pts;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiiss",&m,&p,&q,&nc,&A,&pts)) return NULL;
   /* printf("receiving %d characters in py2c_pieri_homotopies\n",nc); */
   /* printf("the last character is %c\n",A[nc-1]); */
   fail = run_Pieri_homotopies(m,p,q,nc,&r,A,pts);

   return Py_BuildValue("i",r);
}

static PyObject *py2c_schubert_osculating_planes
 ( PyObject *self, PyObject *args )
{
   int fail,m,p,q,nc;
   char *pts;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiis",&m,&p,&q,&nc,&pts)) return NULL;

   {
      const int d = m+p;
      const int n = m*p + q*d;
      const int size = n*m*d*30;
      char planes[size];

      /* printf("passing %d characters through ...\n",nc);
      printf("size = %d\n",size); */
      fail = real_osculating_planes(m,p,q,&nc,pts,planes);

      return Py_BuildValue("s",planes);
   }
}

static PyObject *py2c_schubert_pieri_system ( PyObject *self, PyObject *args )
{
   int fail,m,p,q,nc,r;
   char *A;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiisi",&m,&p,&q,&nc,&A,&r)) return NULL;
   fail = Pieri_polynomial_system(m,p,q,nc,A,r);

   return Py_BuildValue("i",fail);
}

/* wrapping functions in mapcon.h starts from here */

static PyObject *py2c_mapcon_solve_system ( PyObject *self, PyObject *args )
{
   int fail,puretopdim;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&puretopdim)) return NULL;   
   fail = mapcon_solve_system(puretopdim);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_mapcon_write_maps ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = mapcon_write_maps();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_mapcon_clear_maps ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = mapcon_clear_maps();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_mapcon_top_dimension ( PyObject *self, PyObject *args )
{
   int fail,topdim;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = mapcon_top_dimension(&topdim);
              
   return Py_BuildValue("i",topdim);
}

static PyObject *py2c_mapcon_number_of_maps ( PyObject *self, PyObject *args )
{
   int fail,dim,nbmaps;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;   
   fail = mapcon_number_of_maps(dim,&nbmaps);
              
   return Py_BuildValue("i",nbmaps);
}

static PyObject *py2c_mapcon_degree_of_map ( PyObject *self, PyObject *args )
{
   int fail,dim,ind,deg;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&dim,&ind)) return NULL;   
   fail = mapcon_degree_of_map(dim,ind,&deg);
              
   return Py_BuildValue("i",deg);
}

static PyObject *py2c_mapcon_coefficients_of_map
 ( PyObject *self, PyObject *args )
{
   int fail,dim,ind,nbvar,i;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&dim,&ind,&nbvar)) return NULL;   

   double *cff;
   cff = (double*)calloc(2*nbvar,sizeof(double));

   fail = mapcon_coefficients_of_map(dim,ind,nbvar,cff);

   double re_cff,im_cff;
   PyObject *result, *item;
   result = PyList_New(nbvar);
   for(i=0; i<nbvar; i++)
   {
       re_cff = cff[2*i];
       im_cff = cff[2*i+1];
       item = PyComplex_FromDoubles(re_cff,im_cff);
       PyList_SET_ITEM(result,i,item);
   }
   free(cff);

   return result;
}

static PyObject *py2c_mapcon_exponents_of_map
 ( PyObject *self, PyObject *args )
{
   int fail,dim,ind,nbvar,i;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&dim,&ind,&nbvar)) return NULL;   

   int *exp;
   exp = (int*)calloc(dim*nbvar,sizeof(int));

   fail = mapcon_exponents_of_map(dim,ind,nbvar,exp);

   PyObject *result, *item;
   result = PyList_New(dim*nbvar);
   for(i=0; i<dim*nbvar; i++)
   {
       item = PyInt_FromLong(exp[i]);
       PyList_SET_ITEM(result,i,item);
   }
   free(exp);
      
   return result;
}

/* wrapping functions in next_track starts below */

static PyObject *py2c_initialize_standard_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,fixed;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&fixed)) return NULL;
   fail = initialize_standard_homotopy(fixed);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_initialize_dobldobl_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,fixed;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&fixed)) return NULL;
   fail = initialize_dobldobl_homotopy(fixed);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_initialize_quaddobl_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,fixed;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&fixed)) return NULL;
   fail = initialize_quaddobl_homotopy(fixed);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_initialize_multprec_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,deci,fixed;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&fixed,&deci)) return NULL;
   fail = initialize_multprec_homotopy(fixed,deci);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_initialize_standard_solution
 ( PyObject *self, PyObject *args )
{
   int fail,indsol;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&indsol)) return NULL;
   fail = initialize_standard_solution(indsol);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_initialize_dobldobl_solution
 ( PyObject *self, PyObject *args )
{
   int fail,indsol;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&indsol)) return NULL;
   fail = initialize_dobldobl_solution(indsol);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_initialize_quaddobl_solution
 ( PyObject *self, PyObject *args )
{
   int fail,indsol;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&indsol)) return NULL;
   fail = initialize_quaddobl_solution(indsol);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_initialize_multprec_solution
 ( PyObject *self, PyObject *args )
{
   int fail,indsol;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&indsol)) return NULL;
   fail = initialize_multprec_solution(indsol);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_next_standard_solution
 ( PyObject *self, PyObject *args )
{
   int fail,indsol;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&indsol)) return NULL;
   fail = next_standard_solution(indsol);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_next_dobldobl_solution
 ( PyObject *self, PyObject *args )
{
   int fail,indsol;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&indsol)) return NULL;
   fail = next_dobldobl_solution(indsol);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_next_quaddobl_solution
 ( PyObject *self, PyObject *args )
{
   int fail,indsol;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&indsol)) return NULL;
   fail = next_quaddobl_solution(indsol);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_next_multprec_solution
 ( PyObject *self, PyObject *args )
{
   int fail,indsol;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&indsol)) return NULL;
   fail = next_multprec_solution(indsol);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_standard_tracker
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = clear_standard_tracker();
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_dobldobl_tracker
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = clear_dobldobl_tracker();
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_quaddobl_tracker
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = clear_quaddobl_tracker();
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_multprec_tracker
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = clear_multprec_tracker();
   return Py_BuildValue("i",fail);
}

static PyMethodDef phcpy2c_methods[] = 
{
   {"py2c_PHCpack_version_string", py2c_PHCpack_version_string,
     METH_VARARGS, "returns the version string of PHCpack"},
   {"py2c_set_seed", py2c_set_seed,
     METH_VARARGS, "sets the given seed for the random number generator"},

   {"py2c_read_standard_target_system", py2c_read_standard_target_system,
     METH_VARARGS, "reads the target system in standard doubles"},
   {"py2c_read_standard_target_system_from_file",
     py2c_read_standard_target_system_from_file,
     METH_VARARGS, "reads the target system in standard doubles from file"},
   {"py2c_read_standard_start_system", py2c_read_standard_start_system,
     METH_VARARGS, "reads the start system in standard doubles"},
   {"py2c_read_standard_start_system_from_file",
     py2c_read_standard_start_system_from_file,
     METH_VARARGS, "reads the start system in standard doubles from file"},
   {"py2c_read_dobldobl_target_system", py2c_read_dobldobl_target_system,
     METH_VARARGS, "reads the target system in double doubles"},
   {"py2c_read_dobldobl_target_system_from_file",
     py2c_read_dobldobl_target_system_from_file,
     METH_VARARGS, "reads the target system in double doubles from file"},
   {"py2c_read_dobldobl_start_system", py2c_read_dobldobl_start_system,
     METH_VARARGS, "reads the start system in double doubles"},
   {"py2c_read_dobldobl_start_system_from_file",
     py2c_read_dobldobl_start_system_from_file,
     METH_VARARGS, "reads the start system in double doubles from file"},
   {"py2c_read_quaddobl_target_system", py2c_read_quaddobl_target_system,
     METH_VARARGS, "reads the target system in double doubles"},
   {"py2c_read_quaddobl_target_system_from_file",
     py2c_read_quaddobl_target_system_from_file,
     METH_VARARGS, "reads the target system in double doubles from file"},
   {"py2c_read_quaddobl_start_system", py2c_read_quaddobl_start_system,
     METH_VARARGS, "reads the start system in double doubles"},
   {"py2c_read_quaddobl_start_system_from_file",
     py2c_read_quaddobl_start_system_from_file,
     METH_VARARGS, "reads the start system in double doubles from file"},
   {"py2c_define_output_file", py2c_define_output_file,
     METH_VARARGS, "defines the output file"},
   {"py2c_write_standard_target_system", py2c_write_standard_target_system,
     METH_VARARGS, "writes the target system in standard doubles to file"},
   {"py2c_write_standard_start_system", py2c_write_standard_start_system,
     METH_VARARGS, "writes the start system in standard doubles to file"},
   {"py2c_write_start_solutions", py2c_write_start_solutions,
    METH_VARARGS, "writes the start solution to file"},
   {"py2c_copy_target_system_to_container",
     py2c_copy_target_system_to_container,
    METH_VARARGS, "copies target system to container"},
   {"py2c_copy_dobldobl_target_system_to_container",
     py2c_copy_dobldobl_target_system_to_container,
    METH_VARARGS, "copies double double target system to container"},
   {"py2c_copy_quaddobl_target_system_to_container",
     py2c_copy_quaddobl_target_system_to_container,
    METH_VARARGS, "copies quad double target system to container"},
   {"py2c_copy_multprec_target_system_to_container",
     py2c_copy_multprec_target_system_to_container,
    METH_VARARGS, "copies multiprecision target system to container"},
   {"py2c_copy_container_to_target_system",
     py2c_copy_container_to_target_system,
    METH_VARARGS, "copies system in container to target system"},
   {"py2c_copy_dobldobl_container_to_target_system",
     py2c_copy_dobldobl_container_to_target_system,
    METH_VARARGS,"copies system in double double container to target system"},
   {"py2c_copy_quaddobl_container_to_target_system",
     py2c_copy_quaddobl_container_to_target_system,
    METH_VARARGS,"copies system in quad double container to target system"},
   {"py2c_copy_multprec_container_to_target_system",
     py2c_copy_multprec_container_to_target_system,
    METH_VARARGS,"copies system in multiprecision container to target system"},
   {"py2c_copy_start_system_to_container",
     py2c_copy_start_system_to_container,
    METH_VARARGS, "copies start system to container"},
   {"py2c_copy_dobldobl_start_system_to_container",
     py2c_copy_dobldobl_start_system_to_container,
    METH_VARARGS, "copies double double start system to container"},
   {"py2c_copy_quaddobl_start_system_to_container",
     py2c_copy_quaddobl_start_system_to_container,
    METH_VARARGS, "copies quad double start system to container"},
   {"py2c_copy_multprec_start_system_to_container",
     py2c_copy_multprec_start_system_to_container,
    METH_VARARGS, "copies multiprecision start system to container"},
   {"py2c_copy_container_to_start_system",
     py2c_copy_container_to_start_system,
    METH_VARARGS, "copies system in container to start system"},
   {"py2c_copy_dobldobl_container_to_start_system",
     py2c_copy_dobldobl_container_to_start_system,
    METH_VARARGS,"copies system in double double container to start system"},
   {"py2c_copy_quaddobl_container_to_start_system",
     py2c_copy_quaddobl_container_to_start_system,
    METH_VARARGS,"copies system in quad double container to start system"},
   {"py2c_copy_multprec_container_to_start_system",
     py2c_copy_multprec_container_to_start_system,
    METH_VARARGS,"copies system in multiprecision container to start system"},
   {"py2c_create_homotopy", py2c_create_homotopy,
    METH_VARARGS, "creates a homotopy between the target and start system"},
   {"py2c_create_homotopy_with_gamma", py2c_create_homotopy_with_gamma,
    METH_VARARGS, "makes a homotopy with given gamma constant"},
   {"py2c_create_dobldobl_homotopy", py2c_create_dobldobl_homotopy,
    METH_VARARGS,
    "creates a double double homotopy between the target and start system"},
   {"py2c_create_dobldobl_homotopy_with_gamma",
     py2c_create_dobldobl_homotopy_with_gamma,
    METH_VARARGS, "makes a double double homotopy with given gamma constant"},
   {"py2c_create_quaddobl_homotopy", py2c_create_quaddobl_homotopy,
    METH_VARARGS,
    "creates a quad double homotopy between the target and start system"},
   {"py2c_create_quaddobl_homotopy_with_gamma",
     py2c_create_quaddobl_homotopy_with_gamma,
    METH_VARARGS, "makes a quad double homotopy with given gamma constant"},
   {"py2c_create_multprec_homotopy", py2c_create_multprec_homotopy,
    METH_VARARGS,
    "creates a multiprecision homotopy between the target and start system"},
   {"py2c_create_multprec_homotopy_with_gamma",
     py2c_create_multprec_homotopy_with_gamma,
    METH_VARARGS, "makes a multiprecision homotopy with given gamma constant"},
   {"py2c_clear_homotopy", py2c_clear_homotopy,
    METH_VARARGS, "clears the homotopy between the target and start system"},
   {"py2c_clear_dobldobl_homotopy", py2c_clear_dobldobl_homotopy,
    METH_VARARGS,
    "clears the double double homotopy between the target and start system"},
   {"py2c_clear_quaddobl_homotopy", py2c_clear_quaddobl_homotopy,
    METH_VARARGS,
    "clears the quad double homotopy between the target and start system"},
   {"py2c_clear_multprec_homotopy", py2c_clear_multprec_homotopy,
    METH_VARARGS,
    "clears the multiprecision homotopy between the target and start system"},
   {"py2c_tune_continuation_parameters", py2c_tune_continuation_parameters,
    METH_VARARGS, "interactive tuning of the continuation parameters"},
   {"py2c_show_continuation_parameters", py2c_show_continuation_parameters,
    METH_VARARGS, "shows the current values of the continuation parameters"},
   {"py2c_autotune_continuation_parameters",
     py2c_autotune_continuation_parameters,
    METH_VARARGS, 
    "tune continuation parameters with level of difficulty and precision"},
   {"py2c_determine_output_during_continuation", 
     py2c_determine_output_during_continuation,
    METH_VARARGS, "determines the output"},
   {"py2c_solve_by_standard_homotopy_continuation",
     py2c_solve_by_standard_homotopy_continuation,
    METH_VARARGS, "homotopy continuation in standard double precision"},
   {"py2c_solve_by_dobldobl_homotopy_continuation",
     py2c_solve_by_dobldobl_homotopy_continuation,
    METH_VARARGS, "homotopy continuation in double double precision"},
   {"py2c_solve_by_quaddobl_homotopy_continuation",
     py2c_solve_by_quaddobl_homotopy_continuation,
    METH_VARARGS, "homotopy continuation in quad double precision"},
   {"py2c_solve_by_multprec_homotopy_continuation",
     py2c_solve_by_multprec_homotopy_continuation,
    METH_VARARGS, "homotopy continuation in multiprecision"},
   {"py2c_copy_target_solutions_to_container",
     py2c_copy_target_solutions_to_container,
    METH_VARARGS, "copies target solutions to container"},
   {"py2c_copy_dobldobl_target_solutions_to_container",
     py2c_copy_dobldobl_target_solutions_to_container,
    METH_VARARGS, "copies double double target solutions to container"},
   {"py2c_copy_quaddobl_target_solutions_to_container",
     py2c_copy_quaddobl_target_solutions_to_container,
    METH_VARARGS, "copies quad double target solutions to container"},
   {"py2c_copy_multprec_target_solutions_to_container",
     py2c_copy_multprec_target_solutions_to_container,
    METH_VARARGS, "copies multprecicision target solutions to container"},
   {"py2c_copy_container_to_target_solutions",
     py2c_copy_container_to_target_solutions,
    METH_VARARGS, "copies container to target solutions"},
   {"py2c_copy_dobldobl_container_to_target_solutions",
     py2c_copy_dobldobl_container_to_target_solutions,
    METH_VARARGS, "copies double double container to target solutions"},
   {"py2c_copy_quaddobl_container_to_target_solutions",
     py2c_copy_quaddobl_container_to_target_solutions,
    METH_VARARGS, "copies quad double container to target solutions"},
   {"py2c_copy_multprec_container_to_target_solutions",
     py2c_copy_multprec_container_to_target_solutions,
    METH_VARARGS, "copies multiprecision container to target solutions"},
   {"py2c_copy_start_solutions_to_container",
     py2c_copy_start_solutions_to_container,
    METH_VARARGS, "copies start solutions to container"},
   {"py2c_copy_dobldobl_start_solutions_to_container",
     py2c_copy_dobldobl_start_solutions_to_container,
    METH_VARARGS, "copies double double start solutions to container"},
   {"py2c_copy_quaddobl_start_solutions_to_container",
     py2c_copy_quaddobl_start_solutions_to_container,
    METH_VARARGS, "copies quad double start solutions to container"},
   {"py2c_copy_multprec_start_solutions_to_container",
     py2c_copy_multprec_start_solutions_to_container,
    METH_VARARGS, "copies multiprecision start solutions to container"},
   {"py2c_copy_container_to_start_solutions",
     py2c_copy_container_to_start_solutions,
    METH_VARARGS, "copies container to start solutions"},
   {"py2c_copy_dobldobl_container_to_start_solutions",
     py2c_copy_dobldobl_container_to_start_solutions,
    METH_VARARGS, "copies double double container to start solutions"},
   {"py2c_copy_quaddobl_container_to_start_solutions",
     py2c_copy_quaddobl_container_to_start_solutions,
    METH_VARARGS, "copies quad double container to start solutions"},
   {"py2c_copy_multprec_container_to_start_solutions",
     py2c_copy_multprec_container_to_start_solutions,
    METH_VARARGS, "copies multiprecision container to start solutions"},
   {"py2c_solve_system", py2c_solve_system,
    METH_VARARGS, "calls the blackbox solver on a polynomial"},
   {"py2c_solve_Laurent_system", py2c_solve_Laurent_system,
    METH_VARARGS, "calls the blackbox solver on a Laurent system"},
   {"py2c_mixed_volume", py2c_mixed_volume,
    METH_VARARGS, "returns the mixed volume of system in the container"},
   {"py2c_standard_deflate", py2c_standard_deflate,
    METH_VARARGS,
    "applies deflation in standard double precision with default settings"},
   {"py2c_dobldobl_deflate", py2c_dobldobl_deflate,
    METH_VARARGS,
    "applies deflation in double double precision with default settings"},
   {"py2c_quaddobl_deflate", py2c_quaddobl_deflate,
    METH_VARARGS,
    "applies deflation in quad double precision with default settings"},
   {"py2c_standard_Newton_step", py2c_standard_Newton_step,
    METH_VARARGS, "does one Newton step on container data"},
   {"py2c_dobldobl_Newton_step", py2c_dobldobl_Newton_step,
    METH_VARARGS, "does one Newton step on double double container data"},
   {"py2c_quaddobl_Newton_step", py2c_quaddobl_Newton_step,
    METH_VARARGS, "does one Newton step on quad double container data"},
   {"py2c_multprec_Newton_step", py2c_multprec_Newton_step,
    METH_VARARGS, "does one Newton step on multiprecision container data"},
   {"py2c_standard_Newton_Laurent_step", py2c_standard_Newton_Laurent_step,
    METH_VARARGS, 
    "does one Newton step on standard double Laurent systems"},
   {"py2c_dobldobl_Newton_Laurent_step", py2c_dobldobl_Newton_Laurent_step,
    METH_VARARGS,
    "does one Newton step on double double Laurent systems"},
   {"py2c_quaddobl_Newton_Laurent_step", py2c_quaddobl_Newton_Laurent_step,
    METH_VARARGS, 
    "does one Newton step on quad double Laurent systems"},
   {"py2c_multprec_Newton_Laurent_step", py2c_multprec_Newton_Laurent_step,
    METH_VARARGS,
    "does one Newton step on multiprecision Laurent systems"},
   {"py2c_usolve_standard", py2c_usolve_standard, METH_VARARGS,
    "solve polynomial in one variable in standard double precision"},
   {"py2c_usolve_dobldobl", py2c_usolve_dobldobl, METH_VARARGS,
    "solve polynomial in one variable in double double precision"},
   {"py2c_usolve_quaddobl", py2c_usolve_quaddobl, METH_VARARGS,
    "solve polynomial in one variable in quad double precision"},
   {"py2c_usolve_multprec", py2c_usolve_multprec, METH_VARARGS,
    "solve polynomial in one variable in arbitrary multiprecision"},
   {"py2c_giftwrap_planar", py2c_giftwrap_planar, METH_VARARGS,
    "giftwrapping to compute the convex hull for points in the plane"},
   {"py2c_giftwrap_convex_hull", py2c_giftwrap_convex_hull, METH_VARARGS,
    "giftwrapping to compute the convex hull for points 3- or 4-space"},
   {"py2c_giftwrap_number_of_facets", py2c_giftwrap_number_of_facets,
    METH_VARARGS,
    "returns the number of facets of the convex hull in  3- or 4-space"},
   {"py2c_giftwrap_retrieve_facet", py2c_giftwrap_retrieve_facet,
    METH_VARARGS,
    "returns the string representation of a facet in  3- or 4-space"},
   {"py2c_giftwrap_clear_3d_facets", py2c_giftwrap_clear_3d_facets,
    METH_VARARGS,
    "deallocates list of facets of convex hull stored in 3-space"},
   {"py2c_giftwrap_clear_4d_facets", py2c_giftwrap_clear_4d_facets,
    METH_VARARGS,
    "deallocates list of facets of convex hull stored in 4-space"},
   {"py2c_syscon_read_system", py2c_syscon_read_system,
    METH_VARARGS, "reads and puts the system in container"},
   {"py2c_syscon_read_Laurent_system", py2c_syscon_read_Laurent_system,
    METH_VARARGS, "reads and puts the Laurent system in container"},
   {"py2c_syscon_read_dobldobl_system", py2c_syscon_read_dobldobl_system,
    METH_VARARGS, "reads system with double doubles in container"},
   {"py2c_syscon_read_dobldobl_Laurent_system",
     py2c_syscon_read_dobldobl_Laurent_system, METH_VARARGS, 
    "reads Laurent system with double doubles in container"},
   {"py2c_syscon_read_quaddobl_system", py2c_syscon_read_quaddobl_system,
    METH_VARARGS, "reads system with quad doubles in container"},
   {"py2c_syscon_read_quaddobl_Laurent_system",
     py2c_syscon_read_quaddobl_Laurent_system, METH_VARARGS,
    "reads Laurent system with quad doubles in container"},
   {"py2c_syscon_read_multprec_system", py2c_syscon_read_multprec_system,
    METH_VARARGS, "reads system with multiprecision numbers in container"},
   {"py2c_syscon_read_multprec_Laurent_system",
     py2c_syscon_read_multprec_Laurent_system, METH_VARARGS,
    "reads Laurent system with multiprecision numbers in container"},
   {"py2c_syscon_random_system", py2c_syscon_random_system,
    METH_VARARGS, "stores a random system in the container"},
   {"py2c_syscon_write_system", py2c_syscon_write_system,
    METH_VARARGS, "writes Laurent system in container on screen"},
   {"py2c_syscon_write_Laurent_system", py2c_syscon_write_Laurent_system,
    METH_VARARGS, "writes Laurent system in container on screen"},
   {"py2c_syscon_write_dobldobl_system", py2c_syscon_write_dobldobl_system,
    METH_VARARGS, "writes system in double double container on screen"},
   {"py2c_syscon_write_dobldobl_Laurent_system",
     py2c_syscon_write_dobldobl_Laurent_system, METH_VARARGS,
    "writes Laurent system in double double container on screen"},
   {"py2c_syscon_write_quaddobl_system", py2c_syscon_write_quaddobl_system,
    METH_VARARGS, "writes system in quad double container on screen"},
   {"py2c_syscon_write_quaddobl_Laurent_system",
     py2c_syscon_write_quaddobl_Laurent_system, METH_VARARGS,
    "writes Laurent system in quad double container on screen"},
   {"py2c_syscon_write_multprec_system", py2c_syscon_write_multprec_system,
    METH_VARARGS, "writes system in multiprecision container on screen"},
   {"py2c_syscon_write_multprec_Laurent_system",
     py2c_syscon_write_multprec_Laurent_system, METH_VARARGS,
    "writes Laurent system in multiprecision container on screen"},
   {"py2c_syscon_clear_system", py2c_syscon_clear_system,
    METH_VARARGS, "clears the content of the systems container"},
   {"py2c_syscon_clear_Laurent_system", py2c_syscon_clear_Laurent_system,
    METH_VARARGS, "clears the content of the Laurent systems container"},
   {"py2c_syscon_clear_dobldobl_system", py2c_syscon_clear_dobldobl_system,
    METH_VARARGS, "clears the double double polynomial systems container"},
   {"py2c_syscon_clear_dobldobl_Laurent_system",
     py2c_syscon_clear_dobldobl_Laurent_system, METH_VARARGS,
    "clears the double double Laurent polynomial systems container"},
   {"py2c_syscon_clear_quaddobl_system", py2c_syscon_clear_quaddobl_system,
    METH_VARARGS, "clears the quad double polynomial systems container"},
   {"py2c_syscon_clear_quaddobl_Laurent_system",
     py2c_syscon_clear_quaddobl_Laurent_system, METH_VARARGS, 
    "clears the quad double Laurent polynomial systems container"},
   {"py2c_syscon_clear_multprec_system", py2c_syscon_clear_multprec_system,
    METH_VARARGS, "clears the multiprecision polynomial systems container"},
   {"py2c_syscon_clear_multprec_Laurent_system",
     py2c_syscon_clear_multprec_Laurent_system, METH_VARARGS,
    "clears the multiprecision Laurent polynomial systems container"},
   {"py2c_syscon_number_of_symbols", py2c_syscon_number_of_symbols,
    METH_VARARGS, "returns the number of symbols in the symbol table"},
   {"py2c_syscon_write_symbols", py2c_syscon_write_symbols,
    METH_VARARGS, "write symbols in the symbol table"},
   {"py2c_syscon_string_of_symbols", py2c_syscon_string_of_symbols,
    METH_VARARGS, "returns a string of symbols in the symbol table"},
   {"py2c_syscon_remove_symbol_name", py2c_syscon_remove_symbol_name,
    METH_VARARGS, "remove a symbol with given name from the symbol table"},
   {"py2c_syscon_clear_symbol_table", py2c_syscon_clear_symbol_table,
    METH_VARARGS, "clears the symbol table"},
   {"py2c_solcon_read_solutions", py2c_solcon_read_solutions,
    METH_VARARGS,
    "prompts for a file, reads solutions and puts them in container"},
   {"py2c_solcon_read_dobldobl_solutions", py2c_solcon_read_dobldobl_solutions,
    METH_VARARGS,
    "prompts for a file, reads double double solutions into the container"},
   {"py2c_solcon_read_quaddobl_solutions", py2c_solcon_read_quaddobl_solutions,
    METH_VARARGS,
    "prompts for a file, reads quad double solutions into the container"},
   {"py2c_solcon_read_multprec_solutions", py2c_solcon_read_multprec_solutions,
    METH_VARARGS,
    "prompts for a file, reads multiprecision solutions into the container"},
   {"py2c_solcon_write_solutions", py2c_solcon_write_solutions,
    METH_VARARGS, "writes the solutions in the container to screen"},
   {"py2c_solcon_write_dobldobl_solutions",
     py2c_solcon_write_dobldobl_solutions,
    METH_VARARGS,
    "writes the double double solutions in the container to screen"},
   {"py2c_solcon_write_quaddobl_solutions",
     py2c_solcon_write_quaddobl_solutions,
    METH_VARARGS,
    "writes the quad double solutions in the container to screen"},
   {"py2c_solcon_write_multprec_solutions",
     py2c_solcon_write_multprec_solutions,
    METH_VARARGS,
    "writes the multiprecision solutions in the container to screen"},
   {"py2c_solcon_clear_solutions", py2c_solcon_clear_solutions,
    METH_VARARGS, "clears the content of the solution container"},
   {"py2c_solcon_clear_dobldobl_solutions",
     py2c_solcon_clear_dobldobl_solutions,
    METH_VARARGS, "clears the quad double solution container"},
   {"py2c_solcon_clear_quaddobl_solutions",
     py2c_solcon_clear_quaddobl_solutions,
    METH_VARARGS, "clears the quad double solution container"},
   {"py2c_solcon_clear_multprec_solutions",
     py2c_solcon_clear_multprec_solutions,
    METH_VARARGS, "clears the multiprecision solution container"},
   {"py2c_solcon_open_solution_input_file",
    py2c_solcon_open_solution_input_file,
    METH_VARARGS, "prompts the user for the input file and opens it"},
   {"py2c_syscon_number_of_polynomials", py2c_syscon_number_of_polynomials,
    METH_VARARGS, 
    "returns the number of polynomials in the container"},
   {"py2c_syscon_number_of_dobldobl_polynomials",
     py2c_syscon_number_of_dobldobl_polynomials,
    METH_VARARGS, 
    "returns the number of polynomials with double double coefficients"},
   {"py2c_syscon_number_of_quaddobl_polynomials",
     py2c_syscon_number_of_quaddobl_polynomials,
    METH_VARARGS, 
    "returns the number of polynomials with quad double coefficients"},
   {"py2c_syscon_number_of_multprec_polynomials",
     py2c_syscon_number_of_multprec_polynomials,
    METH_VARARGS, 
    "returns the number of polynomials with multiprecision coefficients"},
   {"py2c_syscon_number_of_Laurentials", py2c_syscon_number_of_Laurentials,
    METH_VARARGS, 
    "returns the number of Laurent polynomials in the container"},
   {"py2c_syscon_number_of_dobldobl_Laurentials",
     py2c_syscon_number_of_dobldobl_Laurentials, METH_VARARGS, 
    "returns the number of double double Laurent polynomials"},
   {"py2c_syscon_number_of_quaddobl_Laurentials",
     py2c_syscon_number_of_quaddobl_Laurentials, METH_VARARGS, 
    "returns the number of quad double Laurent polynomials"},
   {"py2c_syscon_number_of_multprec_Laurentials",
     py2c_syscon_number_of_multprec_Laurentials, METH_VARARGS, 
    "returns the number of multiprecision Laurent polynomials"},
   {"py2c_syscon_initialize_number", py2c_syscon_initialize_number,
    METH_VARARGS, "initializes the container with the number of polynomials"},
   {"py2c_syscon_initialize_number_of_dobldobl_polynomials",
     py2c_syscon_initialize_number_of_dobldobl_polynomials,
    METH_VARARGS, "initializes the container of double double polynomials"},
   {"py2c_syscon_initialize_number_of_quaddobl_polynomials",
     py2c_syscon_initialize_number_of_quaddobl_polynomials,
    METH_VARARGS, "initializes the container of quad double polynomials"},
   {"py2c_syscon_initialize_number_of_multprec_polynomials",
     py2c_syscon_initialize_number_of_multprec_polynomials,
    METH_VARARGS, "initializes the container of multiprecision polynomials"},
   {"py2c_syscon_initialize_number_of_Laurentials",
     py2c_syscon_initialize_number_of_Laurentials, METH_VARARGS,
    "initializes the container with the number of Laurentials"},
   {"py2c_syscon_initialize_number_of_dobldobl_Laurentials",
     py2c_syscon_initialize_number_of_dobldobl_Laurentials, METH_VARARGS,
    "initializes the number of double double Laurentials"},
   {"py2c_syscon_initialize_number_of_quaddobl_Laurentials",
     py2c_syscon_initialize_number_of_quaddobl_Laurentials, METH_VARARGS,
    "initializes the number of quad double Laurentials"},
   {"py2c_syscon_initialize_number_of_multprec_Laurentials",
     py2c_syscon_initialize_number_of_multprec_Laurentials, METH_VARARGS,
    "initializes the number of multiprecisionLaurentials"},
   {"py2c_syscon_degree_of_polynomial", py2c_syscon_degree_of_polynomial,
    METH_VARARGS,
    "returns the degree of the k-th polynomial in the standard container"},
   {"py2c_syscon_degree_of_dobldobl_polynomial",
     py2c_syscon_degree_of_dobldobl_polynomial,
    METH_VARARGS,
    "returns the degree of the k-th polynomial with double double precision"},
   {"py2c_syscon_degree_of_quaddobl_polynomial",
     py2c_syscon_degree_of_quaddobl_polynomial,
    METH_VARARGS,
    "returns the degree of the k-th polynomial with quad double precision"},
   {"py2c_syscon_degree_of_multprec_polynomial",
     py2c_syscon_degree_of_multprec_polynomial,
    METH_VARARGS,
    "returns the degree of the k-th polynomial with multiprecision"},
   {"py2c_syscon_number_of_terms", py2c_syscon_number_of_terms,
    METH_VARARGS, "returns in nt the number of terms in the i-th polynomial"},
   {"py2c_syscon_number_of_Laurent_terms", py2c_syscon_number_of_Laurent_terms,
    METH_VARARGS, "returns in nt the number of terms in the i-th Laurential"},
   {"py2c_syscon_retrieve_term", py2c_syscon_retrieve_term,
    METH_VARARGS, "retrieves the j-th term of the i-th polynomial"},
   {"py2c_syscon_store_polynomial", py2c_syscon_store_polynomial,
    METH_VARARGS, "stores the k-th polynomial in the systems container"},
   {"py2c_syscon_store_dobldobl_polynomial",
     py2c_syscon_store_dobldobl_polynomial,
    METH_VARARGS,
    "stores the k-th polynomial with double double coefficients"},
   {"py2c_syscon_store_quaddobl_polynomial",
     py2c_syscon_store_quaddobl_polynomial,
    METH_VARARGS,
    "stores the k-th polynomial with quad double coefficients"},
   {"py2c_syscon_store_multprec_polynomial",
     py2c_syscon_store_multprec_polynomial,
    METH_VARARGS,
    "stores the k-th polynomial with multiprecision coefficients"},
   {"py2c_syscon_load_polynomial", py2c_syscon_load_polynomial,
    METH_VARARGS, "gets the k-th polynomial from the systems container"},
   {"py2c_syscon_load_dobldobl_polynomial",
     py2c_syscon_load_dobldobl_polynomial,
    METH_VARARGS,
    "gets the k-th polynomial from the double double systems container"},
   {"py2c_syscon_load_quaddobl_polynomial",
     py2c_syscon_load_quaddobl_polynomial,
    METH_VARARGS,
    "gets the k-th polynomial from the quad double systems container"},
   {"py2c_syscon_load_multprec_polynomial",
     py2c_syscon_load_multprec_polynomial,
    METH_VARARGS,
    "gets the k-th polynomial from the multiprecision systems container"},
   {"py2c_syscon_store_Laurential",
     py2c_syscon_store_Laurential, METH_VARARGS,
    "defines the k-th polynomial in the standard Laurent systems container"},
   {"py2c_syscon_store_dobldobl_Laurential",
     py2c_syscon_store_dobldobl_Laurential, METH_VARARGS,
    "defines the k-th polynomial in the dobldobl Laurent systems container"},
   {"py2c_syscon_store_quaddobl_Laurential",
     py2c_syscon_store_quaddobl_Laurential, METH_VARARGS,
    "defines the k-th polynomial in the quaddobl Laurent systems container"},
   {"py2c_syscon_store_multprec_Laurential",
     py2c_syscon_store_multprec_Laurential, METH_VARARGS,
    "defines the k-th polynomial in the multprec Laurent systems container"},
   {"py2c_syscon_load_standard_Laurential",
     py2c_syscon_load_standard_Laurential, METH_VARARGS,
    "gets the k-th polynomial from the standard Laurent systems container"},
   {"py2c_syscon_load_dobldobl_Laurential",
     py2c_syscon_load_dobldobl_Laurential, METH_VARARGS,
    "gets the k-th polynomial from the quaddobl Laurent systems container"},
   {"py2c_syscon_load_quaddobl_Laurential",
     py2c_syscon_load_quaddobl_Laurential, METH_VARARGS,
    "gets the k-th polynomial from the quaddobl Laurent systems container"},
   {"py2c_syscon_load_multprec_Laurential",
     py2c_syscon_load_multprec_Laurential, METH_VARARGS,
    "gets the k-th polynomial from the multprec Laurent systems container"},
   {"py2c_syscon_total_degree", py2c_syscon_total_degree,
    METH_VARARGS,"returns the total degree of the system in the container"},
   {"py2c_syscon_standard_drop_variable_by_index",
     py2c_syscon_standard_drop_variable_by_index,
    METH_VARARGS,
    "drops a variable with given index from the standard double system"},
   {"py2c_syscon_standard_drop_variable_by_name",
     py2c_syscon_standard_drop_variable_by_name,
    METH_VARARGS,
    "drops a variable with given name from the standard double system"},
   {"py2c_syscon_dobldobl_drop_variable_by_index",
     py2c_syscon_dobldobl_drop_variable_by_index,
    METH_VARARGS,
    "drops a variable with given index from the double double system"},
   {"py2c_syscon_dobldobl_drop_variable_by_name",
     py2c_syscon_dobldobl_drop_variable_by_name,
    METH_VARARGS,
    "drops a variable with given name from the double double system"},
   {"py2c_syscon_quaddobl_drop_variable_by_index",
     py2c_syscon_quaddobl_drop_variable_by_index,
    METH_VARARGS,
    "drops a variable with given index from the quad double system"},
   {"py2c_syscon_quaddobl_drop_variable_by_name",
     py2c_syscon_quaddobl_drop_variable_by_name,
    METH_VARARGS,
    "drops a variable with given name from the quad double system"},
   {"py2c_solcon_length_solution_string", py2c_solcon_length_solution_string,
    METH_VARARGS, "returns #characters in the k-th solution in the container"},
   {"py2c_solcon_length_dobldobl_solution_string",
     py2c_solcon_length_dobldobl_solution_string,
    METH_VARARGS,
    "returns #characters in the k-th double double solution in the container"},
   {"py2c_solcon_length_quaddobl_solution_string",
     py2c_solcon_length_quaddobl_solution_string,
    METH_VARARGS,
    "returns #characters in the k-th quad double solution in the container"},
   {"py2c_solcon_length_multprec_solution_string",
     py2c_solcon_length_multprec_solution_string,
    METH_VARARGS,
    "returns #characters in the k-th multiprecision solution in the container"},
   {"py2c_solcon_write_solution_string", py2c_solcon_write_solution_string,
    METH_VARARGS, "writes the k-th solution to a string"},
   {"py2c_solcon_write_dobldobl_solution_string",
     py2c_solcon_write_dobldobl_solution_string,
    METH_VARARGS, "writes the k-th double double solution to a string"},
   {"py2c_solcon_write_quaddobl_solution_string",
     py2c_solcon_write_quaddobl_solution_string,
    METH_VARARGS, "writes the k-th quad double solution to a string"},
   {"py2c_solcon_write_multprec_solution_string",
     py2c_solcon_write_multprec_solution_string,
    METH_VARARGS, "writes the k-th multiprecision solution to a string"},
   {"py2c_solcon_append_solution_string", py2c_solcon_append_solution_string,
    METH_VARARGS, "appends a solution string to the container"},
   {"py2c_solcon_append_dobldobl_solution_string",
     py2c_solcon_append_dobldobl_solution_string,
    METH_VARARGS, "appends a solution string to the double double container"},
   {"py2c_solcon_append_quaddobl_solution_string",
     py2c_solcon_append_quaddobl_solution_string,
    METH_VARARGS, "appends a solution string to the quad double container"},
   {"py2c_solcon_append_multprec_solution_string",
     py2c_solcon_append_multprec_solution_string,
    METH_VARARGS, "appends a solution string to the multiprecision container"},
   {"py2c_solcon_number_of_solutions", py2c_solcon_number_of_solutions,
    METH_VARARGS, "returns the size of the solutions container"},
   {"py2c_solcon_number_of_dobldobl_solutions",
     py2c_solcon_number_of_dobldobl_solutions,
    METH_VARARGS, 
    "returns the size of the double double solutions container"},
   {"py2c_solcon_number_of_quaddobl_solutions",
     py2c_solcon_number_of_quaddobl_solutions,
    METH_VARARGS, 
    "returns the size of the quad double solutions container"},
   {"py2c_solcon_number_of_multprec_solutions",
     py2c_solcon_number_of_multprec_solutions,
    METH_VARARGS, 
    "returns the size of the multiprecision solutions container"},
   {"py2c_solcon_standard_drop_coordinate_by_index",
     py2c_solcon_standard_drop_coordinate_by_index,
    METH_VARARGS,
    "drops a coordinate with given index from the standard double solutions"},
   {"py2c_solcon_standard_drop_coordinate_by_name",
     py2c_solcon_standard_drop_coordinate_by_name,
    METH_VARARGS,
    "drops a coordinate with given name from the standard double solutions"},
   {"py2c_solcon_dobldobl_drop_coordinate_by_index",
     py2c_solcon_dobldobl_drop_coordinate_by_index,
    METH_VARARGS,
    "drops a coordinate with given index from the double double solutions"},
   {"py2c_solcon_dobldobl_drop_coordinate_by_name",
     py2c_solcon_dobldobl_drop_coordinate_by_name,
    METH_VARARGS,
    "drops a coordinate with given name from the double double solutions"},
   {"py2c_solcon_quaddobl_drop_coordinate_by_index",
     py2c_solcon_quaddobl_drop_coordinate_by_index,
    METH_VARARGS,
    "drops a coordinate with given index from the quad double solutions"},
   {"py2c_solcon_quaddobl_drop_coordinate_by_name",
     py2c_solcon_quaddobl_drop_coordinate_by_name,
    METH_VARARGS,
    "drops a coordinate with given name from the quad double solutions"},
   {"py2c_product_supporting_set_structure",
     py2c_product_supporting_set_structure,
    METH_VARARGS, "makes a supporting set structure for system in container"},
   {"py2c_product_write_set_structure", py2c_product_write_set_structure,
    METH_VARARGS, "writes the constructed supporting set structure"},
   {"py2c_product_set_structure_string", py2c_product_set_structure_string,
    METH_VARARGS, "return the set structure in a string"},
   {"py2c_product_parse_set_structure", py2c_product_parse_set_structure,
    METH_VARARGS, "parses given string into a set structure"},
   {"py2c_product_is_set_structure_supporting",
     py2c_product_is_set_structure_supporting,
    METH_VARARGS, "verifies if a set structure supports system in container"},
   {"py2c_product_linear_product_root_count",
     py2c_product_linear_product_root_count,
    METH_VARARGS, "returns root count based on the supporting set structure"},
   {"py2c_product_random_linear_product_system",
     py2c_product_random_linear_product_system,
    METH_VARARGS, "puts random linear-product system in the container"},
   {"py2c_product_solve_linear_product_system",
     py2c_product_solve_linear_product_system,
    METH_VARARGS, "puts solutions of linear-product system in the container"},
   {"py2c_product_clear_set_structure", py2c_product_clear_set_structure,
    METH_VARARGS, "clears the constructed supporting set structure"},
   {"py2c_product_m_homogeneous_Bezout_number",
     py2c_product_m_homogeneous_Bezout_number,
    METH_VARARGS, "returns an m-homogeneous Bezout number and a partition"},
   {"py2c_product_m_partition_Bezout_number",
     py2c_product_m_partition_Bezout_number,
    METH_VARARGS, "returns an m-homogeneous Bezout number for given partition"},
   {"py2c_product_m_homogeneous_start_system",
     py2c_product_m_homogeneous_start_system,
    METH_VARARGS, "makes an m-homogeneous start system for given partition"},
   {"py2c_celcon_set_type_of_mixture",
     py2c_celcon_set_type_of_mixture,METH_VARARGS,
    "initializes the cells container with a type of mixture"},
   {"py2c_celcon_type_of_mixture",
     py2c_celcon_type_of_mixture,METH_VARARGS,
    "returns the type of mixture stored in the cells container"},
   {"py2c_celcon_append_lifted_point",
     py2c_celcon_append_lifted_point,METH_VARARGS,
    "appends a lifted point to the cells container"},
   {"py2c_celcon_number_of_cells", py2c_celcon_number_of_cells,
    METH_VARARGS, "returns the number of cells in the cell container"},
   {"py2c_celcon_create_random_coefficient_system",
     py2c_celcon_create_random_coefficient_system,
    METH_VARARGS, "takes cell data to make a random coefficient system"},
   {"py2c_celcon_dobldobl_random_coefficient_system",
     py2c_celcon_dobldobl_random_coefficient_system,
    METH_VARARGS, "use cell data for a random dobldobl coefficient system"},
   {"py2c_celcon_quaddobl_random_coefficient_system",
     py2c_celcon_quaddobl_random_coefficient_system,
    METH_VARARGS, "use cell data for a random quaddobl coefficient system"},
   {"py2c_celcon_copy_into_systems_container",
     py2c_celcon_copy_into_systems_container, METH_VARARGS,
    "copy random coefficient system into the standard double system container"},
   {"py2c_celcon_copy_into_dobldobl_systems_container",
     py2c_celcon_copy_into_dobldobl_systems_container, METH_VARARGS,
    "copy random coefficient system into the double double system container"},
   {"py2c_celcon_copy_into_quaddobl_systems_container",
     py2c_celcon_copy_into_quaddobl_systems_container, METH_VARARGS,
    "copy random coefficient system into the quad double system container"},
   {"py2c_celcon_create_polyhedral_homotopy",
     py2c_celcon_create_polyhedral_homotopy, METH_VARARGS, 
    "polyhedral homotopy to solve a standard random coefficient system"},
   {"py2c_celcon_dobldobl_polyhedral_homotopy",
     py2c_celcon_dobldobl_polyhedral_homotopy, METH_VARARGS, 
    "polyhedral homotopy to solve a dobldobl random coefficient system"},
   {"py2c_celcon_quaddobl_polyhedral_homotopy",
     py2c_celcon_quaddobl_polyhedral_homotopy, METH_VARARGS, 
    "polyhedral homotopy to solve a quaddobl random coefficient system"},
   {"py2c_celcon_solve_start_system",
     py2c_celcon_solve_start_system, METH_VARARGS, 
    "solve start system for a given mixed cell number in standard precision"},
   {"py2c_celcon_solve_dobldobl_start_system",
     py2c_celcon_solve_dobldobl_start_system, METH_VARARGS, 
    "solve start system for a given mixed cell number with double doubles"},
   {"py2c_celcon_solve_quaddobl_start_system",
     py2c_celcon_solve_quaddobl_start_system, METH_VARARGS, 
    "solve start system for a given mixed cell number with quad doubles"},
   {"py2c_celcon_track_solution_path",
     py2c_celcon_track_solution_path, METH_VARARGS, 
    "tracks a solution path for one solution of a cell in standard precision"},
   {"py2c_celcon_track_dobldobl_solution_path",
     py2c_celcon_track_dobldobl_solution_path, METH_VARARGS, 
    "tracks a solution path for one solution of a cell with double doubles"},
   {"py2c_celcon_track_quaddobl_solution_path",
     py2c_celcon_track_quaddobl_solution_path, METH_VARARGS, 
    "tracks a solution path for one solution of a cell with quad doubles"},
   {"py2c_celcon_copy_target_solution_to_container",
     py2c_celcon_copy_target_solution_to_container, METH_VARARGS, 
    "copies one end solution of a cell to standard double solution container"},
   {"py2c_celcon_copy_target_dobldobl_solution_to_container",
     py2c_celcon_copy_target_dobldobl_solution_to_container, METH_VARARGS, 
    "copies one end solution of a cell to double double solution container"},
   {"py2c_celcon_copy_target_quaddobl_solution_to_container",
     py2c_celcon_copy_target_quaddobl_solution_to_container, METH_VARARGS, 
    "copies one end solution of a cell to quad double solution container"},
   {"py2c_celcon_permute_system", py2c_celcon_permute_system,
    METH_VARARGS, "permutes system in standard double systems container"},
   {"py2c_celcon_permute_dobldobl_system",
     py2c_celcon_permute_dobldobl_system,
    METH_VARARGS, "permutes system in double double systems container"},
   {"py2c_celcon_permute_quaddobl_system",
     py2c_celcon_permute_quaddobl_system,
    METH_VARARGS, "permutes system in quad double systems container"},
   {"py2c_celcon_clear_container", py2c_celcon_clear_container,
    METH_VARARGS, "clears the cell container"},
   {"py2c_embed_system", py2c_embed_system,
    METH_VARARGS, "replaces system in container with its embedding"},
   {"py2c_standard_cascade_homotopy", py2c_standard_cascade_homotopy,
    METH_VARARGS, "standard double cascade homotopy to go one dimension down"},
   {"py2c_dobldobl_cascade_homotopy", py2c_dobldobl_cascade_homotopy,
    METH_VARARGS, "double double cascade homotopy to go one dimension down"},
   {"py2c_quaddobl_cascade_homotopy", py2c_quaddobl_cascade_homotopy,
    METH_VARARGS, "quad double cascade homotopy to go one dimension down"},
   {"py2c_factor_set_to_mute", py2c_factor_set_to_mute,
    METH_VARARGS, "set the state to monodromy permutations to silent"},
   {"py2c_factor_define_output_file_with_string",
     py2c_factor_define_output_file_with_string,
    METH_VARARGS, "defines the output file with a given string"},
   {"py2c_factor_assign_labels", py2c_factor_assign_labels,
    METH_VARARGS, "assigns labels to the solutions in the container"},
   {"py2c_factor_initialize_sampler", py2c_factor_initialize_sampler,
    METH_VARARGS, "initialize sampler with the dimension of the set "},
   {"py2c_factor_initialize_monodromy", py2c_factor_initialize_monodromy,
    METH_VARARGS, "initializes state of monodromy permutations "},
   {"py2c_factor_store_solutions", py2c_factor_store_solutions,
    METH_VARARGS, "solutions go from container to monodromy permutations "},
   {"py2c_factor_restore_solutions", py2c_factor_restore_solutions,
    METH_VARARGS, "solutions go from monodromy permutations to container"},
   {"py2c_factor_track_paths", py2c_factor_track_paths,
    METH_VARARGS, "tracks as many paths as defined by the witness set"},
   {"py2c_factor_swap_slices", py2c_factor_swap_slices,
    METH_VARARGS, "swaps the current slides with new slices"},
   {"py2c_factor_new_slices", py2c_factor_new_slices,
    METH_VARARGS, "generates k new slices in n-space"},
   {"py2c_factor_set_trace_slice", py2c_factor_set_trace_slice,
    METH_VARARGS, "sets the coefficient of the slice for the linear trace"},
   {"py2c_factor_store_gammas", py2c_factor_store_gammas,
    METH_VARARGS, "stores n randomly generated complex coefficients"},
   {"py2c_factor_permutation_after_loop", py2c_factor_permutation_after_loop,
    METH_VARARGS, "returns permutation computed after one loop"},
   {"py2c_factor_update_decomposition", py2c_factor_update_decomposition,
    METH_VARARGS, "updates the decomposition with a permutation"},
   {"py2c_factor_number_of_components", py2c_factor_number_of_components,
    METH_VARARGS, "returns the current number of components"},
   {"py2c_factor_witness_points_of_component",
     py2c_factor_witness_points_of_component,
    METH_VARARGS, "returns the labels of witness points in a component"},
   {"py2c_factor_trace_sum_difference", py2c_factor_trace_sum_difference,
    METH_VARARGS, "returns the difference of predicted and actual trace sum"},
   {"py2c_witness_set_of_hypersurface", py2c_witness_set_of_hypersurface,
    METH_VARARGS, "makes a witness set of a hypersurface"},
   {"py2c_create_diagonal_homotopy", py2c_create_diagonal_homotopy,
    METH_VARARGS, "makes diagonal homotopy to intersect two witness sets"},
   {"py2c_start_diagonal_cascade_solutions",
     py2c_start_diagonal_cascade_solutions,
    METH_VARARGS, "makes solutions to start an extrinsic cascade homotopy"},
   {"py2c_extrinsic_top_diagonal_dimension",
     py2c_extrinsic_top_diagonal_dimension,
    METH_VARARGS, "returns dimension at start of extrinsic cascade"},
   {"py2c_collapse_diagonal", py2c_collapse_diagonal,
    METH_VARARGS, "eliminates extrinsic diagonal from system and solutions"},
   {"py2c_schubert_pieri_count", py2c_schubert_pieri_count,
    METH_VARARGS, "the combinatorial Pieri root count"},
   {"py2c_schubert_resolve_conditions", py2c_schubert_resolve_conditions,
    METH_VARARGS, "resolve general Schubert intersection conditions"},
   {"py2c_schubert_littlewood_richardson_homotopies",
     py2c_schubert_littlewood_richardson_homotopies,
    METH_VARARGS, "Littlewood-Richardson homotopies for Schubert problems"},
   {"py2c_schubert_localization_poset", py2c_schubert_localization_poset,
    METH_VARARGS, "localization poset for the Pieri root count"},
   {"py2c_schubert_pieri_homotopies", py2c_schubert_pieri_homotopies,
    METH_VARARGS, "runs the Pieri homotopies on random input data"},
   {"py2c_schubert_osculating_planes", py2c_schubert_osculating_planes,
    METH_VARARGS, "returns real planes osculating a rational normal curve"},
   {"py2c_schubert_pieri_system", py2c_schubert_pieri_system,
    METH_VARARGS, "makes the polynomial system defined by a Pieri problem"},
   {"py2c_mapcon_solve_system", py2c_mapcon_solve_system,
    METH_VARARGS, "computes monomial maps solving a binomial system"},
   {"py2c_mapcon_write_maps", py2c_mapcon_write_maps,
    METH_VARARGS, "writes monomial maps stored in container"},
   {"py2c_mapcon_clear_maps", py2c_mapcon_clear_maps,
    METH_VARARGS, "destroys all monomial maps stored in container"},
   {"py2c_mapcon_top_dimension", py2c_mapcon_top_dimension,
    METH_VARARGS, "returns top dimension of the monomial maps"},
   {"py2c_mapcon_number_of_maps", py2c_mapcon_number_of_maps,
    METH_VARARGS, "returns the number of maps for a given dimension"},
   {"py2c_mapcon_degree_of_map", py2c_mapcon_degree_of_map,
    METH_VARARGS, "returns the degree of map for a given dimension and index"},
   {"py2c_mapcon_coefficients_of_map", py2c_mapcon_coefficients_of_map,
    METH_VARARGS, "returns the coefficients of map for a dimension and index"},
   {"py2c_mapcon_exponents_of_map", py2c_mapcon_exponents_of_map,
    METH_VARARGS, "returns the exponents of map for a dimension and index"},
   {"py2c_initialize_standard_homotopy", py2c_initialize_standard_homotopy,
    METH_VARARGS, "initializes homotopy for tracking with standard doubles"},
   {"py2c_initialize_dobldobl_homotopy", py2c_initialize_dobldobl_homotopy,
    METH_VARARGS, "initializes homotopy for tracking with double doubles"},
   {"py2c_initialize_quaddobl_homotopy", py2c_initialize_quaddobl_homotopy,
    METH_VARARGS, "initializes homotopy for tracking with quad doubles"},
   {"py2c_initialize_multprec_homotopy", py2c_initialize_multprec_homotopy,
    METH_VARARGS, "initializes homotopy for tracking in multiprecision"},
   {"py2c_initialize_standard_solution", py2c_initialize_standard_solution, 
    METH_VARARGS, "initializes standard double start solution for tracking"},
   {"py2c_initialize_dobldobl_solution", py2c_initialize_dobldobl_solution,
    METH_VARARGS, "initializes double double start solution for tracking"},
   {"py2c_initialize_quaddobl_solution", py2c_initialize_quaddobl_solution,
    METH_VARARGS, "initializes quad double start solution for tracking"},
   {"py2c_initialize_multprec_solution", py2c_initialize_multprec_solution,
    METH_VARARGS, "initializes multiprecision start solution for tracking"},
   {"py2c_next_standard_solution", py2c_next_standard_solution, 
    METH_VARARGS, "one predictor-corrector step with standard doubles"},
   {"py2c_next_dobldobl_solution", py2c_next_dobldobl_solution,
    METH_VARARGS, "one predictor-corrector step with double doubles"},
   {"py2c_next_quaddobl_solution", py2c_next_quaddobl_solution,
    METH_VARARGS, "one predictor-corrector step with quad doubles"},
   {"py2c_next_multprec_solution", py2c_next_multprec_solution,
    METH_VARARGS, "one predictor-corrector step in multiprecision"},
   {"py2c_clear_standard_tracker", py2c_clear_standard_tracker,
    METH_VARARGS, "deallocates and resets tracker in standard doubles"},
   {"py2c_clear_dobldobl_tracker", py2c_clear_dobldobl_tracker,
    METH_VARARGS, "deallocates and resets tracker in double doubles"},
   {"py2c_clear_quaddobl_tracker", py2c_clear_quaddobl_tracker,
    METH_VARARGS, "deallocates and resets tracker in quad doubles"},
   {"py2c_clear_multprec_tracker", py2c_clear_multprec_tracker,
    METH_VARARGS, "deallocates and resets tracker in multiprecision"},
   {NULL, NULL, 0, NULL} 
};

/* This is the initialization routine which will be called by the 
 * Python run-time when the library is imported in order to retrieve 
 * a pointer to the above method address table.
 * Note that therefore this routine must be visible in the dynamic library
 * either through the use of a ".def" file or by a compiler instruction 
 * such as "declspec(export)" */

PyMODINIT_FUNC initphcpy2c(void)
{
   Py_InitModule("phcpy2c", phcpy2c_methods);
}
