/* This file contains the definitions of the prototypes in phcpy2c.h. */

#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "phcpack.h"
#include "unisolvers.h"
#include "giftwrappers.h"
#include "schubert.h"
#include "syscon.h"
#include "tabform.h"
#include "syspool.h"
#include "solcon.h"
#include "product.h"
#include "lists_and_strings.h"
#include "celcon.h"
#include "intcelcon.h"
#include "scalers.h"
#include "reducers.h"
#include "numbtrop.h"
#include "sweep.h"
#include "multiplicity.h"
#include "witset.h"
#include "witsols.h"
#include "mapcon.h"
#include "series.h"
#include "padcon.h"
#include "jump_track.h"
#include "next_track.h"
#include "structmember.h"

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern void adafinal( void );
#endif

int g_ada_initialized = 0;
int g_ada_finalized = 0;

void initialize ( void )
{
   if(!g_ada_initialized)
   {
      adainit();
      g_ada_initialized = 1;
      g_ada_finalized = 0;
   }
}

void finalize ( void )
{
   if(!g_ada_finalized)
   {
      adafinal();
      g_ada_finalized = 1;
      g_ada_initialized = 0;
   }
}

static PyObject *py2c_corecount
 ( PyObject *self, PyObject *args )
{
   if(!PyArg_ParseTuple(args,"")) return NULL;   

   int number_of_cores = sysconf(_SC_NPROCESSORS_ONLN);

   return Py_BuildValue("i", number_of_cores);
}

/* The wrapping of functions in phcpack.h starts from here. */

static PyObject *py2c_PHCpack_version_string
 ( PyObject *self, PyObject *args )
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

static PyObject *py2c_set_seed
 ( PyObject *self, PyObject *args )
{
   int fail,seed;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&seed)) return NULL;

   fail = set_seed(seed);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_get_seed
 ( PyObject *self, PyObject *args )
{
   int fail,seed;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;

   fail = get_seed(&seed);
              
   return Py_BuildValue("i",seed);
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

static PyObject *py2c_define_output_file
 ( PyObject *self, PyObject *args )
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

static PyObject *py2c_write_dobldobl_target_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = write_dobldobl_target_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_write_quaddobl_target_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = write_quaddobl_target_system();
   
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

static PyObject *py2c_write_dobldobl_start_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = write_dobldobl_start_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_write_quaddobl_start_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = write_quaddobl_start_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_standard_start_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = read_standard_start_Laurent_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_write_standard_start_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = write_standard_start_Laurent_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_standard_target_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = read_standard_target_Laurent_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_write_standard_target_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = write_standard_target_Laurent_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_dobldobl_start_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = read_dobldobl_start_Laurent_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_write_dobldobl_start_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = write_dobldobl_start_Laurent_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_dobldobl_target_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = read_dobldobl_target_Laurent_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_write_dobldobl_target_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = write_dobldobl_target_Laurent_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_quaddobl_start_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = read_quaddobl_start_Laurent_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_write_quaddobl_start_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = write_quaddobl_start_Laurent_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_read_quaddobl_target_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = read_quaddobl_target_Laurent_system();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_write_quaddobl_target_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = write_quaddobl_target_Laurent_system();
   
   return Py_BuildValue("i",fail);
}

/* The wrapping of copying systems from and to containers starts here. */

static PyObject *py2c_copy_standard_target_system_to_container
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

static PyObject *py2c_copy_standard_container_to_target_system
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

static PyObject *py2c_copy_standard_container_to_start_system 
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

static PyObject *py2c_copy_standard_Laurent_container_to_start_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_standard_Laurent_container_to_start_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_dobldobl_Laurent_container_to_start_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_dobldobl_Laurent_container_to_start_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_quaddobl_Laurent_container_to_start_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_quaddobl_Laurent_container_to_start_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_standard_Laurent_container_to_target_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_standard_Laurent_container_to_target_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_dobldobl_Laurent_container_to_target_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_dobldobl_Laurent_container_to_target_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_quaddobl_Laurent_container_to_target_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_quaddobl_Laurent_container_to_target_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_standard_Laurent_start_system_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_standard_Laurent_start_system_to_container();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_dobldobl_Laurent_start_system_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_dobldobl_Laurent_start_system_to_container();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_quaddobl_Laurent_start_system_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_quaddobl_Laurent_start_system_to_container();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_standard_Laurent_target_system_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_standard_Laurent_target_system_to_container();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_dobldobl_Laurent_target_system_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_dobldobl_Laurent_target_system_to_container();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_quaddobl_Laurent_target_system_to_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = copy_quaddobl_Laurent_target_system_to_container();
              
   return Py_BuildValue("i",fail);
}

/* creation of homotopy and tracking all solution paths */

static PyObject *py2c_create_standard_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = create_homotopy();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_create_standard_homotopy_with_gamma
 ( PyObject *self, PyObject *args )
{
   int fail,pwt;
   double g_re,g_im;

   initialize();
   if(!PyArg_ParseTuple(args,"ddi",&g_re,&g_im,&pwt)) return NULL;   
   fail = create_homotopy_with_given_gamma(g_re,g_im,pwt);
              
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
   int fail,pwt;
   double g_re,g_im;

   initialize();
   if(!PyArg_ParseTuple(args,"ddi",&g_re,&g_im,&pwt)) return NULL;   
   fail = create_dobldobl_homotopy_with_given_gamma(g_re,g_im,pwt);
              
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
   int fail,pwt;
   double g_re,g_im;

   initialize();
   if(!PyArg_ParseTuple(args,"ddi",&g_re,&g_im,&pwt)) return NULL;   
   fail = create_quaddobl_homotopy_with_given_gamma(g_re,g_im,pwt);
              
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
   int fail,pwt;
   double g_re,g_im;

   initialize();
   if(!PyArg_ParseTuple(args,"dd",&g_re,&g_im,&pwt)) return NULL;   
   fail = create_multprec_homotopy_with_given_gamma(g_re,g_im,pwt);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_standard_homotopy
 ( PyObject *self, PyObject *args )
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

static PyObject *py2c_write_start_solutions
 ( PyObject *self, PyObject *args )
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

static PyObject *py2c_get_value_of_continuation_parameter
 ( PyObject *self, PyObject *args )
{
   int fail,idx;
   double val;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&idx)) return NULL;   
   fail = get_value_of_continuation_parameter(idx,&val);
   
   return Py_BuildValue("d",val);
}

static PyObject *py2c_set_value_of_continuation_parameter
 ( PyObject *self, PyObject *args )
{
   int fail,idx;
   double val;

   initialize();
   if(!PyArg_ParseTuple(args,"id",&idx,&val)) return NULL;   
   fail = set_value_of_continuation_parameter(idx,&val);
   
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

static PyObject *py2c_solve_by_standard_Laurent_homotopy_continuation
 ( PyObject *self, PyObject *args )
{
   int fail, nbtasks = 0;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&nbtasks)) return NULL;   
   fail = solve_by_standard_Laurent_homotopy_continuation(nbtasks);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solve_by_dobldobl_Laurent_homotopy_continuation
 ( PyObject *self, PyObject *args )
{
   int fail, nbtasks = 0;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&nbtasks)) return NULL;   
   fail = solve_by_dobldobl_Laurent_homotopy_continuation(nbtasks);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solve_by_quaddobl_Laurent_homotopy_continuation
 ( PyObject *self, PyObject *args )
{
   int fail, nbtasks = 0;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&nbtasks)) return NULL;   
   fail = solve_by_quaddobl_Laurent_homotopy_continuation(nbtasks);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_standard_operations_data
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = clear_data();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_dobldobl_operations_data
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = clear_dobldobl_data();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_quaddobl_operations_data
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = clear_quaddobl_data();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_standard_Laurent_data
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = clear_standard_Laurent_data();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_dobldobl_Laurent_data
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = clear_dobldobl_Laurent_data();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_quaddobl_Laurent_data
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = clear_quaddobl_Laurent_data();

   return Py_BuildValue("i",fail);
}

/* Wrapping of crude path trackers of the jumpstart library starts here. */

static PyObject *py2c_standard_crude_tracker
 ( PyObject *self, PyObject *args )
{
   int fail,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&verbose)) return NULL;   
   fail = standard_crude_tracker(verbose);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_dobldobl_crude_tracker
 ( PyObject *self, PyObject *args )
{
   int fail,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&verbose)) return NULL;   
   fail = dobldobl_crude_tracker(verbose);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_quaddobl_crude_tracker
 ( PyObject *self, PyObject *args )
{
   int fail,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&verbose)) return NULL;   
   fail = quaddobl_crude_tracker(verbose);

   return Py_BuildValue("i",fail);
}

/* The wrapping of copying solutions from and to containers starts here. */

static PyObject *py2c_copy_standard_target_solutions_to_container
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

static PyObject *py2c_copy_standard_container_to_target_solutions
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

static PyObject *py2c_copy_standard_container_to_start_solutions
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

static PyObject *py2c_solve_standard_system
 ( PyObject *self, PyObject *args )
{
   int fail,rc,nrc,silent,vrb,nbtasks,mvfocus = 0;
   char rocos[1024];

   initialize();
   if(!PyArg_ParseTuple(args,"iiii",&silent,&nbtasks,&mvfocus,&vrb))
      return NULL;
   fail = solve_standard_system(&rc,silent,&nrc,rocos,nbtasks,mvfocus,vrb);
   if(silent == 1)
      return Py_BuildValue("i",rc);
   else
      return Py_BuildValue("(i,s)",rc,rocos);
}

static PyObject *py2c_scan_for_symbols
 ( PyObject *self, PyObject *args )
{
   int fail,nbc,dim;
   char *polsys;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nbc,&polsys)) return NULL;

   fail = scan_number_of_variables(nbc,polsys,&dim);

   return Py_BuildValue("i",dim);
}

static PyObject *py2c_solve_dobldobl_system
 ( PyObject *self, PyObject *args )
{
   int fail,rc,nrc,silent,nbtasks,vrb = 0;
   char rocos[1024];

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&silent,&nbtasks,&vrb)) return NULL;
   fail = solve_dobldobl_system(&rc,silent,&nrc,rocos,nbtasks,vrb);
   if(silent == 1)
      return Py_BuildValue("i",rc);
   else
      return Py_BuildValue("(i,s)",rc,rocos);
}

static PyObject *py2c_solve_quaddobl_system
 ( PyObject *self, PyObject *args )
{
   int fail,rc,nrc,silent,vrb,nbtasks = 0;
   char rocos[1024];

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&silent,&nbtasks,&vrb)) return NULL;
   fail = solve_quaddobl_system(&rc,silent,&nrc,rocos,nbtasks,vrb);
   if(silent == 1)
      return Py_BuildValue("i",rc);
   else
      return Py_BuildValue("(i,s)",rc,rocos);
}

static PyObject *py2c_solve_standard_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int silent,fail,rc,nrc,vrb,nbtasks,mvfocus = 0;
   char rocos[1024];

   initialize();
   if(!PyArg_ParseTuple(args,"iiii",&silent,&nbtasks,&mvfocus,&vrb))
      return NULL;
   fail = solve_standard_Laurent_system
            (&rc,silent,&nrc,rocos,nbtasks,mvfocus,vrb);
   if(silent == 1)
      return Py_BuildValue("i",rc);
   else
      return Py_BuildValue("(i,s)",rc,rocos);
}

static PyObject *py2c_solve_dobldobl_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int silent,fail,rc,nrc,vrb,nbtasks = 0;
   char rocos[1024];

   initialize();
   if (!PyArg_ParseTuple(args,"iii",&silent,&nbtasks,&vrb)) return NULL;
   fail = solve_dobldobl_Laurent_system(&rc,silent,&nrc,rocos,nbtasks,vrb);
   if(silent == 1)
      return Py_BuildValue("i",rc);
   else
      return Py_BuildValue("(i,s)",rc,rocos);
}

static PyObject *py2c_solve_quaddobl_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int silent,fail,rc,nrc,nbtasks,vrb = 0;
   char rocos[1024];

   initialize();
   if (!PyArg_ParseTuple(args,"iii",&silent,&nbtasks,&vrb)) return NULL;
   fail = solve_quaddobl_Laurent_system(&rc,silent,&nrc,rocos,nbtasks,vrb);
   if(silent == 1)
      return Py_BuildValue("i",rc);
   else
      return Py_BuildValue("(i,s)",rc,rocos);
}

static PyObject *py2c_set_gamma_constant ( PyObject *self, PyObject *args )
{
   int fail,prc,vrb;
   double regm,imgm;

   initialize();
   if (!PyArg_ParseTuple(args,"ddii",&regm,&imgm,&prc,&vrb)) return NULL;
   fail = set_gamma_constant(regm,imgm,prc,vrb);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_get_gamma_constant ( PyObject *self, PyObject *args )
{
   int fail,prc,vrb;
   double regm,imgm;

   initialize();
   if (!PyArg_ParseTuple(args,"ii",&prc,&vrb)) return NULL;
   fail = get_gamma_constant(&regm,&imgm,prc,vrb);

   return Py_BuildValue("(d,d)",regm,imgm);
}

static PyObject *py2c_mixed_volume
 ( PyObject *self, PyObject *args )
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

static PyObject *py2c_mixed_volume_by_demics 
 ( PyObject *self, PyObject *args )
{
   int stable,fail,mv,smv;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&stable)) return NULL;
   if(stable == 0)
   {
      fail = mixed_volume_by_demics(&mv);
      return Py_BuildValue("i",mv);
   }
   else
   {
      fail = stable_mixed_volume_by_demics(&mv,&smv);
      return Py_BuildValue("(i,i)",mv,smv);
   }
}

static PyObject *py2c_standard_deflate
 ( PyObject *self, PyObject *args )
{
   int mxit,mxdf;
   double terr,tres,trnk;

   initialize();
   if(!PyArg_ParseTuple(args,"iiddd",&mxit,&mxdf,&terr,&tres,&trnk))
      return NULL;
   {
      int fail = standard_deflate(mxit,mxdf,terr,tres,trnk);
      return Py_BuildValue("i",fail);
   }
}

static PyObject *py2c_dobldobl_deflate
 ( PyObject *self, PyObject *args )
{
   int mxit,mxdf;
   double terr,tres,trnk;

   initialize();
   if(!PyArg_ParseTuple(args,"iiddd",&mxit,&mxdf,&terr,&tres,&trnk))
      return NULL;
   {
      int fail = dobldobl_deflate(mxit,mxdf,terr,tres,trnk);
      return Py_BuildValue("i",fail);
   }
}

static PyObject *py2c_quaddobl_deflate
 ( PyObject *self, PyObject *args )
{
   int mxit,mxdf;
   double terr,tres,trnk;

   initialize();
   if(!PyArg_ParseTuple(args,"iiddd",&mxit,&mxdf,&terr,&tres,&trnk))
      return NULL;
   {
      int fail = quaddobl_deflate(mxit,mxdf,terr,tres,trnk);
      return Py_BuildValue("i",fail);
   }
}

static PyObject *py2c_standard_Newton_step
 ( PyObject *self, PyObject *args )
{
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   {
      int fail = standard_Newton_step();
      return Py_BuildValue("i",fail);
   }
}

static PyObject *py2c_dobldobl_Newton_step
 ( PyObject *self, PyObject *args )
{
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   {
      int fail = dobldobl_Newton_step();
      return Py_BuildValue("i",fail);
   }
}

static PyObject *py2c_quaddobl_Newton_step
 ( PyObject *self, PyObject *args )
{
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   {
      int fail = quaddobl_Newton_step();
      return Py_BuildValue("i",fail);
   }
}

static PyObject *py2c_multprec_Newton_step
 ( PyObject *self, PyObject *args )
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

static PyObject *py2c_varbprec_Newton_Laurent_steps
 ( PyObject *self, PyObject *args )
{
   int fail,dim,wanted,maxitr,maxprc,nbstr;
   char *pols;

   initialize();

   if(!PyArg_ParseTuple(args,"iiiiis",
      &dim,&wanted,&maxitr,&maxprc,&nbstr,&pols)) return NULL;

   fail = varbprec_Newton_Laurent_step(dim,wanted,maxitr,maxprc,nbstr,pols);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_standard_condition_report
 ( PyObject *self, PyObject *args )
{
   int fail,maxit,verbose;
   double tolres,tolerr,tolsing;
   char *name;
   int cntfail,cntreal,cntcmplx,cntregu,cntsing,cntclus,nbc;
   int t_err[16];
   int t_rco[16];
   int t_res[16];
   char st_err[256];
   char st_rco[256];
   char st_res[256];

   initialize();
   if(!PyArg_ParseTuple(args,"idddsi",
      &maxit,&tolres,&tolerr,&tolsing,&name,&verbose))
      return NULL;

   nbc = strlen(name);

   if(verbose == 1)
   {
      if(nbc == 0)
         printf("Writing the output to screen.\n");
      else
         printf("Writing the output to file %s.\n", name);
   }
   fail = standard_condition_report
            (maxit,tolres,tolerr,tolsing,nbc,name,
             &cntfail,&cntreal,&cntcmplx,&cntregu,&cntsing,&cntclus,
             t_err,t_rco,t_res,verbose);

   intlist2str(16, t_err, st_err);
   intlist2str(16, t_rco, st_rco);
   intlist2str(16, t_res, st_res);

   return Py_BuildValue("(i, i, i, i, i, i, i, s, s, s)", fail,
                        cntregu, cntsing, cntreal, cntcmplx, cntclus, cntfail,
                        st_err, st_res, st_rco);
}

/* The wrapping of the functions in unisolvers.h starts from here */

static PyObject *py2c_usolve_standard
 ( PyObject *self, PyObject *args )
{
   int fail,max,nit;
   double eps;

   initialize();
   if(!PyArg_ParseTuple(args,"id",&max,&eps)) return NULL;
   fail = solve_with_standard_doubles(max,eps,&nit);

   return Py_BuildValue("i",nit);
}

static PyObject *py2c_usolve_dobldobl
 ( PyObject *self, PyObject *args )
{
   int fail,max,nit;
   double eps;

   initialize();
   if(!PyArg_ParseTuple(args,"id",&max,&eps)) return NULL;
   fail = solve_with_double_doubles(max,eps,&nit);

   return Py_BuildValue("i",nit);
}

static PyObject *py2c_usolve_quaddobl
 ( PyObject *self, PyObject *args )
{
   int fail,max,nit;
   double eps;

   initialize();
   if(!PyArg_ParseTuple(args,"id",&max,&eps)) return NULL;
   fail = solve_with_quad_doubles(max,eps,&nit);

   return Py_BuildValue("i",nit);
}

static PyObject *py2c_usolve_multprec
 ( PyObject *self, PyObject *args )
{
   int fail,dcp,max,nit;
   double eps;

   initialize();
   if(!PyArg_ParseTuple(args,"iid",&dcp,&max,&eps)) return NULL;
   fail = solve_with_multiprecision(dcp,max,eps,&nit);

   return Py_BuildValue("i",nit);
}

/* The wrapping of the functions in giftwrappers.h starts from here. */

static PyObject *py2c_giftwrap_planar
 ( PyObject *self, PyObject *args )
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

static PyObject *py2c_giftwrap_convex_hull
 ( PyObject *self, PyObject *args )
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

static PyObject *py2c_giftwrap_support_size
 ( PyObject *self, PyObject *args )
{
   int idx,nbr;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&idx)) return NULL;   
   nbr = support_size(idx);  

   return Py_BuildValue("i",nbr);
}

static PyObject *py2c_giftwrap_support_string
 ( PyObject *self, PyObject *args )
{
   int nbr,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&nbr)) return NULL;   
   {
      char support[nbr+1];

      fail = support_string(nbr,support);

      return Py_BuildValue("s",support);
   }
}

static PyObject *py2c_giftwrap_clear_support_string
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = clear_support_string();
   
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_giftwrap_initial_form
 ( PyObject *self, PyObject *args )
{
   int fail,dim,nbc;
   char *normal;

   initialize();
   if(!PyArg_ParseTuple(args,"iis",&dim,&nbc,&normal)) return NULL;
   /* printf("the inner normal in wrapper is %s\n", normal); */
   fail = initial_form(dim,nbc,normal);
   
   return Py_BuildValue("i",fail);
}

/* The wrapping of the functions in syscon.h starts from here. */

static PyObject *py2c_syscon_read_standard_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_read_standard_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_read_standard_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_read_standard_Laurent_system();
              
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

static PyObject *py2c_syscon_random_system
 ( PyObject *self, PyObject *args )
{
   int fail,n,m,d,c,neq;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiii",&n,&m,&d,&c,&neq)) return NULL;   
   fail = syscon_random_system(n,m,d,c,neq);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_dobldobl_random_system
 ( PyObject *self, PyObject *args )
{
   int fail,n,m,d,c,neq;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiii",&n,&m,&d,&c,&neq)) return NULL;   
   fail = syscon_dobldobl_random_system(n,m,d,c,neq);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_quaddobl_random_system
 ( PyObject *self, PyObject *args )
{
   int fail,n,m,d,c,neq;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiii",&n,&m,&d,&c,&neq)) return NULL;   
   fail = syscon_quaddobl_random_system(n,m,d,c,neq);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_write_standard_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_write_standard_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_write_standard_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_write_standard_Laurent_system();
              
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

static PyObject *py2c_syscon_clear_standard_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_clear_standard_system();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_clear_standard_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syscon_clear_standard_Laurent_system();
              
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

static PyObject *py2c_syscon_number_of_standard_polynomials
 ( PyObject *self, PyObject *args )
{
   int fail,number;

   initialize();
   if (!PyArg_ParseTuple(args,"")) return NULL;
   fail = syscon_number_of_standard_polynomials(&number);

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

static PyObject *py2c_syscon_number_of_standard_Laurentials
 ( PyObject *self, PyObject *args )
{
   int fail,number;
   static PyObject *a;

   initialize();
   /* if (!PyArg_ParseTuple(args,"i",&number)) return NULL; */
   if (!PyArg_ParseTuple(args,"")) return NULL;
   fail = syscon_number_of_standard_Laurentials(&number);

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

static PyObject *py2c_syscon_initialize_number_of_standard_polynomials
 ( PyObject *self, PyObject *args )
{      
   int fail,dim;

   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;
   initialize();
   fail = syscon_initialize_number_of_standard_polynomials(dim);
                 
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

static PyObject *py2c_syscon_initialize_number_of_standard_Laurentials
 ( PyObject *self, PyObject *args )
{      
   int fail,dim;

   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;
   initialize();
   fail = syscon_initialize_number_of_standard_Laurentials(dim);
                 
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

static PyObject *py2c_syscon_degree_of_standard_polynomial
 ( PyObject *self, PyObject *args )
{
   int fail,equ,deg;

   if(!PyArg_ParseTuple(args,"i",&equ)) return NULL;
   initialize();
   fail = syscon_degree_of_standard_polynomial(equ,&deg);
                 
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

static PyObject *py2c_syscon_number_of_terms
 ( PyObject *self, PyObject *args )
{
   int fail,i,number;

   initialize();
   /* if (!PyArg_ParseTuple(args,"ii",&i,&number)) return NULL; */
   if (!PyArg_ParseTuple(args,"i",&i)) return NULL;
   fail = syscon_number_of_standard_terms(i,&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_syscon_number_of_Laurent_terms
 ( PyObject *self, PyObject *args )
{
   int fail,i,number;

   initialize();
   /* if (!PyArg_ParseTuple(args,"ii",&i,&number)) return NULL; */
   if (!PyArg_ParseTuple(args,"i",&i)) return NULL;
   fail = syscon_number_of_standard_Laurent_terms(i,&number);

   return Py_BuildValue("i",number);
}

static PyObject *py2c_syscon_retrieve_term
 ( PyObject *self, PyObject *args )
{
   int fail, *exp;
   int i,j,n,k;
   double *c;
   static PyObject *a;

   initialize();
   if(!PyArg_ParseTuple(args, "iiiid", &i,&j,&n,&exp,&c)) return NULL;

   exp = (int *)malloc(n * sizeof(int));
   c   = (double *)malloc(2*sizeof(double));

   fail = syscon_retrieve_standard_term(i,j,n,exp,c);
     
   a = Py_BuildValue("i", fail);	 
   
   free(exp);
   free(c);

   return a;           
}

static PyObject *py2c_syscon_store_standard_polynomial
 ( PyObject *self, PyObject *args )
{      
   int fail,n,nc,k;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iiis",&nc,&n,&k,&p)) return NULL;
   fail = syscon_store_standard_polynomial(nc,n,k,p);

   if(fail != 0) printf("Failed to store %s.\n",p);

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
           
   if(fail != 0) printf("Failed to store %s.\n",p);

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

   if(fail != 0) printf("Failed to store %s.\n",p);
                 
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

   if(fail != 0) printf("Failed to store %s.\n",p);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_load_standard_polynomial
 ( PyObject *self, PyObject *args )
{      
   int fail,nc,k,szl;
   
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_standard_size_limit(k,&szl);
   {
      char p[szl];

      fail = syscon_load_standard_polynomial(k,&nc,p);
                 
      return Py_BuildValue("s",p);
   }
}

static PyObject *py2c_syscon_load_dobldobl_polynomial
 ( PyObject *self, PyObject *args )
{      
   int fail,nc,k,szl;
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_dobldobl_size_limit(k,&szl);
   {
      char p[szl];

      fail = syscon_load_dobldobl_polynomial(k,&nc,p);
                 
      return Py_BuildValue("s",p);
   }
}

static PyObject *py2c_syscon_load_quaddobl_polynomial
 ( PyObject *self, PyObject *args )
{      
   int fail,nc,k,szl;
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_quaddobl_size_limit(k,&szl);
   {
      char p[szl];

      fail = syscon_load_quaddobl_polynomial(k,&nc,p);
                 
      return Py_BuildValue("s",p);
   }
}

static PyObject *py2c_syscon_load_multprec_polynomial
 ( PyObject *self, PyObject *args )
{      
   int fail,nc,k,szl;
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_multprec_size_limit(k,&szl);
   {
      char p[szl];

      fail = syscon_load_multprec_polynomial(k,&nc,p);
                 
      return Py_BuildValue("s",p);
   }
}

static PyObject *py2c_syscon_store_standard_Laurential
 ( PyObject *self, PyObject *args )
{      
   int fail,n,nc,k;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iiis",&nc,&n,&k,&p)) return NULL;
   fail = syscon_store_standard_Laurential(nc,n,k,p);
                 
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
   int fail,nc,k,szl;
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_standard_Laurent_size_limit(k,&szl);
   {
      char p[szl];

      fail = syscon_load_standard_Laurential(k,&nc,p);
                 
      return Py_BuildValue("s",p);
   }
}

static PyObject *py2c_syscon_load_dobldobl_Laurential
 ( PyObject *self, PyObject *args )
{      
   int fail,nc,k,szl;
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_dobldobl_Laurent_size_limit(k,&szl);
   {
      char p[szl];

      fail = syscon_load_dobldobl_Laurential(k,&nc,p);
                 
      return Py_BuildValue("s",p);
   }
}

static PyObject *py2c_syscon_load_quaddobl_Laurential
 ( PyObject *self, PyObject *args )
{      
   int fail,nc,k,szl;
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_quaddobl_Laurent_size_limit(k,&szl);
   {
      char p[szl];

      fail = syscon_load_quaddobl_Laurential(k,&nc,p);

      return Py_BuildValue("s",p);
   }
}

static PyObject *py2c_syscon_load_multprec_Laurential
 ( PyObject *self, PyObject *args )
{      
   int fail,nc,k,szl;
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_multprec_Laurent_size_limit(k,&szl);
   {
      char p[szl];

      fail = syscon_load_multprec_Laurential(k,&nc,p);
                 
      return Py_BuildValue("s",p);
   }
}

static PyObject *py2c_syscon_total_degree
 ( PyObject *self, PyObject *args )
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

static PyObject *py2c_syscon_standard_Laurent_drop_variable_by_index
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_standard_Laurent_drop_variable_by_index(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_standard_Laurent_drop_variable_by_name
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *s;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&s)) return NULL;
   fail = syscon_standard_Laurent_drop_variable_by_name(nc,s);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_dobldobl_Laurent_drop_variable_by_index
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_dobldobl_Laurent_drop_variable_by_index(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_dobldobl_Laurent_drop_variable_by_name
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *s;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&s)) return NULL;
   fail = syscon_dobldobl_Laurent_drop_variable_by_name(nc,s);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_quaddobl_Laurent_drop_variable_by_index
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = syscon_quaddobl_Laurent_drop_variable_by_index(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_quaddobl_Laurent_drop_variable_by_name
 ( PyObject *self, PyObject *args )
{
   int nc,fail;
   char *s;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nc,&s)) return NULL;
   fail = syscon_quaddobl_Laurent_drop_variable_by_name(nc,s);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_standard_one_homogenization
 ( PyObject *self, PyObject *args )
{
   int fail,lintype;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&lintype)) return NULL;

   fail = syscon_standard_one_homogenization(lintype);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_dobldobl_one_homogenization
 ( PyObject *self, PyObject *args )
{
   int fail,lintype;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&lintype)) return NULL;

   fail = syscon_dobldobl_one_homogenization(lintype);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_quaddobl_one_homogenization
 ( PyObject *self, PyObject *args )
{
   int fail,lintype;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&lintype)) return NULL;

   fail = syscon_quaddobl_one_homogenization(lintype);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_add_symbol
 ( PyObject *self, PyObject *args )
{
   int fail,nbc;
   char *name;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nbc,&name)) return NULL;

   fail = syscon_add_symbol(nbc,name);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_standard_one_affinization
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if (!PyArg_ParseTuple(args,"")) return NULL;

   fail = syscon_standard_one_affinization();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_dobldobl_one_affinization
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if (!PyArg_ParseTuple(args,"")) return NULL;

   fail = syscon_dobldobl_one_affinization();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syscon_quaddobl_one_affinization
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if (!PyArg_ParseTuple(args,"")) return NULL;

   fail = syscon_quaddobl_one_affinization();

   return Py_BuildValue("i",fail);
}

/* The wrapping of the functions in tabform.h starts from here */

static PyObject *py2c_tabform_store_standard_tableau
 ( PyObject *self, PyObject *args )
{
   int fail,neq,nvr,nc1,nc2,nc3,verbose;
   char *nbt,*cff,*xpc;

   initialize();
   if(!PyArg_ParseTuple(args,"iiisisisi",&neq,&nvr,&nc1,&nbt,&nc2,&cff,
      &nc3,&xpc,&verbose)) return NULL;
   {
      if(verbose > 0)
      {
         printf("the string of terms : %s\n", nbt);
         printf("the string of coefficients : %s\n", cff);
         printf("the string of exponents : %s\n", xpc);
      }
      int ic1 = itemcount(nbt);
      int ic2 = itemcount(xpc);
      int ic3 = itemcount(cff);

      if(verbose > 0)
      {
         printf("number of items in nbt : %d\n", ic1);
         printf("number of items in xpc : %d\n", ic2);
         printf("number of items in cff : %d\n", ic3);
      }
      int nbterms[ic1];
      int exponents[ic2];
      double coefficients[ic3];

      str2intlist(ic1,nbt,nbterms);
      str2intlist(ic2,xpc,exponents);
      str2dbllist(ic3,cff,coefficients);

      if(verbose > 0)
      {
         int idx;

         printf("the number of terms : ");
         for(idx=0; idx<ic1; idx++) printf(" %d", nbterms[idx]);
         printf("\n");
         printf("the coefficients : ");
         for(idx=0; idx<ic3; idx++) printf(" %.15le", coefficients[idx]);
         printf("\n");
         printf("the exponents : ");
         for(idx=0; idx<ic2; idx++) printf(" %d", exponents[idx]);
         printf("\n");
      }
      fail = store_standard_tableau_form
                (neq,nvr,nbterms,coefficients,exponents,verbose);
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_tabform_load_standard_tableau
 ( PyObject *self, PyObject *args )
{
   int fail,neq,nvr,nbt,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&verbose)) return NULL;

   fail = load_standard_tableau_dimensions(&neq,&nvr,&nbt);
   if(verbose > 0)
   {
      printf("the number of equations : %d\n", neq);
      printf("the number of variables : %d\n", nvr);
      printf("total number of terms : %d\n", nbt);
   }
   {
      int nbterms[neq];
      const int nbexp = nvr*nbt;
      const int nbcff = 2*nbt;
      int exponents[nbexp];
      double coefficients[nbcff];
      char strnbt[neq*8];
      char strxps[nbexp*8];
      char strcff[nbcff*26];

      fail = number_of_standard_terms(neq,nbterms,&nbt,verbose);
      fail = standard_tableau_form
               (neq,nvr,nbterms,coefficients,exponents,verbose);

      intlist2str(neq, nbterms, strnbt);
      intlist2str(nbexp, exponents, strxps);
      dbllist2str(nbcff, coefficients, strcff);
      if(verbose > 0)
      {
         int idx;
         printf("the number of terms : ");
         for(idx=0; idx<neq; idx++) printf(" %d", nbterms[idx]);
         printf("\n");
         printf("its string representation : %s\n", strnbt);
         printf("the coefficients : ");
         for(idx=0; idx<nbcff; idx++) printf(" %.15le", coefficients[idx]);
         printf("\n");
         printf("its string representation : %s\n", strcff);
         printf("the exponents : ");
         for(idx=0; idx<nbexp; idx++) printf(" %d", exponents[idx]);
         printf("\n");
         printf("its string representation : %s\n", strxps);
      }

      return Py_BuildValue("(i, i, s, s, s)",
                           neq, nvr, strnbt, strcff, strxps);
   }
}

/* The wrapping of the functions in solcon.h starts from here. */

static PyObject *py2c_solcon_read_standard_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_read_standard_solutions();
              
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

static PyObject *py2c_solcon_read_standard_solutions_from_file
 ( PyObject *self, PyObject *args )
{
   int nbc,fail;
   char *name;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nbc,&name)) return NULL;   
   fail = solcon_read_standard_solutions_from_file(nbc,name);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_read_dobldobl_solutions_from_file
 ( PyObject *self, PyObject *args )
{
   int nbc,fail;
   char *name;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nbc,&name)) return NULL;   
   fail = solcon_read_dobldobl_solutions_from_file(nbc,name);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_read_quaddobl_solutions_from_file
 ( PyObject *self, PyObject *args )
{
   int nbc,fail;
   char *name;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&nbc,&name)) return NULL;   
   fail = solcon_read_quaddobl_solutions_from_file(nbc,name);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_write_standard_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_write_standard_solutions();
              
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

static PyObject *py2c_solcon_clear_standard_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_clear_standard_solutions();
              
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

static PyObject *py2c_solcon_length_standard_solution_string
 ( PyObject *self, PyObject *args )
{
   int fail,i,number;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&i)) return NULL;
   fail = solcon_length_standard_solution_string(i,&number);

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

static PyObject *py2c_solcon_write_standard_solution_string
 ( PyObject *self, PyObject *args )
{      
   int fail;
   int n,k;
   char *p;
   static PyObject* a;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"ii",&n,&k)) return NULL;
   p = (char*)malloc((k+1)*sizeof(char));
   fail = solcon_write_standard_solution_string(n,k,p);
                 
   a = Py_BuildValue("s",p);

   free(p);
   
   return a;
}

static PyObject *py2c_solcon_write_dobldobl_solution_string
 ( PyObject *self, PyObject *args )
{      
   int fail;
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
   int fail;
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
   int fail;
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

static PyObject *py2c_solcon_retrieve_next_standard_initialize
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_retrieve_next_standard_initialize();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_retrieve_next_dobldobl_initialize
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_retrieve_next_dobldobl_initialize();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_retrieve_next_quaddobl_initialize
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_retrieve_next_quaddobl_initialize();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_retrieve_next_multprec_initialize
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_retrieve_next_multprec_initialize();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_move_current_standard_to_next
 ( PyObject *self, PyObject *args )
{
   int fail,cursor;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_move_current_standard_to_next(&cursor);
              
   return Py_BuildValue("i",cursor);
}

static PyObject *py2c_solcon_move_current_dobldobl_to_next
 ( PyObject *self, PyObject *args )
{
   int fail,cursor;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_move_current_dobldobl_to_next(&cursor);
              
   return Py_BuildValue("i",cursor);
}

static PyObject *py2c_solcon_move_current_quaddobl_to_next
 ( PyObject *self, PyObject *args )
{
   int fail,cursor;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_move_current_quaddobl_to_next(&cursor);
              
   return Py_BuildValue("i",cursor);
}

static PyObject *py2c_solcon_move_current_multprec_to_next
 ( PyObject *self, PyObject *args )
{
   int fail,cursor;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_move_current_multprec_to_next(&cursor);
              
   return Py_BuildValue("i",cursor);
}

static PyObject *py2c_solcon_length_current_standard_solution_string
 ( PyObject *self, PyObject *args )
{
   int fail,cursor,len;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_length_current_standard_solution_string(&cursor,&len);
   if(cursor == 0) len = 0;

   return Py_BuildValue("i",len);
}

static PyObject *py2c_solcon_length_current_dobldobl_solution_string
 ( PyObject *self, PyObject *args )
{
   int fail,cursor,len;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_length_current_dobldobl_solution_string(&cursor,&len);
   if(cursor == 0) len = 0;

   return Py_BuildValue("i",len);
}

static PyObject *py2c_solcon_length_current_quaddobl_solution_string
 ( PyObject *self, PyObject *args )
{
   int fail,cursor,len;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_length_current_quaddobl_solution_string(&cursor,&len);
   if(cursor == 0) len = 0;

   return Py_BuildValue("i",len);
}

static PyObject *py2c_solcon_length_current_multprec_solution_string
 ( PyObject *self, PyObject *args )
{
   int fail,cursor,len;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = solcon_length_current_multprec_solution_string(&cursor,&len);
   if(cursor == 0) len = 0;

   return Py_BuildValue("i",len);
}

static PyObject *py2c_solcon_write_current_standard_solution_string
 ( PyObject *self, PyObject *args )
{      
   int fail;
   int n,k;
   char *p;
   static PyObject* a;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&n)) return NULL;
   p = (char*)malloc((n+1)*sizeof(char));
   fail = solcon_write_current_standard_solution_string(&k,n,p);
                 
   a = Py_BuildValue("s",p);

   free(p);
   
   return a;
}

static PyObject *py2c_solcon_write_current_dobldobl_solution_string
 ( PyObject *self, PyObject *args )
{      
   int fail;
   int n,k;
   char *p;
   static PyObject* a;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&n)) return NULL;
   p = (char*)malloc((n+1)*sizeof(char));
   fail = solcon_write_current_dobldobl_solution_string(&k,n,p);
                 
   a = Py_BuildValue("s",p);

   free(p);
   
   return a;
}

static PyObject *py2c_solcon_write_current_quaddobl_solution_string
 ( PyObject *self, PyObject *args )
{      
   int fail;
   int n,k;
   char *p;
   static PyObject* a;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&n)) return NULL;
   p = (char*)malloc((n+1)*sizeof(char));
   fail = solcon_write_current_quaddobl_solution_string(&k,n,p);
                 
   a = Py_BuildValue("s",p);

   free(p);
   
   return a;
}

static PyObject *py2c_solcon_write_current_multprec_solution_string
 ( PyObject *self, PyObject *args )
{      
   int fail;
   int n,k;
   char *p;
   static PyObject* a;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&n)) return NULL;
   p = (char*)malloc((n+1)*sizeof(char));
   fail = solcon_write_current_multprec_solution_string(&k,n,p);
                 
   a = Py_BuildValue("s",p);

   free(p);
   
   return a;
}

static PyObject *py2c_solcon_append_standard_solution_string
 ( PyObject *self, PyObject *args )
{      
   char fail;
   int n,k;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iis",&n,&k,&p)) return NULL;

   /* printf("Calling solcon_append_solution_string ...\n"); */

   fail = solcon_append_standard_solution_string(n,k,p);

   if(fail != 0) printf("Failed to append solution %s.\n",p);

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

   if(fail != 0) printf("Failed to append solution %s.\n",p);
                 
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

   if(fail != 0) printf("Failed to append solution %s.\n",p);
                 
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

   if(fail != 0) printf("Failed to append solution %s.\n",p);
                 
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_number_of_standard_solutions
 ( PyObject *self, PyObject *args )
{      
   int result,n;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   result = solcon_number_of_standard_solutions(&n);
                 
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

static PyObject *py2c_solcon_standard_one_homogenization
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = solcon_standard_one_homogenization();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_dobldobl_one_homogenization
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = solcon_dobldobl_one_homogenization();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_quaddobl_one_homogenization
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = solcon_quaddobl_one_homogenization();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_standard_one_affinization
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = solcon_standard_one_affinization();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_dobldobl_one_affinization
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = solcon_dobldobl_one_affinization();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_solcon_quaddobl_one_affinization
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = solcon_quaddobl_one_affinization();

   return Py_BuildValue("i",fail);
}

/* The wrapping of functions with prototypes in product.h starts here. */

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

static PyObject *py2c_celcon_initialize_supports
 ( PyObject *self, PyObject *args )
{
   int fail,nbr;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&nbr)) return NULL;
   fail = celcon_initialize_supports(nbr);

   return Py_BuildValue("i",fail);
}

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

static PyObject *py2c_celcon_type_of_mixture
 ( PyObject *self, PyObject *args )
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

static PyObject *py2c_celcon_retrieve_lifted_point
 ( PyObject *self, PyObject *args )
{
   int fail,dim,sup,ind,nbc;
   char strpoint[256];

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&dim,&sup,&ind)) return NULL;
   {
      double liftedpoint[dim];

      fail = celcon_get_lifted_point(dim,sup,ind,liftedpoint);
      nbc = dbllist2str(dim,liftedpoint,strpoint);
   }

   return Py_BuildValue("s",strpoint);
}

static PyObject *py2c_celcon_mixed_volume_of_supports
 ( PyObject *self, PyObject *args )
{
   int fail,mv;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_mixed_volume_of_supports(&mv);

   return Py_BuildValue("i",mv);
}

static PyObject *py2c_celcon_number_of_cells
 ( PyObject *self, PyObject *args )
{
   int fail,length;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_number_of_cells(&length);

   return Py_BuildValue("i",length);
}

static PyObject *py2c_celcon_is_stable ( PyObject *self, PyObject *args )
{
   int fail,result;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_is_stable(&result);

   return Py_BuildValue("i",result);
}

static PyObject *py2c_celcon_number_of_original_cells
 ( PyObject *self, PyObject *args )
{
   int fail,result;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_number_of_original_cells(&result);

   return Py_BuildValue("i",result);
}

static PyObject *py2c_celcon_number_of_stable_cells
 ( PyObject *self, PyObject *args )
{
   int fail,result;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_number_of_stable_cells(&result);

   return Py_BuildValue("i",result);
}

static PyObject *py2c_celcon_standard_random_coefficient_system 
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_standard_random_coefficient_system();

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

static PyObject *py2c_celcon_copy_into_standard_systems_container
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_copy_into_standard_systems_container();

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

static PyObject *py2c_celcon_standard_polyhedral_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;
  
   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_standard_polyhedral_homotopy();

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

static PyObject *py2c_celcon_solve_standard_start_system
 ( PyObject *self, PyObject *args )
{
   int fail,k,nb;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = celcon_solve_standard_start_system(k,&nb);

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

static PyObject *py2c_celcon_solve_stable_standard_start_system
 ( PyObject *self, PyObject *args )
{
   int fail,k,nb;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = celcon_solve_stable_standard_start_system(k,&nb);

   return Py_BuildValue("i",nb);
}

static PyObject *py2c_celcon_solve_stable_dobldobl_start_system
 ( PyObject *self, PyObject *args )
{
   int fail,k,nb;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = celcon_solve_stable_dobldobl_start_system(k,&nb);

   return Py_BuildValue("i",nb);
}

static PyObject *py2c_celcon_solve_stable_quaddobl_start_system
 ( PyObject *self, PyObject *args )
{
   int fail,k,nb;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = celcon_solve_stable_quaddobl_start_system(k,&nb);

   return Py_BuildValue("i",nb);
}

static PyObject *py2c_celcon_track_standard_solution_path
 ( PyObject *self, PyObject *args )
{
   int fail,k,i,otp;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&k,&i,&otp)) return NULL;
   fail = celcon_track_standard_solution_path(k,i,otp);

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

static PyObject *py2c_celcon_copy_target_standard_solution_to_container
 ( PyObject *self, PyObject *args )
{
   int fail,k,i;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&k,&i)) return NULL;
   fail = celcon_copy_target_standard_solution_to_container(k,i);

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

static PyObject *py2c_celcon_permute_standard_system
 ( PyObject *self, PyObject *args )
{
   int fail,k,i;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = celcon_permute_standard_system();

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

/* The wrapping of functions with prototypes in intcelcon.h follows. */

static PyObject *py2c_intcelcon_read_mixed_cell_configuration
 ( PyObject *self, PyObject *args )
{
   int fail = 0;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   

   fail = intcelcon_read_mixed_cell_configuration();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_intcelcon_write_mixed_cell_configuration
 ( PyObject *self, PyObject *args )
{
   int fail = 0;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   

   fail = intcelcon_write_mixed_cell_configuration();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_intcelcon_number_of_cells
 ( PyObject *self, PyObject *args )
{
   int fail,length;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = intcelcon_number_of_cells(&length);

   return Py_BuildValue("i",length);
}

static PyObject *py2c_intcelcon_type_of_mixture
 ( PyObject *self, PyObject *args )
{
   int fail,r,nbc;
   int mix[64];
   char strmix[256];

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;

   fail = intcelcon_type_of_mixture(&r,mix);
   nbc = intlist2str(r,mix,strmix);

   return Py_BuildValue("s",strmix);
}

static PyObject *py2c_intcelcon_length_of_supports
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   {
      int dim,nbr,nbc,i;
      int mix[64];

      fail = intcelcon_type_of_mixture(&nbr,mix);

      int lengths[nbr];
      char strlengths[8*nbr];

      fail = intcelcon_length_of_supports(&nbr,lengths);
      nbc = intlist2str(nbr,lengths,strlengths);

      return Py_BuildValue("s",strlengths);
   }
}

static PyObject *py2c_intcelcon_get_lifted_point
 ( PyObject *self, PyObject *args )
{
   int fail,dim,sup,ind,nbc;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&dim,&sup,&ind)) return NULL;
   {
      int liftedpoint[dim];
      char strpoint[8*dim];

      fail = intcelcon_get_lifted_point(dim,sup,ind,liftedpoint);
      nbc = intlist2str(dim,liftedpoint,strpoint);

      return Py_BuildValue("s",strpoint);
   }
}

static PyObject *py2c_intcelcon_get_inner_normal
 ( PyObject *self, PyObject *args )
{
   int fail,dim,idx,nbc;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&dim,&idx)) return NULL;
   {
      int normal[dim];
      char strnormal[8*dim];

      fail = intcelcon_get_inner_normal(dim,idx,normal);
      nbc = intlist2str(dim,normal,strnormal);

      return Py_BuildValue("s",strnormal);
   }
}

static PyObject *py2c_intcelcon_number_of_points_in_cell
 ( PyObject *self, PyObject *args )
{
   int fail,idx,nbr,nbc;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&idx,&nbr)) return NULL;
   {
      int lengths[nbr];
      char strlengths[8*nbr];

      fail = intcelcon_number_of_points_in_cell(idx,lengths);
      nbc = intlist2str(nbr,lengths,strlengths);

      return Py_BuildValue("s",strlengths);
   }
}

static PyObject *py2c_intcelcon_get_point_in_cell 
 ( PyObject *self, PyObject *args )
{
   int fail,dim,idcel,idsup,idpnt,nbc;

   initialize();
   if(!PyArg_ParseTuple(args,"iiii",&dim,&idcel,&idsup,&idpnt)) return NULL;
   {
      int point[dim];
      char strpoint[8*dim];

      fail = intcelcon_get_point_in_cell(dim,idcel,idsup,idpnt,point);
      nbc = intlist2str(dim,point,strpoint);

      return Py_BuildValue("s",strpoint);
   }
}

static PyObject *py2c_intcelcon_mixed_volume
 ( PyObject *self, PyObject *args )
{
   int fail,idx,mv;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&idx)) return NULL;
   fail = intcelcon_mixed_volume(idx,&mv);

   return Py_BuildValue("i",mv);
}

static PyObject *py2c_intcelcon_initialize_supports
 ( PyObject *self, PyObject *args )
{
   int fail,nbr;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&nbr)) return NULL;
   fail = intcelcon_initialize_supports(nbr);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_intcelcon_set_type_of_mixture
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
      fail = intcelcon_set_type_of_mixture(r,mix);
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_intcelcon_append_lifted_point
 ( PyObject *self, PyObject *args )
{
   int fail,dim,ind,cnt;
   char *strpoint;

   initialize();
   if(!PyArg_ParseTuple(args,"iis",&dim,&ind,&strpoint)) return NULL;

   cnt = itemcount(strpoint);
   {
      int point[cnt];

      str2intlist(cnt,strpoint,point);
      fail = intcelcon_append_lifted_point(dim,ind,point);
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_intcelcon_make_subdivision
 ( PyObject *self, PyObject *args )
{
   int fail = 0;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   

   fail = intcelcon_make_subdivision();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_intcelcon_clear_mixed_cell_configuration
 ( PyObject *self, PyObject *args )
{
   int fail = 0;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   

   fail = intcelcon_clear_mixed_cell_configuration();

   return Py_BuildValue("i",fail);
}

/* The wrapping of the functions with prototypes in scalers.h starts here. */

static PyObject *py2c_scale_standard_system
 ( PyObject *self, PyObject *args )
{
   int fail,mode,dim,i;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&mode)) return NULL;
   fail = syscon_number_of_standard_polynomials(&dim);
   if((fail == 0) && (dim > 0))
   {
      double cff[4*dim+2];  

      fail = standard_scale_system(mode,cff);
      if(fail == 0)
      {
         if(mode > 0)
         {
            PyObject *result, *item;
            result = PyList_New(4*dim+1);
            for(i=0; i<4*dim+1; i++)
            {
               item = PyFloat_FromDouble(cff[i]);
               PyList_SET_ITEM(result,i,item);
            }
            return result;
         }
      }
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_scale_dobldobl_system
 ( PyObject *self, PyObject *args )
{
   int fail,mode,dim,i;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&mode)) return NULL;
   fail = syscon_number_of_dobldobl_polynomials(&dim);
   if((fail == 0) && (dim > 0))
   {
      double cff[8*dim+4];  

      fail = dobldobl_scale_system(mode,cff);
      if(fail == 0)
      {
         if(mode > 0)
         {
            PyObject *result, *item;
            result = PyList_New(8*dim+1);
            for(i=0; i<8*dim+1; i++)
            {
               item = PyFloat_FromDouble(cff[i]);
               PyList_SET_ITEM(result,i,item);
            }
            return result;
         }
      }
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_scale_quaddobl_system
 ( PyObject *self, PyObject *args )
{
   int fail,mode,dim,i;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&mode)) return NULL;
   fail = syscon_number_of_quaddobl_polynomials(&dim);
   if((fail == 0) && (dim > 0))
   {
      double cff[16*dim+8];  

      fail = quaddobl_scale_system(mode,cff);
      if(fail == 0)
      {
         if(mode > 0)
         {
            PyObject *result, *item;
            result = PyList_New(16*dim+1);
            for(i=0; i<16*dim+1; i++)
            {
               item = PyFloat_FromDouble(cff[i]);
               PyList_SET_ITEM(result,i,item);
            }
            return result;
         }
      }
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_scale_standard_solutions 
 ( PyObject *self, PyObject *args )
{
   int fail,dim;
   char *cff;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&dim,&cff)) return NULL;
   {
      double scf[dim];

      str2dbllist(dim,cff,scf);
      fail = standard_scale_solutions(dim,10,scf);
   }

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_scale_dobldobl_solutions 
 ( PyObject *self, PyObject *args )
{
   int fail,dim;
   char *cff;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&dim,&cff)) return NULL;
   {
      double scf[dim];

      str2dbllist(dim,cff,scf);
      fail = dobldobl_scale_solutions(dim,10,scf);
   }

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_scale_quaddobl_solutions 
 ( PyObject *self, PyObject *args )
{
   int fail,dim;
   char *cff;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&dim,&cff)) return NULL;
   {
      double scf[dim];

      str2dbllist(dim,cff,scf);
      fail = quaddobl_scale_solutions(dim,10,scf);
   }
   return Py_BuildValue("i",fail);
}

/* The wrapping of the functions with prototypes in reducers.h starts here. */

static PyObject *py2c_linear_reduce_standard_system
 ( PyObject *self, PyObject *args )
{
   int fail,diag;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&diag)) return NULL;

   fail = standard_row_reduce_system(diag);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_linear_reduce_dobldobl_system
 ( PyObject *self, PyObject *args )
{
   int fail,diag;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&diag)) return NULL;

   fail = dobldobl_row_reduce_system(diag);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_linear_reduce_quaddobl_system
 ( PyObject *self, PyObject *args )
{
   int fail,diag;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&diag)) return NULL;

   fail = quaddobl_row_reduce_system(diag);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_nonlinear_reduce_standard_system
 ( PyObject *self, PyObject *args )
{
   int fail,eqmax,spmax,rpmax,eqcnt,spcnt,rpcnt;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&eqmax,&spmax,&rpmax)) return NULL;

   fail = standard_nonlinear_reduce_system
             (eqmax,spmax,rpmax,&eqcnt,&spcnt,&rpcnt);

   return Py_BuildValue("(i,i,i)",eqcnt,spcnt,rpcnt);
}

/* The wrapping of the functions with prototypes in sweep.h starts here. */

static PyObject *py2c_sweep_define_parameters_numerically
 ( PyObject *self, PyObject *args )
{
   int fail,nq,nv,np,cnt;
   char *strpars;

   initialize();
   if(!PyArg_ParseTuple(args,"iiis",&nq,&nv,&np,&strpars)) return NULL;

   cnt = itemcount(strpars);
   {
      int pars[cnt];

      str2intlist(cnt,strpars,pars);
      fail = sweep_define_parameters_numerically(nq,nv,np,pars);
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_sweep_define_parameters_symbolically
 ( PyObject *self, PyObject *args )
{
   int fail,nq,nv,np,nc;
   char *pars;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiis",&nq,&nv,&np,&nc,&pars)) return NULL;

   fail = sweep_define_parameters_symbolically(nq,nv,np,nc,pars);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_sweep_get_number_of_equations
 ( PyObject *self, PyObject *args )
{
   int fail,nb;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   

   fail = sweep_get_number_of_equations(&nb);

   return Py_BuildValue("i",nb);
}

static PyObject *py2c_sweep_get_number_of_variables
 ( PyObject *self, PyObject *args )
{
   int fail,nb;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   

   fail = sweep_get_number_of_variables(&nb);

   return Py_BuildValue("i",nb);
}

static PyObject *py2c_sweep_get_number_of_parameters
 ( PyObject *self, PyObject *args )
{
   int fail,nb;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   

   fail = sweep_get_number_of_parameters(&nb);

   return Py_BuildValue("i",nb);
}

static PyObject *py2c_sweep_get_indices_numerically
 ( PyObject *self, PyObject *args )
{
   int fail,np,nc;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   

   fail = sweep_get_number_of_parameters(&np);
   {
      int idx[np];
      char stridx[np*10];

      fail = sweep_get_indices_numerically(idx);

      nc = intlist2str(np,idx,stridx);

      return Py_BuildValue("s",stridx);
   }
}

static PyObject *py2c_sweep_get_indices_symbolically
 ( PyObject *self, PyObject *args )
{
   int fail,np,nc;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   

   fail = sweep_get_number_of_parameters(&np);
   {
      char buf[np*20]; /* assumes no more than 20 characters per symbol */

      fail = sweep_get_indices_symbolically(&nc,buf);

      buf[nc] = '\0';

      return Py_BuildValue("s",buf);
   }
}

static PyObject *py2c_sweep_clear_definitions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = sweep_clear_definitions();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_sweep_set_standard_start
 ( PyObject *self, PyObject *args )
{
   int fail,m;
   char *strvals;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&m,&strvals)) return NULL;
   {
      const int n = 2*m;
      double vals[n];

      str2dbllist(n,strvals,vals);

      fail = sweep_set_standard_start(n,vals);
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_sweep_set_standard_target
 ( PyObject *self, PyObject *args )
{
   int fail,m;
   char *strvals;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&m,&strvals)) return NULL;
   {
      const int n = 2*m;
      double vals[n];

      str2dbllist(n,strvals,vals);

      fail = sweep_set_standard_target(n,vals);
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_sweep_set_dobldobl_start
 ( PyObject *self, PyObject *args )
{
   int fail,m;
   char *strvals;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&m,&strvals)) return NULL;
   {
      const int n = 4*m;
      double vals[n];

      str2dbllist(n,strvals,vals);

      fail = sweep_set_dobldobl_start(n,vals);
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_sweep_set_dobldobl_target
 ( PyObject *self, PyObject *args )
{
   int fail,m;
   char *strvals;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&m,&strvals)) return NULL;
   {
      const int n = 4*m;
      double vals[n];

      str2dbllist(n,strvals,vals);

      fail = sweep_set_dobldobl_target(n,vals);
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_sweep_set_quaddobl_start
 ( PyObject *self, PyObject *args )
{
   int fail,m;
   char *strvals;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&m,&strvals)) return NULL;
   {
      const int n = 8*m;
      double vals[n];

      str2dbllist(n,strvals,vals);

      fail = sweep_set_quaddobl_start(n,vals);
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_sweep_set_quaddobl_target
 ( PyObject *self, PyObject *args )
{
   int fail,m;
   char *strvals;

   initialize();
   if(!PyArg_ParseTuple(args,"is",&m,&strvals)) return NULL;
   {
      const int n = 8*m;
      double vals[n];

      str2dbllist(n,strvals,vals);

      fail = sweep_set_quaddobl_target(n,vals);
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_sweep_get_standard_start
 ( PyObject *self, PyObject *args )
{
   int fail,n,nc;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&n)) return NULL;
   {
      double cff[n];
      char strcff[n*25];

      fail = sweep_get_standard_start(n,cff);

      nc = dbllist2str(n,cff,strcff);

      return Py_BuildValue("s",strcff);
   }
}

static PyObject *py2c_sweep_get_standard_target
 ( PyObject *self, PyObject *args )
{
   int fail,n,nc;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&n)) return NULL;
   {
      double cff[n];
      char strcff[n*25];

      fail = sweep_get_standard_target(n,cff);

      nc = dbllist2str(n,cff,strcff);

      return Py_BuildValue("s",strcff);
   }
}

static PyObject *py2c_sweep_get_dobldobl_start
 ( PyObject *self, PyObject *args )
{
   int fail,n,nc;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&n)) return NULL;
   {
      double cff[n];
      char strcff[n*25];

      fail = sweep_get_dobldobl_start(n,cff);

      nc = dbllist2str(n,cff,strcff);

      return Py_BuildValue("s",strcff);
   }
}

static PyObject *py2c_sweep_get_dobldobl_target
 ( PyObject *self, PyObject *args )
{
   int fail,n,nc;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&n)) return NULL;
   {
      double cff[n];
      char strcff[n*25];

      fail = sweep_get_dobldobl_target(n,cff);

      nc = dbllist2str(n,cff,strcff);

      return Py_BuildValue("s",strcff);
   }
}

static PyObject *py2c_sweep_get_quaddobl_start
 ( PyObject *self, PyObject *args )
{
   int fail,n,nc;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&n)) return NULL;
   {
      double cff[n];
      char strcff[n*25];

      fail = sweep_get_quaddobl_start(n,cff);

      nc = dbllist2str(n,cff,strcff);

      return Py_BuildValue("s",strcff);
   }
}

static PyObject *py2c_sweep_get_quaddobl_target
 ( PyObject *self, PyObject *args )
{
   int fail,n,nc;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&n)) return NULL;
   {
      double cff[n];
      char strcff[n*25];

      fail = sweep_get_quaddobl_target(n,cff);

      nc = dbllist2str(n,cff,strcff);

      return Py_BuildValue("s",strcff);
   }
}

static PyObject *py2c_sweep_standard_complex_run
 ( PyObject *self, PyObject *args )
{
   int fail,choice;
   double g_re,g_im;

   initialize();
   if(!PyArg_ParseTuple(args,"idd",&choice,&g_re,&g_im)) return NULL;   

   fail = sweep_standard_complex_run(choice,&g_re,&g_im);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_sweep_dobldobl_complex_run
 ( PyObject *self, PyObject *args )
{
   int fail,choice;
   double g_re,g_im;

   initialize();
   if(!PyArg_ParseTuple(args,"idd",&choice,&g_re,&g_im)) return NULL;   
   {
      if(choice < 2)
         fail = sweep_dobldobl_complex_run(choice,&g_re,&g_im);
      else
      {
         double regamma[2];
         double imgamma[2];
         regamma[0] = g_re; regamma[1] = 0.0;
         imgamma[0] = g_im; imgamma[1] = 0.0;
         fail = sweep_dobldobl_complex_run(choice,regamma,imgamma);
      }
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_sweep_quaddobl_complex_run
 ( PyObject *self, PyObject *args )
{
   int fail,choice;
   double g_re,g_im;

   initialize();
   if(!PyArg_ParseTuple(args,"idd",&choice,&g_re,&g_im)) return NULL;   
   {
      if(choice < 2)
         fail = sweep_quaddobl_complex_run(choice,&g_re,&g_im);
      else
      {
         double regamma[4];
         double imgamma[4];
         regamma[0] = g_re; regamma[1] = 0.0;
         regamma[2] = 0.0;  regamma[3] = 0.0;
         imgamma[0] = g_im; imgamma[1] = 0.0;
         imgamma[2] = 0.0;  imgamma[3] = 0.0;
         fail = sweep_quaddobl_complex_run(choice,regamma,imgamma);
      }
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_sweep_standard_real_run
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = sweep_standard_real_run();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_sweep_dobldobl_real_run
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = sweep_dobldobl_real_run();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_sweep_quaddobl_real_run
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = sweep_quaddobl_real_run();

   return Py_BuildValue("i",fail);
}

/* The wrapping for the multiplicity structure starts here. */

static PyObject *py2c_standard_multiplicity_structure
 ( PyObject *self, PyObject *args )
{
   int fail,order,verbose;
   double tol;

   initialize();
   if(!PyArg_ParseTuple(args,"iid",&order,&verbose,&tol)) return NULL;   
   {
      int nbc,mult;
      int hilb[order+1];
      char strhilb[4*(order+1)];

      fail = standard_multiplicity_structure(order,verbose,tol,&mult,hilb);
      nbc = intlist2str(order+1,hilb,strhilb);
   
      return Py_BuildValue("(i,s)",mult,strhilb);
   }
}

static PyObject *py2c_dobldobl_multiplicity_structure
 ( PyObject *self, PyObject *args )
{
   int fail,order,verbose,mult;
   double tol;

   initialize();
   if(!PyArg_ParseTuple(args,"iid",&order,&verbose,&tol)) return NULL;   
   {
      int nbc,mult;
      int hilb[order+1];
      char strhilb[4*(order+1)];
     
      fail = dobldobl_multiplicity_structure(order,verbose,tol,&mult,hilb);
      nbc = intlist2str(order+1,hilb,strhilb);

      return Py_BuildValue("(i,s)",mult,strhilb);
   }
}

static PyObject *py2c_quaddobl_multiplicity_structure
 ( PyObject *self, PyObject *args )
{
   int fail,order,verbose,mult;
   double tol;

   initialize();
   if(!PyArg_ParseTuple(args,"iid",&order,&verbose,&tol)) return NULL;   
   {
      int nbc,mult;
      int hilb[order+1];
      char strhilb[4*(order+1)];

      fail = quaddobl_multiplicity_structure(order,verbose,tol,&mult,hilb);
      nbc = intlist2str(order+1,hilb,strhilb);

      return Py_BuildValue("(i,s)",mult,strhilb);
   }
}

/* The wrapping of the numerical tropisms container starts here. */

static PyObject *py2c_numbtrop_standard_initialize
 ( PyObject *self, PyObject *args )
{
   int fail,nbt,dim,k;
   char *data; /* all numbers come in one long string */

   initialize();
   if(!PyArg_ParseTuple(args,"iis",&nbt,&dim,&data)) return NULL;   
   {
      const int lendata = nbt*(dim+2);
      double numbers[lendata];
      int wnd[nbt];
      double dir[nbt*dim];
      double err[nbt];

      str2dbllist(lendata,data,numbers);

      for(k=0; k<nbt; k++) wnd[k] = (int) numbers[k];
      for(k=0; k<nbt*dim; k++) dir[k] = numbers[nbt+k];
      for(k=0; k<nbt; k++) err[k] = numbers[nbt*(dim+1)+k];

      fail = numbtrop_standard_initialize(nbt,dim,wnd,dir,err);     
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_numbtrop_dobldobl_initialize
 ( PyObject *self, PyObject *args )
{
   int fail,nbt,dim,k;;
   char *data; /* all numbers come in one long string */

   initialize();
   if(!PyArg_ParseTuple(args,"iis",&nbt,&dim,&data)) return NULL;   
   {
      const int lendata = nbt + 2*nbt*(dim+1);
      double numbers[lendata];
      int wnd[nbt];
      double dir[2*nbt*dim];
      double err[2*nbt];

      str2dbllist(lendata,data,numbers);

      for(k=0; k<nbt; k++) wnd[k] = (int) numbers[k];
      for(k=0; k<2*nbt*dim; k++) dir[k] = numbers[nbt+k];
      for(k=0; k<2*nbt; k++) err[k] = numbers[nbt+2*nbt*dim+k];

      fail = numbtrop_dobldobl_initialize(nbt,dim,wnd,dir,err);     
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_numbtrop_quaddobl_initialize
 ( PyObject *self, PyObject *args )
{
   int fail,nbt,dim,k;
   char *data; /* all numbers come in one long string */

   initialize();
   if(!PyArg_ParseTuple(args,"iis",&nbt,&dim,&data)) return NULL;   
   {
      const int lendata = nbt + 4*nbt*(dim+1);
      double numbers[lendata];
      int wnd[nbt];
      double dir[4*nbt*dim];
      double err[4*nbt];

      str2dbllist(lendata,data,numbers);

      for(k=0; k<nbt; k++) wnd[k] = (int) numbers[k];
      for(k=0; k<4*nbt*dim; k++) dir[k] = numbers[nbt+k];
      for(k=0; k<4*nbt; k++) err[k] = numbers[nbt+4*nbt*dim+k];

      fail = numbtrop_quaddobl_initialize(nbt,dim,wnd,dir,err);     
   }

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_numbtrop_standard_retrieve
 ( PyObject *self, PyObject *args )
{
   int fail,nbt,dim,k,nbc;
   char *fltlist;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&nbt,&dim)) return NULL;   
   {
      int wnd[nbt];
      double dir[nbt*dim];
      double err[nbt];
      double data[nbt*(dim+2)];

      fail = numbtrop_standard_retrieve(nbt,dim,wnd,dir,err);

      for(k=0; k<nbt; k++) data[k] = (double) wnd[k];
      for(k=0; k<nbt*dim; k++) data[nbt+k] = dir[k];
      for(k=0; k<nbt; k++) data[nbt*(dim+1)+k] = err[k];

      fltlist = (char*)calloc(25*nbt*(dim+2), sizeof(char));
      nbc = dbllist2str(nbt*(dim+2),data,fltlist);
   }
   return Py_BuildValue("(i,s)",fail,fltlist);
}

static PyObject *py2c_numbtrop_dobldobl_retrieve
 ( PyObject *self, PyObject *args )
{
   int fail,nbt,dim,k,nbc;
   char *fltlist;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&nbt,&dim)) return NULL;   
   {
      int wnd[nbt];
      double dir[2*nbt*dim];
      double err[2*nbt];
      const int lendata = nbt+2*nbt*(dim+1);
      double data[lendata];

      fail = numbtrop_dobldobl_retrieve(nbt,dim,wnd,dir,err);

      for(k=0; k<nbt; k++) data[k] = (double) wnd[k];
      for(k=0; k<2*nbt*dim; k++) data[nbt+k] = dir[k];
      for(k=0; k<2*nbt; k++) data[nbt+2*nbt*dim+k] = err[k];

      fltlist = (char*)calloc(50*nbt*(dim+2), sizeof(char));
      nbc = dbllist2str(lendata,data,fltlist);
   }
   return Py_BuildValue("(i,s)",fail,fltlist);
}

static PyObject *py2c_numbtrop_quaddobl_retrieve
 ( PyObject *self, PyObject *args )
{
   int fail,nbt,dim,k,nbc;
   char *fltlist;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&nbt,&dim)) return NULL;   
   {
      int wnd[nbt];
      double dir[4*nbt*dim];
      double err[4*nbt];
      const int lendata = nbt+4*nbt*(dim+1);
      double data[lendata];

      fail = numbtrop_quaddobl_retrieve(nbt,dim,wnd,dir,err);

      for(k=0; k<nbt; k++) data[k] = (double) wnd[k];
      for(k=0; k<4*nbt*dim; k++) data[nbt+k] = dir[k];
      for(k=0; k<4*nbt; k++) data[nbt+4*nbt*dim+k] = err[k];

      fltlist = (char*)calloc(100*nbt*(dim+2), sizeof(char));
      nbc = dbllist2str(lendata,data,fltlist);
   }
   return Py_BuildValue("(i,s)",fail,fltlist);
}

static PyObject *py2c_numbtrop_standard_size
 ( PyObject *self, PyObject *args )
{
   int fail,nbt;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = numbtrop_standard_size(&nbt);

   return Py_BuildValue("i",nbt);
}

static PyObject *py2c_numbtrop_dobldobl_size
 ( PyObject *self, PyObject *args )
{
   int fail,nbt;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = numbtrop_dobldobl_size(&nbt);

   return Py_BuildValue("i",nbt);
}

static PyObject *py2c_numbtrop_quaddobl_size
 ( PyObject *self, PyObject *args )
{
   int fail,nbt;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = numbtrop_quaddobl_size(&nbt);

   return Py_BuildValue("i",nbt);
}

static PyObject *py2c_numbtrop_standard_dimension
 ( PyObject *self, PyObject *args )
{
   int fail,dim;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = numbtrop_standard_dimension(&dim);

   return Py_BuildValue("i",dim);
}

static PyObject *py2c_numbtrop_dobldobl_dimension
 ( PyObject *self, PyObject *args )
{
   int fail,dim;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = numbtrop_dobldobl_dimension(&dim);

   return Py_BuildValue("i",dim);
}

static PyObject *py2c_numbtrop_quaddobl_dimension
 ( PyObject *self, PyObject *args )
{
   int fail,dim;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = numbtrop_quaddobl_dimension(&dim);

   return Py_BuildValue("i",dim);
}

static PyObject *py2c_numbtrop_store_standard_tropism
 ( PyObject *self, PyObject *args )
{
   int fail,dim,idx,wnd,k;
   char *data; /* all double numbers come in one long string */

   initialize();
   if(!PyArg_ParseTuple(args,"iiis",&dim,&idx,&wnd,&data)) return NULL;   
   {
      const int lendata = dim+1;
      double numbers[lendata];
      double dir[dim];
      double err;

      str2dbllist(lendata,data,numbers);

      for(k=0; k<dim; k++) dir[k] = numbers[k];
      err = numbers[dim];

      fail = numbtrop_store_standard_tropism(dim,idx,wnd,dir,err);     
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_numbtrop_store_dobldobl_tropism
 ( PyObject *self, PyObject *args )
{
   int fail,dim,idx,wnd,k;
   char *data; /* all double numbers come in one long string */

   initialize();
   if(!PyArg_ParseTuple(args,"iiis",&dim,&idx,&wnd,&data)) return NULL;   
   {
      const int lendata = 2*dim+2;
      double numbers[lendata];
      double dir[2*dim];
      double err[2];

      str2dbllist(lendata,data,numbers);

      for(k=0; k<2*dim; k++) dir[k] = numbers[k];
      err[0] = numbers[2*dim];
      err[1] = numbers[2*dim+1];

      fail = numbtrop_store_dobldobl_tropism(dim,idx,wnd,dir,err);     
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_numbtrop_store_quaddobl_tropism
 ( PyObject *self, PyObject *args )
{
   int fail,dim,idx,wnd,k;
   char *data; /* all double numbers come in one long string */

   initialize();
   if(!PyArg_ParseTuple(args,"iiis",&dim,&idx,&wnd,&data)) return NULL;   
   {
      const int lendata = 4*dim+4;
      double numbers[lendata];
      double dir[4*dim];
      double err[4];

      str2dbllist(lendata,data,numbers);

      for(k=0; k<4*dim; k++) dir[k] = numbers[k];
      err[0] = numbers[4*dim];
      err[1] = numbers[4*dim+1];
      err[2] = numbers[4*dim+2];
      err[3] = numbers[4*dim+3];

      fail = numbtrop_store_quaddobl_tropism(dim,idx,wnd,dir,err);     
   }
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_numbtrop_standard_retrieve_tropism
 ( PyObject *self, PyObject *args )
{
   int fail,dim,idx,wnd,k,nbc;
   char *fltlist;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&dim,&idx)) return NULL;   
   {
      double dir[dim];
      double err;
      double direrr[dim+1];

      fail = numbtrop_standard_retrieve_tropism(dim,idx,&wnd,dir,&err);

      fltlist = (char*)calloc(25*(dim+1), sizeof(char));
      for(k=0; k<dim; k++) direrr[k] = dir[k];
      direrr[dim] = err;
      nbc = dbllist2str(dim+1,direrr,fltlist);
   }
   return Py_BuildValue("(i,i,s)",fail,wnd,fltlist);
}

static PyObject *py2c_numbtrop_dobldobl_retrieve_tropism
 ( PyObject *self, PyObject *args )
{
   int fail,dim,idx,wnd,k,nbc;
   char *fltlist;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&dim,&idx)) return NULL;   
   {
      double dir[2*dim];
      double err[2];
      double direrr[2*dim+2];

      fail = numbtrop_dobldobl_retrieve_tropism(dim,idx,&wnd,dir,err);

      fltlist = (char*)calloc(50*(dim+1), sizeof(char));
      for(k=0; k<2*dim; k++) direrr[k] = dir[k];
      direrr[2*dim] = err[0];
      direrr[2*dim+1] = err[1];
      nbc = dbllist2str(2*dim+2,direrr,fltlist);
   }
   return Py_BuildValue("(i,i,s)",fail,wnd,fltlist);
}

static PyObject *py2c_numbtrop_quaddobl_retrieve_tropism
 ( PyObject *self, PyObject *args )
{
   int fail,dim,idx,wnd,k,nbc;
   char *fltlist;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&dim,&idx)) return NULL;   
   {
      double dir[4*dim];
      double err[4];
      double direrr[4*dim+4];

      fail = numbtrop_quaddobl_retrieve_tropism(dim,idx,&wnd,dir,err);

      fltlist = (char*)calloc(100*(dim+1), sizeof(char));
      for(k=0; k<4*dim; k++) direrr[k] = dir[k];
      direrr[4*dim] = err[0];
      direrr[4*dim+1] = err[1];
      direrr[4*dim+2] = err[2];
      direrr[4*dim+3] = err[3];
      nbc = dbllist2str(4*dim+4,direrr,fltlist);
   }
   return Py_BuildValue("(i,i,s)",fail,wnd,fltlist);
}

static PyObject *py2c_numbtrop_standard_clear
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = numbtrop_standard_clear();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_numbtrop_dobldobl_clear
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = numbtrop_dobldobl_clear();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_numbtrop_quaddobl_clear
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = numbtrop_quaddobl_clear();

   return Py_BuildValue("i",fail);
}

/* The wrapping of the functions with prototypes in witset.h starts here. */

static PyObject *py2c_embed_system
 ( PyObject *self, PyObject *args )
{
   int topdim,prc,vrblvl,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&topdim,&prc,&vrblvl)) return NULL;
   fail = embed_system(topdim,prc,vrblvl);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_embed_standard_system
 ( PyObject *self, PyObject *args )
{
   int topdim,vrblvl,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&topdim,&vrblvl)) return NULL;
   fail = embed_standard_system(topdim,vrblvl);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_embed_dobldobl_system
 ( PyObject *self, PyObject *args )
{
   int topdim,vrblvl,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&topdim,&vrblvl)) return NULL;
   fail = embed_dobldobl_system(topdim,vrblvl);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_embed_quaddobl_system
 ( PyObject *self, PyObject *args )
{
   int topdim,vrblvl,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&topdim,&vrblvl)) return NULL;
   fail = embed_quaddobl_system(topdim,vrblvl);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_embed_standard_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int topdim,vrblvl,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&topdim,&vrblvl)) return NULL;
   fail = embed_standard_Laurent_system(topdim,vrblvl);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_embed_dobldobl_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int topdim,vrblvl,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&topdim,&vrblvl)) return NULL;
   fail = embed_dobldobl_Laurent_system(topdim,vrblvl);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_embed_quaddobl_Laurent_system
 ( PyObject *self, PyObject *args )
{
   int topdim,vrblvl,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&topdim,&vrblvl)) return NULL;
   fail = embed_quaddobl_Laurent_system(topdim,vrblvl);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_swap_symbols_for_standard_witness_set
 ( PyObject *self, PyObject *args )
{
   int fail,nvr,dim;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&nvr,&dim)) return NULL;
   fail = swap_symbols_for_standard_witness_set(nvr,dim);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_swap_symbols_for_dobldobl_witness_set
 ( PyObject *self, PyObject *args )
{
   int fail,nvr,dim;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&nvr,&dim)) return NULL;
   fail = swap_symbols_for_dobldobl_witness_set(nvr,dim);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_swap_symbols_for_quaddobl_witness_set
 ( PyObject *self, PyObject *args )
{
   int fail,nvr,dim;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&nvr,&dim)) return NULL;
   fail = swap_symbols_for_quaddobl_witness_set(nvr,dim);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_swap_symbols_for_standard_Laurent_witness_set
 ( PyObject *self, PyObject *args )
{
   int fail,nvr,dim;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&nvr,&dim)) return NULL;
   fail = swap_symbols_for_standard_Laurent_witness_set(nvr,dim);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_swap_symbols_for_dobldobl_Laurent_witness_set
 ( PyObject *self, PyObject *args )
{
   int fail,nvr,dim;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&nvr,&dim)) return NULL;
   fail = swap_symbols_for_dobldobl_Laurent_witness_set(nvr,dim);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_swap_symbols_for_quaddobl_Laurent_witness_set
 ( PyObject *self, PyObject *args )
{
   int fail,nvr,dim;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&nvr,&dim)) return NULL;
   fail = swap_symbols_for_quaddobl_Laurent_witness_set(nvr,dim);

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

static PyObject *py2c_standard_Laurent_cascade_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = create_standard_Laurent_cascade_homotopy();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_dobldobl_Laurent_cascade_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = create_dobldobl_Laurent_cascade_homotopy();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_quaddobl_Laurent_cascade_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = create_quaddobl_Laurent_cascade_homotopy();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_set_standard_to_mute
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = set_standard_state_to_silent();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_set_dobldobl_to_mute
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = set_dobldobl_state_to_silent();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_set_quaddobl_to_mute
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = set_quaddobl_state_to_silent();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_set_standard_to_verbose
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = set_standard_state_to_verbose();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_set_dobldobl_to_verbose
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = set_dobldobl_state_to_verbose();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_set_quaddobl_to_verbose
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = set_quaddobl_state_to_verbose();

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

static PyObject *py2c_factor_standard_assign_labels
 ( PyObject *self, PyObject *args )
{
   int n,nbsols,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&n,&nbsols)) return NULL;
   fail = standard_assign_labels(n,nbsols);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_dobldobl_assign_labels
 ( PyObject *self, PyObject *args )
{
   int n,nbsols,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&n,&nbsols)) return NULL;
   fail = dobldobl_assign_labels(n,nbsols);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_quaddobl_assign_labels
 ( PyObject *self, PyObject *args )
{
   int n,nbsols,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&n,&nbsols)) return NULL;
   fail = quaddobl_assign_labels(n,nbsols);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_initialize_standard_sampler
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = initialize_standard_sampler(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_initialize_dobldobl_sampler
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = initialize_dobldobl_sampler(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_initialize_quaddobl_sampler
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = initialize_quaddobl_sampler(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_initialize_standard_Laurent_sampler
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = initialize_standard_Laurent_sampler(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_initialize_dobldobl_Laurent_sampler
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = initialize_dobldobl_Laurent_sampler(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_initialize_quaddobl_Laurent_sampler
 ( PyObject *self, PyObject *args )
{
   int k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&k)) return NULL;
   fail = initialize_quaddobl_Laurent_sampler(k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_initialize_standard_monodromy
 ( PyObject *self, PyObject *args )
{
   int n,d,k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&n,&d,&k)) return NULL;
   fail = initialize_standard_monodromy(n,d,k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_initialize_dobldobl_monodromy
 ( PyObject *self, PyObject *args )
{
   int n,d,k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&n,&d,&k)) return NULL;
   fail = initialize_dobldobl_monodromy(n,d,k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_initialize_quaddobl_monodromy
 ( PyObject *self, PyObject *args )
{
   int n,d,k,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&n,&d,&k)) return NULL;
   fail = initialize_quaddobl_monodromy(n,d,k);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_standard_trace_grid_diagnostics
 ( PyObject *self, PyObject *args )
{
   int fail;
   double err,dis;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = standard_trace_grid_diagnostics(&err,&dis);

   return Py_BuildValue("(d,d)",err,dis);
}

static PyObject *py2c_factor_dobldobl_trace_grid_diagnostics
 ( PyObject *self, PyObject *args )
{
   int fail;
   double err,dis;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = dobldobl_trace_grid_diagnostics(&err,&dis);

   return Py_BuildValue("(d,d)",err,dis);
}

static PyObject *py2c_factor_quaddobl_trace_grid_diagnostics
 ( PyObject *self, PyObject *args )
{
   int fail;
   double err,dis;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = quaddobl_trace_grid_diagnostics(&err,&dis);

   return Py_BuildValue("(d,d)",err,dis);
}

static PyObject *py2c_factor_store_standard_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = store_standard_solutions();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_store_dobldobl_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = store_dobldobl_solutions();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_store_quaddobl_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = store_quaddobl_solutions();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_restore_standard_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = restore_standard_solutions();

   return Py_BuildValue("i",fail);
}


static PyObject *py2c_factor_restore_dobldobl_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = restore_dobldobl_solutions();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_restore_quaddobl_solutions
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = restore_quaddobl_solutions();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_standard_track_paths
 ( PyObject *self, PyObject *args )
{
   int islaurent,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&islaurent)) return NULL;
   fail = standard_track_paths(islaurent);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_dobldobl_track_paths
 ( PyObject *self, PyObject *args )
{
   int islaurent,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&islaurent)) return NULL;
   fail = dobldobl_track_paths(islaurent);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_quaddobl_track_paths
 ( PyObject *self, PyObject *args )
{
   int islaurent,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&islaurent)) return NULL;
   fail = quaddobl_track_paths(islaurent);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_swap_standard_slices
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = swap_standard_slices();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_swap_dobldobl_slices
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = swap_dobldobl_slices();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_swap_quaddobl_slices
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = swap_quaddobl_slices();

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_new_standard_slices
 ( PyObject *self, PyObject *args )
{
   int k,n,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&k,&n)) return NULL;
   fail = new_standard_slices(k,n);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_new_dobldobl_slices
 ( PyObject *self, PyObject *args )
{
   int k,n,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&k,&n)) return NULL;
   fail = new_dobldobl_slices(k,n);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_new_quaddobl_slices
 ( PyObject *self, PyObject *args )
{
   int k,n,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&k,&n)) return NULL;
   fail = new_quaddobl_slices(k,n);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_set_standard_trace_slice
 ( PyObject *self, PyObject *args )
{
   int first,fail;
   double r[2];

   initialize();
   if(!PyArg_ParseTuple(args,"i",&first)) return NULL;

   r[1] = 0.0;
   if(first == 1)                  /* determine constant coefficient */
      r[0] = -1.0;
   else
      r[0] = +1.0;

   fail = assign_standard_coefficient_of_slice(0,0,r);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_set_dobldobl_trace_slice
 ( PyObject *self, PyObject *args )
{
   int first,fail;
   double r[4];

   initialize();
   if(!PyArg_ParseTuple(args,"i",&first)) return NULL;

   r[1] = 0.0; r[2] = 0.0; r[3] = 0.0;
   if(first == 1)                  /* determine constant coefficient */
      r[0] = -1.0;
   else
      r[0] = +1.0;

   fail = assign_dobldobl_coefficient_of_slice(0,0,r);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_set_quaddobl_trace_slice
 ( PyObject *self, PyObject *args )
{
   int first,fail;
   double r[8];

   initialize();
   if(!PyArg_ParseTuple(args,"i",&first)) return NULL;

   r[1] = 0.0; r[2] = 0.0; r[3] = 0.0; r[4] = 0.0;
   r[5] = 0.0; r[6] = 0.0; r[7] = 0.0;
   if(first == 1)                  /* determine constant coefficient */
      r[0] = -1.0;
   else
      r[0] = +1.0;

   fail = assign_quaddobl_coefficient_of_slice(0,0,r);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_store_standard_gammas
 ( PyObject *self, PyObject *args )
{
   int n,i,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&n)) return NULL;

   {
      double re_gamma[n];
      double im_gamma[n];
    
      for(i=0; i<n; i++)
         random_complex(&re_gamma[i],&im_gamma[i]);
      fail = store_standard_gamma(n,re_gamma,im_gamma);
   }

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_store_dobldobl_gammas
 ( PyObject *self, PyObject *args )
{
   int n,i,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&n)) return NULL;

   {
      double re_gamma[2*n];
      double im_gamma[2*n];
    
      for(i=0; i<n; i++)
         random_dobldobl_complex(&re_gamma[2*i],&im_gamma[2*i]);
      fail = store_dobldobl_gamma(n,re_gamma,im_gamma);
   }

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_store_quaddobl_gammas
 ( PyObject *self, PyObject *args )
{
   int n,i,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&n)) return NULL;

   {
      double re_gamma[4*n];
      double im_gamma[4*n];
    
      for(i=0; i<n; i++)
         random_quaddobl_complex(&re_gamma[4*i],&im_gamma[4*i]);
      fail = store_quaddobl_gamma(n,re_gamma,im_gamma);
   }

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_factor_permutation_after_standard_loop
 ( PyObject *self, PyObject *args )
{
   int d,fail,nb;
   char *result;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&d)) return NULL;
   {
      int i,permutation[d];
      char s[d*10];

      fail = permutation_after_standard_loop(d,permutation);
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

static PyObject *py2c_factor_permutation_after_dobldobl_loop
 ( PyObject *self, PyObject *args )
{
   int d,fail,nb;
   char *result;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&d)) return NULL;
   {
      int i,permutation[d];
      char s[d*10];

      fail = permutation_after_dobldobl_loop(d,permutation);
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

static PyObject *py2c_factor_permutation_after_quaddobl_loop
 ( PyObject *self, PyObject *args )
{
   int d,fail,nb;
   char *result;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&d)) return NULL;
   {
      int i,permutation[d];
      char s[d*10];

      fail = permutation_after_quaddobl_loop(d,permutation);
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

static PyObject *py2c_factor_update_standard_decomposition
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

      fail = update_standard_decomposition(d,perm,nf,&done);
   /* printf("number of factors : %d -> %d\n",nf[0],nf[1]); */
   }
   return Py_BuildValue("i",done);
}

static PyObject *py2c_factor_update_dobldobl_decomposition
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

      fail = update_dobldobl_decomposition(d,perm,nf,&done);
   /* printf("number of factors : %d -> %d\n",nf[0],nf[1]); */
   }
   return Py_BuildValue("i",done);
}

static PyObject *py2c_factor_update_quaddobl_decomposition
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

      fail = update_quaddobl_decomposition(d,perm,nf,&done);
   /* printf("number of factors : %d -> %d\n",nf[0],nf[1]); */
   }
   return Py_BuildValue("i",done);
}

static PyObject *py2c_factor_number_of_standard_components
 ( PyObject *self, PyObject *args )
{
   int fail,nf;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;

   fail = number_of_standard_factors(&nf);

   return Py_BuildValue("i",nf);
}

static PyObject *py2c_factor_number_of_dobldobl_components
 ( PyObject *self, PyObject *args )
{
   int fail,nf;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;

   fail = number_of_dobldobl_factors(&nf);

   return Py_BuildValue("i",nf);
}

static PyObject *py2c_factor_number_of_quaddobl_components
 ( PyObject *self, PyObject *args )
{
   int fail,nf;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;

   fail = number_of_quaddobl_factors(&nf);

   return Py_BuildValue("i",nf);
}

static PyObject *py2c_factor_witness_points_of_standard_component
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

      fail = witness_points_of_standard_factor(k,&deg,w);

      nb = list2str(deg,w,s);
      result = (char*)calloc(nb,sizeof(char));
      for(i=0; i<nb; i++) result[i] = s[i];
   }
   return Py_BuildValue("s",result);
}

static PyObject *py2c_factor_witness_points_of_dobldobl_component
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

      fail = witness_points_of_dobldobl_factor(k,&deg,w);

      nb = list2str(deg,w,s);
      result = (char*)calloc(nb,sizeof(char));
      for(i=0; i<nb; i++) result[i] = s[i];
   }
   return Py_BuildValue("s",result);
}

static PyObject *py2c_factor_witness_points_of_quaddobl_component
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

      fail = witness_points_of_quaddobl_factor(k,&deg,w);

      nb = list2str(deg,w,s);
      result = (char*)calloc(nb,sizeof(char));
      for(i=0; i<nb; i++) result[i] = s[i];
   }
   return Py_BuildValue("s",result);
}

static PyObject *py2c_factor_standard_trace_sum_difference
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

      fail = standard_trace_sum_difference(nb,witset,&tsd);
   /* printf("trace sum difference : %.3e\n",tsd); */
   }
   return Py_BuildValue("d",tsd);
}

static PyObject *py2c_factor_dobldobl_trace_sum_difference
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

      fail = dobldobl_trace_sum_difference(nb,witset,&tsd);
   /* printf("trace sum difference : %.3e\n",tsd); */
   }
   return Py_BuildValue("d",tsd);
}

static PyObject *py2c_factor_quaddobl_trace_sum_difference
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

      fail = quaddobl_trace_sum_difference(nb,witset,&tsd);
   /* printf("trace sum difference : %.3e\n",tsd); */
   }
   return Py_BuildValue("d",tsd);
}

static PyObject *py2c_witset_standard_membertest
 ( PyObject *self, PyObject *args )
{
   int v,n,d,m,fail,onp,ins,nbtasks=0;
   double r,h;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiiidds",&v,&nbtasks,&n,&d,&m,&r,&h,&p))
      return NULL;
   {
      const int dim = 2*n;
      double tpt[dim];

      str2dbllist(dim,p,tpt);
      fail = standard_homotopy_membership_test(v,n,d,r,h,tpt,&onp,&ins,nbtasks);
   }
   return Py_BuildValue("(i,i,i)",fail,onp,ins);
}

static PyObject *py2c_witset_dobldobl_membertest
 ( PyObject *self, PyObject *args )
{
   int v,n,d,m,fail,onp,ins,nbtasks=0;
   double r,h;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiiidds",&v,&nbtasks,&n,&d,&m,&r,&h,&p))
      return NULL;
   {
      const int dim = 4*n;
      double tpt[dim];

      str2dbllist(dim,p,tpt);
      fail = dobldobl_homotopy_membership_test(v,n,d,r,h,tpt,&onp,&ins,nbtasks);
   }
   return Py_BuildValue("(i,i,i)",fail,onp,ins);
}

static PyObject *py2c_witset_quaddobl_membertest
 ( PyObject *self, PyObject *args )
{
   int v,n,d,m,fail,onp,ins,nbtasks=0;
   double r,h;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiiidds",&v,&nbtasks,&n,&d,&m,&r,&h,&p))
      return NULL;
   {
      const int dim = 8*n;
      double tpt[dim];

      str2dbllist(dim,p,tpt);
      fail = quaddobl_homotopy_membership_test(v,n,d,r,h,tpt,&onp,&ins,nbtasks);
   }
   return Py_BuildValue("(i,i,i)",fail,onp,ins);
}


static PyObject *py2c_witset_standard_Laurent_membertest
 ( PyObject *self, PyObject *args )
{
   int v,n,d,m,fail,onp,ins,nbtasks=0;
   double r,h;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiiidds",&v,&nbtasks,&n,&d,&m,&r,&h,&p))
      return NULL;
   {
      const int dim = 2*n;
      double tpt[dim];

      str2dbllist(dim,p,tpt);
      fail = standard_Laurent_homotopy_membership_test
                (v,n,d,r,h,tpt,&onp,&ins,nbtasks);
   }
   return Py_BuildValue("(i,i,i)",fail,onp,ins);
}

static PyObject *py2c_witset_dobldobl_Laurent_membertest
 ( PyObject *self, PyObject *args )
{
   int v,n,d,m,fail,onp,ins,nbtasks=0;
   double r,h;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiiidds",&v,&nbtasks,&n,&d,&m,&r,&h,&p))
      return NULL;
   {
      const int dim = 4*n;
      double tpt[dim];

      str2dbllist(dim,p,tpt);
      fail = dobldobl_Laurent_homotopy_membership_test
                (v,n,d,r,h,tpt,&onp,&ins,nbtasks);
   }
   return Py_BuildValue("(i,i,i)",fail,onp,ins);
}

static PyObject *py2c_witset_quaddobl_Laurent_membertest
 ( PyObject *self, PyObject *args )
{
   int v,n,d,m,fail,onp,ins,nbtasks=0;
   double r,h;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiiidds",&v,&nbtasks,&n,&d,&m,&r,&h,&p))
      return NULL;
   {
      const int dim = 8*n;
      double tpt[dim];

      str2dbllist(dim,p,tpt);
      fail = quaddobl_Laurent_homotopy_membership_test
                (v,n,d,r,h,tpt,&onp,&ins,nbtasks);
   }
   return Py_BuildValue("(i,i,i)",fail,onp,ins);
}

static PyObject *py2c_witset_standard_ismember
 ( PyObject *self, PyObject *args )
{
   int v,n,d,m,fail,onp,ins,nbtasks=0;
   double r,h;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiiidds",&v,&nbtasks,&n,&d,&m,&r,&h,&p))
      return NULL;
   fail = standard_homotopy_ismember(v,n,d,m,p,r,h,&onp,&ins,nbtasks);

   return Py_BuildValue("(i,i,i)",fail,onp,ins);
}

static PyObject *py2c_witset_dobldobl_ismember
 ( PyObject *self, PyObject *args )
{
   int v,n,d,m,fail,onp,ins,nbtasks=0;
   double r,h;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiiidds",&v,&nbtasks,&n,&d,&m,&r,&h,&p))
      return NULL;
   fail = dobldobl_homotopy_ismember(v,n,d,m,p,r,h,&onp,&ins,nbtasks);

   return Py_BuildValue("(i,i,i)",fail,onp,ins);
}

static PyObject *py2c_witset_quaddobl_ismember
 ( PyObject *self, PyObject *args )
{
   int v,n,d,m,fail,onp,ins,nbtasks=0;
   double r,h;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiiidds",&v,&nbtasks,&n,&d,&m,&r,&h,&p))
      return NULL;
   fail = quaddobl_homotopy_ismember(v,n,d,m,p,r,h,&onp,&ins,nbtasks);

   return Py_BuildValue("(i,i,i)",fail,onp,ins);
}

static PyObject *py2c_witset_standard_Laurent_ismember
 ( PyObject *self, PyObject *args )
{
   int v,n,d,m,fail,onp,ins,nbtasks=0;
   double r,h;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiiidds",&v,&nbtasks,&n,&d,&m,&r,&h,&p))
      return NULL;
   fail = standard_Laurent_homotopy_ismember(v,n,d,m,p,r,h,&onp,&ins,nbtasks);

   return Py_BuildValue("(i,i,i)",fail,onp,ins);
}

static PyObject *py2c_witset_dobldobl_Laurent_ismember
 ( PyObject *self, PyObject *args )
{
   int v,n,d,m,fail,onp,ins,nbtasks=0;
   double r,h;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiiidds",&v,&nbtasks,&n,&d,&m,&r,&h,&p))
      return NULL;
   fail = dobldobl_Laurent_homotopy_ismember(v,n,d,m,p,r,h,&onp,&ins,nbtasks);

   return Py_BuildValue("(i,i,i)",fail,onp,ins);
}

static PyObject *py2c_witset_quaddobl_Laurent_ismember
 ( PyObject *self, PyObject *args )
{
   int v,n,d,m,fail,onp,ins,nbtasks=0;
   double r,h;
   char *p;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiiidds",&v,&nbtasks,&n,&d,&m,&r,&h,&p))
      return NULL;
   fail = quaddobl_Laurent_homotopy_ismember(v,n,d,m,p,r,h,&onp,&ins,nbtasks);

   return Py_BuildValue("(i,i,i)",fail,onp,ins);
}

static PyObject *py2c_standard_witset_of_hypersurface
 ( PyObject *self, PyObject *args )
{
   int fail,nv,nc;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iis",&nv,&nc,&p)) return NULL;
   fail = standard_witset_of_hypersurface(nv,nc,p);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_dobldobl_witset_of_hypersurface
 ( PyObject *self, PyObject *args )
{
   int fail,nv,nc;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iis",&nv,&nc,&p)) return NULL;
   fail = dobldobl_witset_of_hypersurface(nv,nc,p);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_quaddobl_witset_of_hypersurface
 ( PyObject *self, PyObject *args )
{
   int fail,nv,nc;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iis",&nv,&nc,&p)) return NULL;
   fail = quaddobl_witset_of_hypersurface(nv,nc,p);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_standard_witset_of_Laurent_hypersurface
 ( PyObject *self, PyObject *args )
{
   int fail,nv,nc;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iis",&nv,&nc,&p)) return NULL;
   fail = standard_witset_of_Laurent_hypersurface(nv,nc,p);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_dobldobl_witset_of_Laurent_hypersurface
 ( PyObject *self, PyObject *args )
{
   int fail,nv,nc;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iis",&nv,&nc,&p)) return NULL;
   fail = dobldobl_witset_of_Laurent_hypersurface(nv,nc,p);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_quaddobl_witset_of_Laurent_hypersurface
 ( PyObject *self, PyObject *args )
{
   int fail,nv,nc;
   char *p;   
                 
   initialize();
   if(!PyArg_ParseTuple(args,"iis",&nv,&nc,&p)) return NULL;
   fail = quaddobl_witset_of_Laurent_hypersurface(nv,nc,p);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_standard_diagonal_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,a,b;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&a,&b)) return NULL;
   fail = standard_diagonal_homotopy(a,b);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_dobldobl_diagonal_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,a,b;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&a,&b)) return NULL;
   fail = dobldobl_diagonal_homotopy(a,b);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_quaddobl_diagonal_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,a,b;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&a,&b)) return NULL;
   fail = quaddobl_diagonal_homotopy(a,b);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_standard_diagonal_cascade_solutions
 ( PyObject *self, PyObject *args )
{
   int fail,a,b;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&a,&b)) return NULL;
   fail = standard_diagonal_cascade_solutions(a,b);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_dobldobl_diagonal_cascade_solutions
 ( PyObject *self, PyObject *args )
{
   int fail,a,b;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&a,&b)) return NULL;
   fail = dobldobl_diagonal_cascade_solutions(a,b);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_quaddobl_diagonal_cascade_solutions
 ( PyObject *self, PyObject *args )
{
   int fail,a,b;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&a,&b)) return NULL;
   fail = quaddobl_diagonal_cascade_solutions(a,b);

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

static PyObject *py2c_diagonal_symbols_doubler 
 ( PyObject *self, PyObject *args )
{
   int fail,n,d,nc;
   char *s;

   initialize();
   if(!PyArg_ParseTuple(args,"iiis",&n,&d,&nc,&s)) return NULL;
   fail = diagonal_symbols_doubler(n,d,nc,s);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_standard_collapse_diagonal
 ( PyObject *self, PyObject *args )
{
   int fail,k,d;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&k,&d)) return NULL;
   fail = standard_collapse_diagonal(k,d);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_dobldobl_collapse_diagonal
 ( PyObject *self, PyObject *args )
{
   int fail,k,d;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&k,&d)) return NULL;
   fail = dobldobl_collapse_diagonal(k,d);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_quaddobl_collapse_diagonal
 ( PyObject *self, PyObject *args )
{
   int fail,k,d;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&k,&d)) return NULL;
   fail = quaddobl_collapse_diagonal(k,d);

   return Py_BuildValue("i",fail);
}

/* The wrapping of functions with prototypes in witsols.h starts here. */

static PyObject *py2c_standard_polysys_solve
 ( PyObject *self, PyObject *args )
{
   int fail,nbtasks,topdim,filter,factor,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiii",
      &nbtasks,&topdim,&filter,&factor,&verbose)) return NULL;
   fail = standard_polysys_solve(nbtasks,topdim,filter,factor,verbose);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_standard_laursys_solve
 ( PyObject *self, PyObject *args )
{
   int fail,nbtasks,topdim,filter,factor,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiii",
      &nbtasks,&topdim,&filter,&factor,&verbose)) return NULL;
   fail = standard_laursys_solve(nbtasks,topdim,filter,factor,verbose);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_dobldobl_polysys_solve
 ( PyObject *self, PyObject *args )
{
   int fail,nbtasks,topdim,filter,factor,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiii",
      &nbtasks,&topdim,&filter,&factor,&verbose)) return NULL;
   fail = dobldobl_polysys_solve(nbtasks,topdim,filter,factor,verbose);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_dobldobl_laursys_solve
 ( PyObject *self, PyObject *args )
{
   int fail,nbtasks,topdim,filter,factor,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiii",
      &nbtasks,&topdim,&filter,&factor,&verbose)) return NULL;
   fail = dobldobl_laursys_solve(nbtasks,topdim,filter,factor,verbose);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_quaddobl_polysys_solve
 ( PyObject *self, PyObject *args )
{
   int fail,nbtasks,topdim,filter,factor,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiii",
      &nbtasks,&topdim,&filter,&factor,&verbose)) return NULL;
   fail = quaddobl_polysys_solve(nbtasks,topdim,filter,factor,verbose);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_quaddobl_laursys_solve
 ( PyObject *self, PyObject *args )
{
   int fail,nbtasks,topdim,filter,factor,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiii",
      &nbtasks,&topdim,&filter,&factor,&verbose)) return NULL;
   fail = quaddobl_laursys_solve(nbtasks,topdim,filter,factor,verbose);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_standard_polysys_witset
 ( PyObject *self, PyObject *args )
{
   int fail,dim;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;

   fail = copy_standard_polysys_witset(dim);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_standard_laursys_witset
 ( PyObject *self, PyObject *args )
{
   int fail,dim;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;

   fail = copy_standard_laursys_witset(dim);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_dobldobl_polysys_witset
 ( PyObject *self, PyObject *args )
{
   int fail,dim;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;

   fail = copy_dobldobl_polysys_witset(dim);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_dobldobl_laursys_witset
 ( PyObject *self, PyObject *args )
{
   int fail,dim;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;

   fail = copy_dobldobl_laursys_witset(dim);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_quaddobl_polysys_witset
 ( PyObject *self, PyObject *args )
{
   int fail,dim;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;

   fail = copy_quaddobl_polysys_witset(dim);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_copy_quaddobl_laursys_witset
 ( PyObject *self, PyObject *args )
{
   int fail,dim;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;

   fail = copy_quaddobl_laursys_witset(dim);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_standard_witsols
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = clear_standard_witsols();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_dobldobl_witsols
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = clear_dobldobl_witsols();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_clear_quaddobl_witsols
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = clear_quaddobl_witsols();
              
   return Py_BuildValue("i",fail);
}

/* The wrapping of functions with prototypes in schubert.h starts here. */

static PyObject *py2c_schubert_pieri_count
 ( PyObject *self, PyObject *args )
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

static PyObject *py2c_schubert_standard_littlewood_richardson_homotopies
 ( PyObject *self, PyObject *args )
{
   int i,n,k,nbc,nc,fail,r,vrb,vrf,mrp,sqr,szn;
   char *cond;
   char *name;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiisiiiiis",
       &n,&k,&nbc,&nc,&cond,&vrb,&vrf,&mrp,&sqr,&szn,&name))
      return NULL;
/*
   printf("the value of szn : %d\n", szn);
   printf("the value of vrb : %d\n", vrb);
   printf("the value of vrf : %d\n", vrf);
   printf("the value of mrp : %d\n", mrp);
   printf("the value of sqr : %d\n", sqr);
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
      fail = standard_Littlewood_Richardson_homotopies
               (n,k,nbc,cds,vrb,vrf,mrp,sqr,szn,name,&r,fg);
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

static PyObject *py2c_schubert_dobldobl_littlewood_richardson_homotopies
 ( PyObject *self, PyObject *args )
{
   int i,n,k,nbc,nc,fail,r,vrb,vrf,mrp,sqr,szn;
   char *cond;
   char *name;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiisiiiiis",
       &n,&k,&nbc,&nc,&cond,&vrb,&vrf,&mrp,&sqr,&szn,&name))
      return NULL;
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
      const int fgsize = 4*(nbc-2)*n*n;
      double fg[fgsize];
      char stfg[fgsize*24+2];
      fail = dobldobl_Littlewood_Richardson_homotopies
               (n,k,nbc,cds,vrb,vrf,mrp,sqr,szn,name,&r,fg);
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

static PyObject *py2c_schubert_quaddobl_littlewood_richardson_homotopies
 ( PyObject *self, PyObject *args )
{
   int i,n,k,nbc,nc,fail,r,vrb,vrf,mrp,sqr,szn;
   char *cond;
   char *name;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiisiiiiis",
       &n,&k,&nbc,&nc,&cond,&vrb,&vrf,&mrp,&sqr,&szn,&name))
      return NULL;
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
      const int fgsize = 8*(nbc-2)*n*n;
      double fg[fgsize];
      char stfg[fgsize*24+2];
      fail = quaddobl_Littlewood_Richardson_homotopies
               (n,k,nbc,cds,vrb,vrf,mrp,sqr,szn,name,&r,fg);
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

static PyObject *py2c_schubert_pieri_system
 ( PyObject *self, PyObject *args )
{
   int fail,m,p,q,nc,r;
   char *A;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiisi",&m,&p,&q,&nc,&A,&r)) return NULL;
   fail = Pieri_polynomial_system(m,p,q,nc,A,r);

   return Py_BuildValue("i",fail);
}

/* The wrapping of the functions with prototypes in mapcon.h starts here. */

static PyObject *py2c_mapcon_solve_system
 ( PyObject *self, PyObject *args )
{
   int fail,puretopdim;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&puretopdim)) return NULL;   
   fail = mapcon_solve_system(puretopdim);
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_mapcon_write_maps
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = mapcon_write_maps();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_mapcon_clear_maps
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = mapcon_clear_maps();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_mapcon_top_dimension
 ( PyObject *self, PyObject *args )
{
   int fail,topdim;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = mapcon_top_dimension(&topdim);
              
   return Py_BuildValue("i",topdim);
}

static PyObject *py2c_mapcon_number_of_maps
 ( PyObject *self, PyObject *args )
{
   int fail,dim,nbmaps;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&dim)) return NULL;   
   fail = mapcon_number_of_maps(dim,&nbmaps);
              
   return Py_BuildValue("i",nbmaps);
}

static PyObject *py2c_mapcon_degree_of_map
 ( PyObject *self, PyObject *args )
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

/* The wrapping of functions with prototypes in series.h starts below. */

static PyObject *py2c_standard_Newton_series ( PyObject *self, PyObject *args )
{
   int idx,maxdeg,nbr,vrb,fail;
   initialize();
   if(!PyArg_ParseTuple(args,"iiii",&idx,&maxdeg,&nbr,&vrb)) return NULL;   
   fail = standard_Newton_series(idx,maxdeg,nbr,vrb);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_dobldobl_Newton_series ( PyObject *self, PyObject *args )
{
   int idx,maxdeg,nbr,vrb,fail;
   initialize();
   if(!PyArg_ParseTuple(args,"iiii",&idx,&maxdeg,&nbr,&vrb)) return NULL;   
   fail = dobldobl_Newton_series(idx,maxdeg,nbr,vrb);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_quaddobl_Newton_series ( PyObject *self, PyObject *args )
{
   int idx,maxdeg,nbr,vrb,fail;
   initialize();
   if(!PyArg_ParseTuple(args,"iiii",&idx,&maxdeg,&nbr,&vrb)) return NULL;   
   fail = quaddobl_Newton_series(idx,maxdeg,nbr,vrb);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_standard_Newton_power_series
 ( PyObject *self, PyObject *args )
{
   int idx,maxdeg,nbr,vrb,fail;
   initialize();
   if(!PyArg_ParseTuple(args,"iiii",&idx,&maxdeg,&nbr,&vrb)) return NULL;   
   fail = standard_Newton_power_series(idx,maxdeg,nbr,vrb);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_dobldobl_Newton_power_series
 ( PyObject *self, PyObject *args )
{
   int idx,maxdeg,nbr,vrb,fail;
   initialize();
   if(!PyArg_ParseTuple(args,"iiii",&idx,&maxdeg,&nbr,&vrb)) return NULL;   
   fail = dobldobl_Newton_power_series(idx,maxdeg,nbr,vrb);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_quaddobl_Newton_power_series
 ( PyObject *self, PyObject *args )
{
   int idx,maxdeg,nbr,vrb,fail;
   initialize();
   if(!PyArg_ParseTuple(args,"iii",&idx,&maxdeg,&nbr,&vrb)) return NULL;   
   fail = quaddobl_Newton_power_series(idx,maxdeg,nbr,vrb);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_standard_Pade_approximant
 ( PyObject *self, PyObject *args )
{
   int idx,numdeg,dendeg,nbr,vrb,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiii",&idx,&numdeg,&dendeg,&nbr,&vrb))
      return NULL;   

   fail = standard_Pade_approximant(idx,numdeg,dendeg,nbr,vrb);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_dobldobl_Pade_approximant
 ( PyObject *self, PyObject *args )
{
   int idx,numdeg,dendeg,nbr,vrb,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiii",&idx,&numdeg,&dendeg,&nbr,&vrb))
      return NULL;   

   fail = dobldobl_Pade_approximant(idx,numdeg,dendeg,nbr,vrb);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_quaddobl_Pade_approximant
 ( PyObject *self, PyObject *args )
{
   int idx,numdeg,dendeg,nbr,vrb,fail;

   initialize();
   if(!PyArg_ParseTuple(args,"iiiii",&idx,&numdeg,&dendeg,&nbr,&vrb))
      return NULL;   

   fail = quaddobl_Pade_approximant(idx,numdeg,dendeg,nbr,vrb);

   return Py_BuildValue("i",fail);
}

/* The wrapping of Pade continuation starts here. */

static PyObject *py2c_padcon_set_default_parameters
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = padcon_set_default_parameters();
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_clear_parameters
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = padcon_clear_parameters();
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_get_homotopy_continuation_parameter
 ( PyObject *self, PyObject *args )
{
   int fail,idx;
   double val[2];

   initialize();
   if(!PyArg_ParseTuple(args,"i",&idx)) return NULL;
   fail = padcon_get_homotopy_continuation_parameter(idx,val);

   if(idx == 1)
      return Py_BuildValue("(d,d)", val[0], val[1]);
   else if((idx == 2) || (idx == 3) || (idx == 11) || (idx == 12))
   {
      int parval = (int) val[0];
      return Py_BuildValue("i",parval);
   }
   else
      return Py_BuildValue("d",val[0]);
}

static PyObject *py2c_padcon_set_homotopy_continuation_gamma
 ( PyObject *self, PyObject *args )
{
   int fail;
   double val[2];

   initialize();
   if(!PyArg_ParseTuple(args,"dd",&val[0],&val[1])) return NULL;
   fail = padcon_set_homotopy_continuation_parameter(1,val);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_set_homotopy_continuation_parameter
 ( PyObject *self, PyObject *args )
{
   int fail,idx;
   double val;

   initialize();
   if(!PyArg_ParseTuple(args,"id",&idx,&val)) return NULL;
   fail = padcon_set_homotopy_continuation_parameter(idx,&val);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_reset_homotopy_continuation_parameters
 ( PyObject *self, PyObject *args )
{
   int fail,prc;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&prc)) return NULL;
   fail = padcon_reset_homotopy_continuation_parameters(prc);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_standard_track
 ( PyObject *self, PyObject *args )
{
   int nbc,fail,localfile,verbose,mhom,nvr,*idz;
   char *name;
   char *stridz;

   initialize();
   if(!PyArg_ParseTuple(args,"isiiiis",
      &nbc,&name,&localfile,&verbose,&mhom,&nvr,&stridz)) return NULL;

   if(mhom == 1)
      padcon_standard_projective_transformation();
   else if(mhom > 1)
   {
      int ic = itemcount(stridz);
      idz = (int*)calloc(ic,sizeof(int));
      str2intlist(ic,stridz,idz);
      if(verbose > 0)
      {
         printf("mhom : %d", mhom);
         printf(" partition :");
         for(int k=0; k<nvr; k++) printf(" %d", idz[k]);
         printf("\n");
      }
      padcon_standard_multi_projective_transformation(nvr,mhom,idz);
   }
   fail = padcon_standard_track(nbc,name,localfile,verbose,mhom,nvr,idz);
   if(mhom == 1)
      fail = solcon_standard_one_affinization();
   else if(mhom > 1)
      fail = solcon_standard_multi_affinization(nvr,mhom,idz);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_dobldobl_track
 ( PyObject *self, PyObject *args )
{
   int nbc,fail,localfile,verbose,mhom,nvr,*idz;
   char *name;
   char *stridz;

   initialize();
   if(!PyArg_ParseTuple(args,"isiiiis",
      &nbc,&name,&localfile,&verbose,&mhom,&nvr,&stridz)) return NULL;

   if(mhom == 1)
      padcon_dobldobl_projective_transformation();
   else if(mhom > 1)
   {
      int ic = itemcount(stridz);
      idz = (int*)calloc(ic,sizeof(int));
      str2intlist(ic,stridz,idz);
      if(verbose > 0)
      {
         printf("mhom : %d", mhom);
         printf(" partition :");
         for(int k=0; k<nvr; k++) printf(" %d", idz[k]);
         printf("\n");
      }
      padcon_dobldobl_multi_projective_transformation(nvr,mhom,idz);
   }
   fail = padcon_dobldobl_track(nbc,name,localfile,verbose,mhom,nvr,idz);
   if(mhom == 1)
      fail = solcon_dobldobl_one_affinization();
   else if(mhom > 1)
      fail = solcon_dobldobl_multi_affinization(nvr,mhom,idz);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_quaddobl_track
 ( PyObject *self, PyObject *args )
{
   int nbc,fail,localfile,verbose,mhom,nvr,*idz;
   char *name;
   char *stridz;

   initialize();
   if(!PyArg_ParseTuple(args,"isiiiis",
      &nbc,&name,&localfile,&verbose,&mhom,&nvr,&stridz)) return NULL;

   if(mhom == 1)
      padcon_quaddobl_projective_transformation();
   else if(mhom > 1)
   {
      int ic = itemcount(stridz);
      idz = (int*)calloc(ic,sizeof(int));
      str2intlist(ic,stridz,idz);
      if(verbose > 0)
      {
         printf("mhom : %d", mhom);
         printf(" partition :");
         for(int k=0; k<nvr; k++) printf(" %d", idz[k]);
         printf("\n");
      }
      padcon_quaddobl_multi_projective_transformation(nvr,mhom,idz);
   }
   fail = padcon_quaddobl_track(nbc,name,localfile,verbose,mhom,nvr,idz);
   if(mhom == 1)
      fail = solcon_quaddobl_one_affinization();
   else if(mhom > 1)
      fail = solcon_quaddobl_multi_affinization(nvr,mhom,idz);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_standard_initialize_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,verbose,homogeneous;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&verbose,&homogeneous)) return NULL;
   fail = padcon_standard_initialize_homotopy(verbose,homogeneous);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_dobldobl_initialize_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,verbose,homogeneous;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&verbose,&homogeneous)) return NULL;
   fail = padcon_dobldobl_initialize_homotopy(verbose,homogeneous);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_quaddobl_initialize_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,verbose,homogeneous;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&verbose,&homogeneous)) return NULL;
   fail = padcon_quaddobl_initialize_homotopy(verbose,homogeneous);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_standard_initialize_parameter_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,index,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&index,&verbose)) return NULL;
   fail = padcon_standard_initialize_parameter_homotopy(index,verbose);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_dobldobl_initialize_parameter_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,index,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&index,&verbose)) return NULL;
   fail = padcon_dobldobl_initialize_parameter_homotopy(index,verbose);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_quaddobl_initialize_parameter_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,index,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&index,&verbose)) return NULL;
   fail = padcon_quaddobl_initialize_parameter_homotopy(index,verbose);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_initialize_standard_solution
 ( PyObject *self, PyObject *args )
{
   int fail,index,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&index,&verbose)) return NULL;
   fail = padcon_initialize_standard_solution(index,verbose);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_initialize_dobldobl_solution
 ( PyObject *self, PyObject *args )
{
   int fail,index,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&index,&verbose)) return NULL;
   fail = padcon_initialize_dobldobl_solution(index,verbose);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_initialize_quaddobl_solution
 ( PyObject *self, PyObject *args )
{
   int fail,index,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&index,&verbose)) return NULL;
   fail = padcon_initialize_quaddobl_solution(index,verbose);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_standard_predict_correct
 ( PyObject *self, PyObject *args )
{
   int fail,precorfail,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&verbose)) return NULL;
   fail = padcon_standard_predict_correct(&precorfail,verbose);
   return Py_BuildValue("i",precorfail);
}

static PyObject *py2c_padcon_dobldobl_predict_correct
 ( PyObject *self, PyObject *args )
{
   int fail,precorfail,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&verbose)) return NULL;
   fail = padcon_dobldobl_predict_correct(&precorfail,verbose);
   return Py_BuildValue("i",precorfail);
}

static PyObject *py2c_padcon_quaddobl_predict_correct
 ( PyObject *self, PyObject *args )
{
   int fail,precorfail,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&verbose)) return NULL;
   fail = padcon_quaddobl_predict_correct(&precorfail,verbose);
   return Py_BuildValue("i",precorfail);
}

static PyObject *py2c_padcon_get_standard_solution
 ( PyObject *self, PyObject *args )
{
   int fail,index,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&index,&verbose)) return NULL;
   fail = padcon_get_standard_solution(index,verbose);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_get_dobldobl_solution
 ( PyObject *self, PyObject *args )
{
   int fail,index,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&index,&verbose)) return NULL;
   fail = padcon_get_dobldobl_solution(index,verbose);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_get_quaddobl_solution
 ( PyObject *self, PyObject *args )
{
   int fail,index,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"ii",&index,&verbose)) return NULL;
   fail = padcon_get_quaddobl_solution(index,verbose);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_standard_pole_radius
 ( PyObject *self, PyObject *args )
{
   int fail;
   double frp;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_standard_pole_radius(&frp);

   return Py_BuildValue("d",frp);
}

static PyObject *py2c_padcon_dobldobl_pole_radius
 ( PyObject *self, PyObject *args )
{
   int fail;
   double frp;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_dobldobl_pole_radius(&frp);

   return Py_BuildValue("d",frp);
}

static PyObject *py2c_padcon_quaddobl_pole_radius
 ( PyObject *self, PyObject *args )
{
   int fail;
   double frp;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_quaddobl_pole_radius(&frp);

   return Py_BuildValue("d",frp);
}

static PyObject *py2c_padcon_standard_closest_pole
 ( PyObject *self, PyObject *args )
{
   int fail;
   double repole,impole;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_standard_closest_pole(&repole,&impole);

   return Py_BuildValue("(d,d)",repole,impole);
}

static PyObject *py2c_padcon_dobldobl_closest_pole
 ( PyObject *self, PyObject *args )
{
   int fail;
   double repole,impole;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_dobldobl_closest_pole(&repole,&impole);

   return Py_BuildValue("(d,d)",repole,impole);
}

static PyObject *py2c_padcon_quaddobl_closest_pole
 ( PyObject *self, PyObject *args )
{
   int fail;
   double repole,impole;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_quaddobl_closest_pole(&repole,&impole);

   return Py_BuildValue("(d,d)",repole,impole);
}

static PyObject *py2c_padcon_standard_t_value
 ( PyObject *self, PyObject *args )
{
   int fail;
   double tval;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_standard_t_value(&tval);

   return Py_BuildValue("d",tval);
}

static PyObject *py2c_padcon_dobldobl_t_value
 ( PyObject *self, PyObject *args )
{
   int fail;
   double tval;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_dobldobl_t_value(&tval);

   return Py_BuildValue("d",tval);
}

static PyObject *py2c_padcon_quaddobl_t_value
 ( PyObject *self, PyObject *args )
{
   int fail;
   double tval;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_quaddobl_t_value(&tval);

   return Py_BuildValue("d",tval);
}

static PyObject *py2c_padcon_standard_step_size
 ( PyObject *self, PyObject *args )
{
   int fail;
   double step;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_standard_step_size(&step);

   return Py_BuildValue("d",step);
}

static PyObject *py2c_padcon_dobldobl_step_size
 ( PyObject *self, PyObject *args )
{
   int fail;
   double step;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_dobldobl_step_size(&step);

   return Py_BuildValue("d",step);
}

static PyObject *py2c_padcon_quaddobl_step_size
 ( PyObject *self, PyObject *args )
{
   int fail;
   double step;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_quaddobl_step_size(&step);

   return Py_BuildValue("d",step);
}

static PyObject *py2c_padcon_standard_series_step
 ( PyObject *self, PyObject *args )
{
   int fail;
   double step;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_standard_series_step(&step);

   return Py_BuildValue("d",step);
}

static PyObject *py2c_padcon_dobldobl_series_step
 ( PyObject *self, PyObject *args )
{
   int fail;
   double step;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_dobldobl_series_step(&step);

   return Py_BuildValue("d",step);
}

static PyObject *py2c_padcon_quaddobl_series_step
 ( PyObject *self, PyObject *args )
{
   int fail;
   double step;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_quaddobl_series_step(&step);

   return Py_BuildValue("d",step);
}

static PyObject *py2c_padcon_standard_pole_step
 ( PyObject *self, PyObject *args )
{
   int fail;
   double step;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_standard_pole_step(&step);

   return Py_BuildValue("d",step);
}

static PyObject *py2c_padcon_dobldobl_pole_step
 ( PyObject *self, PyObject *args )
{
   int fail;
   double step;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_dobldobl_pole_step(&step);

   return Py_BuildValue("d",step);
}

static PyObject *py2c_padcon_quaddobl_pole_step
 ( PyObject *self, PyObject *args )
{
   int fail;
   double step;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_quaddobl_pole_step(&step);

   return Py_BuildValue("d",step);
}

static PyObject *py2c_padcon_standard_estimated_distance
 ( PyObject *self, PyObject *args )
{
   int fail;
   double step;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_standard_estimated_distance(&step);

   return Py_BuildValue("d",step);
}

static PyObject *py2c_padcon_dobldobl_estimated_distance
 ( PyObject *self, PyObject *args )
{
   int fail;
   double step;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_dobldobl_estimated_distance(&step);

   return Py_BuildValue("d",step);
}

static PyObject *py2c_padcon_quaddobl_estimated_distance
 ( PyObject *self, PyObject *args )
{
   int fail;
   double step;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_quaddobl_estimated_distance(&step);

   return Py_BuildValue("d",step);
}

static PyObject *py2c_padcon_standard_hessian_step
 ( PyObject *self, PyObject *args )
{
   int fail;
   double step;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_standard_hessian_step(&step);

   return Py_BuildValue("d",step);
}

static PyObject *py2c_padcon_dobldobl_hessian_step
 ( PyObject *self, PyObject *args )
{
   int fail;
   double step;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_dobldobl_hessian_step(&step);

   return Py_BuildValue("d",step);
}

static PyObject *py2c_padcon_quaddobl_hessian_step
 ( PyObject *self, PyObject *args )
{
   int fail;
   double step;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_get_quaddobl_hessian_step(&step);

   return Py_BuildValue("d",step);
}

static PyObject *py2c_padcon_standard_series_coefficient
 ( PyObject *self, PyObject *args )
{
   int fail,lead,idx,vrb;
   double cre,cim;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&lead,&idx,&vrb)) return NULL;
   fail = padcon_get_standard_series_coefficient(lead,idx,vrb,&cre,&cim);

   return Py_BuildValue("(d,d)",cre,cim);
}

static PyObject *py2c_padcon_dobldobl_series_coefficient
 ( PyObject *self, PyObject *args )
{
   int fail,lead,idx,vrb;
   double cre,cim;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&lead,&idx,&vrb)) return NULL;
   fail = padcon_get_dobldobl_series_coefficient(lead,idx,vrb,&cre,&cim);

   return Py_BuildValue("(d,d)",cre,cim);
}

static PyObject *py2c_padcon_quaddobl_series_coefficient
 ( PyObject *self, PyObject *args )
{
   int fail,lead,idx,vrb;
   double cre,cim;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&lead,&idx,&vrb)) return NULL;
   fail = padcon_get_quaddobl_series_coefficient(lead,idx,vrb,&cre,&cim);

   return Py_BuildValue("(d,d)",cre,cim);
}

static PyObject *py2c_padcon_standard_numerator_coefficient
 ( PyObject *self, PyObject *args )
{
   int fail,lead,idx,vrb;
   double cre,cim;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&lead,&idx,&vrb)) return NULL;
   fail = padcon_get_standard_numerator_coefficient(lead,idx,vrb,&cre,&cim);

   return Py_BuildValue("(d,d)",cre,cim);
}

static PyObject *py2c_padcon_dobldobl_numerator_coefficient
 ( PyObject *self, PyObject *args )
{
   int fail,lead,idx,vrb;
   double cre,cim;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&lead,&idx,&vrb)) return NULL;
   fail = padcon_get_dobldobl_numerator_coefficient(lead,idx,vrb,&cre,&cim);

   return Py_BuildValue("(d,d)",cre,cim);
}

static PyObject *py2c_padcon_quaddobl_numerator_coefficient
 ( PyObject *self, PyObject *args )
{
   int fail,lead,idx,vrb;
   double cre,cim;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&lead,&idx,&vrb)) return NULL;
   fail = padcon_get_quaddobl_numerator_coefficient(lead,idx,vrb,&cre,&cim);

   return Py_BuildValue("(d,d)",cre,cim);
}

static PyObject *py2c_padcon_standard_denominator_coefficient
 ( PyObject *self, PyObject *args )
{
   int fail,lead,idx,vrb;
   double cre,cim;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&lead,&idx,&vrb)) return NULL;
   fail = padcon_get_standard_denominator_coefficient(lead,idx,vrb,&cre,&cim);

   return Py_BuildValue("(d,d)",cre,cim);
}

static PyObject *py2c_padcon_dobldobl_denominator_coefficient
 ( PyObject *self, PyObject *args )
{
   int fail,lead,idx,vrb;
   double cre,cim;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&lead,&idx,&vrb)) return NULL;
   fail = padcon_get_dobldobl_denominator_coefficient(lead,idx,vrb,&cre,&cim);

   return Py_BuildValue("(d,d)",cre,cim);
}

static PyObject *py2c_padcon_quaddobl_denominator_coefficient
 ( PyObject *self, PyObject *args )
{
   int fail,lead,idx,vrb;
   double cre,cim;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&lead,&idx,&vrb)) return NULL;
   fail = padcon_get_quaddobl_denominator_coefficient(lead,idx,vrb,&cre,&cim);

   return Py_BuildValue("(d,d)",cre,cim);
}

static PyObject *py2c_padcon_standard_pole ( PyObject *self, PyObject *args )
{
   int fail,lead,idx,vrb;
   double cre,cim;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&lead,&idx,&vrb)) return NULL;
   fail = padcon_get_standard_pole(lead,idx,vrb,&cre,&cim);

   return Py_BuildValue("(d,d)",cre,cim);
}

static PyObject *py2c_padcon_dobldobl_pole ( PyObject *self, PyObject *args )
{
   int fail,lead,idx,vrb;
   double cre,cim;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&lead,&idx,&vrb)) return NULL;
   fail = padcon_get_dobldobl_pole(lead,idx,vrb,&cre,&cim);

   return Py_BuildValue("(d,d)",cre,cim);
}

static PyObject *py2c_padcon_quaddobl_pole ( PyObject *self, PyObject *args )
{
   int fail,lead,idx,vrb;
   double cre,cim;

   initialize();
   if(!PyArg_ParseTuple(args,"iii",&lead,&idx,&vrb)) return NULL;
   fail = padcon_get_quaddobl_pole(lead,idx,vrb,&cre,&cim);

   return Py_BuildValue("(d,d)",cre,cim);
}

static PyObject *py2c_padcon_clear_standard_data
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_clear_standard_data();
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_clear_dobldobl_data
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_clear_dobldobl_data();
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_padcon_clear_quaddobl_data
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = padcon_clear_quaddobl_data();
   return Py_BuildValue("i",fail);
}

/* The wrapping of functions with prototypes in syspool.h starts below. */

static PyObject *py2c_syspool_standard_init ( PyObject *self, PyObject *args )
{
   int fail,nbr;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&nbr)) return NULL;   
   fail = syspool_standard_initialize(nbr);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syspool_dobldobl_init ( PyObject *self, PyObject *args )
{
   int fail,nbr;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&nbr)) return NULL;   
   fail = syspool_dobldobl_initialize(nbr);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syspool_quaddobl_init ( PyObject *self, PyObject *args )
{
   int fail,nbr;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&nbr)) return NULL;   
   fail = syspool_quaddobl_initialize(nbr);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syspool_standard_size ( PyObject *self, PyObject *args )
{
   int fail,nbr;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syspool_standard_size(&nbr);
              
   return Py_BuildValue("i",nbr);
}

static PyObject *py2c_syspool_dobldobl_size ( PyObject *self, PyObject *args )
{
   int fail,nbr;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syspool_dobldobl_size(&nbr);
              
   return Py_BuildValue("i",nbr);
}

static PyObject *py2c_syspool_quaddobl_size ( PyObject *self, PyObject *args )
{
   int fail,nbr;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syspool_quaddobl_size(&nbr);
              
   return Py_BuildValue("i",nbr);
}

static PyObject *py2c_syspool_standard_create
 ( PyObject *self, PyObject *args )
{
   int fail,idx;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&idx)) return NULL;
   fail = syspool_standard_create(idx);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syspool_dobldobl_create
 ( PyObject *self, PyObject *args )
{
   int fail,idx;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&idx)) return NULL;
   fail = syspool_dobldobl_create(idx);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syspool_quaddobl_create
 ( PyObject *self, PyObject *args )
{
   int fail,idx;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&idx)) return NULL;
   fail = syspool_quaddobl_create(idx);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syspool_copy_to_standard_container
 ( PyObject *self, PyObject *args )
{
   int fail,idx;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&idx)) return NULL;
   fail = syspool_copy_to_standard_container(idx);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syspool_copy_to_dobldobl_container
 ( PyObject *self, PyObject *args )
{
   int fail,idx;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&idx)) return NULL;
   fail = syspool_copy_to_dobldobl_container(idx);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syspool_copy_to_quaddobl_container
 ( PyObject *self, PyObject *args )
{
   int fail,idx;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&idx)) return NULL;
   fail = syspool_copy_to_quaddobl_container(idx);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syspool_standard_clear ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syspool_standard_clear();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syspool_dobldobl_clear ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syspool_dobldobl_clear();
              
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_syspool_quaddobl_clear ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;   
   fail = syspool_quaddobl_clear();
              
   return Py_BuildValue("i",fail);
}

/* The wrapping of functions with prototypes in next_track.h starts below. */

static PyObject *py2c_initialize_standard_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,fixed;
   double regamma,imgamma;

   initialize();
   if(!PyArg_ParseTuple(args,"idd",&fixed,&regamma,&imgamma)) return NULL;
   fail = initialize_standard_homotopy(fixed,regamma,imgamma);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_initialize_dobldobl_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,fixed;
   double regamma,imgamma;

   initialize();
   if(!PyArg_ParseTuple(args,"idd",&fixed,&regamma,&imgamma)) return NULL;
   fail = initialize_dobldobl_homotopy(fixed,regamma,imgamma);
   return Py_BuildValue("i",fail);
}

static PyObject *py2c_initialize_quaddobl_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,fixed;
   double regamma,imgamma;

   initialize();
   if(!PyArg_ParseTuple(args,"idd",&fixed,&regamma,&imgamma)) return NULL;
   fail = initialize_quaddobl_homotopy(fixed,regamma,imgamma);
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

static PyObject *py2c_initialize_varbprec_homotopy
 ( PyObject *self, PyObject *args )
{
   int fail,fixed,nctgt,ncstr;
   char *tgt;
   char *str;

   initialize();
   if(!PyArg_ParseTuple(args,"iisis",&fixed,&nctgt,&tgt,&ncstr,&str))
      return NULL;
   fail = initialize_varbprec_homotopy(fixed,nctgt,tgt,ncstr,str);
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

static PyObject *py2c_initialize_varbprec_solution
 ( PyObject *self, PyObject *args )
{
   int fail,nv,nc;
   char *sol;

   initialize();
   if(!PyArg_ParseTuple(args,"iis",&nv,&nc,&sol)) return NULL;
   fail = initialize_varbprec_solution(nv,nc,sol);

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

static PyObject *py2c_next_varbprec_solution
 ( PyObject *self, PyObject *args )
{
   int fail,want,mxpr,mxit,verb,nc;
   char *sol;

   initialize();
   if(!PyArg_ParseTuple(args,"iiii",&want,&mxpr,&mxit,&verb)) return NULL;
   sol = next_varbprec_solution(want,mxpr,mxit,verb,&nc,&fail);

   return Py_BuildValue("(i,s)",fail,sol);
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

static PyObject *py2c_clear_varbprec_tracker
 ( PyObject *self, PyObject *args )
{
   int fail;

   initialize();
   if(!PyArg_ParseTuple(args,"")) return NULL;
   fail = clear_varbprec_tracker();
   return Py_BuildValue("i",fail);
}

/* The wrapping of Newton's method and path trackers with the evaluation
 * done by algorithmic differentiation in lib2path.h, starts here. */

// For the separate compilation of lib2path, extern declarations are needed.

extern int standard_ade_newton ( int verbose );
extern int standard_ade_onepath
 ( int verbose, double regamma, double imgamma );
extern int standard_ade_manypaths
 ( int verbose, double regamma, double imgamma );

extern int dobldobl_ade_newton ( int verbose );
extern int dobldobl_ade_onepath
 ( int verbose, double regamma, double imgamma );
extern int dobldobl_ade_manypaths
 ( int verbose, double regamma, double imgamma );

extern int quaddobl_ade_newton ( int verbose );
extern int quaddobl_ade_onepath
 ( int verbose, double regamma, double imgamma );
extern int quaddobl_ade_manypaths
 ( int verbose, double regamma, double imgamma );

static PyObject *py2c_ade_newton_d ( PyObject *self, PyObject *args )
{
   int fail,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&verbose)) return NULL;
   fail = standard_ade_newton(verbose);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_ade_newton_dd ( PyObject *self, PyObject *args )
{
   int fail,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&verbose)) return NULL;
   fail = dobldobl_ade_newton(verbose);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_ade_newton_qd ( PyObject *self, PyObject *args )
{
   int fail,verbose;

   initialize();
   if(!PyArg_ParseTuple(args,"i",&verbose)) return NULL;
   fail = quaddobl_ade_newton(verbose);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_ade_onepath_d ( PyObject *self, PyObject *args )
{
   int fail,verbose;
   double reg,img;

   initialize();
   if(!PyArg_ParseTuple(args,"idd",&verbose,&reg,&img)) return NULL;
   fail = standard_ade_onepath(verbose,reg,img);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_ade_onepath_dd ( PyObject *self, PyObject *args )
{
   int fail,verbose;
   double reg,img;

   initialize();
   if(!PyArg_ParseTuple(args,"idd",&verbose,&reg,&img)) return NULL;
   fail = dobldobl_ade_onepath(verbose,reg,img);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_ade_onepath_qd ( PyObject *self, PyObject *args )
{
   int fail,verbose;
   double reg,img;

   initialize();
   if(!PyArg_ParseTuple(args,"idd",&verbose,&reg,&img)) return NULL;
   fail = quaddobl_ade_onepath(verbose,reg,img);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_ade_manypaths_d ( PyObject *self, PyObject *args )
{
   int fail,verbose;
   double reg,img;

   initialize();
   if(!PyArg_ParseTuple(args,"idd",&verbose,&reg,&img)) return NULL;
   fail = standard_ade_manypaths(verbose,reg,img);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_ade_manypaths_dd ( PyObject *self, PyObject *args )
{
   int fail,verbose;
   double reg,img;

   initialize();
   if(!PyArg_ParseTuple(args,"idd",&verbose,&reg,&img)) return NULL;
   fail = dobldobl_ade_manypaths(verbose,reg,img);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_ade_manypaths_qd ( PyObject *self, PyObject *args )
{
   int fail,verbose;
   double reg,img;

   initialize();
   if(!PyArg_ParseTuple(args,"idd",&verbose,&reg,&img)) return NULL;
   fail = quaddobl_ade_manypaths(verbose,reg,img);

   return Py_BuildValue("i",fail);
}

/* For the parameter tuning we need the default path parameters,
 * provided by the lib2path function below,
 * see lib2path.h for its specification. */

extern int get_default_path_parameters
 ( int precision, int* max_step, int* n_predictor,
   double* step_increase, double* step_decrease,
   double* max_delta_t, double* max_delta_t_end, double* min_delta_t,
   double* err_max_res, double* err_max_delta_x, double* err_max_first_delta_x,
   int* max_it, double* err_min_round_off,
   int* max_it_refine, double* err_min_round_off_refine );

static PyObject *py2c_get_default_path_parameters
 ( PyObject *self, PyObject *args )
{
   int fail,precision,max_step,n_predictor,max_it,max_it_refine;
   double step_increase,step_decrease;
   double max_delta_t,max_delta_t_end,min_delta_t;
   double err_max_res,err_max_delta_x,err_max_first_delta_x;
   double err_min_round_off,err_min_round_off_refine;
 
   initialize();
   if(!PyArg_ParseTuple(args,"i",&precision)) return NULL;
   fail = get_default_path_parameters(precision,&max_step,&n_predictor,
      &step_increase,&step_decrease,&max_delta_t,&max_delta_t_end,&min_delta_t,
      &err_max_res,&err_max_delta_x,&err_max_first_delta_x,
      &max_it,&err_min_round_off,&max_it_refine,&err_min_round_off_refine);

   return Py_BuildValue("(i,i,d,d,d,d,d,d,d,d,i,d,i,d)",
      max_step,n_predictor,step_increase,
      step_decrease,max_delta_t,max_delta_t_end,min_delta_t,err_max_res,
      err_max_delta_x,err_max_first_delta_x,max_it,err_min_round_off,
      max_it_refine,err_min_round_off_refine);
}

// for separate compilation, extern declarations are needed, as below

extern int standard_ademanypaths_with_parameters
 ( int verbose, double regamma, double imgamma,
   int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine );

extern int dobldobl_ademanypaths_with_parameters
 ( int verbose, double regamma, double imgamma,
   int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine );

extern int quaddobl_ademanypaths_with_parameters
 ( int verbose, double regamma, double imgamma,
   int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine );

static PyObject *py2c_ade_manypaths_d_pars ( PyObject *self, PyObject *args )
{
   int fail,verbose,max_step,n_predictor,max_it,max_it_refine;
   double reg,img,step_increase,step_decrease;
   double max_delta_t,max_delta_t_end,min_delta_t;
   double err_max_res,err_max_delta_x,err_max_first_delta_x;
   double err_min_round_off,err_min_round_off_refine;

   initialize();
   if(!PyArg_ParseTuple(args,"iddiiddddddddidid",&verbose,&reg,&img,
      &max_step,&n_predictor,&step_increase,&step_decrease,
      &max_delta_t,&max_delta_t_end,&min_delta_t,&err_max_res,
      &err_max_delta_x,&err_max_first_delta_x,&max_it,&err_min_round_off,
      &max_it_refine,&err_min_round_off_refine)) return NULL;

   fail = standard_ademanypaths_with_parameters(verbose,reg,img,
      max_step,n_predictor,step_increase,step_decrease,
      max_delta_t,max_delta_t_end,min_delta_t,err_max_res,
      err_max_delta_x,err_max_first_delta_x,max_it,err_min_round_off,
      max_it_refine,err_min_round_off_refine);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_ade_manypaths_dd_pars ( PyObject *self, PyObject *args )
{
   int fail,verbose,max_step,n_predictor,max_it,max_it_refine;
   double reg,img,step_increase,step_decrease;
   double max_delta_t,max_delta_t_end,min_delta_t;
   double err_max_res,err_max_delta_x,err_max_first_delta_x;
   double err_min_round_off,err_min_round_off_refine;
 
   initialize();
   if(!PyArg_ParseTuple(args,"iddiiddddddddidid",&verbose,&reg,&img,
      &max_step,&n_predictor,&step_increase,&step_decrease,
      &max_delta_t,&max_delta_t_end,&min_delta_t,&err_max_res,
      &err_max_delta_x,&err_max_first_delta_x,&max_it,&err_min_round_off,
      &max_it_refine,&err_min_round_off_refine)) return NULL;

   fail = dobldobl_ademanypaths_with_parameters(verbose,reg,img,
      max_step,n_predictor,step_increase,step_decrease,
      max_delta_t,max_delta_t_end,min_delta_t,err_max_res,
      err_max_delta_x,err_max_first_delta_x,max_it,err_min_round_off,
      max_it_refine,err_min_round_off_refine);

   return Py_BuildValue("i",fail);
}

static PyObject *py2c_ade_manypaths_qd_pars ( PyObject *self, PyObject *args )
{
   int fail,verbose,max_step,n_predictor,max_it,max_it_refine;
   double reg,img,step_increase,step_decrease;
   double max_delta_t,max_delta_t_end,min_delta_t;
   double err_max_res,err_max_delta_x,err_max_first_delta_x;
   double err_min_round_off,err_min_round_off_refine;
 
   initialize();
   if(!PyArg_ParseTuple(args,"iddiiddddddddidid",&verbose,&reg,&img,
      &max_step,&n_predictor,&step_increase,&step_decrease,
      &max_delta_t,&max_delta_t_end,&min_delta_t,&err_max_res,
      &err_max_delta_x,&err_max_first_delta_x,&max_it,&err_min_round_off,
      &max_it_refine,&err_min_round_off_refine)) return NULL;

   fail = quaddobl_ademanypaths_with_parameters(verbose,reg,img,
      max_step,n_predictor,step_increase,step_decrease,
      max_delta_t,max_delta_t_end,min_delta_t,err_max_res,
      err_max_delta_x,err_max_first_delta_x,max_it,err_min_round_off,
      max_it_refine,err_min_round_off_refine);

   return Py_BuildValue("i",fail);
}

static PyMethodDef phcpy2c_methods[] = 
{
   {"py2c_corecount", py2c_corecount, METH_VARARGS,
    "Returns the number of cores available for multithreading."},
   {"py2c_PHCpack_version_string", py2c_PHCpack_version_string, METH_VARARGS,
    "Returns the version string of PHCpack.\n The version string is 40 characters long."},
   {"py2c_set_seed", py2c_set_seed, METH_VARARGS, 
    "Takes the value of the integer given on input\n and sets the seed for the random number generators.\n This fixing of the seed enables reproducible runs." },
   {"py2c_get_seed", py2c_get_seed, METH_VARARGS,
    "Returns the current value of the seed.\n Using this value in py2c_set_seed will ensure that the\n results of previous runs can be reproduced."},
   {"py2c_read_standard_target_system", py2c_read_standard_target_system,
     METH_VARARGS, 
    "Prompts the user to enter a target system that will\n be parsed in standard double precision.\n The failure code is returned, which is zero if all went well."},
   {"py2c_read_standard_target_system_from_file",
     py2c_read_standard_target_system_from_file, METH_VARARGS,
    "The two input arguments are a number and a string:\n 1) The number equals the number of characters in the string.\n 2) The string given on input is the name of a file which contains\n a target system to be parsed in standard double precision.\n The failure code is returned, which is zero if all went well."},
   {"py2c_read_standard_start_system", py2c_read_standard_start_system,
     METH_VARARGS,
    "Prompts the user to enter a start system that will\n be parsed in standard double precision.\n The failure code is returned, which is zero if all went well."},
   {"py2c_read_standard_start_system_from_file",
     py2c_read_standard_start_system_from_file, METH_VARARGS,
    "The two input arguments are a number and a string:\n 1) The number equals the number of characters in the string.\n 2) The string given on input is the name of a file which contains\n a start system to be parsed in standard double precision.\n The failure code is returned, which is zero if all went well."},
   {"py2c_read_dobldobl_target_system", py2c_read_dobldobl_target_system,
     METH_VARARGS,
    "Prompts the user to enter a target system that will\n be parsed in double double precision.\n The failure code is returned, which is zero if all went well."},
   {"py2c_read_dobldobl_target_system_from_file",
     py2c_read_dobldobl_target_system_from_file, METH_VARARGS,
   "The two input arguments are a number and a string:\n 1) The number equals the number of characters in the string.\n 2) The string given on input is the name of a file which contains\n a target system to be parsed in double double precision.\n The failure code is returned, which is zero if all went well."},
   {"py2c_read_dobldobl_start_system", py2c_read_dobldobl_start_system,
     METH_VARARGS,
    "Prompts the user to enter a start system that will\n be parsed in double double precision.\n The failure code is returned, which is zero if all went well."},
   {"py2c_read_dobldobl_start_system_from_file",
     py2c_read_dobldobl_start_system_from_file, METH_VARARGS, 
    "The two input arguments are a number and a string:\n 1) The number equals the number of characters in the string.\n 2) The string given on input is the name of a file which contains\n a start system to be parsed in double double precision.\n The failure code is returned, which is zero if all went well."},
   {"py2c_read_quaddobl_target_system", py2c_read_quaddobl_target_system,
     METH_VARARGS,
    "Prompts the user to enter a target system that will\n be parsed in quad double precision.\n The failure code is returned, which is zero if all went well."},
   {"py2c_read_quaddobl_target_system_from_file",
     py2c_read_quaddobl_target_system_from_file, METH_VARARGS,
    "The two input arguments are a number and a string:\n 1) The number equals the number of characters in the string.\n 2) The string given on input is the name of a file which contains\n a target system to be parsed in quad double precision.\n The failure code is returned, which is zero if all went well."},
   {"py2c_read_quaddobl_start_system", py2c_read_quaddobl_start_system,
     METH_VARARGS,
    "Prompts the user to enter a start system that will\n be parsed in quad double precision.\n The failure code is returned, which is zero if all went well."},
   {"py2c_read_quaddobl_start_system_from_file",
     py2c_read_quaddobl_start_system_from_file, METH_VARARGS,
    "The two input arguments are a number and a string:\n 1) The number equals the number of characters in the string.\n 2) The string given on input is the name of a file which contains\n a start system to be parsed in quad double precision.\n The failure code is returned, which is zero if all went well."},
   {"py2c_define_output_file", py2c_define_output_file, METH_VARARGS,
    "Prompts the user to define the output file.\n On return is the failure code, which is zero if all went well."},
   {"py2c_write_standard_target_system",
     py2c_write_standard_target_system, METH_VARARGS, 
    "Writes the target system as stored in standard double precision\n to screen or to the defined output file."},
   {"py2c_write_dobldobl_target_system",
     py2c_write_dobldobl_target_system, METH_VARARGS, 
    "Writes the target system as stored in double double precision\n to screen or to the defined output file."},
   {"py2c_write_quaddobl_target_system",
     py2c_write_quaddobl_target_system, METH_VARARGS, 
    "Writes the target system as stored in quad double precision\n to screen or to the defined output file."},
   {"py2c_write_standard_start_system",
     py2c_write_standard_start_system, METH_VARARGS,
    "Writes the start system as stored in standard double precision\n to screen or to the defined output file."},
   {"py2c_write_dobldobl_start_system",
     py2c_write_dobldobl_start_system, METH_VARARGS,
    "Writes the start system as stored in double double precision\n to screen or to the defined output file."},
   {"py2c_write_quaddobl_start_system",
     py2c_write_quaddobl_start_system, METH_VARARGS,
    "Writes the start system as stored in quad double precision\n to screen or to the defined output file."},
   {"py2c_read_standard_start_Laurent_system",
     py2c_read_standard_start_Laurent_system, METH_VARARGS,
    "Prompts the user for a file name and reads the start system from file,\n in standard double precision.\n If available on file, also its solutions will be read and stored."}, 
   {"py2c_write_standard_start_Laurent_system",
     py2c_write_standard_start_Laurent_system, METH_VARARGS,
    "Writes the start Laurent system in standard double precision."},
   {"py2c_read_standard_target_Laurent_system",
     py2c_read_standard_target_Laurent_system, METH_VARARGS,
    "Prompts the user for a file name and reads the target system from file,\n in standard double precision.\n If available on file, also its solutions will be read and stored."},
   {"py2c_write_standard_target_Laurent_system",
     py2c_write_standard_target_Laurent_system, METH_VARARGS,
    "Writes the target Laurent system in standard double precision."},
   {"py2c_read_dobldobl_start_Laurent_system",
     py2c_read_dobldobl_start_Laurent_system, METH_VARARGS,
    "Prompts the user for a file name and reads the start system from file,\n in double double precision.\n If available on file, also its solutions will be read and stored."},
   {"py2c_write_dobldobl_start_Laurent_system",
     py2c_write_dobldobl_start_Laurent_system, METH_VARARGS,
    "Writes the start Laurent system in double double precision."},
   {"py2c_read_dobldobl_target_Laurent_system",
     py2c_read_dobldobl_target_Laurent_system, METH_VARARGS,
    "Prompts the user for a file name and reads the target system from file,\n in double double precision.\n If available on file, also its solutions will be read and stored."},
   {"py2c_write_dobldobl_target_Laurent_system",
     py2c_write_dobldobl_target_Laurent_system, METH_VARARGS,
    "Writes the target Laurent system in double double precision."},
   {"py2c_read_quaddobl_start_Laurent_system",
     py2c_read_quaddobl_start_Laurent_system, METH_VARARGS,
    "Prompts the user for a file name and reads the start system from file,\n in quad double precision.\n If available on file, also its solutions will be read and stored."},
   {"py2c_write_quaddobl_start_Laurent_system",
     py2c_write_quaddobl_start_Laurent_system, METH_VARARGS,
    "Writes the start Laurent system in quad double precision."},
   {"py2c_read_quaddobl_target_Laurent_system",
     py2c_read_quaddobl_target_Laurent_system, METH_VARARGS,
    "Prompts the user for a file name and reads the target system from file,\n in quad double precision.\n If available on file, also its solutions will be read and stored."},
   {"py2c_write_quaddobl_target_Laurent_system",
     py2c_write_quaddobl_target_Laurent_system, METH_VARARGS,
    "Writes the target Laurent system in quad double precision."},
   {"py2c_write_start_solutions", py2c_write_start_solutions, METH_VARARGS,
    "Writes the start solutions in standard double precision either to\n the screen (standard output) or to the defined output file.\n On return is the failure code, which is zero if all is well."},
   {"py2c_copy_standard_target_system_to_container",
     py2c_copy_standard_target_system_to_container, METH_VARARGS,
    "Copies the target system to the container for systems\n with coefficients in standard double precision."},
   {"py2c_copy_dobldobl_target_system_to_container",
     py2c_copy_dobldobl_target_system_to_container, METH_VARARGS,
    "Copies the target system to the container for systems\n with coefficients in double double precision."},
   {"py2c_copy_quaddobl_target_system_to_container",
     py2c_copy_quaddobl_target_system_to_container, METH_VARARGS,
    "Copies the target system to the container for systems\n with coefficients in quad double precision."},
   {"py2c_copy_multprec_target_system_to_container",
     py2c_copy_multprec_target_system_to_container, METH_VARARGS,
    "copies multiprecision target system to container"},
   {"py2c_copy_standard_container_to_target_system",
     py2c_copy_standard_container_to_target_system, METH_VARARGS,
    "Copies the system in the container for systems with coefficients\n in standard double precision to the target system."},
   {"py2c_copy_dobldobl_container_to_target_system",
     py2c_copy_dobldobl_container_to_target_system, METH_VARARGS,
    "Copies the system in the container for systems with coefficients\n in double double precision to the target system."},
   {"py2c_copy_quaddobl_container_to_target_system",
     py2c_copy_quaddobl_container_to_target_system, METH_VARARGS,
    "Copies the system in the container for systems with coefficients\n in quad double precision to the target system."},
   {"py2c_copy_multprec_container_to_target_system",
     py2c_copy_multprec_container_to_target_system, METH_VARARGS,
    "Copies the system in the container for systems with coefficients\n in arbitrary multiprecision to the target system."},
   {"py2c_copy_start_system_to_container",
     py2c_copy_start_system_to_container, METH_VARARGS,
    "Copies the start system to the container for systems\n with coefficients in standard double precision."},
   {"py2c_copy_dobldobl_start_system_to_container",
     py2c_copy_dobldobl_start_system_to_container, METH_VARARGS,
    "Copies the start system to the container for systems\n with coefficients in double double precision."},
   {"py2c_copy_quaddobl_start_system_to_container",
     py2c_copy_quaddobl_start_system_to_container, METH_VARARGS, 
    "Copies the start system to the container for systems\n with coefficients in quad double precision."},
   {"py2c_copy_multprec_start_system_to_container",
     py2c_copy_multprec_start_system_to_container, METH_VARARGS, 
    "Copies the start system to the container for systems\n with coefficients in arbitrary multiprecision."},
   {"py2c_copy_standard_container_to_start_system",
     py2c_copy_standard_container_to_start_system, METH_VARARGS, 
    "Copies the system in the container for systems with coefficients\n in standard double precision to the start system."},
   {"py2c_copy_dobldobl_container_to_start_system",
     py2c_copy_dobldobl_container_to_start_system, METH_VARARGS,
    "Copies the system in the container for systems with coefficients\n in double double precision to the start system."},
   {"py2c_copy_quaddobl_container_to_start_system",
     py2c_copy_quaddobl_container_to_start_system, METH_VARARGS,
    "Copies the system in the container for systems with coefficients\n in quad double precision to the start system."},
   {"py2c_copy_multprec_container_to_start_system",
     py2c_copy_multprec_container_to_start_system, METH_VARARGS,
    "Copies the system in the container for systems with coefficients\n in arbitrary multiprecision to the start system."},
   {"py2c_copy_standard_Laurent_container_to_start_system",
     py2c_copy_standard_Laurent_container_to_start_system, METH_VARARGS,
    "Copies the Laurent system in standard double precision\n from the container to the start system."},
   {"py2c_copy_dobldobl_Laurent_container_to_start_system",
     py2c_copy_dobldobl_Laurent_container_to_start_system, METH_VARARGS,
    "Copies the Laurent system in double double precision\n from the container to the start system."},
   {"py2c_copy_quaddobl_Laurent_container_to_start_system",
     py2c_copy_quaddobl_Laurent_container_to_start_system, METH_VARARGS,
    "Copies the Laurent system in quad double precision\n from the container to the start system."},
   {"py2c_copy_standard_Laurent_container_to_target_system",
     py2c_copy_standard_Laurent_container_to_target_system, METH_VARARGS,
    "Copies the Laurent system in standard double precision\n from the container to the target system."},
   {"py2c_copy_dobldobl_Laurent_container_to_target_system",
     py2c_copy_dobldobl_Laurent_container_to_target_system, METH_VARARGS,
    "Copies the Laurent system in double double precision\n from the container to the target system."},
   {"py2c_copy_quaddobl_Laurent_container_to_target_system",
     py2c_copy_quaddobl_Laurent_container_to_target_system, METH_VARARGS,
    "Copies the Laurent system in quad double precision\n from the container to the target system."},
   {"py2c_copy_standard_Laurent_start_system_to_container",
     py2c_copy_standard_Laurent_start_system_to_container, METH_VARARGS,
    "Copies the start Laurent system in standard double precision\n to the systems container for Laurent systems."},
   {"py2c_copy_dobldobl_Laurent_start_system_to_container",
     py2c_copy_dobldobl_Laurent_start_system_to_container, METH_VARARGS,
    "Copies the start Laurent system in double double precision\n to the systems container for Laurent systems."},
   {"py2c_copy_quaddobl_Laurent_start_system_to_container",
     py2c_copy_quaddobl_Laurent_start_system_to_container, METH_VARARGS,
    "Copies the start Laurent system in quad double precision\n to the systems container for Laurent systems."},
   {"py2c_copy_standard_Laurent_target_system_to_container",
     py2c_copy_standard_Laurent_target_system_to_container, METH_VARARGS,
    "Copies the target Laurent system in standard double precision\n to the systems container for Laurent systems."},
   {"py2c_copy_dobldobl_Laurent_target_system_to_container",
     py2c_copy_dobldobl_Laurent_target_system_to_container, METH_VARARGS,
    "Copies the target Laurent system in double double precision\n to the systems container for Laurent systems."},
   {"py2c_copy_quaddobl_Laurent_target_system_to_container",
     py2c_copy_quaddobl_Laurent_target_system_to_container, METH_VARARGS,
    "Copies the target Laurent system in quad double precision\n to the systems container for Laurent systems."},
   {"py2c_create_standard_homotopy", py2c_create_standard_homotopy,
     METH_VARARGS,
    "Initializes the data for a homotopy in standard double precision.\n The failure code is returned, which is zero when all goes well."},
   {"py2c_create_standard_homotopy_with_gamma",
     py2c_create_standard_homotopy_with_gamma, METH_VARARGS,
    "Initializes the data for a homotopy in standard double precision.\n On input are two doubles and one positive integer:\n (1) the real and imaginary part of the gamma constant;\n (2) the power of t in the homotopy.\n The failure code is returned, which is zero when all goes well."},
   {"py2c_create_dobldobl_homotopy", py2c_create_dobldobl_homotopy,
     METH_VARARGS,
    "Initializes the data for a homotopy in double double precision.\n The failure code is returned, which is zero when all goes well."},
   {"py2c_create_dobldobl_homotopy_with_gamma",
     py2c_create_dobldobl_homotopy_with_gamma, METH_VARARGS,
    "Initializes the data for a homotopy in double double precision.\n On input are two doubles and one positive integer:\n (1) the real and imaginary part of the gamma constant;\n (2) the power of t in the homotopy.\n The failure code is returned, which is zero when all goes well."},
   {"py2c_create_quaddobl_homotopy", py2c_create_quaddobl_homotopy,
     METH_VARARGS,
    "Initializes the data for a homotopy in quad double precision.\n The failure code is returned, which is zero when all goes well."},
   {"py2c_create_quaddobl_homotopy_with_gamma",
     py2c_create_quaddobl_homotopy_with_gamma, METH_VARARGS,
    "Initializes the data for a homotopy in quad double precision.\n On input are two doubles and one positive integer:\n (1) the real and imaginary part of the gamma constant;\n (2) the power of t in the homotopy.\n The failure code is returned, which is zero when all goes well."},
   {"py2c_create_multprec_homotopy", py2c_create_multprec_homotopy,
     METH_VARARGS,
    "Initializes the data for a homotopy in arbitrary multiprecision.\n The failure code is returned, which is zero when all goes well."},
   {"py2c_create_multprec_homotopy_with_gamma",
     py2c_create_multprec_homotopy_with_gamma, METH_VARARGS,
    "Initializes the data for a homotopy in arbitrary multiprecision.\n On input are two doubles and one positive integer:\n (1) the real and imaginary part of the gamma constant;\n (2) the power of t in the homotopy.\n The failure code is returned, which is zero when all goes well."},
   {"py2c_clear_standard_homotopy", py2c_clear_standard_homotopy, METH_VARARGS,
    "Deallocation of the homotopy stored in standard double precision.\n On return is the failure code, which equals zero if all is well."},
   {"py2c_clear_dobldobl_homotopy", py2c_clear_dobldobl_homotopy, METH_VARARGS,
    "Deallocation of the homotopy stored in double double precision.\n On return is the failure code, which equals zero if all is well."},
   {"py2c_clear_quaddobl_homotopy", py2c_clear_quaddobl_homotopy, METH_VARARGS,
    "Deallocation of the homotopy stored in quad double precision.\n On return is the failure code, which equals zero if all is well."},
   {"py2c_clear_multprec_homotopy", py2c_clear_multprec_homotopy, METH_VARARGS,
    "Deallocation of the homotopy stored in arbitrary multiprecision.\n On return is the failure code, which equals zero if all is well."},
   {"py2c_tune_continuation_parameters", py2c_tune_continuation_parameters,
     METH_VARARGS,
    "Interactive procedure to tune the continuation parameters."},
   {"py2c_show_continuation_parameters", py2c_show_continuation_parameters,
     METH_VARARGS,
    "Shows the current values of the continuation parameters."},
   {"py2c_autotune_continuation_parameters",
     py2c_autotune_continuation_parameters, METH_VARARGS, 
    "Tunes the values of the continuation parameters.\n On input are two integers:\n 1) the difficulty level of the solution paths; and\n 2) the number of decimal places in the precision."},
   {"py2c_get_value_of_continuation_parameter",
     py2c_get_value_of_continuation_parameter, METH_VARARGS,
   "Returns the value of a continuation parameter.\n On input is the index of this continuation parameter, an integer ranging from 1 to 34.\n On return is a double with the value of the corresponding parameter."},
   {"py2c_set_value_of_continuation_parameter",
     py2c_set_value_of_continuation_parameter, METH_VARARGS,
   "Sets the value of a continuation parameter.\n On input is the index of this continuation parameter, an integer ranging from 1 to 34;\n and the new value for the continuation parameter.\n On return is a double with the value of the corresponding parameter."},
   {"py2c_determine_output_during_continuation", 
     py2c_determine_output_during_continuation, METH_VARARGS, 
    "Interactive procedure to determine the level of output during the path tracking."},
   {"py2c_solve_by_standard_homotopy_continuation",
     py2c_solve_by_standard_homotopy_continuation, METH_VARARGS,
    "Tracks the paths defined by the homotopy in standard double precision.\n On input is one integer: the number of tasks for path tracking.\n If that input number is zero, then no multitasking is applied.\n On return is the failure code, which is zero when all went well."},
   {"py2c_solve_by_dobldobl_homotopy_continuation",
     py2c_solve_by_dobldobl_homotopy_continuation, METH_VARARGS, 
    "Tracks the paths defined by the homotopy in double double precision.\n On input is one integer: the number of tasks for path tracking.\n If that input number is zero, then no multitasking is applied.\n On return is the failure code, which is zero when all went well."},
   {"py2c_solve_by_quaddobl_homotopy_continuation",
     py2c_solve_by_quaddobl_homotopy_continuation, METH_VARARGS,
    "Tracks the paths defined by the homotopy in quad double precision.\n On input is one integer: the number of tasks for path tracking.\n If that input number is zero, then no multitasking is applied.\n On return is the failure code, which is zero when all went well."},
   {"py2c_solve_by_multprec_homotopy_continuation",
     py2c_solve_by_multprec_homotopy_continuation, METH_VARARGS, 
    "Tracks the paths defined by the homotopy in arbitrary multiprecision.\n On input is one integer: the number of decimal places in the precision.\n On return is the failure code, which is zero when all went well."},
   {"py2c_solve_by_standard_Laurent_homotopy_continuation",
     py2c_solve_by_standard_Laurent_homotopy_continuation, METH_VARARGS,
    "Tracks the paths defined by the homotopy in standard double precision\n to solve a Laurent system stored in the systems container,\n starting at the solutions of a stored Laurent start system.\n On input is one integer: the number of tasks for path tracking.\n If that input number is zero, then no multitasking is applied.\n On return is the failure code, which is zero when all went well."},
   {"py2c_solve_by_dobldobl_Laurent_homotopy_continuation",
     py2c_solve_by_dobldobl_Laurent_homotopy_continuation, METH_VARARGS,
    "Tracks the paths defined by the homotopy in double double precision\n to solve a Laurent system stored in the systems container,\n starting at the solutions of a stored Laurent start system.\n On input is one integer: the number of tasks for path tracking.\n If that input number is zero, then no multitasking is applied.\n On return is the failure code, which is zero when all went well."},
   {"py2c_solve_by_quaddobl_Laurent_homotopy_continuation",
     py2c_solve_by_quaddobl_Laurent_homotopy_continuation, METH_VARARGS,
    "Tracks the paths defined by the homotopy in quad double precision\n to solve a Laurent system stored in the systems container,\n starting at the solutions of a stored Laurent start system.\n On input is one integer: the number of tasks for path tracking.\n If that input number is zero, then no multitasking is applied.\n On return is the failure code, which is zero when all went well."},
   {"py2c_clear_standard_operations_data",
     py2c_clear_standard_operations_data, METH_VARARGS,
    "Deallocates the data used by solve_by_standard_homotopy_continuation."},
   {"py2c_clear_dobldobl_operations_data",
     py2c_clear_dobldobl_operations_data, METH_VARARGS,
    "Deallocates the data used by solve_by_dobldobl_homotopy_continuation."},
   {"py2c_clear_quaddobl_operations_data",
     py2c_clear_quaddobl_operations_data, METH_VARARGS,
    "Deallocates the data used by solve_by_quaddobl_homotopy_continuation."},
   {"py2c_clear_standard_Laurent_data",
     py2c_clear_standard_Laurent_data, METH_VARARGS,
    "Deallocates data used to solve Laurent systems by homotopy continuation\n in standard double precision."},
   {"py2c_clear_dobldobl_Laurent_data",
     py2c_clear_dobldobl_Laurent_data, METH_VARARGS,
    "Deallocates data used to solve Laurent systems by homotopy continuation\n in double double precision."},
   {"py2c_clear_quaddobl_Laurent_data",
     py2c_clear_quaddobl_Laurent_data, METH_VARARGS,
    "Deallocates data used to solve Laurent systems by homotopy continuation\n in quad double precision."},
   {"py2c_standard_crude_tracker", py2c_standard_crude_tracker, METH_VARARGS,
    "A crude tracker appends the end point of a path directly to\n the solutions container, without refinement or postprocessing.\n Tracking happens in standard double precision.\n On entry is the verbose parameter which is 1 or 0.\n If 1, then the solution vectors are written to screen, otherwise\n the crude tracker stays mute.\n On return is the failure code, which is zero when all went well.\n The requirement is that\n the target system, start system, and start solutions in standard\n double precision have been initialized in the containers."},
   {"py2c_dobldobl_crude_tracker", py2c_dobldobl_crude_tracker, METH_VARARGS,
    "A crude tracker appends the end point of a path directly to\n the solutions container, without refinement or postprocessing.\n Tracking happens in double double precision.\n On entry is the verbose parameter which is 1 or 0.\n If 1, then the solution vectors are written to screen, otherwise\n the crude tracker stays mute.\n On return is the failure code, which is zero when all went well.\n The requirement is that\n the target system, start system, and start solutions in double\n double precision have been initialized in the containers."},
   {"py2c_quaddobl_crude_tracker", py2c_quaddobl_crude_tracker, METH_VARARGS,
    "A crude tracker appends the end point of a path directly to\n the solutions container, without refinement or postprocessing.\n Tracking happens in quad double precision.\n On entry is the verbose parameter which is 1 or 0.\n If 1, then the solution vectors are written to screen, otherwise\n the crude tracker stays mute.\n On return is the failure code, which is zero when all went well.\n The requirement is that\n the target system, start system, and start solutions in quad\n double precision have been initialized in the containers."},
   {"py2c_copy_standard_target_solutions_to_container",
     py2c_copy_standard_target_solutions_to_container, METH_VARARGS,
    "Copies the target solutions in standard double precision to the\n container for solutions in standard double precision."},
   {"py2c_copy_dobldobl_target_solutions_to_container",
     py2c_copy_dobldobl_target_solutions_to_container, METH_VARARGS,
    "Copies the target solutions in double double precision to the\n container for solutions in double double precision."},
   {"py2c_copy_quaddobl_target_solutions_to_container",
     py2c_copy_quaddobl_target_solutions_to_container, METH_VARARGS,
    "Copies the target solutions in quad double precision to the\n container for solutions in quad double precision."},
   {"py2c_copy_multprec_target_solutions_to_container",
     py2c_copy_multprec_target_solutions_to_container, METH_VARARGS,
    "Copies the target solutions in arbitrary multiprecision to the\n container for solutions in arbitrary multiprecision."},
   {"py2c_copy_standard_container_to_target_solutions",
     py2c_copy_standard_container_to_target_solutions, METH_VARARGS,
    "Copies the solutions in standard double precision from the\n container to the target solutions in standard double precision."},
   {"py2c_copy_dobldobl_container_to_target_solutions",
     py2c_copy_dobldobl_container_to_target_solutions, METH_VARARGS,
    "Copies the solutions in double double precision from the\n container to the target solutions in double double precision."},
   {"py2c_copy_quaddobl_container_to_target_solutions",
     py2c_copy_quaddobl_container_to_target_solutions, METH_VARARGS,
    "Copies the solutions in quad double precision from the\n container to the target solutions in quad double precision."},
   {"py2c_copy_multprec_container_to_target_solutions",
     py2c_copy_multprec_container_to_target_solutions, METH_VARARGS,
    "Copies the solutions in arbitrary multiprecision from the\n container to the target solutions in arbitrary multiprecision."},
   {"py2c_copy_start_solutions_to_container",
     py2c_copy_start_solutions_to_container, METH_VARARGS,
    "Copies the start solutions in standard double precision to the\n container for solutions in standard double precision."},
   {"py2c_copy_dobldobl_start_solutions_to_container",
     py2c_copy_dobldobl_start_solutions_to_container, METH_VARARGS,
    "Copies the start solutions in double double precision to the\n container for solutions in double double precision."},
   {"py2c_copy_quaddobl_start_solutions_to_container",
     py2c_copy_quaddobl_start_solutions_to_container, METH_VARARGS,
    "Copies the start solutions in quad double precision to the\n container for solutions in quad double precision."},
   {"py2c_copy_multprec_start_solutions_to_container",
     py2c_copy_multprec_start_solutions_to_container, METH_VARARGS,
    "Copies the start solutions in arbitrary multiprecision to the\n container for solutions in arbitrary multiprecision."},
   {"py2c_copy_standard_container_to_start_solutions",
     py2c_copy_standard_container_to_start_solutions, METH_VARARGS,
    "Copies the solutions in standard double precision from the\n container to the start solutions in standard double precision."},
   {"py2c_copy_dobldobl_container_to_start_solutions",
     py2c_copy_dobldobl_container_to_start_solutions, METH_VARARGS, 
    "Copies the solutions in double double precision from the\n container to the start solutions in double double precision."},
   {"py2c_copy_quaddobl_container_to_start_solutions",
     py2c_copy_quaddobl_container_to_start_solutions, METH_VARARGS,
    "Copies the solutions in quad double precision from the\n container to the start solutions in quad double precision."},
   {"py2c_copy_multprec_container_to_start_solutions",
     py2c_copy_multprec_container_to_start_solutions, METH_VARARGS,
    "Copies the solutions in arbitrary multiprecision from the\n container to the start solutions in arbitrary multiprecision."},
   {"py2c_solve_standard_system", py2c_solve_standard_system, METH_VARARGS,
    "Calls the blackbox solver on the system stored in the container for\n systems with coefficients in standard double precision.\n Four integers are expected on input:\n 1) a boolean flag silent: if 1, then no intermediate output about\n the root counts is printed, if 0, then the solver is verbose;\n 2) the number of tasks.\n If that number is zero, then no multitasking is applied;\n 3) a flag to focus on mixed volumes and polyhedral homotopies,\n if 1, then no degree bounds will be computed; and\n 4) the verbose level.\n On return, the container for solutions in standard double precision\n contains the solutions to the system in the standard systems container."},
   {"py2c_scan_for_symbols", py2c_scan_for_symbols, METH_VARARGS,
    "Given on input are two arguments: a number and a string.\n The string holds the string representation of a polynomial system,\n where each polynomial is terminated by a semi colon.\n The first argument on input is the number of characters in the string.\n On return is the number of symbols used as variables in the system.\n This function helps to determine whether a system is square or not."},
   {"py2c_solve_dobldobl_system", py2c_solve_dobldobl_system, METH_VARARGS,
    "Calls the blackbox solver on the system stored in the container for\n systems with coefficients in double double precision.\n Three integers are expected on input:\n 1) a boolean flag silent: if 1, then no intermediate output about\n the root counts is printed, if 0, then the solver is verbose;\n 2) the number of tasks.\n If that number is zero, then no multitasking is applied; and\n 3) the verbose level.\n On return, the container for solutions in double double precision\n contains the solutions to the system in the dobldobl systems container."},
   {"py2c_solve_quaddobl_system", py2c_solve_quaddobl_system, METH_VARARGS,
    "Calls the blackbox solver on the system stored in the container for\n systems with coefficients in quad double precision.\n Three integers are expected on input: 1) a boolean flag silent: if 1, then no intermediate output about\n the root counts is printed, if 0, then the solver is verbose;\n 2) the number of tasks.\n If that number is zero, then no multitasking is applied; and\n 3) the verbose level. \n On return, the container for solutions in quad double precision\n contains the solutions to the system in the quaddobl systems container."},
   {"py2c_solve_standard_Laurent_system",
     py2c_solve_standard_Laurent_system, METH_VARARGS,
    "Calls the blackbox solver on the system stored in the container for\n Laurent systems with coefficients in standard double precision.\n Three integers are expected on input:\n 1) a boolean flag silent: if 1, then no intermediate output about\n the root counts is printed, if 0, then the solver is verbose; \n 2) the number of tasks: if 0, then no multitasking is applied,\n otherwise as many tasks as the number will run;\n 3) a flag to focus on mixed volumes and polyhedral homotopies,\n if 1, then no degree bounds will be computed, as when\n the system is genuinely Laurent with negative exponents; and\n 4) the verbose level.\n On return, the container for solutions in standard double precision\n contains the solutions to the system in the standard Laurent systems\n container."},
   {"py2c_solve_dobldobl_Laurent_system",
     py2c_solve_dobldobl_Laurent_system, METH_VARARGS,
    "Calls the blackbox solver on the system stored in the container for\n Laurent systems with coefficients in double double precision.\n Four integers are expected on input:\n 1) a boolean flag silent: if 1, then no intermediate output about\n the root counts is printed, if 0, then the solver is verbose; \n 2) the number of tasks: if 0, then no multitasking is applied,\n otherwise as many tasks as the number will run; and\n 3) 4) the verbose level.\n On return, the container for solutions in double double precision\n contains the solutions to the system in the double double Laurent systems\n container."},
   {"py2c_solve_quaddobl_Laurent_system",
     py2c_solve_quaddobl_Laurent_system, METH_VARARGS,
    "Calls the blackbox solver on the system stored in the container for\n Laurent systems with coefficients in quad double precision.\n Three integers are expected on input:\n 1) a boolean flag silent: if 1, then no intermediate output about\n the root counts is printed, if 0, then the solver is verbose; \n 2) the number of tasks: if 0, then no multitasking is applied,\n otherwise as many tasks as the number will run; and\n 3) the verbose level.\n On return, the container for solutions in quad double precision\n contains the solutions to the system in the quad double Laurent systems\n container."},
   {"py2c_set_gamma_constant", py2c_set_gamma_constant, METH_VARARGS,
    "Stores the gamma constant for later retrieval.\n Four parameters are expected on input, two doubles and two integers.\n The two doubles are the real and imaginary parts of the gamma.\n The two integers are the precision, 1, 2, or 4, respectively for\n double, double double, or quad double; and the verbose level."},
   {"py2c_get_gamma_constant", py2c_get_gamma_constant, METH_VARARGS,
    "Returns the gamma constant used by the solve functions.\n Two integer parameters are expected on input:\n (1) for the precision, 1, 2, or 4, respectively for double,\n double double, or quad double precision; and\n (2) the verbose level.\n The function returns a tuple of two doubles,\n for the real and imaginary part of the gamma constant."},
   {"py2c_mixed_volume", py2c_mixed_volume, METH_VARARGS,
    "Computes the mixed volume, and the stable mixed volume as well if\n the input parameter equals 1.  On return is the mixed volume, or\n a tuple with the mixed volume and the stable mixed volume.\n A regular mixed-cell configuration is in the cells container."},
   {"py2c_mixed_volume_by_demics", py2c_mixed_volume_by_demics, METH_VARARGS,
    "Calls DEMiCs to compute the mixed volume of the system in the\n standard systems container.  If the standard systems container\n is empty, then the system in the standard Laurent systems\n container is taken as input.\n The integer in mv on return equals the mixed volume.\n The regular mixed-cell configuration is in the cells container.\n The above is for the case if the input parameter equals 0.\n If the input parameter equals 1, then on return is a tuple,\n which contains the mixed volume and the stable mixed volume."},
   {"py2c_standard_deflate", py2c_standard_deflate, METH_VARARGS,
    "Applies deflation in standard double precision to the system and\n the solutions stored in the containers.\n There are five input parameters, two integers and three floats:\n (1) maxitr : the maximum number of iterations per root,\n (2) maxdef : the maximum number of deflations per root,\n (3) tolerr : tolerance on the forward error on each root,\n (4) tolres : tolerance on the backward error on each root,\n (5) tolres : tolerance on the numerical rank of the Jacobian matrices.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_dobldobl_deflate", py2c_dobldobl_deflate, METH_VARARGS,
    "Applies deflation in double double precision to the system and\n the solutions stored in the containers.\n There are five input parameters, two integers and three floats:\n (1) maxitr : the maximum number of iterations per root,\n (2) maxdef : the maximum number of deflations per root,\n (3) tolerr : tolerance on the forward error on each root,\n (4) tolres : tolerance on the backward error on each root,\n (5) tolres : tolerance on the numerical rank of the Jacobian matrices.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_quaddobl_deflate", py2c_quaddobl_deflate, METH_VARARGS,
    "Applies deflation in quad double precision to the system and\n the solutions stored in the containers.\n (1) maxitr : the maximum number of iterations per root,\n (2) maxdef : the maximum number of deflations per root,\n (3) tolerr : tolerance on the forward error on each root,\n (4) tolres : tolerance on the backward error on each root,\n (5) tolres : tolerance on the numerical rank of the Jacobian matrices.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_standard_Newton_step", py2c_standard_Newton_step, METH_VARARGS,
    "Applies one Newton step in standard double precision to the system in\n the standard systems container and to the solutions in the container.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_dobldobl_Newton_step", py2c_dobldobl_Newton_step, METH_VARARGS,
    "Applies one Newton step in double double precision to the system in\n the standard systems container and to the solutions in the container.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_quaddobl_Newton_step", py2c_quaddobl_Newton_step, METH_VARARGS,
    "Applies one Newton step in quad double precision to the system in\n the standard systems container and to the solutions in the container.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_multprec_Newton_step", py2c_multprec_Newton_step, METH_VARARGS,
    "Applies one Newton step in arbitrary multiprecision to the system in\n the multprec systems container and to the solutions in the container.\n On input is an integer, the number of decimal places in the precision.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_standard_Newton_Laurent_step", py2c_standard_Newton_Laurent_step,
     METH_VARARGS, 
    "Applies one Newton step in standard double precision to the Laurent\n system in the standard Laurent systems container and to the solutions\n in the container.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_dobldobl_Newton_Laurent_step", py2c_dobldobl_Newton_Laurent_step,
     METH_VARARGS,
    "Applies one Newton step in double double precision to the Laurent\n system in the standard Laurent systems container and to the solutions\n in the container.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_quaddobl_Newton_Laurent_step", py2c_quaddobl_Newton_Laurent_step,
     METH_VARARGS, 
    "Applies one Newton step in quad double precision to the Laurent\n system in the standard Laurent systems container and to the solutions\n in the container.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_multprec_Newton_Laurent_step", py2c_multprec_Newton_Laurent_step,
     METH_VARARGS,
    "Applies one Newton step in arbitrary multiprecision to the Laurent\n system in the multprec Laurent systems container and to the solutions\n in the container.\n On input is an integer: the number of decimal places in the precision.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_varbprec_Newton_Laurent_steps", py2c_varbprec_Newton_Laurent_steps,
     METH_VARARGS,
    "Applies Newton's method in variable precision.\n There are six input parameters:\n 1) the dimension: the number of variables and equations;\n 2) the accuracy, expressed as the correct number of decimal places;\n 3) the maximum number of iterations in Newton's method;\n 4) an upper bound on the number of decimal places in the precision;\n 5) a string, with the representation of the polynomials in the system.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_standard_condition_report", py2c_standard_condition_report,
     METH_VARARGS,
    "For the system and solutions in the containers in double precision,\n computes a condition report.  On input are the following:\n 1) maximum number of Newton iterations per solution;\n 2) tolerance on the residual;\n 3) tolerance on the forward error;\n 4) tolerance on the inverse condition number for singularities;\n 5) a string with the name of the output file,\n this string may be empty if no output to file is needed;\n 6) a verbose flag, either 1 or 0.\n On return are the counts of number of solutions that are\n regular, singular, real, complex, clustered, or failures;\n along with the frequency tables for the forward errors,\n residuals and estimates for the inverse condition numbers."},
   {"py2c_usolve_standard", py2c_usolve_standard, METH_VARARGS,
    "Applies the method of Weierstrass to compute all roots of a\n polynomial in one variable with standard double precision arithmetic.\n On input are two numbers:\n 1) the maximum number of iterations in the method of Weierstrass; and\n 2) the epsilon requirement on the accuracy of the roots.\n Before calling this function, the polynomial should be stored in\n the standard systems container.  After the call of this function,\n the standard solutions container contains the roots of the polynomial.\n On return is the number of iterations done by the solver."},
   {"py2c_usolve_dobldobl", py2c_usolve_dobldobl, METH_VARARGS,
    "Applies the method of Weierstrass to compute all roots of a\n polynomial in one variable with double double precision arithmetic.\n On input are two numbers:\n 1) the maximum number of iterations in the method of Weierstrass; and\n 2) the epsilon requirement on the accuracy of the roots.\n Before calling this function, the polynomial should be stored in\n the dobldobl systems container.  After the call of this function,\n the dobldobl solutions container contains the roots of the polynomial.\n On return is the number of iterations done by the solver."},
   {"py2c_usolve_quaddobl", py2c_usolve_quaddobl, METH_VARARGS,
    "Applies the method of Weierstrass to compute all roots of a\n polynomial in one variable with quad double precision arithmetic.\n On input are two numbers:\n 1) the maximum number of iterations in the method of Weierstrass; and\n 2) the epsilon requirement on the accuracy of the roots.\n Before calling this function, the polynomial should be stored in\n the quaddobl systems container.  After the call of this function,\n the quaddobl solutions container contains the roots of the polynomial.\n On return is the number of iterations done by the solver."},
   {"py2c_usolve_multprec", py2c_usolve_multprec, METH_VARARGS,
    "Applies the method of Weierstrass to compute all roots of a\n polynomial in one variable with arbitrary multiprecision arithmetic.\n On input are three numbers:\n 1) the number of decimal places in the working precision;\n 2) the maximum number of iterations in the method of Weierstrass; and\n 3) the epsilon requirement on the accuracy of the roots.\n Before calling this function, the polynomial should be stored in\n the multprec systems container.  After the call of this function,\n the multprec solutions container contains the roots of the polynomial.\n On return is the number of iterations done by the solver."},
   {"py2c_giftwrap_planar", py2c_giftwrap_planar, METH_VARARGS,
    "Applies the giftwrapping algorithm to a planar point configuration.\n On input are an integer and a string:\n 1) the number of points in the list;\n 2) the string representation of a Python list of tuples.\n On return is the string representation of the vertex points,\n sorted so that each two consecutive points define an edge."},
   {"py2c_giftwrap_convex_hull", py2c_giftwrap_convex_hull, METH_VARARGS,
    "Applies the giftwrapping algorithm to a point configuration.\n On input are an integer and a string:\n 1) the number of points in the list;\n 2) the string representation of a Python list of tuples.\n When the function returns, the internal data structures\n to store the convex hull are defined.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_giftwrap_number_of_facets", py2c_giftwrap_number_of_facets,
    METH_VARARGS,
    "Returns the number of facets of the given dimension.\n On input is an integer, the dimension of the facet."},
   {"py2c_giftwrap_retrieve_facet", py2c_giftwrap_retrieve_facet, METH_VARARGS,
    "Returns the string representation of a facet.\n On input are two integer numbers:\n 1) the dimension of the facet;\n 2) the index of the facet."},
   {"py2c_giftwrap_clear_3d_facets", py2c_giftwrap_clear_3d_facets,
    METH_VARARGS,
    "Deallocates list of facets of convex hull stored in 3-space."},
   {"py2c_giftwrap_clear_4d_facets", py2c_giftwrap_clear_4d_facets,
    METH_VARARGS,
    "Deallocates list of facets of convex hull stored in 4-space."},
   {"py2c_giftwrap_support_size", py2c_giftwrap_support_size, METH_VARARGS,
    "Returns the number of characters in the string representation of\n the support of the k-th Laurent polynomial in the container, where k is given on input."},
   {"py2c_giftwrap_support_string", py2c_giftwrap_support_string,
    METH_VARARGS,
    "Returns the string representation of the support of a Laurent polynomial."},
   {"py2c_giftwrap_clear_support_string", py2c_giftwrap_clear_support_string,
    METH_VARARGS,
    "Deallocates the string representation of the support set\n that was stored internally by the call py2c_giftwrap_support_size."},
   {"py2c_giftwrap_initial_form", py2c_giftwrap_initial_form, METH_VARARGS,
    "Replaces the system in the Laurent systems container by its initial form.\n There are three input parameters:\n 1) the dimension, number of coordinates in the inner normal;\n 2) the number of characters in the string representation for the normal;\n 3) the string representation of the inner normal.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syscon_read_standard_system",
     py2c_syscon_read_standard_system, METH_VARARGS,
    "Interactive procedure to read a polynomial system with coefficients\n in standard double precision.\n The system will be placed in the standard systems container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_read_standard_Laurent_system",
     py2c_syscon_read_standard_Laurent_system, METH_VARARGS,
    "Interactive procedure to read a Laurent polynomial system with\n coefficients in standard double precision.\n The system will be placed in the standard Laurent systems container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_read_dobldobl_system", py2c_syscon_read_dobldobl_system,
    METH_VARARGS,
    "Interactive procedure to read a polynomial system with coefficients\n in double double precision.\n The system will be placed in the dobldobl systems container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_read_dobldobl_Laurent_system",
     py2c_syscon_read_dobldobl_Laurent_system, METH_VARARGS, 
    "Interactive procedure to read a Laurent polynomial system with\n coefficients in double double precision.\n The system will be placed in the dobldobl Laurent systems container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_read_quaddobl_system", py2c_syscon_read_quaddobl_system,
    METH_VARARGS, 
    "Interactive procedure to read a polynomial system with coefficients\n in quad double precision.\n The system will be placed in the quaddobl systems container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_read_quaddobl_Laurent_system",
     py2c_syscon_read_quaddobl_Laurent_system, METH_VARARGS,
    "Interactive procedure to read a Laurent polynomial system with\n coefficients in quad double precision.\n The system will be placed in the quaddobl Laurent systems container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_read_multprec_system", py2c_syscon_read_multprec_system,
     METH_VARARGS,
    "Interactive procedure to read a polynomial system with coefficients\n in arbitrary multiprecision.  The one input parameter is an integer,\n the number of decimal places in the working precision.\n The system will be placed in the multprec systems container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_read_multprec_Laurent_system",
     py2c_syscon_read_multprec_Laurent_system, METH_VARARGS,
    "Interactive procedure to read a Laurent polynomial system with\n coefficients in arbitrary multiprecision.  The one input parameter is\n an integer, the number of decimal places in the working precision.\n The system will be placed in the multprec Laurent systems container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_random_system", py2c_syscon_random_system, METH_VARARGS,
    "Places in the systems container a random polynomial system\n with coefficients in standard double precision.\n There are five integers as input parameters:\n 1) n, the number of polynomials and variables;\n 2) m, the number of monomials per equation;\n 3) d, the largest degree of each monomial;\n 4) c, the type of coefficient: 0 if on the complex unit circle,\n 1, if all coefficients are one, 2, if all coefficients are\n random floats in [-1,+1];\n 5) neq, the number of polynomials in the system."}, 
   {"py2c_syscon_dobldobl_random_system",
     py2c_syscon_dobldobl_random_system, METH_VARARGS,
    "Places in the systems container a random polynomial system\n with coefficients in double double precision.\n There are five integers as input parameters:\n 1) n, the number of polynomials and variables;\n 2) m, the number of monomials per equation;\n 3) d, the largest degree of each monomial;\n 4) c, the type of coefficient: 0 if on the complex unit circle,\n 1, if all coefficients are one, 2, if all coefficients are\n random floats in [-1,+1];\n 5) neq, the number of polynomials in the system."}, 
   {"py2c_syscon_quaddobl_random_system",
     py2c_syscon_quaddobl_random_system, METH_VARARGS,
    "Places in the systems container a random polynomial system\n with coefficients in quad double precision.\n There are five integers as input parameters:\n 1) n, the number of polynomials and variables;\n 2) m, the number of monomials per equation;\n 3) d, the largest degree of each monomial;\n 4) c, the type of coefficient: 0 if on the complex unit circle,\n 1, if all coefficients are one, 2, if all coefficients are\n random floats in [-1,+1];\n 5) neq, the number of polynomials in the system."}, 
   {"py2c_syscon_write_standard_system",
     py2c_syscon_write_standard_system, METH_VARARGS,
    "Writes the polynomial system with standard double precision coefficients\n that is stored in the container."},
   {"py2c_syscon_write_standard_Laurent_system",
     py2c_syscon_write_standard_Laurent_system, METH_VARARGS,
    "Writes the Laurent polynomial system with standard double precision\n coefficients that is stored in the container."},
   {"py2c_syscon_write_dobldobl_system",
     py2c_syscon_write_dobldobl_system, METH_VARARGS,
    "Writes the polynomial system with double double precision coefficients\n that is stored in the container."},
   {"py2c_syscon_write_dobldobl_Laurent_system",
     py2c_syscon_write_dobldobl_Laurent_system, METH_VARARGS,
    "Writes the Laurent polynomial system with double double precision\n coefficients that is stored in the container."},
   {"py2c_syscon_write_quaddobl_system",
     py2c_syscon_write_quaddobl_system, METH_VARARGS,
    "Writes the polynomial system with quad double precision coefficients\n that is stored in the container."},
   {"py2c_syscon_write_quaddobl_Laurent_system",
     py2c_syscon_write_quaddobl_Laurent_system, METH_VARARGS,
    "Writes the Laurent polynomial system with quad double precision\n coefficients that is stored in the container."},
   {"py2c_syscon_write_multprec_system",
     py2c_syscon_write_multprec_system, METH_VARARGS, 
    "Writes the polynomial system with arbitrary multiprecision coefficients\n that is stored in the container."},
   {"py2c_syscon_write_multprec_Laurent_system",
     py2c_syscon_write_multprec_Laurent_system, METH_VARARGS,
    "Writes the Laurent polynomial system with arbitrary multiprecision\n coefficients that is stored in the container."},
   {"py2c_syscon_clear_standard_system",
     py2c_syscon_clear_standard_system, METH_VARARGS, 
    "Deallocates the container for polynomial systems\n with coefficients in standard double precision."},
   {"py2c_syscon_clear_standard_Laurent_system",
     py2c_syscon_clear_standard_Laurent_system, METH_VARARGS,
    "Deallocates the container for Laurent polynomial systems\n with coefficients in standard double precision."},
   {"py2c_syscon_clear_dobldobl_system",
     py2c_syscon_clear_dobldobl_system, METH_VARARGS, 
    "Deallocates the container for polynomial systems\n with coefficients in double double precision."},
   {"py2c_syscon_clear_dobldobl_Laurent_system",
     py2c_syscon_clear_dobldobl_Laurent_system, METH_VARARGS,
    "Deallocates the container for Laurent polynomial systems\n with coefficients in double double precision."},
   {"py2c_syscon_clear_quaddobl_system",
     py2c_syscon_clear_quaddobl_system, METH_VARARGS, 
    "Deallocates the container for polynomial systems\n with coefficients in quad double precision."},
   {"py2c_syscon_clear_quaddobl_Laurent_system",
     py2c_syscon_clear_quaddobl_Laurent_system, METH_VARARGS, 
    "Deallocates the container for Laurent polynomial systems\n with coefficients in quad double precision."},
   {"py2c_syscon_clear_multprec_system",
     py2c_syscon_clear_multprec_system, METH_VARARGS,
    "Deallocates the container for polynomial systems\n with coefficients in arbitrary multiprecision."},
   {"py2c_syscon_clear_multprec_Laurent_system",
     py2c_syscon_clear_multprec_Laurent_system, METH_VARARGS,
    "Deallocates the container for Laurent polynomial systems\n with coefficients in arbitrary multiprecision."},
   {"py2c_syscon_number_of_symbols",
     py2c_syscon_number_of_symbols, METH_VARARGS,
    "Returns the number of symbols in the symbol table."},
   {"py2c_syscon_write_symbols", py2c_syscon_write_symbols, METH_VARARGS, 
    "Writes the symbols in the symbol table to screen.\n Returns the failure code, which equals zero if all went well."},
   {"py2c_syscon_string_of_symbols",
     py2c_syscon_string_of_symbols, METH_VARARGS,
    "Returns a string that contains the symbols in the symbol table.\n The symbols are separate from each other by one space."},
   {"py2c_syscon_remove_symbol_name", py2c_syscon_remove_symbol_name,
     METH_VARARGS,
    "Removes a symbol, given by name, from the symbol table.\n On input are two arguments:\n 1) an integer, as the number of characters in the name;\n 2) a string of characters with the name of the symbol.\n The failure code is returned, which equals zero when all went well."},
   {"py2c_syscon_clear_symbol_table", py2c_syscon_clear_symbol_table,
     METH_VARARGS, "Clears the symbol table."},
   {"py2c_solcon_read_standard_solutions",
     py2c_solcon_read_standard_solutions, METH_VARARGS,
    "Interactive function to read the solutions into the container,\n in standard double precision.\n Returns the failure code, which is zero when all went well."},
   {"py2c_solcon_read_dobldobl_solutions", py2c_solcon_read_dobldobl_solutions,
     METH_VARARGS,
    "Interactive function to read the solutions into the container,\n in double double precision.\n Returns the failure code, which is zero when all went well."},
   {"py2c_solcon_read_quaddobl_solutions", py2c_solcon_read_quaddobl_solutions,
     METH_VARARGS,
    "Interactive function to read the solutions into the container,\n in quad double precision.\n Returns the failure code, which is zero when all went well."},
   {"py2c_solcon_read_multprec_solutions", py2c_solcon_read_multprec_solutions,
     METH_VARARGS,
    "Interactive function to read the solutions into the container,\n in arbitrary multiprecision.\n Returns the failure code, which is zero when all went well."},
   {"py2c_solcon_read_standard_solutions_from_file",
     py2c_solcon_read_standard_solutions_from_file, METH_VARARGS,
    "The two input arguments are a number and a string:\n 1) The number equals the number of characters in the string.\n 2) The string given on input is the name of a file which contains\n a solution list to be parsed in standard double precision.\n Solutions are read from file and stored in the container for\n double precision solutions.\n The failure code is returned, which is zero if all went well."},
   {"py2c_solcon_read_dobldobl_solutions_from_file",
     py2c_solcon_read_dobldobl_solutions_from_file, METH_VARARGS,
    "The two input arguments are a number and a string:\n 1) The number equals the number of characters in the string.\n 2) The string given on input is the name of a file which contains\n a solution list to be parsed in double double precision.\n Solutions are read from file and stored in the container for\n double double precision solutions.\n The failure code is returned, which is zero if all went well."},
   {"py2c_solcon_read_quaddobl_solutions_from_file",
     py2c_solcon_read_quaddobl_solutions_from_file, METH_VARARGS,
    "The two input arguments are a number and a string:\n 1) The number equals the number of characters in the string.\n 2) The string given on input is the name of a file which contains\n a solution list to be parsed in quad double precision.\n Solutions are read from file and stored in the container for\n quad double precision solutions.\n The failure code is returned, which is zero if all went well."},
   {"py2c_solcon_write_standard_solutions",
     py2c_solcon_write_standard_solutions, METH_VARARGS,
    "Writes the solutions in standard double precision to screen.\n Returns the failure code, which equals zero when all went well."},
   {"py2c_solcon_write_dobldobl_solutions",
     py2c_solcon_write_dobldobl_solutions, METH_VARARGS,
    "Writes the solutions in double double precision to screen.\n Returns the failure code, which equals zero when all went well."},
   {"py2c_solcon_write_quaddobl_solutions",
     py2c_solcon_write_quaddobl_solutions, METH_VARARGS,
    "Writes the solutions in quad double precision to screen.\n Returns the failure code, which equals zero when all went well."},
   {"py2c_solcon_write_multprec_solutions",
     py2c_solcon_write_multprec_solutions, METH_VARARGS,
    "Writes the solutions in arbitrary multiprecision to screen.\n Returns the failure code, which equals zero when all went well."},
   {"py2c_solcon_clear_standard_solutions",
     py2c_solcon_clear_standard_solutions, METH_VARARGS,
    "Deallocates the container for solutions in standard double precision.\n Returns the failure code, which equals zero when all went well."},
   {"py2c_solcon_clear_dobldobl_solutions",
     py2c_solcon_clear_dobldobl_solutions, METH_VARARGS,
    "Deallocates the container for solutions in double double precision.\n Returns the failure code, which equals zero when all went well."},
   {"py2c_solcon_clear_quaddobl_solutions",
     py2c_solcon_clear_quaddobl_solutions, METH_VARARGS,
    "Deallocates the container for solutions in quad double precision.\n Returns the failure code, which equals zero when all went well."},
   {"py2c_solcon_clear_multprec_solutions",
     py2c_solcon_clear_multprec_solutions, METH_VARARGS,
    "Deallocates the container for solutions in arbitrary multiprecision.\n Returns the failure code, which equals zero when all went well."},
   {"py2c_solcon_open_solution_input_file",
     py2c_solcon_open_solution_input_file, METH_VARARGS,
    "Prompts the user for the name of the input file for the solutions and\n opens the input file.  All subsequent reading happens from this input.\n Returns the failure code, which equals zero when all went well."},
   {"py2c_syscon_number_of_standard_polynomials",
     py2c_syscon_number_of_standard_polynomials,
     METH_VARARGS, 
    "Returns the number of polynomials with coefficients in standard\n double precision as stored in the systems container."},
   {"py2c_syscon_number_of_dobldobl_polynomials",
     py2c_syscon_number_of_dobldobl_polynomials, METH_VARARGS, 
    "Returns the number of polynomials with coefficients in double\n double precision as stored in the systems container."},
   {"py2c_syscon_number_of_quaddobl_polynomials",
     py2c_syscon_number_of_quaddobl_polynomials, METH_VARARGS, 
    "Returns the number of polynomials with coefficients in quad\n double precision as stored in the systems container."},
   {"py2c_syscon_number_of_multprec_polynomials",
     py2c_syscon_number_of_multprec_polynomials, METH_VARARGS, 
    "Returns the number of polynomials with coefficients in arbitrary\n multiprecision as stored in the systems container."},
   {"py2c_syscon_number_of_standard_Laurentials",
     py2c_syscon_number_of_standard_Laurentials, METH_VARARGS, 
    "Returns the number of Laurent polynomials with coefficients in\n standard double precision as stored in the systems container."},
   {"py2c_syscon_number_of_dobldobl_Laurentials",
     py2c_syscon_number_of_dobldobl_Laurentials, METH_VARARGS, 
    "Returns the number of Laurent polynomials with coefficients in\n double double precision as stored in the systems container."},
   {"py2c_syscon_number_of_quaddobl_Laurentials",
     py2c_syscon_number_of_quaddobl_Laurentials, METH_VARARGS, 
    "Returns the number of Laurent polynomials with coefficients in\n quad double precision as stored in the systems container."},
   {"py2c_syscon_number_of_multprec_Laurentials",
     py2c_syscon_number_of_multprec_Laurentials, METH_VARARGS, 
    "Returns the number of Laurent polynomials with coefficients in\n arbitrary multiprecision as stored in the systems container."},
   {"py2c_syscon_initialize_number_of_standard_polynomials",
     py2c_syscon_initialize_number_of_standard_polynomials, METH_VARARGS,
    "Initializes the container for polynomials with coefficients in\n standard double precision.  The input argument is an integer,\n the number of polynomials in the container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_initialize_number_of_dobldobl_polynomials",
     py2c_syscon_initialize_number_of_dobldobl_polynomials, METH_VARARGS,
    "Initializes the container for polynomials with coefficients in\n double double precision.  The input argument is an integer,\n the number of polynomials in the container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_initialize_number_of_quaddobl_polynomials",
     py2c_syscon_initialize_number_of_quaddobl_polynomials, METH_VARARGS,
    "Initializes the container for polynomials with coefficients in\n quad double precision.  The input argument is an integer,\n the number of polynomials in the container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_initialize_number_of_multprec_polynomials",
     py2c_syscon_initialize_number_of_multprec_polynomials, METH_VARARGS,
    "Initializes the container for polynomials with coefficients in\n arbitrary multiprecision.  The input argument is an integer,\n the number of polynomials in the container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_initialize_number_of_standard_Laurentials",
     py2c_syscon_initialize_number_of_standard_Laurentials, METH_VARARGS,
    "Initializes the container for Laurent polynomials with coefficients\n in standard double precision.  The input argument is an integer,\n the number of polynomials in the container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_initialize_number_of_dobldobl_Laurentials",
     py2c_syscon_initialize_number_of_dobldobl_Laurentials, METH_VARARGS,
    "Initializes the container for Laurent polynomials with coefficients\n in double double precision.  The input argument is an integer,\n the number of polynomials in the container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_initialize_number_of_quaddobl_Laurentials",
     py2c_syscon_initialize_number_of_quaddobl_Laurentials, METH_VARARGS,
    "Initializes the container for Laurent polynomials with coefficients\n in quad double precision.  The input argument is an integer,\n the number of polynomials in the container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_initialize_number_of_multprec_Laurentials",
     py2c_syscon_initialize_number_of_multprec_Laurentials, METH_VARARGS,
    "Initializes the container for Laurent polynomials with coefficients\n in arbitrary multiprecision.  The input argument is an integer,\n the number of polynomials in the container.\n The failure code is returned, which equals zero if all went well."},
   {"py2c_syscon_degree_of_standard_polynomial",
     py2c_syscon_degree_of_standard_polynomial, METH_VARARGS,
    "Returns the degree of the k-th polynomial in the container for\n polynomials with coefficients in standard double precision.\n The index k of the polynomial is the one input argument."},
   {"py2c_syscon_degree_of_dobldobl_polynomial",
     py2c_syscon_degree_of_dobldobl_polynomial, METH_VARARGS,
    "Returns the degree of the k-th polynomial in the container for\n polynomials with coefficients in double double precision.\n The index k of the polynomial is the one input argument."},
   {"py2c_syscon_degree_of_quaddobl_polynomial",
     py2c_syscon_degree_of_quaddobl_polynomial, METH_VARARGS,
    "Returns the degree of the k-th polynomial in the container for\n polynomials with coefficients in quad double precision.\n The index k of the polynomial is the one input argument."},
   {"py2c_syscon_degree_of_multprec_polynomial",
     py2c_syscon_degree_of_multprec_polynomial, METH_VARARGS,
    "Returns the degree of the k-th polynomial in the container for\n polynomials with coefficients in arbitrary multiprecision.\n The index k of the polynomial is the one input argument."},
   {"py2c_syscon_number_of_terms", 
     py2c_syscon_number_of_terms, METH_VARARGS,
    "Returns the number of terms in the k-th polynomial stored in the\n container for systems with coefficients in standard double precision.\n The input parameter k is the index of the polynomial k."},
   {"py2c_syscon_number_of_Laurent_terms",
     py2c_syscon_number_of_Laurent_terms, METH_VARARGS,
    "Returns the number of terms in the k-th Laurent polynomial stored\n in the container for Laurent polynomials systems with coefficients\n in standard double precision.\n The input parameter k is the index of the polynomial k."},
   {"py2c_syscon_retrieve_term",
     py2c_syscon_retrieve_term, METH_VARARGS,
    "Retrieves one term of a polynomial with coefficients in standard\n double precision, that is stored in the systems container.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syscon_store_standard_polynomial",
     py2c_syscon_store_standard_polynomial, METH_VARARGS, 
    "Defines the k-th polynomial in the systems container for polynomials\n with coefficients in standard double precision.\n As a precondition for this function, the container must be initialized\n for sufficiently many polynomials, in any case >= k.\n There are four input parameters, three integers and one string:\n 1) nc, the number of characters in the string p;\n 2) n, the number of variables in the multivariate polynomial;\n 3) k, the index of the polynomial in the system;\n 4) p, a valid string representation for a polynomial.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_syscon_store_dobldobl_polynomial",
     py2c_syscon_store_dobldobl_polynomial, METH_VARARGS,
    "Defines the k-th polynomial in the systems container for polynomials\n with coefficients in double double precision.\n As a precondition for this function, the container must be initialized\n for sufficiently many polynomials, in any case >= k.\n There are four input parameters, three integers and one string:\n 1) nc, the number of characters in the string p;\n 2) n, the number of variables in the multivariate polynomial;\n 3) k, the index of the polynomial in the system;\n 4) p, a valid string representation for a polynomial.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_syscon_store_quaddobl_polynomial",
     py2c_syscon_store_quaddobl_polynomial, METH_VARARGS,
    "Defines the k-th polynomial in the systems container for polynomials\n with coefficients in quad double precision.\n As a precondition for this function, the container must be initialized\n for sufficiently many polynomials, in any case >= k.\n There are four input parameters, three integers and one string:\n 1) nc, the number of characters in the string p;\n 2) n, the number of variables in the multivariate polynomial;\n 3) k, the index of the polynomial in the system;\n 4) p, a valid string representation for a polynomial.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_syscon_store_multprec_polynomial",
     py2c_syscon_store_multprec_polynomial, METH_VARARGS,
    "Defines the k-th polynomial in the systems container for polynomials\n with coefficients in arbitrary multiprecision.\n As a precondition for this function, the container must be initialized\n for sufficiently many polynomials, in any case >= k.\n There are five input parameters, four integers and one string:\n 1) nc, the number of characters in the string p;\n 2) n, the number of variables in the multivariate polynomial;\n 3) k, the index of the polynomial in the system;\n 4) dp, the number of decimal places to parse the coefficients;\n 5) p, a valid string representation for a polynomial.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_syscon_load_standard_polynomial",
     py2c_syscon_load_standard_polynomial, METH_VARARGS,
    "Returns the k-th polynomial in the systems container\n with standard double complex coefficients as a string.\n The value for k is in the one integer parameter of this function."},
   {"py2c_syscon_load_dobldobl_polynomial",
     py2c_syscon_load_dobldobl_polynomial, METH_VARARGS,
    "Returns the k-th polynomial in the systems container\n with double double complex coefficients as a string.\n The value for k is in the one integer parameter of this function."},
   {"py2c_syscon_load_quaddobl_polynomial",
     py2c_syscon_load_quaddobl_polynomial, METH_VARARGS,
    "Returns the k-th polynomial in the systems container\n with quad double complex coefficients as a string.\n The value for k is in the one integer parameter of this function."},
   {"py2c_syscon_load_multprec_polynomial",
     py2c_syscon_load_multprec_polynomial, METH_VARARGS,
    "Returns the k-th polynomial in the systems container\n with arbitrary multiprecision complex coefficients as a string.\n The value for k is in the one integer parameter of this function."},
   {"py2c_syscon_store_standard_Laurential",
     py2c_syscon_store_standard_Laurential, METH_VARARGS,
    "Defines the k-th polynomial in the systems container for Laurent\n polynomials with coefficients in standard double precision.\n As a precondition for this function, the container must be initialized\n for sufficiently many polynomials, in any case >= k.\n There are four input parameters, three integers and one string:\n 1) nc, the number of characters in the string p;\n 2) n, the number of variables in the multivariate polynomial;\n 3) k, the index of the polynomial in the system;\n 4) p, a valid string representation for a polynomial.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_syscon_store_dobldobl_Laurential",
     py2c_syscon_store_dobldobl_Laurential, METH_VARARGS,
    "Defines the k-th polynomial in the systems container for Laurent\n polynomials with coefficients in double double precision.\n As a precondition for this function, the container must be initialized\n for sufficiently many polynomials, in any case >= k.\n There are four input parameters, three integers and one string:\n 1) nc, the number of characters in the string p;\n 2) n, the number of variables in the multivariate polynomial;\n 3) k, the index of the polynomial in the system;\n 4) p, a valid string representation for a polynomial.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_syscon_store_quaddobl_Laurential",
     py2c_syscon_store_quaddobl_Laurential, METH_VARARGS,
    "Defines the k-th polynomial in the systems container for Laurent\n polynomials with coefficients in quad double precision.\n As a precondition for this function, the container must be initialized\n for sufficiently many polynomials, in any case >= k.\n There are four input parameters, three integers and one string:\n 1) nc, the number of characters in the string p;\n 2) n, the number of variables in the multivariate polynomial;\n 3) k, the index of the polynomial in the system;\n 4) p, a valid string representation for a polynomial.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_syscon_store_multprec_Laurential",
     py2c_syscon_store_multprec_Laurential, METH_VARARGS,
    "Defines the k-th polynomial in the systems container for Laurent\n polynomials with coefficients in arbitrary multiprecision.\n As a precondition for this function, the container must be initialized\n for sufficiently many polynomials, in any case >= k.\n There are five input parameters, four integers and one string:\n 1) nc, the number of characters in the string p;\n 2) n, the number of variables in the multivariate polynomial;\n 3) k, the index of the polynomial in the system;\n 4) dp, the number of decimal places to parse the coefficients;\n 5) p, a valid string representation for a polynomial.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_syscon_load_standard_Laurential",
     py2c_syscon_load_standard_Laurential, METH_VARARGS,
    "Returns the k-th polynomial in the Laurent systems container\n with standard double complex coefficients as a string.\n The value for k is in the one integer parameter of this function."},
   {"py2c_syscon_load_dobldobl_Laurential",
     py2c_syscon_load_dobldobl_Laurential, METH_VARARGS,
    "Returns the k-th polynomial in the Laurent systems container\n with double double complex coefficients as a string.\n The value for k is in the one integer parameter of this function."},
   {"py2c_syscon_load_quaddobl_Laurential",
     py2c_syscon_load_quaddobl_Laurential, METH_VARARGS,
    "Returns the k-th polynomial in the Laurent systems container\n with quad double complex coefficients as a string.\n The value for k is in the one integer parameter of this function."},
   {"py2c_syscon_load_multprec_Laurential",
     py2c_syscon_load_multprec_Laurential, METH_VARARGS,
    "Returns the k-th polynomial in the Laurent systems container\n with arbitrary multiprecision complex coefficients as a string.\n The value for k is in the one integer parameter of this function."},
   {"py2c_syscon_total_degree", py2c_syscon_total_degree, METH_VARARGS,
    "Returns in d the total degree of the system with coefficients in\n standard double precision, as stored in the container."},
   {"py2c_syscon_standard_drop_variable_by_index",
     py2c_syscon_standard_drop_variable_by_index, METH_VARARGS,
    "Replaces the system in the standard double precision container\n with the same system that has its k-th variable dropped.\n The index k of the variable is given as an input parameter.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syscon_standard_drop_variable_by_name",
     py2c_syscon_standard_drop_variable_by_name, METH_VARARGS,
    "Replaces the system in the standard double precision container\n with the same system that have that variable dropped\n corresponding to the name in the string s of nc characters long.\n The function has two input parameters, an integer and a string:\n 1) nc, the number of characters in the string with the name;\n 2) s, a string that holds the name of the variable.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syscon_dobldobl_drop_variable_by_index",
     py2c_syscon_dobldobl_drop_variable_by_index, METH_VARARGS,
    "Replaces the system in the double double precision container\n with the same system that has its k-th variable dropped.\n The index k of the variable is given as an input parameter.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syscon_dobldobl_drop_variable_by_name",
     py2c_syscon_dobldobl_drop_variable_by_name, METH_VARARGS,
    "Replaces the system in the double double precision container\n with the same system that have that variable dropped\n corresponding to the name in the string s of nc characters long.\n The function has two input parameters, an integer and a string:\n 1) nc, the number of characters in the string with the name;\n 2) s, a string that holds the name of the variable.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syscon_quaddobl_drop_variable_by_index",
     py2c_syscon_quaddobl_drop_variable_by_index, METH_VARARGS,
    "Replaces the system in the quad double precision container\n with the same system that has its k-th variable dropped.\n The index k of the variable is given as an input parameter.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syscon_quaddobl_drop_variable_by_name",
     py2c_syscon_quaddobl_drop_variable_by_name, METH_VARARGS,
    "Replaces the system in the quad double precision container\n with the same system that have that variable dropped\n corresponding to the name in the string s of nc characters long.\n The function has two input parameters, an integer and a string:\n 1) nc, the number of characters in the string with the name;\n 2) s, a string that holds the name of the variable.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syscon_standard_Laurent_drop_variable_by_index",
     py2c_syscon_standard_Laurent_drop_variable_by_index, METH_VARARGS,
    "Replaces the Laurent system in the standard double precision container\n with the same Laurent system that has its k-th variable dropped.\n The index k of the variable is given as an input parameter.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syscon_standard_Laurent_drop_variable_by_name",
     py2c_syscon_standard_Laurent_drop_variable_by_name, METH_VARARGS,
    "Replaces the Laurent system in the standard double precision container\n with the same Laurent system that have that variable dropped\n corresponding to the name in the string s of nc characters long.\n The function has two input parameters, an integer and a string:\n 1) nc, the number of characters in the string with the name;\n 2) s, a string that holds the name of the variable.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syscon_dobldobl_Laurent_drop_variable_by_index",
     py2c_syscon_dobldobl_Laurent_drop_variable_by_index, METH_VARARGS,
    "Replaces the Laurent system in the double double precision container\n with the same Laurent system that has its k-th variable dropped.\n The index k of the variable is given as an input parameter.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syscon_dobldobl_Laurent_drop_variable_by_name",
     py2c_syscon_dobldobl_Laurent_drop_variable_by_name, METH_VARARGS,
    "Replaces the Laurent system in the double double precision container\n with the same Laurent system that have that variable dropped\n corresponding to the name in the string s of nc characters long.\n The function has two input parameters, an integer and a string:\n 1) nc, the number of characters in the string with the name;\n 2) s, a string that holds the name of the variable.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syscon_quaddobl_Laurent_drop_variable_by_index",
     py2c_syscon_quaddobl_Laurent_drop_variable_by_index, METH_VARARGS,
    "Replaces the Laurent system in the quad double precision container\n with the same Laurent system that has its k-th variable dropped.\n The index k of the variable is given as an input parameter.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syscon_quaddobl_Laurent_drop_variable_by_name",
     py2c_syscon_quaddobl_Laurent_drop_variable_by_name, METH_VARARGS,
    "Replaces the Laurent system in the quad double precision container\n with the same Laurent system that have that variable dropped\n corresponding to the name in the string s of nc characters long.\n The function has two input parameters, an integer and a string:\n 1) nc, the number of characters in the string with the name;\n 2) s, a string that holds the name of the variable.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syscon_standard_one_homogenization",
     py2c_syscon_standard_one_homogenization, METH_VARARGS,
    "Replaces the system in the standard double precision container\n with its transformation in 1-homogeneous coordinates.\n There is one integer on input.\n If 0, then a random linear equation is added,\n otherwise, the linear equation z0 - 1 = 0 is added,\n where z0 is the extra homogeneous coordinate."},
   {"py2c_syscon_dobldobl_one_homogenization",
     py2c_syscon_dobldobl_one_homogenization, METH_VARARGS,
    "Replaces the system in the double double precision container\n with its transformation in 1-homogeneous coordinates.\n There is one integer on input.\n If 0, then a random linear equation is added,\n otherwise, the linear equation z0 - 1 = 0 is added,\n where z0 is the extra homogeneous coordinate."},
   {"py2c_syscon_quaddobl_one_homogenization",
     py2c_syscon_quaddobl_one_homogenization, METH_VARARGS,
    "Replaces the system in the quad double precision container\n with its transformation in 1-homogeneous coordinates.\n There is one integer on input.\n If 0, then a random linear equation is added,\n otherwise, the linear equation z0 - 1 = 0 is added,\n where z0 is the extra homogeneous coordinate."},
   {"py2c_syscon_add_symbol",
     py2c_syscon_add_symbol, METH_VARARGS,
    "Adds a symbol to the table, with name given in the string,\n where the number of characters in the name equals the first\n integer argument.  The second input parameter is the string.\n This symbol represents the last variable added in the homogeneous\n coordinate transformation."},
   {"py2c_syscon_standard_one_affinization",
     py2c_syscon_standard_one_affinization, METH_VARARGS,
    "Replaces the system in the standard double precision container\n by its transformation to affine coordinates, substituting the\n value of the last coordinate by one and removing the last equation."},
   {"py2c_syscon_dobldobl_one_affinization",
     py2c_syscon_dobldobl_one_affinization, METH_VARARGS,
    "Replaces the system in the double double precision container\n by its transformation to affine coordinates, substituting the\n value of the last coordinate by one and removing the last equation."},
   {"py2c_syscon_quaddobl_one_affinization",
     py2c_syscon_quaddobl_one_affinization, METH_VARARGS,
    "Replaces the system in the quad double precision container\n by its transformation to affine coordinates, substituting the\n value of the last coordinate by one and removing the last equation."},
   {"py2c_tabform_store_standard_tableau",
     py2c_tabform_store_standard_tableau, METH_VARARGS,
    "On input is the tableau form of a polynomial system, given by\n 1) the number of equations as an integer,\n 2) the number of equations as an integer,\n 3) the number of characters in the 4-th string input,\n 4) the number of terms in each polynomial, given as a string,\n the string representation of a list of integers,\n 5) the number of characters in the 6-th string input,\n 6) the coefficients of all terms, given as a string,\n the string representation of a list of doubles,\n each pair of consecutive doubles represents a complex coefficient,\n 7) the number of characters in the 7-th string input,\n 8) the exponents of all terms, given as a string,\n the string representation of a list of integers.\n The tableau form is parsed and the container for systems with\n standard double precision coefficients is initialized."},
   {"py2c_tabform_load_standard_tableau",
     py2c_tabform_load_standard_tableau, METH_VARARGS,
    "Returns a 5-tuple with the tableau form of the system with\n standard double precision coefficients in the container.\n The five items in the returned tuple are\n 1) the number of equations as an integer,\n 2) the number of equations as an integer,\n 3) the number of terms in each polynomial, given as a string,\n the string representation of a list of integers,\n 4) the coefficients of all terms, given as a string,\n the string representation of a list of doubles,\n each pair of consecutive doubles represents a complex coefficient,\n 5) the exponents of all terms, given as a string,\n the string representation of a list of integers."},
   {"py2c_solcon_length_standard_solution_string",
     py2c_solcon_length_standard_solution_string,
     METH_VARARGS,
    "On input is the index k to a solution in standard double precision,\n stored in the container.  On return is the length of the string\n representation for that k-th solution in the container."},
   {"py2c_solcon_length_dobldobl_solution_string",
     py2c_solcon_length_dobldobl_solution_string, METH_VARARGS,
    "On input is the index k to a solution in double double precision,\n stored in the container.  On return is the length of the string\n representation for that k-th solution in the container."},
   {"py2c_solcon_length_quaddobl_solution_string",
     py2c_solcon_length_quaddobl_solution_string, METH_VARARGS,
    "On input is the index k to a solution in quad double precision,\n stored in the container.  On return is the length of the string\n representation for that k-th solution in the container."},
   {"py2c_solcon_length_multprec_solution_string",
     py2c_solcon_length_multprec_solution_string, METH_VARARGS,
    "On input is the index k to a solution in arbitrary multiprecision,\n stored in the container.  On return is the length of the string\n representation for that k-th solution in the container."},
   {"py2c_solcon_write_standard_solution_string",
     py2c_solcon_write_standard_solution_string,
     METH_VARARGS,
    "Returns the string representation for the k-th solution stored\n in standard double precision in the container.\n On input are two integers:\n 1) the index to the solution; and\n 2) the number of characters in the string representation\n for that solution."},
   {"py2c_solcon_write_dobldobl_solution_string",
     py2c_solcon_write_dobldobl_solution_string, METH_VARARGS,
    "Returns the string representation for the k-th solution stored\n in double double precision in the container.\n On input are two integers:\n 1) the index to the solution; and\n 2) the number of characters in the string representation\n for that solution."},
   {"py2c_solcon_write_quaddobl_solution_string",
     py2c_solcon_write_quaddobl_solution_string, METH_VARARGS, 
    "Returns the string representation for the k-th solution stored\n in quad double precision in the container.\n On input are two integers:\n 1) the index to the solution; and\n 2) the number of characters in the string representation\n for that solution."},
   {"py2c_solcon_write_multprec_solution_string",
     py2c_solcon_write_multprec_solution_string, METH_VARARGS,
    "Returns the string representation for the k-th solution stored\n in arbitrary multiprecision in the container.\n On input are two integers:\n 1) the index to the solution; and\n 2) the number of characters in the string representation\n for that solution."},
   {"py2c_solcon_retrieve_next_standard_initialize",
     py2c_solcon_retrieve_next_standard_initialize, METH_VARARGS,
    "Resets the pointer to the current standard solution in the container\n to the first solution in the list.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_solcon_retrieve_next_dobldobl_initialize",
     py2c_solcon_retrieve_next_dobldobl_initialize, METH_VARARGS,
    "Resets the pointer to the current dobldobl solution in the container\n to the first solution in the list.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_solcon_retrieve_next_quaddobl_initialize",
     py2c_solcon_retrieve_next_quaddobl_initialize, METH_VARARGS,
    "Resets the pointer to the current quaddobl solution in the container\n to the first solution in the list.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_solcon_retrieve_next_multprec_initialize",
     py2c_solcon_retrieve_next_multprec_initialize, METH_VARARGS,
    "Resets the pointer to the current multprec solution in the container\n to the first solution in the list.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_solcon_move_current_standard_to_next",
     py2c_solcon_move_current_standard_to_next, METH_VARARGS,
    "Moves the pointer to the current solution in standard double precision\n to the next solution and returns the value of the cursor.\n If cursor on return is zero, then either the pointer was null\n or there is no next solution."},
   {"py2c_solcon_move_current_dobldobl_to_next",
     py2c_solcon_move_current_dobldobl_to_next, METH_VARARGS,
    "Moves the pointer to the current solution in double double precision\n to the next solution and returns the value of the cursor.\n If cursor on return is zero, then either the pointer was null\n or there is no next solution."},
   {"py2c_solcon_move_current_quaddobl_to_next",
     py2c_solcon_move_current_quaddobl_to_next, METH_VARARGS,
    "Moves the pointer to the current solution in quad double precision\n to the next solution and returns the value of the cursor.\n If cursor on return is zero, then either the pointer was null\n or there is no next solution."},
   {"py2c_solcon_move_current_multprec_to_next",
     py2c_solcon_move_current_multprec_to_next, METH_VARARGS,
    "Moves the pointer to the current solution in arbitrary multiprecision\n to the next solution and returns the value of the cursor.\n If cursor on return is zero, then either the pointer was null\n or there is no next solution."},
   {"py2c_solcon_length_current_standard_solution_string",
     py2c_solcon_length_current_standard_solution_string, METH_VARARGS,
    "Returns the number of characters in the string representation\n of the current standard double solution in the container,\n at the place indicated by the value of the cursor.\n If this value equals zero, then there is no current solution,\n and then the length on return equals zero."},
   {"py2c_solcon_length_current_dobldobl_solution_string",
     py2c_solcon_length_current_dobldobl_solution_string, METH_VARARGS,
    "Returns the number of characters in the string representation\n of the current double double solution in the container,\n at the place indicated by the value of the cursor.\n If this value equals zero, then there is no current solution,\n and then the length on return equals zero."},
   {"py2c_solcon_length_current_quaddobl_solution_string",
     py2c_solcon_length_current_quaddobl_solution_string, METH_VARARGS,
    "Returns the number of characters in the string representation\n of the current quad double solution in the container,\n at the place indicated by the value of the cursor.\n If this value equals zero, then there is no current solution,\n and then the length on return equals zero."},
   {"py2c_solcon_length_current_multprec_solution_string",
     py2c_solcon_length_current_multprec_solution_string, METH_VARARGS,
    "Returns the number of characters in the string representation\n of the current arbitrary multiprecision solution in the container,\n at the place indicated by the value of the cursor.\n If this value equals zero, then there is no current solution,\n and then the length on return equals zero."},
   {"py2c_solcon_write_current_standard_solution_string",
     py2c_solcon_write_current_standard_solution_string, METH_VARARGS,
    "Writes the current standard double solution in the solution container\n to the string s of n+1 characters.\n The last character is the end of string symbol.\n The value of n is given as the one input parameter to this function.\n On return is the string that contains the string representation of\n the current solution in standard double precision in the container."},
   {"py2c_solcon_write_current_dobldobl_solution_string",
     py2c_solcon_write_current_dobldobl_solution_string, METH_VARARGS,
    "Writes the current double double solution in the solution container\n to the string s of n+1 characters.\n The last character is the end of string symbol.\n The value of n is given as the one input parameter to this function.\n On return is the string that contains the string representation of\n the current solution in standard double precision in the container."},
   {"py2c_solcon_write_current_quaddobl_solution_string",
     py2c_solcon_write_current_quaddobl_solution_string, METH_VARARGS,
    "Writes the current quad double solution in the solution container\n to the string s of n+1 characters.\n The last character is the end of string symbol.\n The value of n is given as the one input parameter to this function.\n On return is the string that contains the string representation of\n the current solution in standard double precision in the container."},
   {"py2c_solcon_write_current_multprec_solution_string",
     py2c_solcon_write_current_multprec_solution_string, METH_VARARGS,
    "Writes the current arbitrary multiprecision solution in the solution container\n to the string s of n+1 characters.\n The last character is the end of string symbol.\n The value of n is given as the one input parameter to this function.\n On return is the string that contains the string representation of\n the current solution in standard double precision in the container."},
   {"py2c_solcon_append_standard_solution_string",
     py2c_solcon_append_standard_solution_string,
     METH_VARARGS, 
    "Appends a solution in standard double precision to the list\n of solutions already stored in the container.\n There are three input parameters:\n 1) the number of variables;\n 2) the number of characters in the string;\n 3) the string representing the solution to append to the list.\n Returns the failure code, which equals zero if all went well."},
   {"py2c_solcon_append_dobldobl_solution_string",
     py2c_solcon_append_dobldobl_solution_string, METH_VARARGS,
    "Appends a solution in double double precision to the list\n of solutions already stored in the container.\n There are three input parameters:\n 1) the number of variables;\n 2) the number of characters in the string;\n 3) the string representing the solution to append to the list.\n Returns the failure code, which equals zero if all went well."},
   {"py2c_solcon_append_quaddobl_solution_string",
     py2c_solcon_append_quaddobl_solution_string, METH_VARARGS,
    "Appends a solution in quad double precision to the list\n of solutions already stored in the container.\n There are three input parameters:\n 1) the number of variables;\n 2) the number of characters in the string;\n 3) the string representing the solution to append to the list.\n Returns the failure code, which equals zero if all went well."},
   {"py2c_solcon_append_multprec_solution_string",
     py2c_solcon_append_multprec_solution_string, METH_VARARGS,
    "Appends a solution in arbitrary multiprecision to the list\n of solutions already stored in the container.\n There are three input parameters:\n 1) the number of variables;\n 2) the number of characters in the string;\n 3) the string representing the solution to append to the list.\n Returns the failure code, which equals zero if all went well."},
   {"py2c_solcon_number_of_standard_solutions",
     py2c_solcon_number_of_standard_solutions, METH_VARARGS,
    "Returns the number of solutions in standard double precision,\n as stored in the container."},
   {"py2c_solcon_number_of_dobldobl_solutions",
     py2c_solcon_number_of_dobldobl_solutions, METH_VARARGS, 
    "Returns the number of solutions in double double precision,\n as stored in the container."},
   {"py2c_solcon_number_of_quaddobl_solutions",
     py2c_solcon_number_of_quaddobl_solutions, METH_VARARGS, 
    "Returns the number of solutions in quad double precision,\n as stored in the container."},
   {"py2c_solcon_number_of_multprec_solutions",
     py2c_solcon_number_of_multprec_solutions, METH_VARARGS, 
    "Returns the number of solutions in arbitrary multiprecision,\n as stored in the container."},
   {"py2c_solcon_standard_drop_coordinate_by_index",
     py2c_solcon_standard_drop_coordinate_by_index, METH_VARARGS,
    "Replaces the solutions in the standard double precision container\n with the same solutions that have their k-th coordinate dropped.\n There is one input parameter: the index k of the coordinate.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_solcon_standard_drop_coordinate_by_name",
     py2c_solcon_standard_drop_coordinate_by_name, METH_VARARGS,
    "Replaces the solutions in the standard double precision container\n with the same solutions that have their coordinate dropped\n corresponding to the name in the string s of nc characters long.\n There are two input parameters, an integer and a string:\n 1) nc, the number of characters in the string with the name;\n 2) s, the string with the name of the variable.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_solcon_dobldobl_drop_coordinate_by_index",
     py2c_solcon_dobldobl_drop_coordinate_by_index, METH_VARARGS,
    "Replaces the solutions in the double double precision container\n with the same solutions that have their k-th coordinate dropped.\n There is one input parameter: the index k of the coordinate.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_solcon_dobldobl_drop_coordinate_by_name",
     py2c_solcon_dobldobl_drop_coordinate_by_name, METH_VARARGS,
    "Replaces the solutions in the double double precision container\n with the same solutions that have their coordinate dropped\n corresponding to the name in the string s of nc characters long.\n There are two input parameters, an integer and a string:\n 1) nc, the number of characters in the string with the name;\n 2) s, the string with the name of the variable.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_solcon_quaddobl_drop_coordinate_by_index",
     py2c_solcon_quaddobl_drop_coordinate_by_index, METH_VARARGS,
    "Replaces the solutions in the quad double precision container\n with the same solutions that have their k-th coordinate dropped.\n There is one input parameter: the index k of the coordinate.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_solcon_quaddobl_drop_coordinate_by_name",
     py2c_solcon_quaddobl_drop_coordinate_by_name, METH_VARARGS,
    "Replaces the solutions in the quad double precision container\n with the same solutions that have their coordinate dropped\n corresponding to the name in the string s of nc characters long.\n There are two input parameters, an integer and a string:\n 1) nc, the number of characters in the string with the name;\n 2) s, the string with the name of the variable.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_solcon_standard_one_homogenization",
     py2c_solcon_standard_one_homogenization, METH_VARARGS,
    "Add one extra coordinate one to every solution in the container\n for solutions in standard double precision."},
   {"py2c_solcon_dobldobl_one_homogenization",
     py2c_solcon_dobldobl_one_homogenization, METH_VARARGS,
    "Add one extra coordinate one to every solution in the container\n for solutions in double double precision."},
   {"py2c_solcon_quaddobl_one_homogenization",
     py2c_solcon_quaddobl_one_homogenization, METH_VARARGS,
    "Add one extra coordinate one to every solution in the container\n for solutions in double double precision."},
   {"py2c_solcon_standard_one_affinization",
     py2c_solcon_standard_one_affinization, METH_VARARGS,
    "Divides every coordinate by the last coordinate of every solution\n in the container for solutions in standard double precision."},
   {"py2c_solcon_dobldobl_one_affinization",
     py2c_solcon_dobldobl_one_affinization, METH_VARARGS,
    "Divides every coordinate by the last coordinate of every solution\n in the container for solutions in double double precision."},
   {"py2c_solcon_quaddobl_one_affinization",
     py2c_solcon_quaddobl_one_affinization, METH_VARARGS,
    "Divides every coordinate by the last coordinate of every solution\n in the container for solutions in quad double precision."},
   {"py2c_product_supporting_set_structure",
     py2c_product_supporting_set_structure, METH_VARARGS,
    "Builds a supporting set structure for the system stored in the\n container with coefficients in standard double precision."},
   {"py2c_product_write_set_structure", py2c_product_write_set_structure,
     METH_VARARGS, 
    "Writes the supporting set structure to screen."},
   {"py2c_product_set_structure_string", py2c_product_set_structure_string,
     METH_VARARGS,
    "Returns the string representation of the set structure."},
   {"py2c_product_parse_set_structure", py2c_product_parse_set_structure,
     METH_VARARGS, 
    "Parses a given string into a set structure.\n On input are two parameters, one integer and one string:\n 1) the number of characters in the given string; and\n 2) the characters in the string.\n On return is the failure code, if zero, then the string\n has been parsed into a valid set structure."},
   {"py2c_product_is_set_structure_supporting",
     py2c_product_is_set_structure_supporting, METH_VARARGS,
    "Checks whether the stored set structure is supporting\n for the system in the standard systems container.\n Returns an integer which represents true (1) or false (0)."},
   {"py2c_product_linear_product_root_count",
     py2c_product_linear_product_root_count, METH_VARARGS,
    "Returns the linear-product root count, computed from\n the supporting set structure."},
   {"py2c_product_random_linear_product_system",
     py2c_product_random_linear_product_system, METH_VARARGS, 
    "Builds a random linear-product system based on the\n stored set structure.   On return is the failure code,\n which equals zero if all went well."},
   {"py2c_product_solve_linear_product_system",
     py2c_product_solve_linear_product_system, METH_VARARGS,
    "Computes all solutions to the linear-product system\n and stores the solutions in the container for solutions\n in standard double precision.  On return is the failure\n code, which equals zero if all went well."},
   {"py2c_product_clear_set_structure", py2c_product_clear_set_structure,
     METH_VARARGS,
    "Deallocates the set structure."},
   {"py2c_product_m_homogeneous_Bezout_number",
     py2c_product_m_homogeneous_Bezout_number, METH_VARARGS,
    "For the system in the standard systems container,\n a heuristic partition of the set of variables may\n lead to a Bezout number that is smaller than the total degree.\n On return is the m-homogeneous Bezout number for the\n string representation of the partition that is returned\n as the second argument in the tuple."},
   {"py2c_product_m_partition_Bezout_number",
     py2c_product_m_partition_Bezout_number, METH_VARARGS,
    "Given a partition of the set of variables, computes\n the m-homogeneous Bezout number for the system in\n the standard systems container.\n On input are two arguments:\n 1) the number of characters in the string (second argument); and\n 2) the string representation for a partition of the variables.\n On return is the m-homogeneous Bezout number."},
   {"py2c_product_m_homogeneous_start_system",
     py2c_product_m_homogeneous_start_system, METH_VARARGS,
    "Given a partition of the set of variables, constructs\n an m-homogeneous Bezout number for the system in\n the standard systems container.\n On input are two arguments:\n 1) the number of characters in the string (second argument); and\n 2) the string representation for a partition of the variables.\n On return is the m-homogeneous Bezout number."},
   {"py2c_celcon_initialize_supports",
     py2c_celcon_initialize_supports, METH_VARARGS,
    "Initializes the cell container with the number of distinct supports,\n this number is given as the one input parameter.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_celcon_set_type_of_mixture",
     py2c_celcon_set_type_of_mixture, METH_VARARGS,
    "Defines the type of mixture of the support sets.\n On input are two parameters, an integer and a string:\n 1) the integer equals the number of distinct supports;\n 2) the string is a string representation of a Python list of integers,\n there are as many integers as the value of the first parameter.\n Each integer is a positive number, equal to the number of occurrences\n of each support set."},
   {"py2c_celcon_type_of_mixture",
     py2c_celcon_type_of_mixture, METH_VARARGS,
    "Returns the string representation of the type of mixture of the support sets.\n This string is the string representation of a Python list of integers."},
   {"py2c_celcon_append_lifted_point",
     py2c_celcon_append_lifted_point, METH_VARARGS,
    "Appends a lifted point to the cells container.\n There are three input parameters:\n 1) the dimension of the point;\n 2) the index of the support to where to append to; and\n 3) the string representation of the lifted point.\n Returns the failure code, which equals zero when all went well."},
   {"py2c_celcon_retrieve_lifted_point",
     py2c_celcon_retrieve_lifted_point, METH_VARARGS,
    "Returns a string representation of a lifted point.\n On input are three integer numbers:\n 1) the number of coordinates in the lifted point;\n 2) the index to the support set; and\n 3) the index to the point in that support set."},
   {"py2c_celcon_mixed_volume_of_supports",
     py2c_celcon_mixed_volume_of_supports, METH_VARARGS,
    "Returns the mixed volume of the supports stored in the cell container."},
   {"py2c_celcon_number_of_cells", py2c_celcon_number_of_cells,
    METH_VARARGS, "returns the number of cells in the cell container"},
   {"py2c_celcon_is_stable", py2c_celcon_is_stable,
    METH_VARARGS, "returns 1 if stable mixed cells were stored, 0 otherwise"},
   {"py2c_celcon_number_of_original_cells",
     py2c_celcon_number_of_original_cells,
    METH_VARARGS, "returns the number of original cells in the cell container"},
   {"py2c_celcon_number_of_stable_cells",
     py2c_celcon_number_of_stable_cells,
    METH_VARARGS, "returns the number of stable cells in the cell container"},
   {"py2c_celcon_standard_random_coefficient_system",
     py2c_celcon_standard_random_coefficient_system, METH_VARARGS,
    "Based on the lifted supports stored in the container,\n a random coefficient system with coefficients in standard double\n precision is stored in the cell container."},
   {"py2c_celcon_dobldobl_random_coefficient_system",
     py2c_celcon_dobldobl_random_coefficient_system, METH_VARARGS,
    "Based on the lifted supports stored in the container,\n a random coefficient system with coefficients in double double\n precision is stored in the cell container."},
   {"py2c_celcon_quaddobl_random_coefficient_system",
     py2c_celcon_quaddobl_random_coefficient_system, METH_VARARGS,
    "Based on the lifted supports stored in the container,\n a random coefficient system with coefficients in quad double\n precision is stored in the cell container."},
   {"py2c_celcon_copy_into_standard_systems_container",
     py2c_celcon_copy_into_standard_systems_container, METH_VARARGS,
    "The random coefficient system in standard double precision is copied\n from the cell container to the container for systems with\n coefficients in standard double precision."},
   {"py2c_celcon_copy_into_dobldobl_systems_container",
     py2c_celcon_copy_into_dobldobl_systems_container, METH_VARARGS,
    "The random coefficient system in double double precision is copied\n from the cell container to the container for systems with\n coefficients in double double precision."},
   {"py2c_celcon_copy_into_quaddobl_systems_container",
     py2c_celcon_copy_into_quaddobl_systems_container, METH_VARARGS,
    "The random coefficient system in quad double precision is copied\n from the cell container to the container for systems with\n coefficients in quad double precision."},
   {"py2c_celcon_standard_polyhedral_homotopy",
     py2c_celcon_standard_polyhedral_homotopy, METH_VARARGS, 
    "Based on the lifting and the random coefficient system,\n the polyhedral homotopy to solve the random coefficient system\n in standard double precision is constructed.\n This function also initializes the internal data structures to store\n the solutions of start and target systems.\n The lifted supports and the random coefficient system are defined.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_celcon_dobldobl_polyhedral_homotopy",
     py2c_celcon_dobldobl_polyhedral_homotopy, METH_VARARGS, 
    "Based on the lifting and the random coefficient system,\n the polyhedral homotopy to solve the random coefficient system\n in double double precision is constructed.\n This function also initializes the internal data structures to store\n the solutions of start and target systems.\n The lifted supports and the random coefficient system are defined.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_celcon_quaddobl_polyhedral_homotopy",
     py2c_celcon_quaddobl_polyhedral_homotopy, METH_VARARGS, 
    "Based on the lifting and the random coefficient system,\n the polyhedral homotopy to solve the random coefficient system\n in quad double precision is constructed.\n This function also initializes the internal data structures to store\n the solutions of start and target systems.\n The lifted supports and the random coefficient system are defined.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_celcon_solve_standard_start_system",
     py2c_celcon_solve_standard_start_system, METH_VARARGS, 
    "Solves the start system corresponding to the k-th mixed cell,\n using standard double precision arithmetic.\n The precondition for this function is that the creation of\n the polyhedral homotopy in standard double precision ended well.\n On return is the number of solution found, which must equal\n the mixed volume of the k-th mixed cell."},
   {"py2c_celcon_solve_dobldobl_start_system",
     py2c_celcon_solve_dobldobl_start_system, METH_VARARGS, 
    "Solves the start system corresponding to the k-th mixed cell,\n using double double precision arithmetic.\n The precondition for this function is that the creation of\n the polyhedral homotopy in double double precision ended well.\n On return is the number of solution found, which must equal\n the mixed volume of the k-th mixed cell."},
   {"py2c_celcon_solve_quaddobl_start_system",
     py2c_celcon_solve_quaddobl_start_system, METH_VARARGS, 
    "Solves the start system corresponding to the k-th mixed cell,\n using quad double precision arithmetic.\n The precondition for this function is that the creation of\n the polyhedral homotopy in quad double precision ended well.\n On return is the number of solution found, which must equal\n the mixed volume of the k-th mixed cell."},
   {"py2c_celcon_solve_stable_standard_start_system",
     py2c_celcon_solve_stable_standard_start_system, METH_VARARGS, 
    "Solves the start system corresponding to the k-th stable mixed cell,\n using standard double precision arithmetic.\n The precondition for this function is that the creation of\n the polyhedral homotopy in standard double precision ended well.\n On return is the number of solution found, which must equal\n the mixed volume of the k-th stable mixed cell."},
   {"py2c_celcon_solve_stable_dobldobl_start_system",
     py2c_celcon_solve_stable_dobldobl_start_system, METH_VARARGS, 
    "Solves the start system corresponding to the k-th stable mixed cell,\n using double double precision arithmetic.\n The precondition for this function is that the creation of\n the polyhedral homotopy in double double precision ended well.\n On return is the number of solution found, which must equal\n the mixed volume of the k-th stable mixed cell."},
   {"py2c_celcon_solve_stable_quaddobl_start_system",
     py2c_celcon_solve_stable_quaddobl_start_system, METH_VARARGS, 
    "Solves the start system corresponding to the k-th stable mixed cell,\n using quad double precision arithmetic.\n The precondition for this function is that the creation of\n the polyhedral homotopy in quad double precision ended well.\n On return is the number of solution found, which must equal\n the mixed volume of the k-th stable mixed cell."},
   {"py2c_celcon_track_standard_solution_path",
     py2c_celcon_track_standard_solution_path, METH_VARARGS, 
    "Tracks a solution path starting at the i-th solution of the k-th cell,\n using standard double precision arithmetic.\n The precondition for this function is that the start system defined\n by the k-th mixed cell is solved in standard double precision.\n There are three input parameters:\n 1) k, the index to a mixed cell in the cell container;\n 2) i, the index to a solution path defined by that mixed cell;\n 3) otp, the level for intermediate output during path tracking.\n A target solution corresponding to the k-th cell is added on return."},
   {"py2c_celcon_track_dobldobl_solution_path",
     py2c_celcon_track_dobldobl_solution_path, METH_VARARGS, 
    "Tracks a solution path starting at the i-th solution of the k-th cell,\n using double double precision arithmetic.\n The precondition for this function is that the start system defined\n by the k-th mixed cell is solved in double double precision.\n There are three input parameters:\n 1) k, the index to a mixed cell in the cell container;\n 2) i, the index to a solution path defined by that mixed cell;\n 3) otp, the level for intermediate output during path tracking.\n A target solution corresponding to the k-th cell is added on return."},
   {"py2c_celcon_track_quaddobl_solution_path",
     py2c_celcon_track_quaddobl_solution_path, METH_VARARGS, 
    "Tracks a solution path starting at the i-th solution of the k-th cell,\n using quad double precision arithmetic.\n The precondition for this function is that the start system defined\n by the k-th mixed cell is solved in quad double precision.\n There are three input parameters:\n 1) k, the index to a mixed cell in the cell container;\n 2) i, the index to a solution path defined by that mixed cell;\n 3) otp, the level for intermediate output during path tracking.\n A target solution corresponding to the k-th cell is added on return."},
   {"py2c_celcon_copy_target_standard_solution_to_container",
     py2c_celcon_copy_target_standard_solution_to_container, METH_VARARGS, 
    "Copies the i-th target solution corresponding to the k-th mixed cell\n to the container for solutions in standard double precision.\n There are two input parameters for this function:\n 1) k, the index to the mixed cell;\n 2) i, the index to the i-th solution path defined by the cell.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_celcon_copy_target_dobldobl_solution_to_container",
     py2c_celcon_copy_target_dobldobl_solution_to_container, METH_VARARGS, 
    "Copies the i-th target solution corresponding to the k-th mixed cell\n to the container for solutions in double double precision.\n There are two input parameters for this function:\n 1) k, the index to the mixed cell;\n 2) i, the index to the i-th solution path defined by the cell.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_celcon_copy_target_quaddobl_solution_to_container",
     py2c_celcon_copy_target_quaddobl_solution_to_container, METH_VARARGS, 
    "Copies the i-th target solution corresponding to the k-th mixed cell\n to the container for solutions in quad double precision.\n There are two input parameters for this function:\n 1) k, the index to the mixed cell;\n 2) i, the index to the i-th solution path defined by the cell.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_celcon_permute_standard_system",
     py2c_celcon_permute_standard_system, METH_VARARGS,
    "Permutes the systems in the container for polynomial and Laurent systems\n with standard double coefficients corresponding to the permutation\n used to compute the mixed-cell configuration.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_celcon_permute_dobldobl_system",
     py2c_celcon_permute_dobldobl_system, METH_VARARGS,
    "Permutes the systems in the container for polynomial and Laurent systems\n with double double coefficients corresponding to the permutation\n used to compute the mixed-cell configuration.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_celcon_permute_quaddobl_system",
     py2c_celcon_permute_quaddobl_system, METH_VARARGS,
    "Permutes the systems in the container for polynomial and Laurent systems\n with quad double coefficients corresponding to the permutation\n used to compute the mixed-cell configuration.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_celcon_clear_container", py2c_celcon_clear_container, METH_VARARGS,
    "Deallocates the data in the cell container."},
   {"py2c_intcelcon_read_mixed_cell_configuration",
     py2c_intcelcon_read_mixed_cell_configuration, METH_VARARGS,
    "Reads a mixed-cell configuration"},
   {"py2c_intcelcon_write_mixed_cell_configuration",
     py2c_intcelcon_write_mixed_cell_configuration, METH_VARARGS,
    "Writes the mixed-cell configuration to screen."},
   {"py2c_intcelcon_number_of_cells",
     py2c_intcelcon_number_of_cells, METH_VARARGS,
    "Returns the number of cells in the mixed subdivision by integer lifting."},
   {"py2c_intcelcon_type_of_mixture",
     py2c_intcelcon_type_of_mixture, METH_VARARGS,
    "Returns the type of mixture for the integer cells container."},
   {"py2c_intcelcon_length_of_supports",
     py2c_intcelcon_length_of_supports, METH_VARARGS,
    "Returns the string representation of a list of lengths of each support."},
   {"py2c_intcelcon_append_lifted_point",
     py2c_intcelcon_append_lifted_point, METH_VARARGS,
    "Appends a lifted point to the cells container.\n There are three input parameters:\n 1) the dimension of the point;\n 2) the index of the support to where to append to; and\n 3) the string representation of the lifted point.\n Returns the failure code, which equals zero when all went well."},
   {"py2c_intcelcon_get_lifted_point",
     py2c_intcelcon_get_lifted_point, METH_VARARGS,
    "Returns the string representation of the coordinates of a lifted point."},
   {"py2c_intcelcon_get_inner_normal",
     py2c_intcelcon_get_inner_normal, METH_VARARGS,
    "Given on input the dimension of the lifted points and the\n index of the mixed cell of interest, returns the string\n representation of the inner normal of the mixed cell."},
   {"py2c_intcelcon_number_of_points_in_cell",
     py2c_intcelcon_number_of_points_in_cell, METH_VARARGS,
    "Given are two integer numbers: the index to a cell\n (starting the count at one) and the number of different supports.\n On return is the string representation of the number of points\n which span each component of the mixed cell."},
   {"py2c_intcelcon_get_point_in_cell",
     py2c_intcelcon_get_point_in_cell, METH_VARARGS,
    "Returns the string representation of the n coordinates of\n the k-th point from the j-th list of the i-th cell.\n On input are the four integers: n, i, j, k, respectively\n the length of the lifted vectors in the supports,\n the index to a cell in the container,\n the index to a support of the i-th cell, and\n the index to a point in the j-th support of the i-th cell."},
   {"py2c_intcelcon_mixed_volume",
     py2c_intcelcon_mixed_volume, METH_VARARGS,
    "Returns the mixed volume of a mixed cell."},
   {"py2c_intcelcon_initialize_supports",
     py2c_intcelcon_initialize_supports, METH_VARARGS,
    "Initializes the supports with an integer valued lifting."},
   {"py2c_intcelcon_set_type_of_mixture",
     py2c_intcelcon_set_type_of_mixture, METH_VARARGS,
    "Defines the type of mixture of the tuple of supports."},
   {"py2c_intcelcon_make_subdivision",
     py2c_intcelcon_make_subdivision, METH_VARARGS,
    "Computes the cells in the regular subdivision induced by an integer\n valued lifting function."},
   {"py2c_intcelcon_clear_mixed_cell_configuration",
     py2c_intcelcon_clear_mixed_cell_configuration, METH_VARARGS,
    "Deallocates the data in the integer cell container."},
   {"py2c_scale_standard_system", py2c_scale_standard_system, METH_VARARGS,
    "Applies scaling to the system in the standard systems container,\n with standard double precision arithmetic.  The system in the standard\n systems container is replaced by the scaled system.\n On entry is one integer, which should be either 0, 1, or 2:\n 0 for only scaling of the equations,\n 1 variable scaling without variability reduction,\n 2 variable scaling with variability reduction.\n On return is a tuple with the scaling coefficients (if mode > 0)\n and the estimated inverse condition number of the scaling problem."},
   {"py2c_scale_dobldobl_system", py2c_scale_dobldobl_system, METH_VARARGS,
    "Applies scaling to the system in the dobldobl systems container,\n with double double precision arithmetic.  The system in the dobldobl\n systems container is replaced by the scaled system.\n On entry is one integer, which should be either 0, 1, or 2:\n 0 for only scaling of the equations,\n 1 variable scaling without variability reduction,\n 2 variable scaling with variability reduction.\n On return is a tuple with the scaling coefficients (if mode > 0)\n and the estimated inverse condition number of the scaling problem."},
   {"py2c_scale_quaddobl_system", py2c_scale_quaddobl_system, METH_VARARGS,
    "Applies scaling to the system in the quaddobl systems container,\n with quad double precision arithmetic.  The system in the quaddobl\n systems container is replaced by the scaled system.\n On entry is one integer, which should be either 0, 1, or 2:\n 0 for only scaling of the equations,\n 1 variable scaling without variability reduction,\n 2 variable scaling with variability reduction.\n On return is a tuple with the scaling coefficients (if mode > 0)\n and the estimated inverse condition number of the scaling problem."},
   {"py2c_scale_standard_solutions",
     py2c_scale_standard_solutions, METH_VARARGS,
    "Replaces the solutions in the standard solutions container with\n the scaled solutions, scaled with standard double precision arithmetic,\n using the given scaling coefficients.\n On entry are two parameters: an integer and a string.\n The integer contains the number of elements in the list\n of scaling coefficients (doubles) stored in the string.\n The format of the string is the Python string representation\n of a list of doubles, i.e.: starting with [ and ending with ]."},
   {"py2c_scale_dobldobl_solutions",
     py2c_scale_dobldobl_solutions, METH_VARARGS,
    "Replaces the solutions in the dobldobl solutions container with\n the scaled solutions, scaled with double double precision arithmetic,\n using the given scaling coefficients.\n On entry are two parameters: an integer and a string.\n The integer contains the number of elements in the list\n of scaling coefficients (doubles) stored in the string.\n The format of the string is the Python string representation\n of a list of doubles, i.e.: starting with [ and ending with ]."},
   {"py2c_scale_quaddobl_solutions",
     py2c_scale_quaddobl_solutions, METH_VARARGS,
    "Replaces the solutions in the quaddobl solutions container with\n the scaled solutions, scaled with quad double precision arithmetic,\n using the given scaling coefficients.\n On entry are two parameters: an integer and a string.\n The integer contains the number of elements in the list\n of scaling coefficients (doubles) stored in the string.\n The format of the string is the Python string representation\n of a list of doubles, i.e.: starting with [ and ending with ]."},
   {"py2c_linear_reduce_standard_system",
     py2c_linear_reduce_standard_system, METH_VARARGS,
    "Applies linear reduction on the coefficient matrix of the system\n in the container for standard double precision.\n There is one integer parameter: whether to diagonalize or not."},
   {"py2c_linear_reduce_dobldobl_system",
     py2c_linear_reduce_dobldobl_system, METH_VARARGS,
    "Applies linear reduction on the coefficient matrix of the system\n in the container for double double precision.\n There is one integer parameter: whether to diagonalize or not."},
   {"py2c_linear_reduce_quaddobl_system",
     py2c_linear_reduce_quaddobl_system, METH_VARARGS,
    "Applies linear reduction on the coefficient matrix of the system\n in the container for quad double precision.\n There is one integer parameter: whether to diagonalize or not."},
   {"py2c_nonlinear_reduce_standard_system",
     py2c_nonlinear_reduce_standard_system, METH_VARARGS,
    "Applies nonlinear reduction on the system in the container\n for standard double precision.\n Three integer numbers are expected on input:\n (1) the maximum number of equal degree replacements,\n (2) the maximum number of computed S-polynomials,\n (3) the maximum number of computed R-polynomials.\n The system in the standard container is replace by the reduced system.\n Three numbers are returned:\n (1) the number of equal degree replacements,\n (2) the number of computed S-polynomials,\n (3) the number of computed R-polynomials."},
   {"py2c_sweep_define_parameters_numerically",
     py2c_sweep_define_parameters_numerically, METH_VARARGS,
    "Defines the indices to the variables that serve as parameters\n numerically, that is: via integer indices.\n On entry are three integer numbers and a string.\n The string is a string representation of a Python list of integers,\n The three integers are the number of equations, the number of variables,\n and the number of parameters.  The number of variables m includes the\n number of parameters.  Then there should be as many as m indices in\n the list of integers to define which of the variables are parameters."},
   {"py2c_sweep_define_parameters_symbolically",
     py2c_sweep_define_parameters_symbolically, METH_VARARGS,
    "Defines the indices to the variables that serve as parameters\n symbolically, that is, as names of variables.\n For this to work, the symbol table must be initialized.\n On entry are four integer numbers and a string.\n The four integers are the number of equations, the number of variables,\n the number of parameters (the number of variables m includes the\n number of parameters), and the number of characters in the string.\n The string contains the names of the parameters, separated by one comma.\n For this to work, the symbol table must be initialized, e.g.:\n via the reading of a polynomial system."},
   {"py2c_sweep_get_number_of_equations",
     py2c_sweep_get_number_of_equations, METH_VARARGS,
    "Returns the number of equations in the sweep homotopy."},
   {"py2c_sweep_get_number_of_variables",
     py2c_sweep_get_number_of_variables, METH_VARARGS,
    "Returns the number of variables in the sweep homotopy."},
   {"py2c_sweep_get_number_of_parameters",
     py2c_sweep_get_number_of_parameters, METH_VARARGS,
    "Returns the number of parameters in the sweep homotopy."},
   {"py2c_sweep_get_indices_numerically",
     py2c_sweep_get_indices_numerically, METH_VARARGS,
    "Returns the indices of the variables that are parameters,\n as the string representation of a Python list of integers."},
   {"py2c_sweep_get_indices_symbolically",
     py2c_sweep_get_indices_symbolically, METH_VARARGS,
    "Returns a string with the names of the parameters,\n each separated by one space."},
   {"py2c_sweep_clear_definitions",
     py2c_sweep_clear_definitions, METH_VARARGS,
    "Clears the definitions in the sweep homotopy."},
   {"py2c_sweep_set_standard_start",
     py2c_sweep_set_standard_start, METH_VARARGS,
    "Sets the start values for the m parameters in standard double precision,\n giving on input an integer m and 2*m doubles, with the consecutive\n real and imaginary parts for the start values of all m parameters.\n The doubles are given in a string representation of a Python\n list of doubles."},
   {"py2c_sweep_set_standard_target",
     py2c_sweep_set_standard_target, METH_VARARGS,
    "Sets the target values for the m parameters in standard double precision,\n giving on input an integer m and 2*m doubles, with the consecutive\n real and imaginary parts for the target values of all m parameters."},
   {"py2c_sweep_set_dobldobl_start",
     py2c_sweep_set_dobldobl_start, METH_VARARGS,
    "Sets the start values for the m parameters in double double precision,\n giving on input an integer m and 4*m doubles, with the consecutive\n real and imaginary parts for the start values of all m parameters."},
   {"py2c_sweep_set_dobldobl_target",
     py2c_sweep_set_dobldobl_target, METH_VARARGS,
    "Sets the target values for the m parameters in double double precision,\n giving on input an integer m and 4*m doubles, with the consecutive\n real and imaginary parts for the target values of all m parameters."},
   {"py2c_sweep_set_quaddobl_start",
     py2c_sweep_set_quaddobl_start, METH_VARARGS,
    "Sets the start values for the m parameters in quad double precision,\n giving on input an integer m and 8*m doubles, with the consecutive\n real and imaginary parts for the start values of all m parameters."},
   {"py2c_sweep_set_quaddobl_target",
     py2c_sweep_set_quaddobl_target, METH_VARARGS,
    "Sets the target values for the m parameters in quad double precision,\n giving on input an integer m and 8*m doubles, with the consecutive\n real and imaginary parts for the target values of all m parameters."},
   {"py2c_sweep_get_standard_start",
     py2c_sweep_get_standard_start, METH_VARARGS,
    "Gets the start values for the parameters in standard double precision,\n giving on input the number n of doubles that need to be returned.\n On return will be n doubles, for the consecutive real and imaginary\n parts for the start values of all parameters,\n stored in the string representation of a Python list of doubles."},
   {"py2c_sweep_get_standard_target",
     py2c_sweep_get_standard_target, METH_VARARGS,
    "Gets the target values for the parameters in standard double precision,\n giving on input the number n of doubles that need to be returned.\n On return will be n doubles, for the consecutive real and imaginary\n parts for the target values of all parameters,\n stored in the string representation of a Python list of doubles."},
   {"py2c_sweep_get_dobldobl_start",
     py2c_sweep_get_dobldobl_start, METH_VARARGS,
    "Gets the start values for the parameters in double double precision,\n giving on input the number n of doubles that need to be returned.\n On return will be n doubles, for the consecutive real and imaginary\n parts for the start values of all parameters,\n stored in the string representation of a Python list of doubles."},
   {"py2c_sweep_get_dobldobl_target",
     py2c_sweep_get_dobldobl_target, METH_VARARGS,
    "Gets the target values for the parameters in double double precision,\n giving on input the number n of doubles that need to be returned.\n On return will be n doubles, for the consecutive real and imaginary\n parts for the target values of all parameters,\n stored in the string representation of a Python list of doubles."},
   {"py2c_sweep_get_quaddobl_start",
     py2c_sweep_get_quaddobl_start, METH_VARARGS,
    "Gets the start values for the parameters in quad double precision,\n giving on input the number n of doubles that need to be returned.\n On return will be n doubles, for the consecutive real and imaginary\n parts for the start values of all parameters,\n stored in the string representation of a Python list of doubles."},
   {"py2c_sweep_get_quaddobl_target",
     py2c_sweep_get_quaddobl_target, METH_VARARGS,
    "Returns the target values for the parameters in quad double precision,\n giving on input the number n of doubles that need to be returned.\n On return will be n doubles, for the consecutive real and imaginary\n parts for the target values of all parameters,\n stored in the string representation of a Python list of doubles."},
   {"py2c_sweep_standard_complex_run",
     py2c_sweep_standard_complex_run, METH_VARARGS,
    "Starts the trackers in a complex convex parameter homotopy,\n in standard double precision, where the indices to the parameters,\n start and target values are already defined.  Moreover, the containers\n of systems and solutions in standard double precision have been\n initialized with a parametric systems and start solutions.\n The first input parameter is 0, 1, or 2, for respectively\n a randomly generated gamma (0), or no gamma (1), or a user given\n gamma with real and imaginary parts given in 2 pointers to doubles."},
   {"py2c_sweep_dobldobl_complex_run",
     py2c_sweep_dobldobl_complex_run, METH_VARARGS,
    "Starts the trackers in a complex convex parameter homotopy,\n in double double precision, where the indices to the parameters,\n start and target values are already defined.  Moreover, the containers\n of systems and solutions in double double precision have been\n initialized with a parametric systems and start solutions.\n The first input parameter is 0, 1, or 2, for respectively\n a randomly generated gamma (0), or no gamma (1), or a user given\n gamma with real and imaginary parts given in 2 pointers to doubles."},
   {"py2c_sweep_quaddobl_complex_run",
     py2c_sweep_quaddobl_complex_run, METH_VARARGS,
    "Starts the trackers in a complex convex parameter homotopy,\n in quad double precision, where the indices to the parameters,\n start and target values are already defined.  Moreover, the containers\n of systems and solutions in quad double precision have been\n initialized with a parametric systems and start solutions.\n The first input parameter is 0, 1, or 2, for respectively\n a randomly generated gamma (0), or no gamma (1), or a user given\n gamma with real and imaginary parts given in 2 pointers to doubles."},
   {"py2c_sweep_standard_real_run",
     py2c_sweep_standard_real_run, METH_VARARGS, 
    "There are no input arguments to this routine.\n Starts a sweep with a natural parameter in a family of n equations\n in n+1 variables, where the last variable is the artificial parameter s\n that moves the one natural parameter from a start to target value.\n The last equation is of the form (1-s)*(A - v[0]) + s*(A - v[1]),\n where A is the natural parameter, going from the start value v[0]\n to the target value v[1].\n This family must be stored in the systems container in standard double\n precision and the corresponding start solutions in the standard solutions\n container, where every solution has the value v[0] for the A variable.\n The sweep stops when s reaches the value v[1], or when a singularity\n is encountered on the path."},
   {"py2c_sweep_dobldobl_real_run",
     py2c_sweep_dobldobl_real_run, METH_VARARGS, 
    "There are no input arguments to this routine.\n Starts a sweep with a natural parameter in a family of n equations\n in n+1 variables, where the last variable is the artificial parameter s\n that moves the one natural parameter from a start to target value.\n The last equation is of the form (1-s)*(A - v[0]) + s*(A - v[1]),\n where A is the natural parameter, going from the start value v[0]\n to the target value v[1].\n This family must be stored in the systems container in double double\n precision and the corresponding start solutions in the dobldobl solutions\n container, where every solution has the value v[0] for the A variable.\n The sweep stops when s reaches the value v[1], or when a singularity\n is encountered on the path."},
   {"py2c_sweep_quaddobl_real_run",
     py2c_sweep_quaddobl_real_run, METH_VARARGS, 
    "There are no input arguments to this routine.\n Starts a sweep with a natural parameter in a family of n equations\n in n+1 variables, where the last variable is the artificial parameter s\n that moves the one natural parameter from a start to target value.\n The last equation is of the form (1-s)*(A - v[0]) + s*(A - v[1]),\n where A is the natural parameter, going from the start value v[0]\n to the target value v[1].\n This family must be stored in the systems container in quad double\n precision and the corresponding start solutions in the quaddobl solutions\n container, where every solution has the value v[0] for the A variable.\n The sweep stops when s reaches the value v[1], or when a singularity\n is encountered on the path."},
   {"py2c_standard_multiplicity_structure",
     py2c_standard_multiplicity_structure, METH_VARARGS,
    "Computes the multiplicity structure in standard double precision.\n Required is the presence of a polynomial system in the standard\n systems container and a solution in the standard solutions container.\n The input parameters are two integers and one double:\n order : the maximum differentiation order,\n verbose : 1 for verbose, 0 for silent, and\n tol : tolerance on the numerical rank.\n On return is a tuple: the multiplicity and the values\n of the Hilbert function."},
   {"py2c_dobldobl_multiplicity_structure",
     py2c_dobldobl_multiplicity_structure, METH_VARARGS,
    "Computes the multiplicity structure in double double precision.\n Required is the presence of a polynomial system in the dobldobl\n systems container and a solution in the dobldobl solutions container.\n The input parameters are two integers and one double:\n order : the maximum differentiation order,\n verbose : 1 for verbose, 0 for silent, and\n tol : tolerance on the numerical rank.\n On return is a tuple: the multiplicity and the values\n of the Hilbert function."},
   {"py2c_quaddobl_multiplicity_structure",
     py2c_quaddobl_multiplicity_structure, METH_VARARGS,
    "Computes the multiplicity structure in quad double precision.\n Required is the presence of a polynomial system in the quaddobl\n systems container and a solution in the quaddobl solutions container.\n The input parameters are two integers and one double:\n order : the maximum differentiation order,\n verbose : 1 for verbose, 0 for silent, and\n tol : tolerance on the numerical rank.\n On return is a tuple: the multiplicity and the values\n of the Hilbert function."},
   {"py2c_numbtrop_standard_initialize",
     py2c_numbtrop_standard_initialize, METH_VARARGS,
    "Initializes the numerical tropisms container,\n in standard double precision.  The input parameters are\n nbt : number of tropisms;\n dim : length_of_each tropism;\n wnd : winding numbers, as many as nbt;\n dir : nbt*dim doubles with the coordinates of the tropisms;\n err : errors on the tropisms, as many doubles as the value of nbt.\n The numbers in wnd, dir, and err must be given in one string,\n as the string representation of a list of doubles.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_numbtrop_dobldobl_initialize",
     py2c_numbtrop_dobldobl_initialize, METH_VARARGS,
    "Initializes the numerical tropisms container,\n in double double precision.  The input parameters are\n nbt : number of tropisms;\n dim : length_of_each tropism;\n wnd : winding numbers, as many as nbt;\n dir : 2*nbt*dim doubles with the coordinates of the tropisms;\n err : errors on the tropisms, as many doubles as the value of 2*nbt.\n The numbers in wnd, dir, and err must be given in one string,\n as the string representation of a list of doubles.\n On return is the the failure code, which equals zero if all went well."},
   {"py2c_numbtrop_quaddobl_initialize",
     py2c_numbtrop_quaddobl_initialize, METH_VARARGS,
    "Initializes the numerical tropisms container,\n in quad double precision.  The input parameters are\n nbt : number of tropisms;\n dim : length_of_each tropism;\n wnd : winding numbers, as many as nbt;\n dir : 4*nbt*dim doubles with the coordinates of the tropisms;\n err : errors on the tropisms, as many doubles as the value of 4*nbt.\n The numbers in wnd, dir, and err must be given in one string,\n as the string representation of a list of doubles.\n On return is the the failure code, which equals zero if all went well."},
   {"py2c_numbtrop_standard_retrieve",
     py2c_numbtrop_standard_retrieve, METH_VARARGS,
    "Retrieves all tropisms stored in standard double precision.\n The input parameters are two integers:\n nbt : number of tropisms;\n dim : length_of_each tropism.\n On return are\n wnd : winding numbers, as many as nbt;\n dir : nbt*dim doubles with the coordinates of the tropisms;\n err : errors on the tropisms, as many doubles as the value of nbt.\n All numbers are returns in one string, as the string representation\n of a list of doubles.\n The failure code, which equals zero if all went well."},
   {"py2c_numbtrop_dobldobl_retrieve",
     py2c_numbtrop_dobldobl_retrieve, METH_VARARGS,
    "Retrieves all tropisms stored in double double precision.\n The input parameters are two integers:\n nbt : number of tropisms;\n dim : length_of_each tropism.\n On return are\n wnd : winding numbers, as many as nbt;\n dir : 2*nbt*dim doubles with the coordinates of the tropisms;\n err : errors on the tropisms, as many doubles as the value of 2*nbt.\n All numbers are returns in one string, as the string representation\n of a list of doubles.\n The failure code, which equals zero if all went well."},
   {"py2c_numbtrop_quaddobl_retrieve",
     py2c_numbtrop_quaddobl_retrieve, METH_VARARGS,
    "Retrieves all tropisms stored in quad double precision.\n The input parameters are two integers:\n nbt : number of tropisms;\n dim : length_of_each tropism.\n On return are\n wnd : winding numbers, as many as nbt;\n dir : 4*nbt*dim doubles with the coordinates of the tropisms;\n err : errors on the tropisms, as many doubles as the value of 4*nbt.\n All numbers are returns in one string, as the string representation\n of a list of doubles.\n The failure code, which equals zero if all went well."},
   {"py2c_numbtrop_standard_size", py2c_numbtrop_standard_size, METH_VARARGS,
    "Returns the number of tropisms, stored in standard double\n precision, in the numerical tropisms container."},
   {"py2c_numbtrop_dobldobl_size", py2c_numbtrop_dobldobl_size, METH_VARARGS,
    "Returns the number of tropisms, stored in double double\n precision, in the numerical tropisms container."},
   {"py2c_numbtrop_quaddobl_size", py2c_numbtrop_quaddobl_size, METH_VARARGS,
    "Returns the number of tropisms, stored in quad double\n precision, in the numerical tropisms container."},
   {"py2c_numbtrop_standard_dimension",
     py2c_numbtrop_standard_dimension, METH_VARARGS,
    "Returns the dimension of the tropisms, stored in standard double\n precision, in the numerical tropisms container."},
   {"py2c_numbtrop_dobldobl_dimension",
     py2c_numbtrop_dobldobl_dimension, METH_VARARGS,
    "Returns the dimension of the tropisms, stored in double double\n precision, in the numerical tropisms container."},
   {"py2c_numbtrop_quaddobl_dimension",
     py2c_numbtrop_quaddobl_dimension, METH_VARARGS,
    "Returns the dimension of the tropisms, stored in quad double\n precision, in the numerical tropisms container."},
   {"py2c_numbtrop_store_standard_tropism",
     py2c_numbtrop_store_standard_tropism, METH_VARARGS,
    "Stores a tropism given in standard double precision.\n The first three input parmeters are integers:\n dim : the length of the tropism vector;\n idx : the index of the tropism, indexing starts at one,\n and ends at nbt, what is returned by standard_size;\n wnd : estimated winding number;\n The other input parameters are of type double:\n dir : coordinates of the tropisms, as many as dim;\n err : the error on the tropism.\n All dim+1 doubles are given in one string,\n the string representation of a list of doubles."},
   {"py2c_numbtrop_store_dobldobl_tropism",
     py2c_numbtrop_store_dobldobl_tropism, METH_VARARGS,
    "Stores a tropism given in double double precision.\n The first three input parameters are integers:\n dim : the length of the tropism vector;\n idx : the index of the tropism, indexing starts at one,\n and ends at nbt, what is returned by dobldobl_size;\n wnd : estimated winding number;\n The other input parameters are of type double:\n dir : coordinates of the tropisms, as many as 2*dim;\n err : the error on the tropism, two doubles.\n All 2*dim+2 doubles are given in one string,\n the string representatin of a list of doubles."},
   {"py2c_numbtrop_store_quaddobl_tropism",
     py2c_numbtrop_store_quaddobl_tropism, METH_VARARGS,
    "Stores a tropism given in quad double precision.\n The first three input parameters are integers:\n dim : the length of the tropism vector;\n idx : the index of the tropism, indexing starts at one,\n and ends at nbt, what is returned by quaddobl_size;\n The other input parameters are of type double:\n wnd : estimated winding number;\n dir : coordinates of the tropisms, as many as 4*dim;\n err : the error on the tropism, four double\n All 4*dim+4 doubles are given in one string,\n the string representatin of a list of doubles."},
   {"py2c_numbtrop_standard_retrieve_tropism",
     py2c_numbtrop_standard_retrieve_tropism, METH_VARARGS,
    "Returns one tropism, stored in standard double precision.\n The input parameters are two integers:\n dim : the length of the tropism vector;\n idx : the index of the tropism, indexing starts at one,\n and ends at nbt, what is returned by numbtrop_standard_size.\n The first parameter on return is an integer:\n wnd : estimated winding number;\n The other output parameters are of type double:\n dir : coordinates of the tropisms, as many as dim;\n err : the error on the tropism.\n All dim+1 doubles are returned in one string,\n the string representation of a list of doubles."},
   {"py2c_numbtrop_dobldobl_retrieve_tropism",
     py2c_numbtrop_dobldobl_retrieve_tropism, METH_VARARGS,
    "Returns one tropism, stored in double double precision.\n The input parameters are two integers:\n dim : the length of the tropism vector;\n idx : the index of the tropism, indexing starts at one,\n and ends at nbt, what is returned by numbtrop_dobldobl_size.\n The first parameter on return is an integer:\n wnd : estimated winding number;\n The other output parameters are of type double:\n dir : coordinates of the tropisms, as many as 2*dim;\n err : the error on the tropism, two doubles.\n All 2*dim+2 doubles are returned in one string,\n the string representation of a list of doubles."},
   {"py2c_numbtrop_quaddobl_retrieve_tropism",
     py2c_numbtrop_quaddobl_retrieve_tropism, METH_VARARGS,
    "Returns one tropism, stored in quad double precision.\n The input parameters are two integers:\n dim : the length of the tropism vector;\n idx : the index of the tropism, indexing starts at one,\n and ends at nbt, what is returned by numbtrop_quaddobl_size.\n The first parameter on return is an integer:\n wnd : estimated winding number;\n The other output parameters are of type double:\n dir : coordinates of the tropisms, as many as 4*dim;\n err : the error on the tropism, four doubles.\n All 4*dim+4 doubles are returned in one string,\n the string representation of a list of doubles."},
   {"py2c_numbtrop_standard_clear", py2c_numbtrop_standard_clear, METH_VARARGS,
    "Deallocates the stored numerically computed tropisms,\n computed in standard double precision."},
   {"py2c_numbtrop_dobldobl_clear", py2c_numbtrop_dobldobl_clear, METH_VARARGS,
    "Deallocates the stored numerically computed tropisms,\n computed in double double precision."},
   {"py2c_numbtrop_quaddobl_clear", py2c_numbtrop_quaddobl_clear, METH_VARARGS,
    "Deallocates the stored numerically computed tropisms,\n computed in quad double precision."},
   {"py2c_embed_system", py2c_embed_system, METH_VARARGS,
    "Replaces the system in the container with its embedding of dimension d.\n The dimension d is given as the first integer parameter on input.\n The second integer parameter indicates the precision, either 0, 1, or 2,\n respectively for double, double double, or quad double precision.\n The third integer parameter is the verbose level.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_embed_standard_system", py2c_embed_standard_system, METH_VARARGS,
    "Replaces the system with coefficients in standard double precision\n in the container with its embedding of dimension d.\n The dimension d is given as an integer parameter on input.\n The second integer parameter is the verbose level.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_embed_dobldobl_system", py2c_embed_dobldobl_system, METH_VARARGS,
    "Replaces the system with coefficients in double double precision\n in the container with its embedding of dimension d.\n The dimension d is given as an integer parameter on input.\n The second integer parameter is the verbose level.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_embed_quaddobl_system", py2c_embed_quaddobl_system, METH_VARARGS,
    "Replaces the system with coefficients in quad double precision\n in the container with its embedding of dimension d.\n The dimension d is given as an integer parameter on input.\n The second integer parameter is the verbose level.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_embed_standard_Laurent_system",
     py2c_embed_standard_Laurent_system, METH_VARARGS,
    "Replaces the Laurent system with coefficients in standard double\n precision in the container with its embedding of dimension d.\n The dimension d is given as an integer parameter on input.\n The second integer parameter is the verbose level.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_embed_dobldobl_Laurent_system",
     py2c_embed_dobldobl_Laurent_system, METH_VARARGS,
    "Replaces the Laurent system with coefficients in double double\n precision in the container with its embedding of dimension d.\n The dimension d is given as an integer parameter on input.\n The second integer parameter is the verbose level.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_embed_quaddobl_Laurent_system",
     py2c_embed_quaddobl_Laurent_system, METH_VARARGS,
    "Replaces the Laurent system with coefficients in quad double\n precision in the container with its embedding of dimension d.\n The dimension d is given as an integer parameter on input.\n The second integer parameter is the verbose level.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_swap_symbols_for_standard_witness_set",
     py2c_swap_symbols_for_standard_witness_set, METH_VARARGS,
    "Permutes the slack variables in the polynomial system with standard\n double precision coefficients and its corresponding solutions in the\n containers so the slack variables appear at the end.  On input are\n two integers: the total number of variables; and\n the number of slack variables, or the dimension of the set.\n This permutation is necessary to consider the system and solutions\n stored in containers as a witness set."},
   {"py2c_swap_symbols_for_dobldobl_witness_set",
     py2c_swap_symbols_for_dobldobl_witness_set, METH_VARARGS,
    "Permutes the slack variables in the polynomial system with double\n double precision coefficients and its corresponding solutions in the\n containers so the slack variables appear at the end.  On input are\n two integers: the total number of variables; and\n the number of slack variables, or the dimension of the set.\n This permutation is necessary to consider the system and solutions\n stored in containers as a witness set."},
   {"py2c_swap_symbols_for_quaddobl_witness_set",
     py2c_swap_symbols_for_quaddobl_witness_set, METH_VARARGS,
    "Permutes the slack variables in the polynomial system with quad\n double precision coefficients and its corresponding solutions in the\n containers so the slack variables appear at the end.  On input are\n two integers: the total number of variables; and\n the number of slack variables, or the dimension of the set.\n This permutation is necessary to consider the system and solutions\n stored in containers as a witness set."},
   {"py2c_swap_symbols_for_standard_Laurent_witness_set",
     py2c_swap_symbols_for_standard_Laurent_witness_set, METH_VARARGS,
    "Permutes the slack variables in the Laurent system with standard\n double precision coefficients and its corresponding solutions in the\n containers so the slack variables appear at the end.  On input are\n two integers: the total number of variables; and\n the number of slack variables, or the dimension of the set.\n This permutation is necessary to consider the system and solutions\n stored in containers as a witness set."},
   {"py2c_swap_symbols_for_dobldobl_Laurent_witness_set",
     py2c_swap_symbols_for_dobldobl_Laurent_witness_set, METH_VARARGS,
    "Permutes the slack variables in the Laurent system with double\n double precision coefficients and its corresponding solutions in the\n containers so the slack variables appear at the end.  On input are\n two integers: the total number of variables; and\n the number of slack variables, or the dimension of the set.\n This permutation is necessary to consider the system and solutions\n stored in containers as a witness set."},
   {"py2c_swap_symbols_for_quaddobl_Laurent_witness_set",
     py2c_swap_symbols_for_quaddobl_Laurent_witness_set, METH_VARARGS,
    "Permutes the slack variables in the Laurent system with quad\n double precision coefficients and its corresponding solutions in the\n containers so the slack variables appear at the end.  On input are\n two integers: the total number of variables; and\n the number of slack variables, or the dimension of the set.\n This permutation is necessary to consider the system and solutions\n stored in containers as a witness set."},
   {"py2c_standard_cascade_homotopy", py2c_standard_cascade_homotopy,
     METH_VARARGS,
    "Creates a homotopy in standard double precision using the stored\n systems to go one level down the cascade, removing one slice.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_dobldobl_cascade_homotopy", py2c_dobldobl_cascade_homotopy,
     METH_VARARGS,
    "Creates a homotopy in double double precision using the stored\n systems to go one level down the cascade, removing one slice.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_quaddobl_cascade_homotopy", py2c_quaddobl_cascade_homotopy,
     METH_VARARGS,
    "Creates a homotopy in quad double precision using the stored\n systems to go one level down the cascade, removing one slice.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_standard_Laurent_cascade_homotopy",
     py2c_standard_Laurent_cascade_homotopy, METH_VARARGS,
    "Creates a homotopy in standard double precision using the stored\n Laurent systems to go one level down the cascade, removing one slice.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_dobldobl_Laurent_cascade_homotopy",
     py2c_dobldobl_Laurent_cascade_homotopy, METH_VARARGS,
    " Creates a homotopy in double double precision using the stored\n Laurent systems to go one level down the cascade, removing one slice.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_quaddobl_Laurent_cascade_homotopy",
     py2c_quaddobl_Laurent_cascade_homotopy, METH_VARARGS,
    "Creates a homotopy in quad double precision using the stored\n Laurent systems to go one level down the cascade, removing one slice.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_factor_set_standard_to_mute",
     py2c_factor_set_standard_to_mute, METH_VARARGS,
    "Sets the state of monodromy permutations in standard double\n precision to silent."},
   {"py2c_factor_set_dobldobl_to_mute",
     py2c_factor_set_dobldobl_to_mute, METH_VARARGS,
    "Sets the state of monodromy permutations in double double\n precision to silent."},
   {"py2c_factor_set_quaddobl_to_mute",
     py2c_factor_set_quaddobl_to_mute, METH_VARARGS,
    "Sets the state of monodromy permutations in quad double\n precision to silent."},
   {"py2c_factor_set_standard_to_verbose",
     py2c_factor_set_standard_to_verbose, METH_VARARGS,
    "Sets the state of monodromy permutations in standard double\n precision to verbose."},
   {"py2c_factor_set_dobldobl_to_verbose",
     py2c_factor_set_dobldobl_to_verbose, METH_VARARGS,
    "Sets the state of monodromy permutations in double double\n precision to verbose."},
   {"py2c_factor_set_quaddobl_to_verbose",
     py2c_factor_set_quaddobl_to_verbose, METH_VARARGS,
    "Sets the state of monodromy permutations in quad double\n precision to verbose."},
   {"py2c_factor_define_output_file_with_string",
     py2c_factor_define_output_file_with_string, METH_VARARGS,
    "Defines the output file for the factorization.\n On input are an integer and a string:\n 1) the integer equals the number of characters in the string; and\n 2) the string contains the name of a file.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_factor_standard_assign_labels",
     py2c_factor_standard_assign_labels, METH_VARARGS,
    "Assigns labels, replacing the multiplicity field of each solution\n in standard double precision stored in the container.\n On entry are two integers:\n 1) n, the number of coordinates of the solutions;\n 2) nbsols, the number of solutions in the container.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_factor_dobldobl_assign_labels",
     py2c_factor_dobldobl_assign_labels, METH_VARARGS,
    "Assigns labels, replacing the multiplicity field of each solution\n in double double precision stored in the container.\n On entry are two integers:\n 1) n, the number of coordinates of the solutions;\n 2) nbsols, the number of solutions in the container.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_factor_quaddobl_assign_labels",
     py2c_factor_quaddobl_assign_labels, METH_VARARGS,
    "Assigns labels, replacing the multiplicity field of each solution\n in quad double precision stored in the container.\n On entry are two integers:\n 1) n, the number of coordinates of the solutions;\n 2) nbsols, the number of solutions in the container.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_factor_initialize_standard_sampler",
     py2c_factor_initialize_standard_sampler, METH_VARARGS,
    "Initializes the sampling machine with a witness set,\n defined by an ordinary polynomial system in standard double precision.\n The embedded system is taken from the polynomial systems container\n and the generic points from the solutions container.\n On entry is the dimension or the number of hyperplanes\n to slide the positive dimensional solution set."},
   {"py2c_factor_initialize_dobldobl_sampler",
     py2c_factor_initialize_dobldobl_sampler, METH_VARARGS,
    "Initializes the sampling machine with a witness set,\n defined by an ordinary polynomial system in double double precision.\n The embedded system is taken from the polynomial systems container\n and the generic points from the solutions container.\n On entry is the dimension or the number of hyperplanes\n to slide the positive dimensional solution set."},
   {"py2c_factor_initialize_quaddobl_sampler",
     py2c_factor_initialize_quaddobl_sampler, METH_VARARGS,
    "Initializes the sampling machine with a witness set,\n defined by an ordinary polynomial system in quad double precision.\n The embedded system is taken from the polynomial systems container\n and the generic points from the solutions container.\n On entry is the dimension or the number of hyperplanes\n to slide the positive dimensional solution set."},
   {"py2c_factor_initialize_standard_Laurent_sampler",
     py2c_factor_initialize_standard_Laurent_sampler, METH_VARARGS,
    "Initializes the sampling machine with a witness set,\n defined by a Laurent polynomial system in standard double precision.\n The embedded system is taken from the Laurent systems container\n and the generic points from the solutions container.\n On entry is the dimension or the number of hyperplanes\n to slide the positive dimensional solution set."},
   {"py2c_factor_initialize_dobldobl_Laurent_sampler",
     py2c_factor_initialize_dobldobl_Laurent_sampler, METH_VARARGS,
    "Initializes the sampling machine with a witness set,\n defined by a Laurent polynomial system in double double precision.\n The embedded system is taken from the Laurent systems container\n and the generic points from the solutions container.\n On entry is the dimension or the number of hyperplanes\n to slide the positive dimensional solution set."},
   {"py2c_factor_initialize_quaddobl_Laurent_sampler",
     py2c_factor_initialize_quaddobl_Laurent_sampler, METH_VARARGS,
    "Initializes the sampling machine with a witness set,\n defined by a Laurent polynomial system in quad double precision.\n The embedded system is taken from the Laurent systems container\n and the generic points from the solutions container.\n On entry is the dimension or the number of hyperplanes\n to slide the positive dimensional solution set."},
   {"py2c_factor_initialize_standard_monodromy",
     py2c_factor_initialize_standard_monodromy, METH_VARARGS,
    "Initializes the internal data structures for n loops,\n to factor a k-dimensional solution component of degree d,\n in standard double precision.\n There are three integers on input, in the following order:\n 1) n, the number of loops;\n 2) d, the degree of the solution set;\n 3) k, the dimensional of the solution set.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_factor_initialize_dobldobl_monodromy",
     py2c_factor_initialize_dobldobl_monodromy, METH_VARARGS,
    "Initializes the internal data structures for n loops,\n to factor a k-dimensional solution component of degree d,\n in double double precision.\n There are three integers on input, in the following order:\n 1) n, the number of loops;\n 2) d, the degree of the solution set;\n 3) k, the dimensional of the solution set.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_factor_initialize_quaddobl_monodromy",
     py2c_factor_initialize_quaddobl_monodromy, METH_VARARGS,
    "Initializes the internal data structures for n loops,\n to factor a k-dimensional solution component of degree d,\n in quad double precision.\n There are three integers on input, in the following order:\n 1) n, the number of loops;\n 2) d, the degree of the solution set;\n 3) k, the dimensional of the solution set.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_factor_standard_trace_grid_diagnostics",
     py2c_factor_standard_trace_grid_diagnostics, METH_VARARGS,
    "Returns a tuple of two doubles with the diagnostics on the\n trace grid computed in standard double precision.\n The first double is the largest error of the samples.\n The second double is the smallest distance between two samples."},
   {"py2c_factor_dobldobl_trace_grid_diagnostics",
     py2c_factor_dobldobl_trace_grid_diagnostics, METH_VARARGS,
    "Returns a tuple of two doubles with the diagnostics on the\n trace grid computed in double double precision.\n The first double is the largest error of the samples.\n The second double is the smallest distance between two samples."},
   {"py2c_factor_quaddobl_trace_grid_diagnostics",
     py2c_factor_quaddobl_trace_grid_diagnostics, METH_VARARGS,
    "Returns a tuple of two doubles with the diagnostics on the\n trace grid computed in quad double precision.\n The first double is the largest error of the samples.\n The second double is the smallest distance between two samples."},
   {"py2c_factor_store_standard_solutions",
     py2c_factor_store_standard_solutions, METH_VARARGS,
    "Stores the solutions in the container, in standard double precision,\n to the data for monodromy loops."},
   {"py2c_factor_store_dobldobl_solutions",
     py2c_factor_store_dobldobl_solutions, METH_VARARGS,
    "Stores the solutions in the container, in double double precision,\n to the data for monodromy loops."},
   {"py2c_factor_store_quaddobl_solutions",
     py2c_factor_store_quaddobl_solutions, METH_VARARGS,
    "Stores the solutions in the container, in quad double precision,\n to the data for monodromy loops."},
   {"py2c_factor_restore_standard_solutions",
     py2c_factor_restore_standard_solutions, METH_VARARGS,
    "Restores the first initialized solutions, in standard double precision,\n from sampler to the container."},
   {"py2c_factor_restore_dobldobl_solutions",
     py2c_factor_restore_dobldobl_solutions, METH_VARARGS,
    "Restores the first initialized solutions, in double double precision,\n from sampler to the container."},
   {"py2c_factor_restore_quaddobl_solutions",
     py2c_factor_restore_quaddobl_solutions, METH_VARARGS,
    "Restores the first initialized solutions, in quad double precision,\n from sampler to the container."},
   {"py2c_factor_standard_track_paths",
     py2c_factor_standard_track_paths, METH_VARARGS,
    "Tracks as many paths as defined by witness set,\n in standard double precision.\n On input is an integer, which must be 1 if the witness set is\n defined by a Laurent polynomial system.\n On return is the failure code, which is zero when all went well."},
   {"py2c_factor_dobldobl_track_paths",
     py2c_factor_dobldobl_track_paths, METH_VARARGS,
    "Tracks as many paths as defined by witness set,\n in double double precision.\n On input is an integer, which must be 1 if the witness set is\n defined by a Laurent polynomial system.\n On return is the failure code, which is zero when all went well."},
   {"py2c_factor_quaddobl_track_paths",
     py2c_factor_quaddobl_track_paths, METH_VARARGS,
    "Tracks as many paths as defined by witness set,\n in quad double precision.\n On input is an integer, which must be 1 if the witness set is\n defined by a Laurent polynomial system.\n On return is the failure code, which is zero when all went well."},
   {"py2c_factor_swap_standard_slices",
     py2c_factor_swap_standard_slices, METH_VARARGS,
    "Swaps the current slices with new slices and takes new solutions\n as start to turn back, in standard double precision.\n On return is the failure code, which is zero when all went well."},
   {"py2c_factor_swap_dobldobl_slices",
     py2c_factor_swap_dobldobl_slices, METH_VARARGS,
    "Swaps the current slices with new slices and takes new solutions\n as start to turn back, in double double precision.\n On return is the failure code, which is zero when all went well."},
   {"py2c_factor_swap_quaddobl_slices",
     py2c_factor_swap_quaddobl_slices, METH_VARARGS,
    "Swaps the current slices with new slices and takes new solutions\n as start to turn back, in quad double precision.\n On return is the failure code, which is zero when all went well."},
   {"py2c_factor_new_standard_slices",
     py2c_factor_new_standard_slices, METH_VARARGS,
    "Generates k random slides in n-space, in standard double precision.\n The k and the n are the two input parameters.\n On return is the failure code, which is zero when all went well."},
   {"py2c_factor_new_dobldobl_slices",
     py2c_factor_new_dobldobl_slices, METH_VARARGS,
    "Generates k random slides in n-space, in double double precision.\n The k and the n are the two input parameters.\n On return is the failure code, which is zero when all went well."},
   {"py2c_factor_new_quaddobl_slices",
     py2c_factor_new_quaddobl_slices, METH_VARARGS,
    "Generates k random slides in n-space, in quad double precision.\n The k and the n are the two input parameters.\n On return is the failure code, which is zero when all went well."},
   {"py2c_factor_set_standard_trace_slice",
     py2c_factor_set_standard_trace_slice, METH_VARARGS,
    "Assigns the constant coefficient of the first slice,\n in standard double precision.\n On entry is a flag to indicate if it was the first time or not.\n On return is the failure code, which is zero if all went well."},
   {"py2c_factor_set_dobldobl_trace_slice",
     py2c_factor_set_dobldobl_trace_slice, METH_VARARGS,
    "Assigns the constant coefficient of the first slice,\n in double double precision.\n On entry is a flag to indicate if it was the first time or not.\n On return is the failure code, which is zero if all went well."},
   {"py2c_factor_set_quaddobl_trace_slice",
     py2c_factor_set_quaddobl_trace_slice, METH_VARARGS,
    "Assigns the constant coefficient of the first slice,\n in quad double precision.\n On entry is a flag to indicate if it was the first time or not.\n On return is the failure code, which is zero if all went well."},
   {"py2c_factor_store_standard_gammas",
     py2c_factor_store_standard_gammas, METH_VARARGS,
    "Stores the gamma constants in standard double precision\n for the sampler in the monodromy loops.\n Generates as many random complex constants as the value on input.\n On return is the failure code, which is zero if all went well."},
   {"py2c_factor_store_dobldobl_gammas",
     py2c_factor_store_dobldobl_gammas, METH_VARARGS,
    "Stores the gamma constants in double double precision\n for the sampler in the monodromy loops.\n Generates as many random complex constants as the value on input.\n On return is the failure code, which is zero if all went well."},
   {"py2c_factor_store_quaddobl_gammas",
     py2c_factor_store_quaddobl_gammas, METH_VARARGS,
    "Stores the gamma constants in quad double precision\n for the sampler in the monodromy loops.\n Generates as many random complex constants as the value on input.\n On return is the failure code, which is zero if all went well."},
   {"py2c_factor_permutation_after_standard_loop",
     py2c_factor_permutation_after_standard_loop, METH_VARARGS,
    "For a set of degree d, computes the permutation using the solutions\n most recently stored, after a loop in standard double precision.\n The number d is the input parameter of this function.\n On return is the string representation of the permutation."},
   {"py2c_factor_permutation_after_dobldobl_loop",
     py2c_factor_permutation_after_dobldobl_loop, METH_VARARGS,
    "For a set of degree d, computes the permutation using the solutions\n most recently stored, after a loop in double double precision.\n The number d is the input parameter of this function.\n On return is the string representation of the permutation."},
   {"py2c_factor_permutation_after_quaddobl_loop",
     py2c_factor_permutation_after_quaddobl_loop, METH_VARARGS,
    "For a set of degree d, computes the permutation using the solutions\n most recently stored, after a loop in quad double precision.\n The number d is the input parameter of this function.\n On return is the string representation of the permutation."},
   {"py2c_factor_update_standard_decomposition",
     py2c_factor_update_standard_decomposition, METH_VARARGS,
    "Updates the decomposition with the given permutation of d elements,\n computed in standard double precision.\n On entry are two integers and one string:\n 1) d, the number of elements in the permutation;\n 2) nc, the number of characters in the string;\n 3) p, the string representation of the permutation.\n Returns one if the current decomposition is certified,\n otherwise returns zero."},
   {"py2c_factor_update_dobldobl_decomposition",
     py2c_factor_update_dobldobl_decomposition, METH_VARARGS,
    "Updates the decomposition with the given permutation of d elements,\n computed in double double precision.\n On entry are two integers and one string:\n 1) d, the number of elements in the permutation;\n 2) nc, the number of characters in the string;\n 3) p, the string representation of the permutation.\n Returns one if the current decomposition is certified,\n otherwise returns zero."},
   {"py2c_factor_update_quaddobl_decomposition",
     py2c_factor_update_quaddobl_decomposition, METH_VARARGS,
    "Updates the decomposition with the given permutation of d elements,\n computed in quad double precision.\n On entry are two integers and one string:\n 1) d, the number of elements in the permutation;\n 2) nc, the number of characters in the string;\n 3) p, the string representation of the permutation.\n Returns one if the current decomposition is certified,\n otherwise returns zero."},
   {"py2c_factor_number_of_standard_components",
     py2c_factor_number_of_standard_components, METH_VARARGS,
    "Returns the number of irreducible factors in the current standard double\n precision decomposition of the witness set."},
   {"py2c_factor_number_of_dobldobl_components",
     py2c_factor_number_of_dobldobl_components, METH_VARARGS,
    "Returns the number of irreducible factors in the current double double\n precision decomposition of the witness set."},
   {"py2c_factor_number_of_quaddobl_components",
     py2c_factor_number_of_quaddobl_components, METH_VARARGS,
    "Returns the number of irreducible factors in the current quad double\n precision decomposition of the witness set."},
   {"py2c_factor_witness_points_of_standard_component",
     py2c_factor_witness_points_of_standard_component, METH_VARARGS,
    "Returns a string which represents an irreducible component,\n computed in standard double precision.\n On entry are two integers:\n 1) the sum of the degrees of all components;\n 2) the index of the component."},
   {"py2c_factor_witness_points_of_dobldobl_component",
     py2c_factor_witness_points_of_dobldobl_component, METH_VARARGS,
    "Returns a string which represents an irreducible component,\n computed in double double precision.\n On entry are two integers:\n 1) the sum of the degrees of all components;\n 2) the index of the component."},
   {"py2c_factor_witness_points_of_quaddobl_component",
     py2c_factor_witness_points_of_quaddobl_component, METH_VARARGS,
    "Returns a string which represents an irreducible component,\n computed in quad double precision.\n On entry are two integers:\n 1) the sum of the degrees of all components;\n 2) the index of the component."},
   {"py2c_factor_standard_trace_sum_difference",
     py2c_factor_standard_trace_sum_difference, METH_VARARGS,
    "Returns the difference between the actual sum at the samples\n defined by the labels to the generic points in the factor,\n and the trace sum, computed in standard double precision.\n On entry are three integer numbers and one string:\n 1) d, the number of points in the witness set;\n 2) k, the dimension of the solution set;\n 3) nc, the number of characters in the string;\n 4) ws, the string representing the labels of the witness set."},
   {"py2c_factor_dobldobl_trace_sum_difference",
     py2c_factor_dobldobl_trace_sum_difference, METH_VARARGS,
    "Returns the difference between the actual sum at the samples\n defined by the labels to the generic points in the factor,\n and the trace sum, computed in double double precision.\n On entry are three integer numbers and one string:\n 1) d, the number of points in the witness set;\n 2) k, the dimension of the solution set;\n 3) nc, the number of characters in the string;\n 4) ws, the string representing the labels of the witness set."},
   {"py2c_factor_quaddobl_trace_sum_difference",
     py2c_factor_quaddobl_trace_sum_difference, METH_VARARGS,
    "Returns the difference between the actual sum at the samples\n defined by the labels to the generic points in the factor,\n and the trace sum, computed in quad double precision.\n On entry are three integer numbers and one string:\n 1) d, the number of points in the witness set;\n 2) k, the dimension of the solution set;\n 3) nc, the number of characters in the string;\n 4) ws, the string representing the labels of the witness set."},
   {"py2c_witset_standard_membertest",
     py2c_witset_standard_membertest, METH_VARARGS,
    "Executes the homotopy membership test for a point to belong to\n a witness set defined by an ordinary polynomial system\n in standard double precision.\n The containers in standard double precision must contain the embedded\n polynomial system and its corresponding solutions for the witness set\n of a positive dimensional solution set.\n On entry are the seven parameters, the first four are integers:\n 1) vrb, an integer flag (0 or 1) for the verbosity of the test,\n 2) nvr, the ambient dimension, number of coordinates of the point,\n 3) dim, the dimension of the witness set,\n 4) nbc, the number of characters in the string representing the point;\n the next two parameters are two doubles:\n 5) restol, tolerance on the residual for the valuation of the point,\n 6) homtol, tolerance on the homotopy membership test for the point;\n and the last parameter is a string:\n 7) tpt, the string representation of the point as a list with as\n many as 2*nvr doubles for the real and imaginary parts of the\n standard double precision coordinates of the test point.\n On return are three 0/1 integers, to be interpreted as booleans:\n 1) fail, the failure code of the procedure,\n 2) onsys, 0 if the evaluation test failed, 1 if success,\n 3) onset, 0 if not a member of the witness set, 1 if a member."},
   {"py2c_witset_dobldobl_membertest",
     py2c_witset_dobldobl_membertest, METH_VARARGS,
    "Executes the homotopy membership test for a point to belong to\n a witness set defined by an ordinary polynomial system\n in double double precision.\n The containers in double double precision must contain the embedded\n polynomial system and its corresponding solutions for the witness set\n of a positive dimensional solution set.\n On entry are the seven parameters, the first four are integers:\n 1) vrb, an integer flag (0 or 1) for the verbosity of the test,\n 2) nvr, the ambient dimension, number of coordinates of the point,\n 3) dim, the dimension of the witness set,\n 4) nbc, the number of characters in the string representing the point;\n the next two parameters are two doubles:\n 5) restol, tolerance on the residual for the valuation of the point,\n 6) homtol, tolerance on the homotopy membership test for the point;\n and the last parameter is a string:\n 7) tpt, the string representation of the point as a list with as\n many as 4*nvr doubles for the real and imaginary parts of the\n double double precision coordinates of the test point.\n On return are three 0/1 integers, to be interpreted as booleans:\n 1) fail, the failure code of the procedure,\n 2) onsys, 0 if the evaluation test failed, 1 if success,\n 3) onset, 0 if not a member of the witness set, 1 if a member."},
   {"py2c_witset_quaddobl_membertest",
     py2c_witset_quaddobl_membertest, METH_VARARGS,
    "Executes the homotopy membership test for a point to belong to\n a witness set defined by an ordinary polynomial system\n in quad double precision.\n The containers in quad double precision must contain the embedded\n polynomial system and its corresponding solutions for the witness set\n of a positive dimensional solution set.\n On entry are the seven parameters, the first four are integers:\n 1) vrb, an integer flag (0 or 1) for the verbosity of the test,\n 2) nvr, the ambient dimension, number of coordinates of the point,\n 3) dim, the dimension of the witness set,\n 4) nbc, the number of characters in the string representing the point;\n the next two parameters are two doubles:\n 5) restol, tolerance on the residual for the valuation of the point,\n 6) homtol, tolerance on the homotopy membership test for the point;\n and the last parameter is a string:\n 7) tpt, the string representation of the point as a list with as\n many as 8*nvr doubles for the real and imaginary parts of the\n quad double precision coordinates of the test point.\n On return are three 0/1 integers, to be interpreted as booleans:\n 1) fail, the failure code of the procedure,\n 2) onsys, 0 if the evaluation test failed, 1 if success,\n 3) onset, 0 if not a member of the witness set, 1 if a member."},
   {"py2c_witset_standard_Laurent_membertest",
     py2c_witset_standard_Laurent_membertest, METH_VARARGS,
    "Executes the homotopy membership test for a point to belong to\n a witness set defined by a Laurent polynomial system\n in standard double precision.\n The containers in standard double precision must contain the embedded\n Laurent system and its corresponding solutions for the witness set\n of a positive dimensional solution set.\n On entry are the seven parameters, the first four are integers:\n 1) vrb, an integer flag (0 or 1) for the verbosity of the test,\n 2) nvr, the ambient dimension, number of coordinates of the point,\n 3) dim, the dimension of the witness set,\n 4) nbc, the number of characters in the string representing the point;\n the next two parameters are two doubles:\n 5) restol, tolerance on the residual for the valuation of the point,\n 6) homtol, tolerance on the homotopy membership test for the point;\n and the last parameter is a string:\n 7) tpt, the string representation of the point as a list with as\n many as 2*nvr doubles for the real and imaginary parts of the\n standard double precision coordinates of the test point.\n On return are three 0/1 integers, to be interpreted as booleans:\n 1) fail, the failure code of the procedure,\n 2) onsys, 0 if the evaluation test failed, 1 if success,\n 3) onset, 0 if not a member of the witness set, 1 if a member."},
   {"py2c_witset_dobldobl_Laurent_membertest",
     py2c_witset_dobldobl_Laurent_membertest, METH_VARARGS,
    "Executes the homotopy membership test for a point to belong to\n a witness set defined by a Laurent polynomial system\n in double double precision.\n The containers in double double precision must contain the embedded\n Laurent system and its corresponding solutions for the witness set\n of a positive dimensional solution set.\n On entry are the seven parameters, the first four are integers:\n 1) vrb, an integer flag (0 or 1) for the verbosity of the test,\n 2) nvr, the ambient dimension, number of coordinates of the point,\n 3) dim, the dimension of the witness set,\n 4) nbc, the number of characters in the string representing the point;\n the next two parameters are two doubles:\n 5) restol, tolerance on the residual for the valuation of the point,\n 6) homtol, tolerance on the homotopy membership test for the point;\n and the last parameter is a string:\n 7) tpt, the string representation of the point as a list with as\n many as 4*nvr doubles for the real and imaginary parts of the\n double double precision coordinates of the test point.\n On return are three 0/1 integers, to be interpreted as booleans:\n 1) fail, the failure code of the procedure,\n 2) onsys, 0 if the evaluation test failed, 1 if success,\n 3) onset, 0 if not a member of the witness set, 1 if a member."},
   {"py2c_witset_quaddobl_Laurent_membertest",
     py2c_witset_quaddobl_Laurent_membertest, METH_VARARGS,
    "Executes the homotopy membership test for a point to belong to\n a witness set defined by a Laurent polynomial system\n in quad double precision.\n The containers in quad double precision must contain the embedded\n Laurent system and its corresponding solutions for the witness set\n of a positive dimensional solution set.\n On entry are the seven parameters, the first four are integers:\n 1) vrb, an integer flag (0 or 1) for the verbosity of the test,\n 2) nvr, the ambient dimension, number of coordinates of the point,\n 3) dim, the dimension of the witness set,\n 4) nbc, the number of characters in the string representing the point;\n the next two parameters are two doubles:\n 5) restol, tolerance on the residual for the valuation of the point,\n 6) homtol, tolerance on the homotopy membership test for the point;\n and the last parameter is a string:\n 7) tpt, the string representation of the point as a list with as\n many as 8*nvr doubles for the real and imaginary parts of the\n quad double precision coordinates of the test point.\n On return are three 0/1 integers, to be interpreted as booleans:\n 1) fail, the failure code of the procedure,\n 2) onsys, 0 if the evaluation test failed, 1 if success,\n 3) onset, 0 if not a member of the witness set, 1 if a member."},
   {"py2c_witset_standard_ismember",
     py2c_witset_standard_ismember, METH_VARARGS,
    "Runs the homotopy membership test for a point to belong to a witness set\n defined by an ordinary polynomial system in standard double precision,\n where the test point is given as a string in PHCpack format.\n The containers in standard double precision must contain the\n embedded system and the corresponding generic points.\n On entry are seven parameters.  The first four are integers:\n 1) vrb, an integer flag (0 or 1) for the verbosity of the test,\n 2) nvr, the ambient dimension, number of coordinates of the point,\n 3) dim, the dimension of the witness set,\n 4) nbc, the number of characters in the string representing the point,\n the test point is represented as a solution string in symbolic format,\n including the symbols for the variables, before the coordinates;\n the next two parameters are two doubles:\n 5) restol, tolerance on the residual for the valuation of the point,\n 6) homtol, tolerance on the homotopy membership test for the point;\n and the last parameter is a string:\n 7) tpt, the string representation of a solution which contains\n the coordinates of the test point in symbolic format.\n On return are three 0/1 integers, to be interpreted as booleans:\n 1) fail, the failure code of the procedure,\n 2) onsys, 0 if the evaluation test failed, 1 if success,\n 3) onset, 0 if not a member of the witness set, 1 if a member."},
   {"py2c_witset_dobldobl_ismember",
     py2c_witset_dobldobl_ismember, METH_VARARGS,
    "Runs the homotopy membership test for a point to belong to a witness set\n defined by an ordinary polynomial system in double double precision,\n where the test point is given as a string in PHCpack format.\n The containers in double double precision must contain the\n embedded system and the corresponding generic points.\n On entry are seven parameters.  The first four are integers:\n 1) vrb, an integer flag (0 or 1) for the verbosity of the test,\n 2) nvr, the ambient dimension, number of coordinates of the point,\n 3) dim, the dimension of the witness set,\n 4) nbc, the number of characters in the string representing the point,\n the test point is represented as a solution string in symbolic format,\n including the symbols for the variables, before the coordinates;\n the next two parameters are two doubles:\n 5) restol, tolerance on the residual for the valuation of the point,\n 6) homtol, tolerance on the homotopy membership test for the point;\n and the last parameter is a string:\n 7) tpt, the string representation of a solution which contains\n the coordinates of the test point in symbolic format.\n On return are three 0/1 integers, to be interpreted as booleans:\n 1) fail, the failure code of the procedure,\n 2) onsys, 0 if the evaluation test failed, 1 if success,\n 3) onset, 0 if not a member of the witness set, 1 if a member."},
   {"py2c_witset_quaddobl_ismember",
     py2c_witset_quaddobl_ismember, METH_VARARGS,
    "Runs the homotopy membership test for a point to belong to a witness set\n defined by an ordinary polynomial system in quad double precision,\n where the test point is given as a string in PHCpack format.\n The containers in quad double precision must contain the\n embedded system and the corresponding generic points.\n On entry are seven parameters.  The first four are integers:\n 1) vrb, an integer flag (0 or 1) for the verbosity of the test,\n 2) nvr, the ambient dimension, number of coordinates of the point,\n 3) dim, the dimension of the witness set,\n 4) nbc, the number of characters in the string representing the point,\n the test point is represented as a solution string in symbolic format,\n including the symbols for the variables, before the coordinates;\n the next two parameters are two doubles:\n 5) restol, tolerance on the residual for the valuation of the point,\n 6) homtol, tolerance on the homotopy membership test for the point;\n and the last parameter is a string:\n 7) tpt, the string representation of a solution which contains\n the coordinates of the test point in symbolic format.\n On return are three 0/1 integers, to be interpreted as booleans:\n 1) fail, the failure code of the procedure,\n 2) onsys, 0 if the evaluation test failed, 1 if success,\n 3) onset, 0 if not a member of the witness set, 1 if a member."},
   {"py2c_witset_standard_Laurent_ismember",
     py2c_witset_standard_Laurent_ismember, METH_VARARGS,
    "Runs the homotopy membership test for a point to belong to a witness set\n defined by a Laurent polynomial system in standard double precision,\n where the test point is given as a string in PHCpack format.\n The containers in standard double precision must contain the\n embedded Laurent system and the corresponding generic points.\n On entry are seven parameters.  The first four are integers:\n 1) vrb, an integer flag (0 or 1) for the verbosity of the test,\n 2) nvr, the ambient dimension, number of coordinates of the point,\n 3) dim, the dimension of the witness set,\n 4) nbc, the number of characters in the string representing the point,\n the test point is represented as a solution string in symbolic format,\n including the symbols for the variables, before the coordinates;\n the next two parameters are two doubles:\n 5) restol, tolerance on the residual for the valuation of the point,\n 6) homtol, tolerance on the homotopy membership test for the point;\n and the last parameter is a string:\n 7) tpt, the string representation of a solution which contains\n the coordinates of the test point in symbolic format.\n On return are three 0/1 integers, to be interpreted as booleans:\n 1) fail, the failure code of the procedure,\n 2) onsys, 0 if the evaluation test failed, 1 if success,\n 3) onset, 0 if not a member of the witness set, 1 if a member."},
   {"py2c_witset_dobldobl_Laurent_ismember",
     py2c_witset_dobldobl_Laurent_ismember, METH_VARARGS,
    "Runs the homotopy membership test for a point to belong to a witness set\n defined by a Laurent polynomial system in double double precision,\n where the test point is given as a string in PHCpack format.\n The containers in double double precision must contain the\n embedded Laurent system and the corresponding generic points.\n On entry are seven parameters.  The first four are integers:\n 1) vrb, an integer flag (0 or 1) for the verbosity of the test,\n 2) nvr, the ambient dimension, number of coordinates of the point,\n 3) dim, the dimension of the witness set,\n 4) nbc, the number of characters in the string representing the point,\n the test point is represented as a solution string in symbolic format,\n including the symbols for the variables, before the coordinates;\n the next two parameters are two doubles:\n 5) restol, tolerance on the residual for the valuation of the point,\n 6) homtol, tolerance on the homotopy membership test for the point;\n and the last parameter is a string:\n 7) tpt, the string representation of a solution which contains\n the coordinates of the test point in symbolic format.\n On return are three 0/1 integers, to be interpreted as booleans:\n 1) fail, the failure code of the procedure,\n 2) onsys, 0 if the evaluation test failed, 1 if success,\n 3) onset, 0 if not a member of the witness set, 1 if a member."},
   {"py2c_witset_quaddobl_Laurent_ismember",
     py2c_witset_quaddobl_Laurent_ismember, METH_VARARGS,
    "Runs the homotopy membership test for a point to belong to a witness set\n defined by a Laurent polynomial system in quad double precision,\n where the test point is given as a string in PHCpack format.\n The containers in quad double precision must contain the\n embedded Laurent system and the corresponding generic points.\n On entry are seven parameters.  The first four are integers:\n 1) vrb, an integer flag (0 or 1) for the verbosity of the test,\n 2) nvr, the ambient dimension, number of coordinates of the point,\n 3) dim, the dimension of the witness set,\n 4) nbc, the number of characters in the string representing the point,\n the test point is represented as a solution string in symbolic format,\n including the symbols for the variables, before the coordinates;\n the next two parameters are two doubles:\n 5) restol, tolerance on the residual for the valuation of the point,\n 6) homtol, tolerance on the homotopy membership test for the point;\n and the last parameter is a string:\n 7) tpt, the string representation of a solution which contains\n the coordinates of the test point in symbolic format.\n On return are three 0/1 integers, to be interpreted as booleans:\n 1) fail, the failure code of the procedure,\n 2) onsys, 0 if the evaluation test failed, 1 if success,\n 3) onset, 0 if not a member of the witness set, 1 if a member."},
   {"py2c_standard_witset_of_hypersurface",
     py2c_standard_witset_of_hypersurface, METH_VARARGS,
    "Given in the string p of nc characters a polynomial in nv variables,\n terminated by a semicolon, the systems and solutions container\n in standard double precision on return contain a witness set for\n the hypersurface defined by the ordinary polynomial in p.\n On entry are two integers and one string, in the following order:\n 1) nv, the number of variables of the polynomials;\n 2) nc, the number of characters in the string p;\n 3) p, string representation of an ordinary polynomial in several\n variables, terminates with a semicolon.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_dobldobl_witset_of_hypersurface",
     py2c_dobldobl_witset_of_hypersurface, METH_VARARGS,
    "Given in the string p of nc characters a polynomial in nv variables,\n terminated by a semicolon, the systems and solutions container\n in double double precision on return contain a witness set for\n the hypersurface defined by the ordinary polynomial in p.\n On entry are two integers and one string, in the following order:\n 1) nv, the number of variables of the polynomials;\n 2) nc, the number of characters in the string p;\n 3) p, string representation of an ordinary polynomial in several\n variables, terminates with a semicolon.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_quaddobl_witset_of_hypersurface",
     py2c_quaddobl_witset_of_hypersurface, METH_VARARGS,
    "Given in the string p of nc characters a polynomial in nv variables,\n terminated by a semicolon, the systems and solutions container\n in quad double precision on return contain a witness set for\n the hypersurface defined by the ordinary polynomial in p.\n On entry are two integers and one string, in the following order:\n 1) nv, the number of variables of the polynomials;\n 2) nc, the number of characters in the string p;\n 3) p, string representation of an ordinary polynomial in several\n variables, terminates with a semicolon.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_standard_witset_of_Laurent_hypersurface",
     py2c_standard_witset_of_Laurent_hypersurface, METH_VARARGS,
    "Given in the string p of nc characters a polynomial in nv variables,\n terminated by a semicolon, the systems and solutions container\n in standard double precision on return contain a witness set for\n the hypersurface defined by the Laurent polynomial in p.\n On entry are two integers and one string, in the following order:\n 1) nv, the number of variables of the polynomials;\n 2) nc, the number of characters in the string p;\n 3) p, string representation of a Laurent polynomial in several\n variables, terminates with a semicolon.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_dobldobl_witset_of_Laurent_hypersurface",
     py2c_dobldobl_witset_of_Laurent_hypersurface, METH_VARARGS,
    "Given in the string p of nc characters a polynomial in nv variables,\n terminated by a semicolon, the systems and solutions container\n in double double precision on return contain a witness set for\n the hypersurface defined by the Laurent polynomial in p.\n On entry are two integers and one string, in the following order:\n 1) nv, the number of variables of the polynomials;\n 2) nc, the number of characters in the string p;\n 3) p, string representation of a Laurent polynomial in several\n variables, terminates with a semicolon.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_quaddobl_witset_of_Laurent_hypersurface",
     py2c_quaddobl_witset_of_Laurent_hypersurface, METH_VARARGS,
    "Given in the string p of nc characters a polynomial in nv variables,\n terminated by a semicolon, the systems and solutions container\n in quad double precision on return contain a witness set for\n the hypersurface defined by the Laurent polynomial in p.\n On entry are two integers and one string, in the following order:\n 1) nv, the number of variables of the polynomials;\n 2) nc, the number of characters in the string p;\n 3) p, string representation of a Laurent polynomial in several\n variables, terminates with a semicolon.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_standard_diagonal_homotopy", py2c_standard_diagonal_homotopy,
     METH_VARARGS,
    "Creates a diagonal homotopy to intersect two solution sets of\n dimensions a and b respectively, where a >= b.\n The two input parameters are values for a and b.\n The systems stored as target and start system in the container,\n in standard double precision, define the witness sets for\n these two solution sets."},
   {"py2c_dobldobl_diagonal_homotopy", py2c_dobldobl_diagonal_homotopy,
     METH_VARARGS,
    "Creates a diagonal homotopy to intersect two solution sets of\n dimensions a and b respectively, where a >= b.\n The two input parameters are values for a and b.\n The systems stored as target and start system in the container,\n in double double precision, define the witness sets for\n these two solution sets."},
   {"py2c_quaddobl_diagonal_homotopy", py2c_quaddobl_diagonal_homotopy,
     METH_VARARGS,
    "Creates a diagonal homotopy to intersect two solution sets of\n dimensions a and b respectively, where a >= b.\n The two input parameters are values for a and b.\n The systems stored as target and start system in the container,\n in quad double precision, define the witness sets for\n these two solution sets."},
   {"py2c_standard_diagonal_cascade_solutions",
     py2c_standard_diagonal_cascade_solutions, METH_VARARGS,
    "Makes the start solutions to start the cascade homotopy to\n intersect two solution sets of dimensions a and b, where a >= b,\n in standard double precision.\n The dimensions a and b are given as input parameters.\n The systems stored as target and start system in the container\n define the witness sets for these two solution sets.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_dobldobl_diagonal_cascade_solutions",
     py2c_dobldobl_diagonal_cascade_solutions, METH_VARARGS,
    "Makes the start solutions to start the cascade homotopy to\n intersect two solution sets of dimensions a and b, where a >= b,\n in double double precision.\n The dimensions a and b are given as input parameters.\n The systems stored as target and start system in the container\n define the witness sets for these two solution sets.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_quaddobl_diagonal_cascade_solutions",
     py2c_quaddobl_diagonal_cascade_solutions, METH_VARARGS,
    "Makes the start solutions to start the cascade homotopy to\n intersect two solution sets of dimensions a and b, where a >= b,\n in quad double precision.\n The dimensions a and b are given as input parameters.\n The systems stored as target and start system in the container\n define the witness sets for these two solution sets.\n On return is the failure code, which equals zero when all went well."},
   {"py2c_extrinsic_top_diagonal_dimension",
     py2c_extrinsic_top_diagonal_dimension, METH_VARARGS,
    "Returns the dimension of the start and target system to\n start the extrinsic cascade to intersect two witness sets,\n respectively of dimensions a and b, with ambient dimensions\n respectively equal to n1 and n2.\n There are four integers as parameters on input: n1, n2, a and b."},
   {"py2c_diagonal_symbols_doubler",
     py2c_diagonal_symbols_doubler, METH_VARARGS, 
    "Doubles the number of symbols in the symbol table to enable the\n writing of the target system to string properly when starting the\n cascade of a diagonal homotopy in extrinsic coordinates.\n On input are three integers, n, d, nc, and one string s.\n On input are n, the ambient dimension = #variables before the embedding,\n d is the number of slack variables, or the dimension of the first set,\n and in s (nc characters) are the symbols for the first witness set.\n This function takes the symbols in s and combines those symbols with\n those in the current symbol table for the second witness set stored\n in the standard systems container.  On return, the symbol table\n contains then all symbols to write the top system in the cascade\n to start the diagonal homotopy."},
   {"py2c_standard_collapse_diagonal",
     py2c_standard_collapse_diagonal, METH_VARARGS,
    "Eliminates the extrinsic diagonal for the system and solutions\n in the containers for standard doubles.  On input are two integers:\n 1) k, the current number of slack variables in the embedding;\n 2) d, the number of slack variables to add to the final embedding.\n The system in the container has its diagonal eliminated and is\n embedded with k+d slack variables.  The solutions corresponding\n to this system are in the solutions container.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_dobldobl_collapse_diagonal",
     py2c_dobldobl_collapse_diagonal, METH_VARARGS,
    "Eliminates the extrinsic diagonal for the system and solutions\n in the containers for double doubles.  On input are two integers:\n 1) k, the current number of slack variables in the embedding;\n 2) d, the number of slack variables to add to the final embedding.\n The system in the container has its diagonal eliminated and is\n embedded with k+d slack variables.  The solutions corresponding\n to this system are in the solutions container.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_quaddobl_collapse_diagonal",
     py2c_quaddobl_collapse_diagonal, METH_VARARGS,
    "Eliminates the extrinsic diagonal for the system and solutions\n in the containers for quad doubles.  On input are two integers:\n 1) k, the current number of slack variables in the embedding;\n 2) d, the number of slack variables to add to the final embedding.\n The system in the container has its diagonal eliminated and is\n embedded with k+d slack variables.  The solutions corresponding\n to this system are in the solutions container.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_standard_polysys_solve",
     py2c_standard_polysys_solve, METH_VARARGS,
    "Runs the cascades of homotopies on the polynomial system in\n the standard systems container.  Runs in standard double precision.\n On input are five integers :\n 1) nbtasks equals the number of tasks for multitasking,\n 2) topdim is the top dimension to start the homotopy cascades,\n 3) filter is a 0 or 1 flag to filter the witness supersets,\n 4) factor is a 0 or 1 flag to factor the witness sets,\n 5) verbose is a flag for intermediate output."},
   {"py2c_standard_laursys_solve",
     py2c_standard_laursys_solve, METH_VARARGS,
    "Runs the cascades of homotopies on the Laurent polynomial system in\n the standard systems container.  Runs in standard double precision.\n On input are five integers :\n 1) nbtasks equals the number of tasks for multitasking,\n 2) topdim is the top dimension to start the homotopy cascades,\n 3) filter is a 0 or 1 flag to filter the witness supersets,\n 4) factor is a 0 or 1 flag to factor the witness sets,\n 5) verbose is a flag for intermediate output."},
   {"py2c_dobldobl_polysys_solve",
     py2c_dobldobl_polysys_solve, METH_VARARGS,
    "Runs the cascades of homotopies on the polynomial system in\n the dobldobl systems container.  Runs in double double precision.\n On input are five integers :\n 1) nbtasks equals the number of tasks for multitasking,\n 2) topdim is the top dimension to start the homotopy cascades,\n 3) filter is a 0 or 1 flag to filter the witness supersets,\n 4) factor is a 0 or 1 flag to factor the witness sets,\n 5) verbose is a flag for intermediate output."},
   {"py2c_dobldobl_laursys_solve",
     py2c_dobldobl_laursys_solve, METH_VARARGS,
    "Runs the cascades of homotopies on the Laurent polynomial system in\n the dobldobl systems container.  Runs in double double precision.\n On input are five integers :\n 1) nbtasks equals the number of tasks for multitasking,\n 2) topdim is the top dimension to start the homotopy cascades,\n 3) filter is a 0 or 1 flag to filter the witness supersets,\n 4) factor is a 0 or 1 flag to factor the witness sets,\n 5) verbose is a flag for intermediate output."},
   {"py2c_quaddobl_polysys_solve",
     py2c_quaddobl_polysys_solve, METH_VARARGS,
    "Runs the cascades of homotopies on the polynomial system in\n the quaddobl systems container.  Runs in quad double precision.\n On input are five integers :\n 1) nbtasks equals the number of tasks for multitasking,\n 2) topdim is the top dimension to start the homotopy cascades,\n 3) filter is a 0 or 1 flag to filter the witness supersets,\n 4) factor is a 0 or 1 flag to factor the witness sets,\n 5) verbose is a flag for intermediate output."},
   {"py2c_quaddobl_laursys_solve",
     py2c_quaddobl_laursys_solve, METH_VARARGS,
    "Runs the cascades of homotopies on the Laurent polynomial system in\n the quaddobl systems container.  Runs in quad double precision.\n On input are five integers :\n 1) nbtasks equals the number of tasks for multitasking,\n 2) topdim is the top dimension to start the homotopy cascades,\n 3) filter is a 0 or 1 flag to filter the witness supersets,\n 4) factor is a 0 or 1 flag to factor the witness sets,\n 5) verbose is a flag for intermediate output."},
   {"py2c_copy_standard_polysys_witset",
     py2c_copy_standard_polysys_witset, METH_VARARGS,
    "There is one integer parameter dim on input,\n which represents the dimension of the witness set.\n Copies the witness set representation for a solution set\n of dimension dim into the systems and solutions container,\n in standard double precision.  REQUIRED :\n 1) py2c_standard_polysys_solve was executed successfully, and\n 2) dim is in the range 0..topdim."},
   {"py2c_copy_standard_laursys_witset",
     py2c_copy_standard_laursys_witset, METH_VARARGS,
    "There is one integer parameter dim on input,\n which represents the dimension of the witness set.\n Copies the witness set representation for a solution set\n of dimension dim into the Laurent systems and solutions container,\n in standard double precision.  REQUIRED :\n 1) py2c_standard_laursys_solve was executed successfully, and\n 2) dim is in the range 0..topdim."},
   {"py2c_copy_dobldobl_polysys_witset",
     py2c_copy_dobldobl_polysys_witset, METH_VARARGS,
    "There is one integer parameter dim on input,\n which represents the dimension of the witness set.\n Copies the witness set representation for a solution set\n of dimension dim into the systems and solutions container,\n in double double precision.  REQUIRED :\n 1) py2c_dobldobl_polysys_solve was executed successfully, and\n 2) dim is in the range 0..topdim."},
   {"py2c_copy_dobldobl_laursys_witset",
     py2c_copy_dobldobl_laursys_witset, METH_VARARGS,
    "There is one integer parameter dim on input,\n which represents the dimension of the witness set.\n Copies the witness set representation for a solution set\n of dimension dim into the Laurent systems and solutions container,\n in double double precision.  REQUIRED :\n 1) py2c_dobldobl_laursys_solve was executed successfully, and\n 2) dim is in the range 0..topdim."},
    {"py2c_copy_quaddobl_polysys_witset",
      py2c_copy_quaddobl_polysys_witset, METH_VARARGS,
     "There is one integer parameter dim on input,\n which represents the dimension of the witness set.\n Copies the witness set representation for a solution set\n of dimension dim into the systems and solutions container,\n in quad double precision.  REQUIRED :\n 1) py2c_quaddobl_polysys_solve was executed successfully, and\n 2) dim is in the range 0..topdim."},
    {"py2c_copy_quaddobl_laursys_witset",
      py2c_copy_quaddobl_laursys_witset, METH_VARARGS,
     "There is one integer parameter dim on input,\n which represents the dimension of the witness set.\n Copies the witness set representation for a solution set\n of dimension dim into the Laurent systems and solutions container,\n in quad double precision.  REQUIRED :\n 1) py2c_quaddobl_laursys_solve was executed successfully, and\n 2) dim is in the range 0..topdim."},
   {"py2c_clear_standard_witsols",
     py2c_clear_standard_witsols, METH_VARARGS,
    "Clears the witness solutions in standard double precision."},
   {"py2c_clear_dobldobl_witsols",
     py2c_clear_dobldobl_witsols, METH_VARARGS,
    "Clears the witness solutions in double double precision."},
   {"py2c_clear_quaddobl_witsols",
     py2c_clear_quaddobl_witsols, METH_VARARGS,
    "Clears the witness solutions in quad double precision."},
   {"py2c_schubert_pieri_count", py2c_schubert_pieri_count, METH_VARARGS,
    "Returns the number of p-plane producing curves of degree q\n that meet m*p + q*(m+p) given general m-planes.\n On input are three integer numbers:\n 1) m, the dimension of the input planes;\n 2) p, the dimension of the output planes; and\n 3) q, the degree of the curve that produces p-planes.\n The dimension of the ambient space of this Pieri problem is m+p."},
   {"py2c_schubert_resolve_conditions", py2c_schubert_resolve_conditions,
     METH_VARARGS,
    "Resolves a general Schubert intersection condition in n-space\n for k-planes subject to conditions defined by brackers.\n On return is the root count, the number of k-planes that satisfy\n the intersection conditions imposed by the brackets for general flags.\n On entry are five integers and one string:\n 1) n, the ambient dimension, where the k-planes live;\n 2) k, the dimension of the solution planes;\n 3) c, the number of intersection conditions;\n 4) nc, the number of characters in the string brackets;\n 5) brackets is a string representation of c brackets, where the numbers\n in each bracket are separated by spaces;\n 6) the flag verbose: when 0, no intermediate output is written,\n when 1, then the resolution is dispayed on screen."},
   {"py2c_schubert_standard_littlewood_richardson_homotopies",
     py2c_schubert_standard_littlewood_richardson_homotopies, METH_VARARGS,
    "Runs the Littlewood-Richardson homotopies to resolve a number of\n general Schubert intersection conditions on k-planes in n-space,\n in standard double precision.\n The polynomial system that was solved is in the container for\n systems with coefficients in standard double precision and the\n corresponding solutions are in the standard solutions container.\n On entry are seven integers and two strings, in the following order:\n 1) n, the ambient dimension, where the k-planes live;\n 2) k, the dimension of the solution planes;\n 3) c,the number of intersection conditions;\n 4) nc, the number of characters in the string brackets;\n 5) brackets is a string representation of c brackets, where the numbers\n in each bracket are separated by spaces;\n 6) the flag verbose: if 0, then no intermediate output is written,\n when 1, then the resolution is dispayed on screen;\n 7) the flag verify: if 0, then diagnostic verification is done,\n if 1, then diagnostic verification is written to file;\n 8) the flag minrep: if 0, then all minors are used in the system,\n if 1, then a minimal representation of the problem is used;\n 9) the flag tosquare: if 0, then Gauss-Newton path trackers run,\n if 1, then the overdetermined systems are squared;\n 10) nbchar, the number of characters in the string filename;\n 11) filename is the name of the output file.\n The function returns a tuple of an integer and a string:\n 0) r is the formal root count as the number of k-planes\n for conditions imposed by the brackets for general flags;\n 1) flags, a string with the coefficients of the general flags."},
   {"py2c_schubert_dobldobl_littlewood_richardson_homotopies",
     py2c_schubert_dobldobl_littlewood_richardson_homotopies, METH_VARARGS,
    "Runs the Littlewood-Richardson homotopies to resolve a number of\n general Schubert intersection conditions on k-planes in n-space,\n in double double precision.\n The polynomial system that was solved is in the container for\n systems with coefficients in double double precision and the\n corresponding solutions are in the dobldobl solutions container.\n On entry are seven integers and two strings, in the following order:\n 1) n, the ambient dimension, where the k-planes live;\n 2) k, the dimension of the solution planes;\n 3) c,the number of intersection conditions;\n 4) nc, the number of characters in the string brackets;\n 5) brackets is a string representation of c brackets, where the numbers\n in each bracket are separated by spaces;\n 6) the flag verbose: if 0, then no intermediate output is written,\n when 1, then the resolution is dispayed on screen;\n 7) the flag verify: if 0, then no diagnostic verification is done,\n if 1, then diagnostic verification is written to file;\n 8) the flag minrep: if 0, then all minors are used in the system,\n if 1, then a minimal representation of the problem is used;\n 9) the flag tosquare: if 0, then Gauss-Newton path trackers run,\n if 1, then the overdetermined systems are squared;\n 10) nbchar, the number of characters in the string filename;\n 11) filename is the name of the output file.\n The function returns a tuple of an integer and a string:\n 0) r is the formal root count as the number of k-planes\n for conditions imposed by the brackets for general flags;\n 1) flags, a string with the coefficients of the general flags."},
   {"py2c_schubert_quaddobl_littlewood_richardson_homotopies",
     py2c_schubert_quaddobl_littlewood_richardson_homotopies, METH_VARARGS,
    "Runs the Littlewood-Richardson homotopies to resolve a number of\n general Schubert intersection conditions on k-planes in n-space,\n in quad double precision.\n The polynomial system that was solved is in the container for\n systems with coefficients in quad double precision and the\n corresponding solutions are in the quaddobl solutions container.\n On entry are seven integers and two strings, in the following order:\n 1) n, the ambient dimension, where the k-planes live;\n 2) k, the dimension of the solution planes;\n 3) c,the number of intersection conditions;\n 4) nc, the number of characters in the string brackets;\n 5) brackets is a string representation of c brackets, where the numbers\n in each bracket are separated by spaces;\n 6) the flag verbose: if 0, then no intermediate output is written,\n if 1, then the resolution is dispayed on screen;\n 7) the flag verify: if 0, then no diagnostic verification is done,\n if 1, then diagnostic verification is written to file;\n 8) the flag minrep: if 0, then all minors are used in the system,\n if 1, then a minimal representation of the problem is used;\n 9) the flag tosquare: if 0, then Gauss-Newton path trackers run,\n if 1, then the overdetermined systems are squared;\n 10) nbchar, the number of characters in the string filename;\n 11) filename is the name of the output file.\n The function returns a tuple of an integer and a string:\n 0) r is the formal root count as the number of k-planes\n for conditions imposed by the brackets for general flags;\n 1) flags, a string with the coefficients of the general flags."},
   {"py2c_schubert_localization_poset", py2c_schubert_localization_poset,
     METH_VARARGS,
    "Returns the string representation of the localization poset for the\n Pieri root count for m, p, and q.  The input parameters are the\n integer values for m, p, and q:\n 1) m, the dimension of the input planes;\n 2) p, the dimension of the output planes;\n 3) q, the degree of the curves that produce p-planes."},
   {"py2c_schubert_pieri_homotopies", py2c_schubert_pieri_homotopies,
    METH_VARARGS,
   "Runs the Pieri homotopies for (m,p,q) dimensions on generic input data.\n On return the systems container for systems with coefficients in standard\n double precision contains the polynomial system solved and in the\n solutions in standard double precision are in the solutions container.\n On entry are four integers and two strings:\n 1) m, the dimension of the input planes;\n 2) p, the dimension of the output planes;\n 3) q, the degree of the solution maps;\n 4) nc, the number of characters in the string A;\n 5) A, the string with m*p + q*(m+p) random complex input m-planes,\n where the real and imaginary parts are separated by a space;\n 6) pts, the string with m*p + q*(m+p) random complex interpolation\n points, only needed if q > 0.\n The function returns the combinatorial Pieri root count,\n which should equal the number of solutions in the container."},
   {"py2c_schubert_osculating_planes", py2c_schubert_osculating_planes,
     METH_VARARGS, 
    "Returns the string representation of n real m-planes in\n d-space osculating a rational normal curve\n at the n points in s, where n = m*p + q*(m+p) and d = m+p.\n On entry are four integers and one string:\n 1) m, the dimension of the input planes;\n 2) p, the dimension of the output planes;\n 3) q, the degree of the solution maps;\n 4) nc, the number of characters in the string pts; and\n 5) pts, the string with m*p + q*(m+p) interpolation points."},
   {"py2c_schubert_pieri_system", py2c_schubert_pieri_system, METH_VARARGS,
    "Fills the container of systems with coefficients in standard\n double precision with a polynomial system that expresses the\n intersection conditions of a general Pieri problem.\n On input are five integers and one string:\n 1) m, the dimension of the input planes;\n 2) p, the dimension of the output planes;\n 3) q, the degree of the solution maps;\n 4) nc, the number of characters in the string A;\n 5) A,  m*p + q*(m+p) random complex input m-planes, where\n the real and imaginary parts are separated by a space;\n 6) a flag is_real: if == 1, then the coefficients of A are real,\n if == 0, then the coefficients of A are complex.\n Returns the failure code, which equals zero if all went well."},
   {"py2c_mapcon_solve_system", py2c_mapcon_solve_system, METH_VARARGS,
    "Solves the binomial system stored in the Laurent systems container.\n There is one input argument, either one or zero.\n If one, then only the pure top dimensional solutions are computed.\n If zero, then all solution sets are computed.\n Returns the failure code, which equals zero if all went well."},
   {"py2c_mapcon_write_maps", py2c_mapcon_write_maps, METH_VARARGS,
    "Writes the maps stored in the container to screen.\n Returns the failure code, which equals zero if all went well."},
   {"py2c_mapcon_clear_maps", py2c_mapcon_clear_maps, METH_VARARGS, 
    "Deallocates the maps stored in the container.\n Returns the failure code, which equals zero if all went well."},
   {"py2c_mapcon_top_dimension", py2c_mapcon_top_dimension, METH_VARARGS,
    "Returns the top dimension of the maps in the container."},
   {"py2c_mapcon_number_of_maps", py2c_mapcon_number_of_maps, METH_VARARGS, 
    "Returns the number of maps in the container."},
   {"py2c_mapcon_degree_of_map", py2c_mapcon_degree_of_map, METH_VARARGS,
    "Given the dimension and index of a map, given as two integers as\n input parameters, returns the degree of that map."},
   {"py2c_mapcon_coefficients_of_map", py2c_mapcon_coefficients_of_map,
     METH_VARARGS,
    "Returns the coefficients of a monomial map stored in the container.\n On entry are three parameters:\n 1) the dimension of the map;\n 2) the index of the map in all maps of that dimension;\n 3) the number of variables.\n On return is a Python list of complex doubles."},
   {"py2c_mapcon_exponents_of_map", py2c_mapcon_exponents_of_map, METH_VARARGS,
    "Returns the exponents of a monomial map stored in the container.\n On entry are three parameters:\n 1) the dimension of the map;\n 2) the index of the map in all maps of that dimension;\n 3) the number of variables.\n On return is a Python list of integers."},
   {"py2c_standard_Newton_series", py2c_standard_Newton_series, METH_VARARGS,
    "Given in the systems container a polynomial system with coefficients\n in standard double precision, and in the solutions container the\n leading coefficients of the power series, this function runs Newton's\n method to compute power series solutions of the system in the container,\n in standard double precision.  There are four integers on input:\n 1) the index of the series parameter;\n 2) the maximal degree of the series;\n 3) the number of Newton steps to be done on each solution;\n 4) a 0/1-flag to indicate whether additional diagnostic output needs\n to be written to screen.\n The solution series are stored in the standard systems pool.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_dobldobl_Newton_series", py2c_dobldobl_Newton_series, METH_VARARGS,
    "Given in the systems container a polynomial system with coefficients\n in standard double precision, and in the solutions container the\n leading coefficients of the power series, this function runs Newton's\n method to compute power series solutions of the system in the container,\n in double double precision.  There are four integers on input:\n 1) the index of the series parameter;\n 2) the maximal degree of the series;\n 3) the number of Newton steps to be done on each solution;\n 4) a 0/1-flag to indicate whether additional diagnostic output needs\n to be written to screen.\n The solution series are stored in the dobldobl systems pool.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_quaddobl_Newton_series", py2c_quaddobl_Newton_series, METH_VARARGS,
    "Given in the systems container a polynomial system with coefficients\n in standard double precision, and in the solutions container the\n leading coefficients of the power series, this function runs Newton's\n method to compute power series solutions of the system in the container,\n in quad double precision.  There are four integers on input:\n 1) the index of the series parameter;\n 2) the maximal degree of the series;\n 3) the number of Newton steps to be done on each solution;\n 4) a 0/1-flag to indicate whether additional diagnostic output needs\n to be written to screen.\n The solution series are stored in the quaddobl systems pool.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_standard_Newton_power_series",
     py2c_standard_Newton_power_series, METH_VARARGS,
    "Given in the systems container a polynomial system with coefficients\n in standard double precision, and in the standard systems pool the\n leading terms of the power series, this function runs Newton's\n method to compute power series solutions of the system in the container,\n in standard double precision.  There are four integers on input:\n 1) the index of the series parameter;\n 2) the maximal degree of the series;\n 3) the number of Newton steps to be done on each solution;\n 4) a 0/1-flag to indicate whether additional diagnostic output needs\n to be written to screen.\n The solution series are stored in the standard systems pool.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_dobldobl_Newton_power_series",
     py2c_dobldobl_Newton_power_series, METH_VARARGS,
    "Given in the systems container a polynomial system with coefficients\n in standard double precision, and in the dobldobl systems pool the\n leading terms of the power series, this function runs Newton's\n method to compute power series solutions of the system in the container,\n in double double precision.  There are four integers on input:\n 1) the index of the series parameter;\n 2) the maximal degree of the series;\n 3) the number of Newton steps to be done on each solution;\n 4) a 0/1-flag to indicate whether additional diagnostic output needs\n to be written to screen.\n The solution series are stored in the dobldobl systems pool.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_quaddobl_Newton_power_series",
     py2c_quaddobl_Newton_power_series, METH_VARARGS,
    "Given in the systems container a polynomial system with coefficients\n in standard double precision, and in the quaddobl systems pool the\n leading terms of the power series, this function runs Newton's\n method to compute power series solutions of the system in the container,\n in quad double precision.  There are four integers on input:\n 1) the index of the series parameter;\n 2) the maximal degree of the series;\n 3) the number of Newton steps to be done on each solution;\n 4) a 0/1-flag to indicate whether additional diagnostic output needs\n to be written to screen.\n The solution series are store in the quaddobl systems pool.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_standard_Pade_approximant",
     py2c_standard_Pade_approximant, METH_VARARGS,
    "Given in the systems container a polynomial system with coefficients\n in standard double precision, and in the solutions container the\n leading coefficients of the power series, this function runs Newton's\n method to compute power series solutions of the system in the container,\n in standard double precision, followed by the construction of the\n Pade approximants, for each solution. There are five integers on input:\n 1) the index of the series parameter;\n 2) the degree of the numerator of the Pade approximant;\n 3) the degree of the denominator of the Pade approximant;\n 4) the number of Newton steps to be done on each solution;\n 5) a 0/1-flag to indicate whether additional diagnostic output needs\n to be written to screen.\n The Pade approximants are stored in the standard systems pool,\n numerators in the odd indexed entries and denominators in the entries\n with even index in each system.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_dobldobl_Pade_approximant",
     py2c_dobldobl_Pade_approximant, METH_VARARGS,
    "Given in the systems container a polynomial system with coefficients\n in double double precision, and in the solutions container the\n leading coefficients of the power series, this function runs Newton's\n method to compute power series solutions of the system in the container,\n in double double precision, followed by the construction of the\n Pade approximants, for each solution. There are five integers on input:\n 1) the index of the series parameter;\n 2) the degree of the numerator of the Pade approximant;\n 3) the degree of the denominator of the Pade approximant;\n 4) the number of Newton steps to be done on each solution;\n 5) a 0/1-flag to indicate whether additional diagnostic output needs\n to be written to screen.\n The Pade approximants are stored in the dobldobl systems pool,\n numerators in the odd indexed entries and denominators in the entries\n with even index in each system.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_quaddobl_Pade_approximant",
     py2c_quaddobl_Pade_approximant, METH_VARARGS,
    "Given in the systems container a polynomial system with coefficients\n in quad double precision, and in the solutions container the\n leading coefficients of the power series, this function runs Newton's\n method to compute power series solutions of the system in the container,\n in quad double precision, followed by the construction of the\n Pade approximants, for each solution. There are five integers on input:\n 1) the index of the series parameter;\n 2) the degree of the numerator of the Pade approximant;\n 3) the degree of the denominator of the Pade approximant;\n 4) the number of Newton steps to be done on each solution;\n 5) a 0/1-flag to indicate whether additional diagnostic output needs\n to be written to screen.\n The Pade approximants are stored in the quaddobl systems pool,\n numerators in the odd indexed entries and denominators in the entries\n with even index in each system.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_padcon_set_default_parameters",
     py2c_padcon_set_default_parameters, METH_VARARGS,
    "Sets the default values of the homotopy continuation parameters."},
   {"py2c_padcon_clear_parameters",
     py2c_padcon_clear_parameters, METH_VARARGS,
    "Deallocates the allocated space for the parameters."},
   {"py2c_padcon_get_homotopy_continuation_parameter",
     py2c_padcon_get_homotopy_continuation_parameter, METH_VARARGS,
    "Returns the value of the k-th continuation parameter,\n if k ranges between 1 and 13.  The integer k is given on entry."},
   {"py2c_padcon_set_homotopy_continuation_gamma",
     py2c_padcon_set_homotopy_continuation_gamma, METH_VARARGS,
    "The gamma constant is the first homotopy continuation parameter.\n The gamma is a complex number and it should be given as two\n doubles, as its real and imaginary part respectively."},
   {"py2c_padcon_set_homotopy_continuation_parameter",
     py2c_padcon_set_homotopy_continuation_parameter, METH_VARARGS,
    "Sets the value of the k-th continuation parameter to the given value.\n The first parameter k is an integer number between 2 and 13.\n The second parameter is the value of the k-th parameter,\n parsed as a floating point number."},
   {"py2c_padcon_reset_homotopy_continuation_parameters",
     py2c_padcon_reset_homotopy_continuation_parameters, METH_VARARGS,
    "Resets the value of the homotopy continuation parameters\n for the step-by-step path trackers.\n The first parameter is an integer number, 0, 1, or 2,\n respectively for double, double double, or quad double precision."},
   {"py2c_padcon_standard_track",
     py2c_padcon_standard_track, METH_VARARGS,
    "For the defined target, start system, and start solutions,\n launches the Pade continuation in standard double precision.\n Seven input parameters are expected:\n 1) the number of characters in the name of the output file;\n 2) a string which defines the name of the output file,\n if the string is empty, then no file is created;\n 3) a flag to indicate whether the output file is the defined output file\n (value 1 of the flag), or whether the file is local (value 0);\n 4) an integer for the verbose flag, if zero, then no extra\n information is written to file or screen;\n 5) an integer for the homogenization, if zero, tracking happens in\n affine space, if one, then tracking happens in 1-projective space,\n if m, for m > 1, then multihomogenization is applied;\n 6) an integer for the number of variables, 0 if the fifth parameter m\n is zero or one;\n 7) a string with the index representation for the partition of the\n set of variables, if the fifth parameter m is larger than one."},
   {"py2c_padcon_dobldobl_track",
     py2c_padcon_dobldobl_track, METH_VARARGS, 
    "For the defined target, start system, and start solutions,\n launches the Pade continuation in double double precision.\n Seven input parameters are expected:\n 1) the number of characters in the name of the output file;\n 2) a string which defines the name of the output file,\n if the string is empty, then no file is created;\n 3) a flag to indicate whether the output file is the defined output file\n (value 1 of the flag), or whether the file is local (value 0);\n 4) an integer for the verbose flag, if zero, then no extra\n information is written to file or screen;\n 5) an integer for the homogenization, if zero, tracking happens in\n affine space, if one, then tracking happens in 1-projective space,\n if m, for m > 1, then multihomogenization is applied;\n 6) an integer for the number of variables, 0 if the fifth parameter m\n is zero or one;\n 7) a string with the index representation for the partition of the\n set of variables, if the fifth parameter m is larger than one."},
   {"py2c_padcon_quaddobl_track",
     py2c_padcon_quaddobl_track, METH_VARARGS,
    "For the defined target, start system, and start solutions,\n launches the Pade continuation in quad double precision.\n Seven input parameters are expected:\n 1) the number of characters in the name of the output file;\n 2) a string which defines the name of the output file,\n if the string is empty, then no file is created;\n 3) a flag to indicate whether the output file is the defined output file\n (value 1 of the flag), or whether the file is local (value 0);\n 4) an integer for the verbose flag, if zero, then no extra\n information is written to file or screen;\n 5) an integer for the homogenization, if zero, tracking happens in\n affine space, if one, then tracking happens in 1-projective space,\n if m, for m > 1, then multihomogenization is applied;\n 6) an integer for the number of variables, 0 if the fifth parameter m\n is zero or one;\n 7) a string with the index representation for the partition of the\n set of variables, if the fifth parameter m is larger than one."},
   {"py2c_padcon_standard_initialize_homotopy",
     py2c_padcon_standard_initialize_homotopy, METH_VARARGS,
    "For the defined target and start system,\n initializes the homotopy in standard double precision,\n for the step-by-step Pade continuation.\n On entry is one parameter, the verbose flag which is zero or one.\n If the verbose flag is 1, then extra output will be written."},
   {"py2c_padcon_dobldobl_initialize_homotopy",
     py2c_padcon_dobldobl_initialize_homotopy, METH_VARARGS,
    "For the defined target and start system,\n initializes the homotopy in double double precision,\n for the step-by-step Pade continuation.\n On entry is one parameter, the verbose flag which is zero or one.\n If the verbose flag is 1, then extra output will be written."},
   {"py2c_padcon_quaddobl_initialize_homotopy",
     py2c_padcon_quaddobl_initialize_homotopy, METH_VARARGS,
    "For the defined target and start system,\n initializes the homotopy in quad double precision,\n for the step-by-step Pade continuation.\n On entry is one parameter, the verbose flag which is zero or one.\n If the verbose flag is 1, then extra output will be written."},
   {"py2c_padcon_standard_initialize_parameter_homotopy",
     py2c_padcon_standard_initialize_parameter_homotopy, METH_VARARGS,
    "On entry are two integers: 1) the index for the continuation\n parameter in the natural homotopy and 2) the verbose flag.\n With the system, defined as target system, and the index\n for the continuation parameter, initializes the homotopy in\n standard double precision for the step-by-step Pade continuation.\n If the verbose flag is 1, then extra output will be written."},
   {"py2c_padcon_dobldobl_initialize_parameter_homotopy",
     py2c_padcon_dobldobl_initialize_parameter_homotopy, METH_VARARGS,
    "On entry are two integers: 1) the index for the continuation\n parameter in the natural homotopy and 2) the verbose flag.\n With the system, defined as target system, and the index\n for the continuation parameter, initializes the homotopy in\n double double precision for the step-by-step Pade continuation.\n If the verbose flag is 1, then extra output will be written."},
   {"py2c_padcon_quaddobl_initialize_parameter_homotopy",
     py2c_padcon_quaddobl_initialize_parameter_homotopy, METH_VARARGS,
    "On entry are two integers: 1) the index for the continuation\n parameter in the natural homotopy and 2) the verbose flag.\n With the system, defined as target system, and the index\n for the continuation parameter, initializes the homotopy in\n quad double precision for the step-by-step Pade continuation.\n If the verbose flag is 1, then extra output will be written."},
   {"py2c_padcon_initialize_standard_solution",
     py2c_padcon_initialize_standard_solution, METH_VARARGS,
    "Takes the solution with a given index in the solutions container in\n standard double precision and initializes the series-Pade tracker.\n On entry are two integers: 1) the index of the position of the solution\n in the container and 2) the verbose flag, which is zero or one.\n If the verbose flag is 1, then extra output will be written."},
   {"py2c_padcon_initialize_dobldobl_solution",
     py2c_padcon_initialize_dobldobl_solution, METH_VARARGS,
    "Takes the solution with a given index in the solutions container in\n double double precision and initializes the series-Pade tracker.\n On entry are two integers: 1) the index of the position of the solution\n in the container and 2) the verbose flag, which is zero or one.\n If the verbose flag is 1, then extra output will be written."},
   {"py2c_padcon_initialize_quaddobl_solution",
     py2c_padcon_initialize_quaddobl_solution, METH_VARARGS,
    "Takes the solution with a given index in the solutions container in\n quad double precision and initializes the series-Pade tracker.\n On entry are two integers: 1) the index of the position of the solution\n in the container and 2) the verbose flag, which is zero or one.\n If the verbose flag is 1, then extra output will be written."},
   {"py2c_padcon_standard_predict_correct",
     py2c_padcon_standard_predict_correct, METH_VARARGS,
    "Executes one predict-correct step on the current solution and\n the defined homotopy in standard double precision.\n On entry is one integer, the verbose flag which is zero or one.\n On return is the failure code of the predict-correct step:\n if zero, then the required accuracies were met,\n otherwise, either the predict or the correct step failed.\n If the verbose flag is 1, then extra output will be written."},
   {"py2c_padcon_dobldobl_predict_correct",
     py2c_padcon_dobldobl_predict_correct, METH_VARARGS,
    "Executes one predict-correct step on the current solution and\n the defined homotopy in double double precision.\n On entry is one integer, the verbose flag which is zero or one.\n On return is the failure code of the predict-correct step:\n if zero, then the required accuracies were met,\n otherwise, either the predict or the correct step failed.\n If the verbose flag is 1, then extra output will be written."},
   {"py2c_padcon_quaddobl_predict_correct",
     py2c_padcon_quaddobl_predict_correct, METH_VARARGS,
    "Executes one predict-correct step on the current solution and\n the defined homotopy in quad double precision.\n On entry is one integer, the verbose flag which is zero or one.\n On return is the failure code of the predict-correct step:\n if zero, then the required accuracies were met,\n otherwise, either the predict or the correct step failed.\n If the verbose flag is 1, then extra output will be written."},
   {"py2c_padcon_get_standard_solution",
     py2c_padcon_get_standard_solution, METH_VARARGS,
    "On entry are two integer parameters: 1) the index of the position of\n the solution and 2) the verbose flag, which is zero or one.\n Retrieves the current solution and places it at the given position\n in the solutions container in standard double precision.\n If the verbose flag is 1, then extra output will be written."},
   {"py2c_padcon_get_dobldobl_solution",
     py2c_padcon_get_dobldobl_solution, METH_VARARGS,
    "On entry are two integer parameters: 1) the index of the position of\n the solution and 2) the verbose flag, which is zero or one.\n Retrieves the current solution and places it at the given position\n in the solutions container in double double precision.\n If the verbose flag is 1, then extra output will be written."},
   {"py2c_padcon_get_quaddobl_solution",
     py2c_padcon_get_quaddobl_solution, METH_VARARGS,
    "On entry are two integer parameters: 1) the index of the position of\n the solution and 2) the verbose flag, which is zero or one.\n Retrieves the current solution and places it at the given position\n in the solutions container in quad double precision.\n If the verbose flag is 1, then extra output will be written."},
   {"py2c_padcon_standard_pole_radius",
     py2c_padcon_standard_pole_radius, METH_VARARGS,
    "Returns the smallest pole radius computed\n by the predictor in standard double precision."},
   {"py2c_padcon_dobldobl_pole_radius",
     py2c_padcon_dobldobl_pole_radius, METH_VARARGS,
    "Returns the smallest pole radius computed\n by the predictor in double double precision.\n The returned number is the high part of the double double number."},
   {"py2c_padcon_quaddobl_pole_radius",
     py2c_padcon_quaddobl_pole_radius, METH_VARARGS,
    "Returns the smallest pole radius computed\n by the predictor in quad double precision.\n The returned number is the highest part of the quad double number."},
   {"py2c_padcon_standard_closest_pole",
     py2c_padcon_standard_closest_pole, METH_VARARGS,
    "Returns the complex number representation of the closest pole,\n computed by the predictor in standard double precision.\n Results are meaningful only if the real part >= 0.0."},
   {"py2c_padcon_dobldobl_closest_pole",
     py2c_padcon_dobldobl_closest_pole, METH_VARARGS,
    "Returns the complex number representation of the closest pole,\n computed by the predictor in double double precision.\n The returned numbers are the high parts of the double doubles.\n Results are meaningful only if the real part >= 0.0."},
   {"py2c_padcon_quaddobl_closest_pole",
     py2c_padcon_quaddobl_closest_pole, METH_VARARGS,
    "Returns the complex number representation of the closest pole,\n computed by the predictor in quad double precision.\n The returned numbers are the highest parts of the quad doubles.\n Results are meaningful only if the real part >= 0.0."},
   {"py2c_padcon_standard_t_value",
     py2c_padcon_standard_t_value, METH_VARARGS,
    "Returns the current t value of the path tracker\n which runs in standard double precision."},
   {"py2c_padcon_dobldobl_t_value",
     py2c_padcon_dobldobl_t_value, METH_VARARGS,
    "Returns the current t value of the path tracker\n which runs in double double precision."},
   {"py2c_padcon_quaddobl_t_value",
     py2c_padcon_quaddobl_t_value, METH_VARARGS,
    "Returns the current t value of the path tracker\n which runs in quad double precision."},
   {"py2c_padcon_standard_step_size",
     py2c_padcon_standard_step_size, METH_VARARGS,
    "Returns the current step size of the path tracker\n which runs in standard double precision."},
   {"py2c_padcon_dobldobl_step_size",
     py2c_padcon_dobldobl_step_size, METH_VARARGS,
    "Returns the current step size of the path tracker\n which runs in double double precision."},
   {"py2c_padcon_quaddobl_step_size",
     py2c_padcon_quaddobl_step_size, METH_VARARGS,
    "Returns the current step size of the path tracker\n which runs in quad double precision."},
   {"py2c_padcon_standard_series_step",
     py2c_padcon_standard_series_step, METH_VARARGS,
    "Returns the current series step size of the path tracker which runs in standard double precision."},
   {"py2c_padcon_dobldobl_series_step",
     py2c_padcon_dobldobl_series_step, METH_VARARGS,
    "Returns the current series step size of the path tracker which runs in double double precision."},
   {"py2c_padcon_quaddobl_series_step",
     py2c_padcon_quaddobl_series_step, METH_VARARGS,
    "Returns the current series step size of the path tracker which runs in quad double precision."},
   {"py2c_padcon_standard_pole_step",
     py2c_padcon_standard_pole_step, METH_VARARGS,
    "Returns the current pole step size of the path tracker which runs in standard double precision."},
   {"py2c_padcon_dobldobl_pole_step",
     py2c_padcon_dobldobl_pole_step, METH_VARARGS,
    "Returns the current pole step size of the path tracker which runs in double double precision."},
   {"py2c_padcon_quaddobl_pole_step",
     py2c_padcon_quaddobl_pole_step, METH_VARARGS,
    "Returns the current pole step size of the path tracker which runs in quad double precision."},
   {"py2c_padcon_standard_estimated_distance",
     py2c_padcon_standard_estimated_distance, METH_VARARGS,
    "Returns the estimated distance to the closest solution by the path tracker which runs in standard double precision."},
   {"py2c_padcon_dobldobl_estimated_distance",
     py2c_padcon_dobldobl_estimated_distance, METH_VARARGS,
    "Returns the estimated distance to the closest solution by the path tracker which runs in double double precision."},
   {"py2c_padcon_quaddobl_estimated_distance",
     py2c_padcon_quaddobl_estimated_distance, METH_VARARGS,
    "Returns the estimated distance to the closest solution by the path tracker which runs in quad double precision."},
   {"py2c_padcon_standard_hessian_step",
     py2c_padcon_standard_hessian_step, METH_VARARGS,
    "Returns the current Hessian step size of the path tracker which runs in standard double precision."},
   {"py2c_padcon_dobldobl_hessian_step",
     py2c_padcon_dobldobl_hessian_step, METH_VARARGS,
    "Returns the current Hessian step size of the path tracker which runs in double double precision."},
   {"py2c_padcon_quaddobl_hessian_step",
     py2c_padcon_quaddobl_hessian_step, METH_VARARGS,
    "Returns the current Hessian step size of the path tracker which runs in quad double precision."},
   {"py2c_padcon_standard_series_coefficient",
     py2c_padcon_standard_series_coefficient, METH_VARARGS,
    "Returns a tuple: the real and imaginary parts of the series\n coefficient of component with leadidx at position idx,\n of the series computed by the predictor in double precision.\n The integers leadidx and idx are two input parameters,\n the third input integer is the verbose flag."},
   {"py2c_padcon_dobldobl_series_coefficient",
     py2c_padcon_dobldobl_series_coefficient, METH_VARARGS,
    "Returns a tuple: the real and imaginary parts of the series\n coefficient of component with leadidx at position idx, of the\n series computed by the predictor in double double precision.\n The doubles are the highest parts of the double doubles.\n The integers leadidx and idx are two input parameters,\n the third input integer is the verbose flag."},
   {"py2c_padcon_quaddobl_series_coefficient",
     py2c_padcon_quaddobl_series_coefficient, METH_VARARGS,
    "Returns a tuple: the real and imaginary parts of the series\n coefficient of component with leadidx at position idx, of the\n series computed by the predictor in quad double precision.\n The doubles are the highest parts of the quad doubles.\n The integers leadidx and idx are two input parameters,\n the third input integer is the verbose flag."},
   {"py2c_padcon_standard_numerator_coefficient",
     py2c_padcon_standard_numerator_coefficient, METH_VARARGS,
    "Returns a tuple: the real and imaginary parts of the\n coefficient of the numerator of the Pade approximant,\n at the component with leadidx at position idx,\n computed by the predictor in double precision.\n The integers leadidx and idx are two input parameters,\n the third input integer is the verbose flag."},
   {"py2c_padcon_dobldobl_numerator_coefficient",
     py2c_padcon_dobldobl_numerator_coefficient, METH_VARARGS,
    "Returns a tuple: the real and imaginary parts of the\n coefficient of the numerator of the Pade approximant,\n at the component with leadidx at position idx,\n computed by the predictor in double double precision.\n The doubles are the highest parts of the double doubles.\n The integers leadidx and idx are two input parameters,\n the third input integer is the verbose flag."},
   {"py2c_padcon_quaddobl_numerator_coefficient",
     py2c_padcon_quaddobl_numerator_coefficient, METH_VARARGS,
    "Returns a tuple: the real and imaginary parts of the series\n coefficient of the numerator of the Pade approximant,\n at the component with leadidx at position idx,\n computed by the predictor in quad double precision.\n The doubles are the highest parts of the quad doubles.\n The integers leadidx and idx are two input parameters,\n the third input integer is the verbose flag."},
   {"py2c_padcon_standard_denominator_coefficient",
     py2c_padcon_standard_denominator_coefficient, METH_VARARGS,
    "Returns a tuple: the real and imaginary parts of the\n coefficient of the denominator of the Pade approximant,\n at the component with leadidx at position idx,\n computed by the predictor in double precision.\n The integers leadidx and idx are two input parameters,\n the third input integer is the verbose flag."},
   {"py2c_padcon_dobldobl_denominator_coefficient",
     py2c_padcon_dobldobl_denominator_coefficient, METH_VARARGS,
    "Returns a tuple: the real and imaginary parts of the\n coefficient of the denominator of the Pade approximant,\n at the component with leadidx at position idx,\n computed by the predictor in double double precision.\n The doubles are the highest parts of the double doubles.\n The integers leadidx and idx are two input parameters,\n the third input integer is the verbose flag."},
   {"py2c_padcon_quaddobl_denominator_coefficient",
     py2c_padcon_quaddobl_denominator_coefficient, METH_VARARGS,
    "Returns a tuple: the real and imaginary parts of the series\n coefficient of the denominator of the Pade approximant,\n at the component with leadidx at position idx,\n computed by the predictor in quad double precision.\n The doubles are the highest parts of the quad doubles.\n The integers leadidx and idx are two input parameters,\n the third input integer is the verbose flag."},
   {"py2c_padcon_standard_pole", py2c_padcon_standard_pole, METH_VARARGS,
    "Returns a tuple: the real and imaginary parts of the pole\n Pade approximant with leadidx at position poleidx,\n computed by the predictor in double precision.\n The integers leadidx and poleidx are two input parameters,\n the third input integer is the verbose flag."},
   {"py2c_padcon_dobldobl_pole", py2c_padcon_dobldobl_pole, METH_VARARGS,
    "Returns a tuple: the real and imaginary parts of the pole\n Pade approximant with leadidx at position poleidx,\n computed by the predictor in double double precision.\n The integers leadidx and poleidx are two input parameters,\n the third input integer is the verbose flag.\n The returned doubles are the highest parts of the double doubles."},
   {"py2c_padcon_quaddobl_pole", py2c_padcon_quaddobl_pole, METH_VARARGS,
    "Returns a tuple: the real and imaginary parts of the pole\n Pade approximant with leadidx at position poleidx,\n computed by the predictor in quad double precision.\n The integers leadidx and poleidx are two input parameters,\n the third input integer is the verbose flag.\n The returned doubles are the highest parts of the quad doubles."},
   {"py2c_padcon_clear_standard_data",
     py2c_padcon_clear_standard_data, METH_VARARGS,
    "Deallocates data for the series-Pade tracker in double precision."},
   {"py2c_padcon_clear_dobldobl_data",
     py2c_padcon_clear_dobldobl_data, METH_VARARGS,
    "Deallocates data for the series-Pade tracker in double double precision."},
   {"py2c_padcon_clear_quaddobl_data",
     py2c_padcon_clear_quaddobl_data, METH_VARARGS,
    "Deallocates data for the series-Pade tracker in quad double precision."},
   {"py2c_syspool_standard_init", py2c_syspool_standard_init, METH_VARARGS,
    "Initializes the pool for systems in standard double precision."},
   {"py2c_syspool_dobldobl_init", py2c_syspool_dobldobl_init, METH_VARARGS,
    "Initializes the pool for systems in double double precision."},
   {"py2c_syspool_quaddobl_init", py2c_syspool_quaddobl_init, METH_VARARGS,
    "Initializes the pool for systems in quad double precision."},
   {"py2c_syspool_standard_size", py2c_syspool_standard_size, METH_VARARGS,
    "Returns the size of the pool for systems in standard double precision."},
   {"py2c_syspool_dobldobl_size", py2c_syspool_dobldobl_size, METH_VARARGS,
    "Returns the size of the pool for systems in double double precision."},
   {"py2c_syspool_quaddobl_size", py2c_syspool_quaddobl_size, METH_VARARGS,
    "Returns the size of the pool for systems in quad double precision."},
   {"py2c_syspool_standard_create",
     py2c_syspool_standard_create, METH_VARARGS,
    "Defines the k-th system in the standard system pool,\n using the system in the standard container."},
   {"py2c_syspool_dobldobl_create",
     py2c_syspool_dobldobl_create, METH_VARARGS,
    "Defines the k-th system in the dobldobl system pool,\n using the system in the dobldobl container."},
   {"py2c_syspool_quaddobl_create",
     py2c_syspool_quaddobl_create, METH_VARARGS,
    "Defines the k-th system in the quaddobl system pool,\n using the system in the quaddobl container."},
   {"py2c_syspool_copy_to_standard_container",
     py2c_syspool_copy_to_standard_container, METH_VARARGS,
    "Copies the k-th system in the pool for systems in standard double\n precision to the standard systems container.\n The value for k is given as an integer input parameter.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syspool_copy_to_dobldobl_container",
     py2c_syspool_copy_to_dobldobl_container, METH_VARARGS,
    "Copies the k-th system in the pool for systems in double double\n precision to the dobldobl systems container.\n The value for k is given as an integer input parameter.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syspool_copy_to_quaddobl_container",
     py2c_syspool_copy_to_quaddobl_container, METH_VARARGS,
    "Copies the k-th system in the pool for systems in quad double\n precision to the quaddobl systems container.\n The value for k is given as an integer input parameter.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_syspool_standard_clear", py2c_syspool_standard_clear, METH_VARARGS,
    "Clears the pool for systems in standard double precision."},
   {"py2c_syspool_dobldobl_clear", py2c_syspool_dobldobl_clear, METH_VARARGS,
    "Clears the pool for systems in double double precision."},
   {"py2c_syspool_quaddobl_clear", py2c_syspool_quaddobl_clear, METH_VARARGS,
    "Clears the pool for systems in quad double precision."},
   {"py2c_initialize_standard_homotopy", py2c_initialize_standard_homotopy,
     METH_VARARGS,
    "Initializes the homotopy to track a path with a generator,\n using standard double precision arithmetic.\n There is one integer number on input to be considered as a boolean,\n as an indicator whether a fixed gamma constant will be used.\n Before calling this routine the target and start system must\n be copied over from the standard systems container.\n The two other input parameters are two doubles: the real and imaginary part\n of the gamma constant.  If the integer parameter equals zero and if the two\n input doubles are not both zero, then the input gamma constant will be used,\n otherwise, if the two input doubles are zero and the first integer parameter\n is zero as well, then a random gamma constant will be generated."},
   {"py2c_initialize_dobldobl_homotopy", py2c_initialize_dobldobl_homotopy,
     METH_VARARGS,
    "Initializes the homotopy to track a path with a generator,\n using double double precision arithmetic.\n There is one integer number on input to be considered as a boolean,\n as an indicator whether a fixed gamma constant will be used.\n Before calling this routine the target and start system must\n be copied over from the dobldobl systems container.\n The two other input parameters are two doubles: the real and imaginary part\n of the gamma constant.  If the integer parameter equals zero and if the two\n input doubles are not both zero, then the input gamma constant will be used,\n otherwise, if the two input doubles are zero and the first integer parameter\n is zero as well, then a random gamma constant will be generated."},
   {"py2c_initialize_quaddobl_homotopy", py2c_initialize_quaddobl_homotopy,
     METH_VARARGS,
    "Initializes the homotopy to track a path with a generator,\n using quad double precision arithmetic.\n There is one integer number on input to be considered as a boolean,\n as an indicator whether a fixed gamma constant will be used.\n Before calling this routine the target and start system must\n be copied over from the quaddobl systems container.\n The two other input parameters are two doubles: the real and imaginary part\n of the gamma constant.  If the integer parameter equals zero and if the two\n input doubles are not both zero, then the input gamma constant will be used,\n otherwise, if the two input doubles are zero and the first integer parameter\n is zero as well, then a random gamma constant will be generated."},
   {"py2c_initialize_multprec_homotopy", py2c_initialize_multprec_homotopy,
     METH_VARARGS,
    "Initializes the homotopy to track a path with a generator,\n using arbitrary multiprecision arithmetic.\n There is are two integer numbers on input:\n 1) one to be considered as a boolean,\n as an indicator whether a fixed gamma constant will be used; and\n 2) the number of decimal places in the working precision.\n Before calling this routine the target and start system must\n be copied over from the multprec systems container."},
   {"py2c_initialize_varbprec_homotopy", py2c_initialize_varbprec_homotopy,
     METH_VARARGS,
    "Initializes the variable precision homotopy with the target and\n start system stored in the strings.\n On entry are three integers and two strings, in the following order:\n 1) fixed_gamma is a flag: if 1, then a fixed value for the gamma constant\n is used, if 0, a random value for gamma will be generated;\n 2) nc_target, the number of characters in the string target;\n 3) target, the string representation of the target system;\n 4) nc_start, the number of characters in the string start;\n 5) start, the string representation of the start system.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_initialize_standard_solution", py2c_initialize_standard_solution, 
     METH_VARARGS,
    "Initializes the path tracker with a generator with a solution\n from the standard solutions container.  The index to the solution\n is given as an integer input parameter.  The counting of the\n indices starts at one, so the first solution has index one."},
   {"py2c_initialize_dobldobl_solution", py2c_initialize_dobldobl_solution,
     METH_VARARGS,
    "Initializes the path tracker with a generator with a solution\n from the dobldobl solutions container.  The index to the solution\n is given as an integer input parameter.  The counting of the\n indices starts at one, so the first solution has index one."},
   {"py2c_initialize_quaddobl_solution", py2c_initialize_quaddobl_solution,
     METH_VARARGS,
    "Initializes the path tracker with a generator with a solution\n from the quaddobl solutions container.  The index to the solution\n is given as an integer input parameter.  The counting of the\n indices starts at one, so the first solution has index one."},
   {"py2c_initialize_multprec_solution", py2c_initialize_multprec_solution,
     METH_VARARGS, 
    "Initializes the path tracker with a generator with a solution\n from the multprec solutions container.  The index to the solution\n is given as an integer input parameter.  The counting of the\n indices starts at one, so the first solution has index one."},
   {"py2c_initialize_varbprec_solution", py2c_initialize_varbprec_solution,
     METH_VARARGS,
    "Uses the string representation of a solution to initialize the\n variable precision path tracker with.\n There are three input parameters, two integers and one string:\n 1) nv, the number of variables in the solution;\n 2) nc, the number of characters in the string sol;\n 3) sol, the string representation of a solution.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_next_standard_solution", py2c_next_standard_solution, 
     METH_VARARGS,
    "Computes the next point on the solution path with standard double\n precision for the given index.  This index is given as an input\n parameter.  The index to the solution path starts its count at one.\n The point itself is stored in the standard solutions container.\n The functions py2c_initialized_standard_tracker and\n py2c_initialize_standard_solution must have been executed earlier.\n The failcode is returned, which equals zero if all is well."},
   {"py2c_next_dobldobl_solution", py2c_next_dobldobl_solution,
     METH_VARARGS,
    "Computes the next point on the solution path with double double\n precision for the given index.  This index is given as an input\n parameter.  The index to the solution path starts its count at one.\n The point itself is stored in the dobldobl solutions container.\n The functions py2c_initialized_dobldobl_tracker and\n py2c_initialize_dobldobl_solution must have been executed earlier.\n The failcode is returned, which equals zero if all is well."},
   {"py2c_next_quaddobl_solution", py2c_next_quaddobl_solution,
     METH_VARARGS,
    "Computes the next point on the solution path with quad double\n precision for the given index.  This index is given as an input\n parameter.  The index to the solution path starts its count at one.\n The point itself is stored in the quaddobl solutions container.\n The functions py2c_initialized_quaddobl_tracker and\n py2c_initialize_quaddobl_solution must have been executed earlier.\n The failcode is returned, which equals zero if all is well."},
   {"py2c_next_multprec_solution", py2c_next_multprec_solution,
     METH_VARARGS,
    "Computes the next point on the solution path with arbitrary\n multiprecision for the given index.  This index is given as an input\n parameter.  The index to the solution path starts its count at one.\n The point itself is stored in the multprec solutions container.\n The functions py2c_initialized_multprec_tracker and\n py2c_initialize_multprec_solution must have been executed earlier.\n The failcode is returned, which equals zero if all is well."},
   {"py2c_next_varbprec_solution", py2c_next_varbprec_solution,
     METH_VARARGS,
    "Computes the next point on a solution path in variable precision.\n There are four integer input parameters:\n 1) the number of correct decimal places in the solution;\n 2) an upper bound on the number of decimal places in the precision;\n 3) the maximum number of Newton iterations;\n 4) a flag zero or one to indicate the verbose level.\n On return is a tuple:\n 0) the failure code, which equals zero if all went well; and\n 1) the string representation of the next solution on the path."},
   {"py2c_clear_standard_tracker", py2c_clear_standard_tracker, METH_VARARGS, 
    "Deallocates data used in the standard double precision tracker\n with a generator."},
   {"py2c_clear_dobldobl_tracker", py2c_clear_dobldobl_tracker, METH_VARARGS, 
    "Deallocates data used in the double double precision tracker\n with a generator."},
   {"py2c_clear_quaddobl_tracker", py2c_clear_quaddobl_tracker, METH_VARARGS,
    "Deallocates data used in the quad double precision tracker\n with a generator."},
   {"py2c_clear_multprec_tracker", py2c_clear_multprec_tracker, METH_VARARGS,
    "Deallocates data used in the arbitrary multiprecision tracker\n with a generator."},
   {"py2c_clear_varbprec_tracker", py2c_clear_varbprec_tracker, METH_VARARGS, 
    "Deallocates data used in the variable precision tracker\n with a generator."},
   {"py2c_ade_newton_d", py2c_ade_newton_d, METH_VARARGS,
    "Runs Newton's method with algorithmic differentation\n in double precision on the data in the systems and solutions container.\n The standard systems container must contain a valid polynomial system\n and the standard solutions container must hold a valid solution.\n On entry is the verbose flag, which equals zero if no output is wanted,\n or 1 if extra information should be written to screen.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_ade_onepath_d", py2c_ade_onepath_d, METH_VARARGS,
    "Tracks one solution path with algorithmic differentation\n in double precision on the data in the systems and solutions container.\n The start and target systems must have been defined\n and the standard solutions container must holds valid solution.\n On entry is the verbose flag, which equals zero if no output is wanted,\n or 1 if extra information should be written to screen.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_ade_manypaths_d", py2c_ade_manypaths_d, METH_VARARGS,
    "Tracks many solution paths with algorithmic differentation\n in double precision on the data in the systems and solutions container.\n The start and target systems must have been defined\n and the standard solutions container holds valid solutions.\n On entry is the verbose flag, which equals zero if no output is wanted,\n or 1 if extra information should be written to screen.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_ade_newton_dd", py2c_ade_newton_dd, METH_VARARGS,
    "Runs Newton's method with algorithmic differentation\n in double double precision on the data in the systems and solutions container.\n The dobldobl systems container must contain a valid polynomial system\n and the dobldobl solutions container must hold a valid solution.\n On entry is the verbose flag, which equals zero if no output is wanted,\n or 1 if extra information should be written to screen.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_ade_onepath_dd", py2c_ade_onepath_dd, METH_VARARGS,
    "Tracks one solution path with algorithmic differentation\n in double double precision on the data in the systems and solutions container.\n The start and target systems must have been defined\n and the dobldobl solutions container must holds valid solution.\n On entry is the verbose flag, which equals zero if no output is wanted,\n or 1 if extra information should be written to screen.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_ade_manypaths_dd", py2c_ade_manypaths_dd, METH_VARARGS,
    "Tracks many solution paths with algorithmic differentation\n in double precision on the data in the systems and solutions container.\n The start and target systems must have been defined\n and the dobldobl solutions container holds valid solutions.\n On entry is the verbose flag, which equals zero if no output is wanted,\n or 1 if extra information should be written to screen.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_ade_newton_qd", py2c_ade_newton_qd, METH_VARARGS,
    "Runs Newton's method with algorithmic differentation\n in quad double precision on the data in the systems and solutions container.\n The quaddobl systems container must contain a valid polynomial system\n and the quaddobl solutions container must hold a valid solution.\n On entry is the verbose flag, which equals zero if no output is wanted,\n or 1 if extra information should be written to screen.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_ade_onepath_qd", py2c_ade_onepath_qd, METH_VARARGS,
    "Tracks one solution path with algorithmic differentation\n in quad double precision on the data in the systems and solutions container.\n The start and target systems must have been defined\n and the quaddobl solutions container must holds valid solution.\n On entry is the verbose flag, which equals zero if no output is wanted,\n or 1 if extra information should be written to screen.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_ade_manypaths_qd", py2c_ade_manypaths_qd, METH_VARARGS,
    "Tracks many solution paths with algorithmic differentation\n in quad double precision on the data in the systems and solutions container.\n The start and target systems must have been defined\n and the quaddobl solutions container holds valid solutions.\n On entry is the verbose flag, which equals zero if no output is wanted,\n or 1 if extra information should be written to screen.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_get_default_path_parameters", 
     py2c_get_default_path_parameters, METH_VARARGS, 
    "Given the working precision (16, 32, or 64), returns the default values\n of the path parameters, for the path trackers with algorithmic differentiation."},
   {"py2c_ade_manypaths_d_pars",
     py2c_ade_manypaths_d_pars, METH_VARARGS,
    "Tracks many solution paths with algorithmic differentation\n in double precision on the data in the systems and solutions container.\n The start and target systems must have been defined\n and the standard solutions container holds valid solutions.\n On entry is the verbose flag, which equals zero if no output is wanted,\n or 1 if extra information should be written to screen.\n Other input parameters are the real and imaginary parts of the gamma constant.\n Then, the 14 values of the path parameters has to be provided.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_ade_manypaths_dd_pars",
     py2c_ade_manypaths_dd_pars, METH_VARARGS,
    "Tracks many solution paths with algorithmic differentation\n in double double precision on the data in the systems and solutions container.\n The start and target systems must have been defined\n and the dobldobl solutions container holds valid solutions.\n On entry is the verbose flag, which equals zero if no output is wanted,\n or 1 if extra information should be written to screen.\n Other input parameters are the real and imaginary parts of the gamma constant.\n Then, the 14 values of the path parameters has to be provided.\n On return is the failure code, which equals zero if all went well."},
   {"py2c_ade_manypaths_qd_pars",
     py2c_ade_manypaths_qd_pars, METH_VARARGS,
    "Tracks many solution paths with algorithmic differentation\n in quad double precision on the data in the systems and solutions container.\n The start and target systems must have been defined\n and the quaddobl solutions container holds valid solutions.\n On entry is the verbose flag, which equals zero if no output is wanted,\n or 1 if extra information should be written to screen.\n Other input parameters are the real and imaginary parts of the gamma constant.\n Then, the 14 values of the path parameters has to be provided.\n On return is the failure code, which equals zero if all went well."},
   {NULL, NULL, 0, NULL} 
};

/* This is the initialization routine which will be called by the 
 * Python run-time when the library is imported in order to retrieve 
 * a pointer to the above method address table.
 * Note that therefore this routine must be visible in the dynamic library
 * either through the use of a ".def" file or by a compiler instruction 
 * such as "declspec(export)" */

PyMODINIT_FUNC initphcpy2c2(void)
{
   Py_InitModule("phcpy2c2", phcpy2c_methods);
}
