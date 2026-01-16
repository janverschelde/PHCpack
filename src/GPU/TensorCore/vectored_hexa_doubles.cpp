/* Defines the functions with prototypes in vectored_hexa_doubles.h. */

#include "hexa_double.h"
#include "hexa_double_functions.h"
#include "splitting_doubles.h"

void quarter_hexa_double
 ( double xhihihihi, double xlohihihi, double xhilohihi, double xlolohihi,
   double xhihilohi, double xlohilohi, double xhilolohi, double xlololohi,
   double xhihihilo, double xlohihilo, double xhilohilo, double xlolohilo,
   double xhihilolo, double xlohilolo, double xhilololo, double xlolololo,
   double *xhihihihi0, double *xhihihihi1,
   double *xhihihihi2, double *xhihihihi3,
   double *xlohihihi0, double *xlohihihi1,
   double *xlohihihi2, double *xlohihihi3,
   double *xhilohihi0, double *xhilohihi1,
   double *xhilohihi2, double *xhilohihi3,
   double *xlolohihi0, double *xlolohihi1,
   double *xlolohihi2, double *xlolohihi3,
   double *xhihilohi0, double *xhihilohi1,
   double *xhihilohi2, double *xhihilohi3,
   double *xlohilohi0, double *xlohilohi1,
   double *xlohilohi2, double *xlohilohi3,
   double *xhilolohi0, double *xhilolohi1,
   double *xhilolohi2, double *xhilolohi3,
   double *xlololohi0, double *xlololohi1,
   double *xlololohi2, double *xlololohi3,
   double *xhihihilo0, double *xhihihilo1,
   double *xhihihilo2, double *xhihihilo3,
   double *xlohihilo0, double *xlohihilo1,
   double *xlohihilo2, double *xlohihilo3,
   double *xhilohilo0, double *xhilohilo1,
   double *xhilohilo2, double *xhilohilo3,
   double *xlolohilo0, double *xlolohilo1,
   double *xlolohilo2, double *xlolohilo3,
   double *xhihilolo0, double *xhihilolo1,
   double *xhihilolo2, double *xhihilolo3,
   double *xlohilolo0, double *xlohilolo1,
   double *xlohilolo2, double *xlohilolo3,
   double *xhilololo0, double *xhilololo1,
   double *xhilololo2, double *xhilololo3,
   double *xlolololo0, double *xlolololo1,
   double *xlolololo2, double *xlolololo3 )
{
   quarter_split(xhihihihi, xhihihihi0, xhihihihi1, xhihihihi2, xhihihihi3);
   quarter_split(xlohihihi, xlohihihi0, xlohihihi1, xlohihihi2, xlohihihi3);
   quarter_split(xhilohihi, xhilohihi0, xhilohihi1, xhilohihi2, xhilohihi3);
   quarter_split(xlolohihi, xlolohihi0, xlolohihi1, xlolohihi2, xlolohihi3);
   quarter_split(xhihilohi, xhihilohi0, xhihilohi1, xhihilohi2, xhihilohi3);
   quarter_split(xlohilohi, xlohilohi0, xlohilohi1, xlohilohi2, xlohilohi3);
   quarter_split(xhilolohi, xhilolohi0, xhilolohi1, xhilolohi2, xhilolohi3);
   quarter_split(xlololohi, xlololohi0, xlololohi1, xlololohi2, xlololohi3);
   quarter_split(xhihihilo, xhihihilo0, xhihihilo1, xhihihilo2, xhihihilo3);
   quarter_split(xlohihilo, xlohihilo0, xlohihilo1, xlohihilo2, xlohihilo3);
   quarter_split(xhilohilo, xhilohilo0, xhilohilo1, xhilohilo2, xhilohilo3);
   quarter_split(xlolohilo, xlolohilo0, xlolohilo1, xlolohilo2, xlolohilo3);
   quarter_split(xhihilolo, xhihilolo0, xhihilolo1, xhihilolo2, xhihilolo3);
   quarter_split(xlohilolo, xlohilolo0, xlohilolo1, xlohilolo2, xlohilolo3);
   quarter_split(xhilololo, xhilololo0, xhilololo1, xhilololo2, xhilololo3);
   quarter_split(xlolololo, xlolololo0, xlolololo1, xlolololo2, xlolololo3);
}

void quarter_hd_vector
 ( int dim,
   double *xhihihihi, double *xlohihihi,
   double *xhilohihi, double *xlolohihi,
   double *xhihilohi, double *xlohilohi,
   double *xhilolohi, double *xlololohi,
   double *xhihihilo, double *xlohihilo,
   double *xhilohilo, double *xlolohilo,
   double *xhihilolo, double *xlohilolo,
   double *xhilololo, double *xlolololo,
   double *xhihihihi0, double *xhihihihi1,
   double *xhihihihi2, double *xhihihihi3,
   double *xlohihihi0, double *xlohihihi1,
   double *xlohihihi2, double *xlohihihi3,
   double *xhilohihi0, double *xhilohihi1,
   double *xhilohihi2, double *xhilohihi3,
   double *xlolohihi0, double *xlolohihi1,
   double *xlolohihi2, double *xlolohihi3,
   double *xhihilohi0, double *xhihilohi1,
   double *xhihilohi2, double *xhihilohi3,
   double *xlohilohi0, double *xlohilohi1,
   double *xlohilohi2, double *xlohilohi3,
   double *xhilolohi0, double *xhilolohi1,
   double *xhilolohi2, double *xhilolohi3,
   double *xlololohi0, double *xlololohi1,
   double *xlololohi2, double *xlololohi3,
   double *xhihihilo0, double *xhihihilo1,
   double *xhihihilo2, double *xhihihilo3,
   double *xlohihilo0, double *xlohihilo1,
   double *xlohihilo2, double *xlohihilo3,
   double *xhilohilo0, double *xhilohilo1,
   double *xhilohilo2, double *xhilohilo3,
   double *xlolohilo0, double *xlolohilo1,
   double *xlolohilo2, double *xlolohilo3,
   double *xhihilolo0, double *xhihilolo1,
   double *xhihilolo2, double *xhihilolo3,
   double *xlohilolo0, double *xlohilolo1,
   double *xlohilolo2, double *xlohilolo3,
   double *xhilololo0, double *xhilololo1,
   double *xhilololo2, double *xhilololo3,
   double *xlolololo0, double *xlolololo1,
   double *xlolololo2, double *xlolololo3 )
{
   for(int i=0; i<dim; i++)
   {
      quarter_hexa_double
         (xhihihihi[i], xlohihihi[i], xhilohihi[i], xlolohihi[i],
          xhihilohi[i], xlohilohi[i], xhilolohi[i], xlololohi[i],
          xhihihilo[i], xlohihilo[i], xhilohilo[i], xlolohilo[i],
          xhihilolo[i], xlohilolo[i], xhilololo[i], xlolololo[i],
          &xhihihihi0[i], &xhihihihi1[i], &xhihihihi2[i], &xhihihihi3[i],
          &xlohihihi0[i], &xlohihihi1[i], &xlohihihi2[i], &xlohihihi3[i],
          &xhilohihi0[i], &xhilohihi1[i], &xhilohihi2[i], &xhilohihi3[i],
          &xlolohihi0[i], &xlolohihi1[i], &xlolohihi2[i], &xlolohihi3[i],
          &xhihilohi0[i], &xhihilohi1[i], &xhihilohi2[i], &xhihilohi3[i],
          &xlohilohi0[i], &xlohilohi1[i], &xlohilohi2[i], &xlohilohi3[i],
          &xhilolohi0[i], &xhilolohi1[i], &xhilolohi2[i], &xhilolohi3[i],
          &xlololohi0[i], &xlololohi1[i], &xlololohi2[i], &xlololohi3[i],
          &xhihihilo0[i], &xhihihilo1[i], &xhihihilo2[i], &xhihihilo3[i],
          &xlohihilo0[i], &xlohihilo1[i], &xlohihilo2[i], &xlohihilo3[i],
          &xhilohilo0[i], &xhilohilo1[i], &xhilohilo2[i], &xhilohilo3[i],
          &xlolohilo0[i], &xlolohilo1[i], &xlolohilo2[i], &xlolohilo3[i],
          &xhihilolo0[i], &xhihilolo1[i], &xhihilolo2[i], &xhihilolo3[i],
          &xlohilolo0[i], &xlohilolo1[i], &xlohilolo2[i], &xlohilolo3[i],
          &xhilololo0[i], &xhilololo1[i], &xhilololo2[i], &xhilololo3[i],
          &xlolololo0[i], &xlolololo1[i], &xlolololo2[i], &xlolololo3[i]);
   }
}

void to_hexa_double
 ( double xhihihihi0, double xhihihihi1, double xhihihihi2, double xhihihihi3,
   double xlohihihi0, double xlohihihi1, double xlohihihi2, double xlohihihi3,
   double xhilohihi0, double xhilohihi1, double xhilohihi2, double xhilohihi3,
   double xlolohihi0, double xlolohihi1, double xlolohihi2, double xlolohihi3,
   double xhihilohi0, double xhihilohi1, double xhihilohi2, double xhihilohi3,
   double xlohilohi0, double xlohilohi1, double xlohilohi2, double xlohilohi3,
   double xhilolohi0, double xhilolohi1, double xhilolohi2, double xhilolohi3,
   double xlololohi0, double xlololohi1, double xlololohi2, double xlololohi3,
   double xhihihilo0, double xhihihilo1, double xhihihilo2, double xhihihilo3,
   double xlohihilo0, double xlohihilo1, double xlohihilo2, double xlohihilo3,
   double xhilohilo0, double xhilohilo1, double xhilohilo2, double xhilohilo3,
   double xlolohilo0, double xlolohilo1, double xlolohilo2, double xlolohilo3,
   double xhihilolo0, double xhihilolo1, double xhihilolo2, double xhihilolo3,
   double xlohilolo0, double xlohilolo1, double xlohilolo2, double xlohilolo3,
   double xhilololo0, double xhilololo1, double xhilololo2, double xhilololo3,
   double xlolololo0, double xlolololo1, double xlolololo2, double xlolololo3,
   double *xhihihihi, double *xlohihihi, double *xhilohihi, double *xlolohihi,
   double *xhihilohi, double *xlohilohi, double *xhilolohi, double *xlololohi,
   double *xhihihilo, double *xlohihilo, double *xhilohilo, double *xlolohilo,
   double *xhihilolo, double *xlohilolo, double *xhilololo, double *xlolololo )
{
   *xhihihihi = xhihihihi0; *xlohihihi = 0.0;
   *xhilohihi = 0.0; *xlolohihi = 0.0;
   *xhihilohi = 0.0; *xlohilohi = 0.0;
   *xhilolohi = 0.0; *xlololohi = 0.0;
   *xhihihilo = 0.0; *xlohihilo = 0.0;
   *xhilohilo = 0.0; *xlolohilo = 0.0;
   *xhihilolo = 0.0; *xlohilolo = 0.0;
   *xhilololo = 0.0; *xlolololo = 0.0;

   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhihihihi1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhihihihi2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhihihihi3);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohihihi0);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohihihi1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohihihi2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohihihi3);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilohihi0);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilohihi1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilohihi2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilohihi3);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlolohihi0);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlolohihi1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlolohihi2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlolohihi3);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhihilohi0);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhihilohi1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhihilohi2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhihilohi3);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohilohi0);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohilohi1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohilohi2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohilohi3);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilolohi0);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilolohi1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilolohi2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilolohi3);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlololohi0);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlololohi1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlololohi2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlololohi3);

   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhihihilo0);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhihihilo1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhihihilo2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhihihilo3);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohihilo0);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohihilo1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohihilo2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohihilo3);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilohilo0);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilohilo1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilohilo2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilohilo3);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlolohilo0);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlolohilo1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlolohilo2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlolohilo3);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhihilolo0);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhihilolo1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhihilolo2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhihilolo3);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohilolo0);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohilolo1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohilolo2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlohilolo3);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilololo0);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilololo1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilololo2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xhilololo3);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlolololo0);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlolololo1);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlolololo2);
   hdf_inc_d(xhihihihi, xlohihihi, xhilohihi, xlolohihi,
             xhihilohi, xlohilohi, xhilolohi, xlololohi,
             xhihihilo, xlohihilo, xhilohilo, xlolohilo,
             xhihilolo, xlohilolo, xhilololo, xlolololo, xlolololo3);
}

void hd_write_vector
 ( int dim,
   double *xhihihihi, double *xlohihihi, double *xhilohihi, double *xlolohihi,
   double *xhihilohi, double *xlohilohi, double *xhilolohi, double *xlololohi,
   double *xhihihilo, double *xlohihilo, double *xhilohilo, double *xlolohilo,
   double *xhihilolo, double *xlohilolo, double *xhilololo, double *xlolololo )
{
   double x[16];

   for(int i=0; i<dim; i++)
   {
      x[0] = xhihihihi[i]; x[1] = xlohihihi[i];
      x[2] = xhilohihi[i]; x[3] = xlolohihi[i];
      x[4] = xhihilohi[i]; x[5] = xlohilohi[i];
      x[6] = xhilolohi[i]; x[7] = xlololohi[i];
      x[8] = xhihihilo[i]; x[9] = xlohihilo[i];
      x[10] = xhilohilo[i]; x[11] = xlolohilo[i];
      x[12] = xhihilolo[i]; x[13] = xlohilolo[i];
      x[14] = xhilololo[i]; x[15] = xlolololo[i];
      hd_write_doubles(x);
   }
}

void hexa_double_product
 ( int dim,
   double *xhihihihi, double *xlohihihi, double *xhilohihi, double *xlolohihi,
   double *xhihilohi, double *xlohilohi, double *xhilolohi, double *xlololohi,
   double *xhihihilo, double *xlohihilo, double *xhilohilo, double *xlolohilo,
   double *xhihilolo, double *xlohilolo, double *xhilololo, double *xlolololo,
   double *yhihihihi, double *ylohihihi, double *yhilohihi, double *ylolohihi,
   double *yhihilohi, double *ylohilohi, double *yhilolohi, double *ylololohi,
   double *yhihihilo, double *ylohihilo, double *yhilohilo, double *ylolohilo,
   double *yhihilolo, double *ylohilolo, double *yhilololo, double *ylolololo,
   double *phihihihi, double *plohihihi, double *philohihi, double *plolohihi,
   double *phihilohi, double *plohilohi, double *philolohi, double *plololohi,
   double *phihihilo, double *plohihilo, double *philohilo, double *plolohilo,
   double *phihilolo, double *plohilolo, double *philololo, double *plolololo )
{
   double acchihihihi,acclohihihi,acchilohihi,acclolohihi;
   double acchihilohi,acclohilohi,acchilolohi,acclololohi;
   double acchihihilo,acclohihilo,acchilohilo,acclolohilo;
   double acchihilolo,acclohilolo,acchilololo,acclolololo;

   *phihihihi = 0.0; *plohihihi = 0.0; *philohihi = 0.0; *plolohihi = 0.0;
   *phihilohi = 0.0; *plohilohi = 0.0; *philolohi = 0.0; *plololohi = 0.0;
   *phihihilo = 0.0; *plohihilo = 0.0; *philohilo = 0.0; *plolohilo = 0.0;
   *phihilolo = 0.0; *plohilolo = 0.0; *philololo = 0.0; *plolololo = 0.0;

   for(int i=0; i<dim; i++)
   {
      hdf_mul(xhihihihi[i], xlohihihi[i], xhilohihi[i], xlolohihi[i],
              xhihilohi[i], xlohilohi[i], xhilolohi[i], xlololohi[i],
              xhihihilo[i], xlohihilo[i], xhilohilo[i], xlolohilo[i],
              xhihilolo[i], xlohilolo[i], xhilololo[i], xlolololo[i],
              yhihihihi[i], ylohihihi[i], yhilohihi[i], ylolohihi[i],
              yhihilohi[i], ylohilohi[i], yhilolohi[i], ylololohi[i],
              yhihihilo[i], ylohihilo[i], yhilohilo[i], ylolohilo[i],
              yhihilolo[i], ylohilolo[i], yhilololo[i], ylolololo[i],
              &acchihihihi, &acclohihihi, &acchilohihi, &acclolohihi,
              &acchihilohi, &acclohilohi, &acchilolohi, &acclololohi,
              &acchihihilo, &acclohihilo, &acchilohilo, &acclolohilo,
              &acchihilolo, &acclohilolo, &acchilololo, &acclolololo);
      hdf_inc(phihihihi, plohihihi, philohihi, plolohihi,
              phihilohi, plohilohi, philolohi, plololohi,
              phihihilo, plohihilo, philohilo, plolohilo,
              phihilolo, plohilolo, philololo, plolololo,
              acchihihihi, acclohihihi, acchilohihi, acclolohihi,
              acchihilohi, acclohilohi, acchilolohi, acclololohi,
              acchihihilo, acclohihilo, acchilohilo, acclolohilo,
              acchihilolo, acclohilolo, acchilololo, acclolololo);
   }
}

void vectored_hd_product
 ( int dim,
   double *x0, double *x1, double *x2, double *x3,
   double *x4, double *x5, double *x6, double *x7,
   double *x8, double *x9, double *x10, double *x11,
   double *x12, double *x13, double *x14, double *x15,
   double *x16, double *x17, double *x18, double *x19,
   double *x20, double *x21, double *x22, double *x23,
   double *x24, double *x25, double *x26, double *x27,
   double *x28, double *x29, double *x30, double *x31,
   double *x32, double *x33, double *x34, double *x35,
   double *x36, double *x37, double *x38, double *x39,
   double *x40, double *x41, double *x42, double *x43,
   double *x44, double *x45, double *x46, double *x47,
   double *x48, double *x49, double *x50, double *x51,
   double *x52, double *x53, double *x54, double *x55,
   double *x56, double *x57, double *x58, double *x59,
   double *x60, double *x61, double *x62, double *x63,
   double *y0, double *y1, double *y2, double *y3,
   double *y4, double *y5, double *y6, double *y7,
   double *y8, double *y9, double *y10, double *y11,
   double *y12, double *y13, double *y14, double *y15,
   double *y16, double *y17, double *y18, double *y19,
   double *y20, double *y21, double *y22, double *y23,
   double *y24, double *y25, double *y26, double *y27,
   double *y28, double *y29, double *y30, double *y31,
   double *y32, double *y33, double *y34, double *y35,
   double *y36, double *y37, double *y38, double *y39,
   double *y40, double *y41, double *y42, double *y43,
   double *y44, double *y45, double *y46, double *y47,
   double *y48, double *y49, double *y50, double *y51,
   double *y52, double *y53, double *y54, double *y55,
   double *y56, double *y57, double *y58, double *y59,
   double *y60, double *y61, double *y62, double *y63,
   double *s0, double *s1, double *s2, double *s3,
   double *s4, double *s5, double *s6, double *s7,
   double *s8, double *s9, double *s10, double *s11,
   double *s12, double *s13, double *s14, double *s15,
   double *s16, double *s17, double *s18, double *s19,
   double *s20, double *s21, double *s22, double *s23,
   double *s24, double *s25, double *s26, double *s27,
   double *s28, double *s29, double *s30, double *s31,
   double *s32, double *s33, double *s34, double *s35,
   double *s36, double *s37, double *s38, double *s39,
   double *s40, double *s41, double *s42, double *s43,
   double *s44, double *s45, double *s46, double *s47,
   double *s48, double *s49, double *s50, double *s51,
   double *s52, double *s53, double *s54, double *s55,
   double *s56, double *s57, double *s58, double *s59,
   double *s60, double *s61, double *s62, double *s63 )
{
   *s0 = 0.0; *s1 = 0.0; *s2 = 0.0; *s3 = 0.0;
   *s4 = 0.0; *s5 = 0.0; *s6 = 0.0; *s7 = 0.0;
   *s8 = 0.0; *s9 = 0.0; *s10 = 0.0; *s11 = 0.0;
   *s12 = 0.0; *s13 = 0.0; *s14 = 0.0; *s15 = 0.0;
   *s16 = 0.0; *s17 = 0.0; *s18 = 0.0; *s19 = 0.0;
   *s20 = 0.0; *s21 = 0.0; *s22 = 0.0; *s23 = 0.0;
   *s24 = 0.0; *s25 = 0.0; *s26 = 0.0; *s27 = 0.0;
   *s28 = 0.0; *s29 = 0.0; *s30 = 0.0; *s31 = 0.0;
   *s32 = 0.0; *s33 = 0.0; *s34 = 0.0; *s35 = 0.0;
   *s36 = 0.0; *s37 = 0.0; *s38 = 0.0; *s39 = 0.0;
   *s40 = 0.0; *s41 = 0.0; *s42 = 0.0; *s43 = 0.0;
   *s44 = 0.0; *s45 = 0.0; *s46 = 0.0; *s47 = 0.0;
   *s48 = 0.0; *s49 = 0.0; *s50 = 0.0; *s51 = 0.0;
   *s52 = 0.0; *s53 = 0.0; *s54 = 0.0; *s55 = 0.0;
   *s56 = 0.0; *s57 = 0.0; *s58 = 0.0; *s59 = 0.0;
   *s60 = 0.0; *s61 = 0.0; *s62 = 0.0; *s63 = 0.0;

   for(int i=0; i<dim; i++)
   {
      *s0 += x0[i]*y0[i];
      *s1 += x0[i]*y1[i] + x1[i]*y0[i];
      *s2 += x0[i]*y2[i] + x1[i]*y1[i] + x2[i]*y0[i];
      *s3 += x0[i]*y3[i] + x1[i]*y2[i] + x2[i]*y1[i] + x3[i]*y0[i];
      *s4 += x0[i]*y4[i] + x1[i]*y3[i] + x2[i]*y2[i] + x3[i]*y1[i]
           + x4[i]*y0[i];
      *s5 += x0[i]*y5[i] + x1[i]*y4[i] + x2[i]*y3[i] + x3[i]*y2[i]
           + x4[i]*y1[i] + x5[i]*y0[i];
      *s6 += x0[i]*y6[i] + x1[i]*y5[i] + x2[i]*y4[i] + x3[i]*y3[i]
           + x4[i]*y2[i] + x5[i]*y1[i] + x6[i]*y0[i];
      *s7 += x0[i]*y7[i] + x1[i]*y6[i] + x2[i]*y5[i] + x3[i]*y4[i]
           + x4[i]*y3[i] + x5[i]*y2[i] + x6[i]*y1[i] + x7[i]*y0[i];
      *s8 += x0[i]*y8[i] + x1[i]*y7[i] + x2[i]*y6[i] + x3[i]*y5[i]
           + x4[i]*y4[i] + x5[i]*y3[i] + x6[i]*y2[i] + x7[i]*y1[i]
           + x8[i]*y0[i];
      *s9 += x0[i]*y9[i] + x1[i]*y8[i] + x2[i]*y7[i] + x3[i]*y6[i]
           + x4[i]*y5[i] + x5[i]*y4[i] + x6[i]*y3[i] + x7[i]*y2[i]
           + x8[i]*y1[i] + x9[i]*y0[i];
      *s10 += x0[i]*y10[i] + x1[i]*y9[i] +  x2[i]*y8[i] + x3[i]*y7[i]
            + x4[i]*y6[i]  + x5[i]*y5[i] +  x6[i]*y4[i] + x7[i]*y3[i]
            + x8[i]*y2[i]  + x9[i]*y1[i] + x10[i]*y0[i];
      *s11 += x0[i]*y11[i] + x1[i]*y10[i] + x2[i]*y9[i] + x3[i]*y8[i]
            + x4[i]*y7[i]  + x5[i]*y6[i]  + x6[i]*y5[i] + x7[i]*y4[i]
            + x8[i]*y3[i]  + x9[i]*y2[i] + x10[i]*y1[i] + x11[i]*y0[i];
      *s12 += x0[i]*y12[i] + x1[i]*y11[i] + x2[i]*y10[i] + x3[i]*y9[i]
            + x4[i]*y8[i] + x5[i]*y7[i] + x6[i]*y6[i] + x7[i]*y5[i]
            + x8[i]*y4[i] + x9[i]*y3[i] + x10[i]*y2[i] + x11[i]*y1[i]
            + x12[i]*y0[i];
      *s13 += x0[i]*y13[i] + x1[i]*y12[i] + x2[i]*y11[i] + x3[i]*y10[i]
            + x4[i]*y9[i] + x5[i]*y8[i] + x6[i]*y7[i] + x7[i]*y6[i]
            + x8[i]*y5[i] + x9[i]*y4[i] + x10[i]*y3[i] + x11[i]*y2[i]
            + x12[i]*y1[i] + x13[i]*y0[i];
      *s14 += x0[i]*y14[i] + x1[i]*y13[i] + x2[i]*y12[i] + x3[i]*y11[i]
            + x4[i]*y10[i] + x5[i]*y9[i] + x6[i]*y8[i] + x7[i]*y7[i]
            + x8[i]*y6[i] + x9[i]*y5[i] + x10[i]*y4[i] + x11[i]*y3[i]
            + x12[i]*y2[i] + x13[i]*y1[i] + x14[i]*y0[i];
      *s15 += x0[i]*y15[i] + x1[i]*y14[i] + x2[i]*y13[i] + x3[i]*y12[i]
            + x4[i]*y11[i] + x5[i]*y10[i] + x6[i]*y9[i] + x7[i]*y8[i]
            + x8[i]*y7[i] + x9[i]*y6[i] + x10[i]*y5[i] + x11[i]*y4[i]
            + x12[i]*y3[i] + x13[i]*y2[i] + x14[i]*y1[i] + x15[i]*y0[i];
      *s16 += x0[i]*y16[i] + x1[i]*y15[i] + x2[i]*y14[i] + x3[i]*y13[i]
            + x4[i]*y12[i] + x5[i]*y11[i] + x6[i]*y10[i] + x7[i]*y9[i]
            + x8[i]*y8[i] + x9[i]*y7[i] + x10[i]*y6[i] + x11[i]*y5[i]
            + x12[i]*y4[i] + x13[i]*y3[i] + x14[i]*y2[i] + x15[i]*y1[i]
            + x16[i]*y0[i];
      *s17 += x0[i]*y17[i] + x1[i]*y16[i] + x2[i]*y15[i] + x3[i]*y14[i]
            + x4[i]*y13[i] + x5[i]*y12[i] + x6[i]*y11[i] + x7[i]*y10[i]
            + x8[i]*y9[i] + x9[i]*y8[i] + x10[i]*y7[i] + x11[i]*y6[i]
            + x12[i]*y5[i] + x13[i]*y4[i] + x14[i]*y3[i] + x15[i]*y2[i]
            + x16[i]*y1[i] + x17[i]*y0[i];
      *s18 += x0[i]*y18[i] + x1[i]*y17[i] + x2[i]*y16[i] + x3[i]*y15[i]
            + x4[i]*y14[i] + x5[i]*y13[i] + x6[i]*y12[i] + x7[i]*y11[i]
            + x8[i]*y10[i] + x9[i]*y9[i] + x10[i]*y8[i] + x11[i]*y7[i]
            + x12[i]*y6[i] + x13[i]*y5[i] + x14[i]*y4[i] + x15[i]*y3[i]
            + x16[i]*y2[i] + x17[i]*y1[i] + x18[i]*y0[i];
      *s19 += x0[i]*y19[i] + x1[i]*y18[i] + x2[i]*y17[i] + x3[i]*y16[i]
            + x4[i]*y15[i] + x5[i]*y14[i] + x6[i]*y13[i] + x7[i]*y12[i]
            + x8[i]*y11[i] + x9[i]*y10[i] + x10[i]*y9[i] + x11[i]*y8[i]
            + x12[i]*y7[i] + x13[i]*y6[i] + x14[i]*y5[i] + x15[i]*y4[i]
            + x16[i]*y3[i] + x17[i]*y2[i] + x18[i]*y1[i] + x19[i]*y0[i];
      *s20 += x0[i]*y20[i] + x1[i]*y19[i] + x2[i]*y18[i] + x3[i]*y17[i]
            + x4[i]*y16[i] + x5[i]*y15[i] + x6[i]*y14[i] + x7[i]*y13[i]
            + x8[i]*y12[i] + x9[i]*y11[i] + x10[i]*y10[i] + x11[i]*y9[i]
            + x12[i]*y8[i] + x13[i]*y7[i] + x14[i]*y6[i] + x15[i]*y5[i]
            + x16[i]*y4[i] + x17[i]*y3[i] + x18[i]*y2[i] + x19[i]*y1[i]
            + x20[i]*y0[i];
      *s21 += x0[i]*y21[i] + x1[i]*y20[i] + x2[i]*y19[i] + x3[i]*y18[i]
            + x4[i]*y17[i] + x5[i]*y16[i] + x6[i]*y15[i] + x7[i]*y14[i]
            + x8[i]*y13[i] + x9[i]*y12[i] + x10[i]*y11[i] + x11[i]*y10[i]
            + x12[i]*y9[i] + x13[i]*y8[i] + x14[i]*y7[i] + x15[i]*y6[i]
            + x16[i]*y5[i] + x17[i]*y4[i] + x18[i]*y3[i] + x19[i]*y2[i]
            + x20[i]*y1[i] + x21[i]*y0[i];
      *s22 += x0[i]*y22[i] + x1[i]*y21[i] + x2[i]*y20[i] + x3[i]*y19[i]
            + x4[i]*y18[i] + x5[i]*y17[i] + x6[i]*y16[i] + x7[i]*y15[i]
            + x8[i]*y14[i] + x9[i]*y13[i] + x10[i]*y12[i] + x11[i]*y11[i]
            + x12[i]*y10[i] + x13[i]*y9[i] + x14[i]*y8[i] + x15[i]*y7[i]
            + x16[i]*y6[i] + x17[i]*y5[i] + x18[i]*y4[i] + x19[i]*y3[i]
            + x20[i]*y2[i] + x21[i]*y1[i] + x22[i]*y0[i];
      *s23 += x0[i]*y23[i] + x1[i]*y22[i] + x2[i]*y21[i] + x3[i]*y20[i]
            + x4[i]*y19[i] + x5[i]*y18[i] + x6[i]*y17[i] + x7[i]*y16[i]
            + x8[i]*y15[i] + x9[i]*y14[i] + x10[i]*y13[i] + x11[i]*y12[i]
            + x12[i]*y11[i] + x13[i]*y10[i] + x14[i]*y9[i] + x15[i]*y8[i]
            + x16[i]*y7[i] + x17[i]*y6[i] + x18[i]*y5[i] + x19[i]*y4[i]
            + x20[i]*y3[i] + x21[i]*y2[i] + x22[i]*y1[i] + x23[i]*y0[i];
      *s24 += x0[i]*y24[i] + x1[i]*y23[i] + x2[i]*y22[i] + x3[i]*y21[i]
            + x4[i]*y20[i] + x5[i]*y19[i] + x6[i]*y18[i] + x7[i]*y17[i]
            + x8[i]*y16[i] + x9[i]*y15[i] + x10[i]*y14[i] + x11[i]*y13[i]
            + x12[i]*y12[i] + x13[i]*y11[i] + x14[i]*y10[i] + x15[i]*y9[i]
            + x16[i]*y8[i] + x17[i]*y7[i] + x18[i]*y6[i] + x19[i]*y5[i]
            + x20[i]*y4[i] + x21[i]*y3[i] + x22[i]*y2[i] + x23[i]*y1[i]
            + x24[i]*y0[i];
      *s25 += x0[i]*y25[i] + x1[i]*y24[i] + x2[i]*y23[i] + x3[i]*y22[i]
            + x4[i]*y21[i] + x5[i]*y20[i] + x6[i]*y19[i] + x7[i]*y18[i]
            + x8[i]*y17[i] + x9[i]*y16[i] + x10[i]*y15[i] + x11[i]*y14[i]
            + x12[i]*y13[i] + x13[i]*y12[i] + x14[i]*y11[i] + x15[i]*y10[i]
            + x16[i]*y9[i] + x17[i]*y8[i] + x18[i]*y7[i] + x19[i]*y6[i]
            + x20[i]*y5[i] + x21[i]*y4[i] + x22[i]*y3[i] + x23[i]*y2[i]
            + x24[i]*y1[i] + x25[i]*y0[i];
      *s26 += x0[i]*y26[i] + x1[i]*y25[i] + x2[i]*y24[i] + x3[i]*y23[i]
            + x4[i]*y22[i] + x5[i]*y21[i] + x6[i]*y20[i] + x7[i]*y19[i]
            + x8[i]*y18[i] + x9[i]*y17[i] + x10[i]*y16[i] + x11[i]*y15[i]
            + x12[i]*y14[i] + x13[i]*y13[i] + x14[i]*y12[i] + x15[i]*y11[i]
            + x16[i]*y10[i] + x17[i]*y9[i] + x18[i]*y8[i] + x19[i]*y7[i]
            + x20[i]*y6[i] + x21[i]*y5[i] + x22[i]*y4[i] + x23[i]*y3[i]
            + x24[i]*y2[i] + x25[i]*y1[i] + x26[i]*y0[i];
      *s27 += x0[i]*y27[i] + x1[i]*y26[i] + x2[i]*y25[i] + x3[i]*y24[i]
            + x4[i]*y23[i] + x5[i]*y22[i] + x6[i]*y21[i] + x7[i]*y20[i]
            + x8[i]*y19[i] + x9[i]*y18[i] + x10[i]*y17[i] + x11[i]*y16[i]
            + x12[i]*y15[i] + x13[i]*y14[i] + x14[i]*y13[i] + x15[i]*y12[i]
            + x16[i]*y11[i] + x17[i]*y10[i] + x18[i]*y9[i] + x19[i]*y8[i]
            + x20[i]*y7[i] + x21[i]*y6[i] + x22[i]*y5[i] + x23[i]*y4[i]
            + x24[i]*y3[i] + x25[i]*y2[i] + x26[i]*y1[i] + x27[i]*y0[i];
      *s28 += x0[i]*y28[i] + x1[i]*y27[i] + x2[i]*y26[i] + x3[i]*y25[i]
            + x4[i]*y24[i] + x5[i]*y23[i] + x6[i]*y22[i] + x7[i]*y21[i]
            + x8[i]*y20[i] + x9[i]*y19[i] + x10[i]*y18[i] + x11[i]*y17[i]
            + x12[i]*y16[i] + x13[i]*y15[i] + x14[i]*y14[i] + x15[i]*y13[i]
            + x16[i]*y12[i] + x17[i]*y11[i] + x18[i]*y10[i] + x19[i]*y9[i]
            + x20[i]*y8[i] + x21[i]*y7[i] + x22[i]*y6[i] + x23[i]*y5[i]
            + x24[i]*y4[i] + x25[i]*y3[i] + x26[i]*y2[i] + x27[i]*y1[i]
            + x28[i]*y0[i];
      *s29 += x0[i]*y29[i] + x1[i]*y28[i] + x2[i]*y27[i] + x3[i]*y26[i]
            + x4[i]*y25[i] + x5[i]*y24[i] + x6[i]*y23[i] + x7[i]*y22[i]
            + x8[i]*y21[i] + x9[i]*y20[i] + x10[i]*y19[i] + x11[i]*y18[i]
            + x12[i]*y17[i] + x13[i]*y16[i] + x14[i]*y15[i] + x15[i]*y14[i]
            + x16[i]*y13[i] + x17[i]*y12[i] + x18[i]*y11[i] + x19[i]*y10[i]
            + x20[i]*y9[i] + x21[i]*y8[i] + x22[i]*y7[i] + x23[i]*y6[i]
            + x24[i]*y5[i] + x25[i]*y4[i] + x26[i]*y3[i] + x27[i]*y2[i]
            + x28[i]*y1[i] + x29[i]*y0[i];
      *s30 += x0[i]*y30[i] + x1[i]*y29[i] + x2[i]*y28[i] + x3[i]*y27[i]
            + x4[i]*y26[i] + x5[i]*y25[i] + x6[i]*y24[i] + x7[i]*y23[i]
            + x8[i]*y22[i] + x9[i]*y21[i] + x10[i]*y20[i] + x11[i]*y19[i]
            + x12[i]*y18[i] + x13[i]*y17[i] + x14[i]*y16[i] + x15[i]*y15[i]
            + x16[i]*y14[i] + x17[i]*y13[i] + x18[i]*y12[i] + x19[i]*y11[i]
            + x20[i]*y10[i] + x21[i]*y9[i] + x22[i]*y8[i] + x23[i]*y7[i]
            + x24[i]*y6[i] + x25[i]*y5[i] + x26[i]*y4[i] + x27[i]*y3[i]
            + x28[i]*y2[i] + x29[i]*y1[i] + x30[i]*y0[i];
      *s31 += x0[i]*y31[i] + x1[i]*y30[i] + x2[i]*y29[i] + x3[i]*y28[i]
            + x4[i]*y27[i] + x5[i]*y26[i] + x6[i]*y25[i] + x7[i]*y24[i]
            + x8[i]*y23[i] + x9[i]*y22[i] + x10[i]*y21[i] + x11[i]*y20[i]
            + x12[i]*y19[i] + x13[i]*y18[i] + x14[i]*y17[i] + x15[i]*y16[i]
            + x16[i]*y15[i] + x17[i]*y14[i] + x18[i]*y13[i] + x19[i]*y12[i]
            + x20[i]*y11[i] + x21[i]*y10[i] + x22[i]*y9[i] + x23[i]*y8[i]
            + x24[i]*y7[i] + x25[i]*y6[i] + x26[i]*y5[i] + x27[i]*y4[i]
            + x28[i]*y3[i] + x29[i]*y2[i] + x30[i]*y1[i] + x31[i]*y0[i];
      *s32 += x0[i]*y32[i] + x1[i]*y31[i] + x2[i]*y30[i] + x3[i]*y29[i]
            + x4[i]*y28[i] + x5[i]*y27[i] + x6[i]*y26[i] + x7[i]*y25[i]
            + x8[i]*y24[i] + x9[i]*y23[i] + x10[i]*y22[i] + x11[i]*y21[i]
            + x12[i]*y20[i] + x13[i]*y19[i] + x14[i]*y18[i] + x15[i]*y17[i]
            + x16[i]*y16[i] + x17[i]*y15[i] + x18[i]*y14[i] + x19[i]*y13[i]
            + x20[i]*y12[i] + x21[i]*y11[i] + x22[i]*y10[i] + x23[i]*y9[i]
            + x24[i]*y8[i] + x25[i]*y7[i] + x26[i]*y6[i] + x27[i]*y5[i]
            + x28[i]*y4[i] + x29[i]*y3[i] + x30[i]*y2[i] + x31[i]*y1[i]
            + x32[i]*y0[i];
      *s33 += x0[i]*y33[i] + x1[i]*y32[i] + x2[i]*y31[i] + x3[i]*y30[i]
            + x4[i]*y29[i] + x5[i]*y28[i] + x6[i]*y27[i] + x7[i]*y26[i]
            + x8[i]*y25[i] + x9[i]*y24[i] + x10[i]*y23[i] + x11[i]*y22[i]
            + x12[i]*y21[i] + x13[i]*y20[i] + x14[i]*y19[i] + x15[i]*y18[i]
            + x16[i]*y17[i] + x17[i]*y16[i] + x18[i]*y15[i] + x19[i]*y14[i]
            + x20[i]*y13[i] + x21[i]*y12[i] + x22[i]*y11[i] + x23[i]*y10[i]
            + x24[i]*y9[i] + x25[i]*y8[i] + x26[i]*y7[i] + x27[i]*y6[i]
            + x28[i]*y5[i] + x29[i]*y4[i] + x30[i]*y3[i] + x31[i]*y2[i]
            + x32[i]*y1[i] + x33[i]*y0[i];
      *s34 += x0[i]*y34[i] + x1[i]*y33[i] + x2[i]*y32[i] + x3[i]*y31[i]
            + x4[i]*y30[i] + x5[i]*y29[i] + x6[i]*y28[i] + x7[i]*y27[i]
            + x8[i]*y26[i] + x9[i]*y25[i] + x10[i]*y24[i] + x11[i]*y23[i]
            + x12[i]*y22[i] + x13[i]*y21[i] + x14[i]*y20[i] + x15[i]*y19[i]
            + x16[i]*y18[i] + x17[i]*y17[i] + x18[i]*y16[i] + x19[i]*y15[i]
            + x20[i]*y14[i] + x21[i]*y13[i] + x22[i]*y12[i] + x23[i]*y11[i]
            + x24[i]*y10[i] + x25[i]*y9[i] + x26[i]*y8[i] + x27[i]*y7[i]
            + x28[i]*y6[i] + x29[i]*y5[i] + x30[i]*y4[i] + x31[i]*y3[i]
            + x32[i]*y2[i] + x33[i]*y1[i] + x34[i]*y0[i];
      *s35 += x0[i]*y35[i] + x1[i]*y34[i] + x2[i]*y33[i] + x3[i]*y32[i]
            + x4[i]*y31[i] + x5[i]*y30[i] + x6[i]*y29[i] + x7[i]*y28[i]
            + x8[i]*y27[i] + x9[i]*y26[i] + x10[i]*y25[i] + x11[i]*y24[i]
            + x12[i]*y23[i] + x13[i]*y22[i] + x14[i]*y21[i] + x15[i]*y20[i]
            + x16[i]*y19[i] + x17[i]*y18[i] + x18[i]*y17[i] + x19[i]*y16[i]
            + x20[i]*y15[i] + x21[i]*y14[i] + x22[i]*y13[i] + x23[i]*y12[i]
            + x24[i]*y11[i] + x25[i]*y10[i] + x26[i]*y9[i] + x27[i]*y8[i]
            + x28[i]*y7[i] + x29[i]*y6[i] + x30[i]*y5[i] + x31[i]*y4[i]
            + x32[i]*y3[i] + x33[i]*y2[i] + x34[i]*y1[i] + x35[i]*y0[i];
      *s36 += x0[i]*y36[i] + x1[i]*y35[i] + x2[i]*y34[i] + x3[i]*y33[i]
            + x4[i]*y32[i] + x5[i]*y31[i] + x6[i]*y30[i] + x7[i]*y29[i]
            + x8[i]*y28[i] + x9[i]*y27[i] + x10[i]*y26[i] + x11[i]*y25[i]
            + x12[i]*y24[i] + x13[i]*y23[i] + x14[i]*y22[i] + x15[i]*y21[i]
            + x16[i]*y20[i] + x17[i]*y19[i] + x18[i]*y18[i] + x19[i]*y17[i]
            + x20[i]*y16[i] + x21[i]*y15[i] + x22[i]*y14[i] + x23[i]*y13[i]
            + x24[i]*y12[i] + x25[i]*y11[i] + x26[i]*y10[i] + x27[i]*y9[i]
            + x28[i]*y8[i] + x29[i]*y7[i] + x30[i]*y6[i] + x31[i]*y5[i]
            + x32[i]*y4[i] + x33[i]*y3[i] + x34[i]*y2[i] + x35[i]*y1[i]
            + x36[i]*y0[i];
      *s37 += x0[i]*y37[i] + x1[i]*y36[i] + x2[i]*y35[i] + x3[i]*y34[i]
            + x4[i]*y33[i] + x5[i]*y32[i] + x6[i]*y31[i] + x7[i]*y30[i]
            + x8[i]*y29[i] + x9[i]*y28[i] + x10[i]*y27[i] + x11[i]*y26[i]
            + x12[i]*y25[i] + x13[i]*y24[i] + x14[i]*y23[i] + x15[i]*y22[i]
            + x16[i]*y21[i] + x17[i]*y20[i] + x18[i]*y19[i] + x19[i]*y18[i]
            + x20[i]*y17[i] + x21[i]*y16[i] + x22[i]*y15[i] + x23[i]*y14[i]
            + x24[i]*y13[i] + x25[i]*y12[i] + x26[i]*y11[i] + x27[i]*y10[i]
            + x28[i]*y9[i] + x29[i]*y8[i] + x30[i]*y7[i] + x31[i]*y6[i]
            + x32[i]*y5[i] + x33[i]*y4[i] + x34[i]*y3[i] + x35[i]*y2[i]
            + x36[i]*y1[i] + x37[i]*y0[i];
      *s38 += x0[i]*y38[i] + x1[i]*y37[i] + x2[i]*y36[i] + x3[i]*y35[i]
            + x4[i]*y34[i] + x5[i]*y33[i] + x6[i]*y32[i] + x7[i]*y31[i]
            + x8[i]*y30[i] + x9[i]*y29[i] + x10[i]*y28[i] + x11[i]*y27[i]
            + x12[i]*y26[i] + x13[i]*y25[i] + x14[i]*y24[i] + x15[i]*y23[i]
            + x16[i]*y22[i] + x17[i]*y21[i] + x18[i]*y20[i] + x19[i]*y19[i]
            + x20[i]*y18[i] + x21[i]*y17[i] + x22[i]*y16[i] + x23[i]*y15[i]
            + x24[i]*y14[i] + x25[i]*y13[i] + x26[i]*y12[i] + x27[i]*y11[i]
            + x28[i]*y10[i] + x29[i]*y9[i] + x30[i]*y8[i] + x31[i]*y7[i]
            + x32[i]*y6[i] + x33[i]*y5[i] + x34[i]*y4[i] + x35[i]*y3[i]
            + x36[i]*y2[i] + x37[i]*y1[i] + x38[i]*y0[i];
      *s39 += x0[i]*y39[i] + x1[i]*y38[i] + x2[i]*y37[i] + x3[i]*y36[i]
            + x4[i]*y35[i] + x5[i]*y34[i] + x6[i]*y33[i] + x7[i]*y32[i]
            + x8[i]*y31[i] + x9[i]*y30[i] + x10[i]*y29[i] + x11[i]*y28[i]
            + x12[i]*y27[i] + x13[i]*y26[i] + x14[i]*y25[i] + x15[i]*y24[i]
            + x16[i]*y23[i] + x17[i]*y22[i] + x18[i]*y21[i] + x19[i]*y20[i]
            + x20[i]*y19[i] + x21[i]*y18[i] + x22[i]*y17[i] + x23[i]*y16[i]
            + x24[i]*y15[i] + x25[i]*y14[i] + x26[i]*y13[i] + x27[i]*y12[i]
            + x28[i]*y11[i] + x29[i]*y10[i] + x30[i]*y9[i] + x31[i]*y8[i]
            + x32[i]*y7[i] + x33[i]*y6[i] + x34[i]*y5[i] + x35[i]*y4[i]
            + x36[i]*y3[i] + x37[i]*y2[i] + x38[i]*y1[i] + x39[i]*y0[i];
      *s40 += x0[i]*y40[i] + x1[i]*y39[i] + x2[i]*y38[i] + x3[i]*y37[i]
            + x4[i]*y36[i] + x5[i]*y35[i] + x6[i]*y34[i] + x7[i]*y33[i]
            + x8[i]*y32[i] + x9[i]*y31[i] + x10[i]*y30[i] + x11[i]*y29[i]
            + x12[i]*y28[i] + x13[i]*y27[i] + x14[i]*y26[i] + x15[i]*y25[i]
            + x16[i]*y24[i] + x17[i]*y23[i] + x18[i]*y22[i] + x19[i]*y21[i]
            + x20[i]*y20[i] + x21[i]*y19[i] + x22[i]*y18[i] + x23[i]*y17[i]
            + x24[i]*y16[i] + x25[i]*y15[i] + x26[i]*y14[i] + x27[i]*y13[i]
            + x28[i]*y12[i] + x29[i]*y11[i] + x30[i]*y10[i] + x31[i]*y9[i]
            + x32[i]*y8[i] + x33[i]*y7[i] + x34[i]*y6[i] + x35[i]*y5[i]
            + x36[i]*y4[i] + x37[i]*y3[i] + x38[i]*y2[i] + x39[i]*y1[i]
            + x40[i]*y0[i];
      *s41 += x0[i]*y41[i] + x1[i]*y40[i] + x2[i]*y39[i] + x3[i]*y38[i]
            + x4[i]*y37[i] + x5[i]*y36[i] + x6[i]*y35[i] + x7[i]*y34[i]
            + x8[i]*y33[i] + x9[i]*y32[i] + x10[i]*y31[i] + x11[i]*y30[i]
            + x12[i]*y29[i] + x13[i]*y28[i] + x14[i]*y27[i] + x15[i]*y26[i]
            + x16[i]*y25[i] + x17[i]*y24[i] + x18[i]*y23[i] + x19[i]*y22[i]
            + x20[i]*y21[i] + x21[i]*y20[i] + x22[i]*y19[i] + x23[i]*y18[i]
            + x24[i]*y17[i] + x25[i]*y16[i] + x26[i]*y15[i] + x27[i]*y14[i]
            + x28[i]*y13[i] + x29[i]*y12[i] + x30[i]*y11[i] + x31[i]*y10[i]
            + x32[i]*y9[i] + x33[i]*y8[i] + x34[i]*y7[i] + x35[i]*y6[i]
            + x36[i]*y5[i] + x37[i]*y4[i] + x38[i]*y3[i] + x39[i]*y2[i]
            + x40[i]*y1[i] + x41[i]*y0[i];
      *s42 += x0[i]*y42[i] + x1[i]*y41[i] + x2[i]*y40[i] + x3[i]*y39[i]
            + x4[i]*y38[i] + x5[i]*y37[i] + x6[i]*y36[i] + x7[i]*y35[i]
            + x8[i]*y34[i] + x9[i]*y33[i] + x10[i]*y32[i] + x11[i]*y31[i]
            + x12[i]*y30[i] + x13[i]*y29[i] + x14[i]*y28[i] + x15[i]*y27[i]
            + x16[i]*y26[i] + x17[i]*y25[i] + x18[i]*y24[i] + x19[i]*y23[i]
            + x20[i]*y22[i] + x21[i]*y21[i] + x22[i]*y20[i] + x23[i]*y19[i]
            + x24[i]*y18[i] + x25[i]*y17[i] + x26[i]*y16[i] + x27[i]*y15[i]
            + x28[i]*y14[i] + x29[i]*y13[i] + x30[i]*y12[i] + x31[i]*y11[i]
            + x32[i]*y10[i] + x33[i]*y9[i] + x34[i]*y8[i] + x35[i]*y7[i]
            + x36[i]*y6[i] + x37[i]*y5[i] + x38[i]*y4[i] + x39[i]*y3[i]
            + x40[i]*y2[i] + x41[i]*y1[i] + x42[i]*y0[i];
      *s43 += x0[i]*y43[i] + x1[i]*y42[i] + x2[i]*y41[i] + x3[i]*y40[i]
            + x4[i]*y39[i] + x5[i]*y38[i] + x6[i]*y37[i] + x7[i]*y36[i]
            + x8[i]*y35[i] + x9[i]*y34[i] + x10[i]*y33[i] + x11[i]*y32[i]
            + x12[i]*y31[i] + x13[i]*y30[i] + x14[i]*y29[i] + x15[i]*y28[i]
            + x16[i]*y27[i] + x17[i]*y26[i] + x18[i]*y25[i] + x19[i]*y24[i]
            + x20[i]*y23[i] + x21[i]*y22[i] + x22[i]*y21[i] + x23[i]*y20[i]
            + x24[i]*y19[i] + x25[i]*y18[i] + x26[i]*y17[i] + x27[i]*y16[i]
            + x28[i]*y15[i] + x29[i]*y14[i] + x30[i]*y13[i] + x31[i]*y12[i]
            + x32[i]*y11[i] + x33[i]*y10[i] + x34[i]*y9[i] + x35[i]*y8[i]
            + x36[i]*y7[i] + x37[i]*y6[i] + x38[i]*y5[i] + x39[i]*y4[i]
            + x40[i]*y3[i] + x41[i]*y2[i] + x42[i]*y1[i] + x43[i]*y0[i];
      *s44 += x0[i]*y44[i] + x1[i]*y43[i] + x2[i]*y42[i] + x3[i]*y41[i]
            + x4[i]*y40[i] + x5[i]*y39[i] + x6[i]*y38[i] + x7[i]*y37[i]
            + x8[i]*y36[i] + x9[i]*y35[i] + x10[i]*y34[i] + x11[i]*y33[i]
            + x12[i]*y32[i] + x13[i]*y31[i] + x14[i]*y30[i] + x15[i]*y29[i]
            + x16[i]*y28[i] + x17[i]*y27[i] + x18[i]*y26[i] + x19[i]*y25[i]
            + x20[i]*y24[i] + x21[i]*y23[i] + x22[i]*y22[i] + x23[i]*y21[i]
            + x24[i]*y20[i] + x25[i]*y19[i] + x26[i]*y18[i] + x27[i]*y17[i]
            + x28[i]*y16[i] + x29[i]*y15[i] + x30[i]*y14[i] + x31[i]*y13[i]
            + x32[i]*y12[i] + x33[i]*y11[i] + x34[i]*y10[i] + x35[i]*y9[i]
            + x36[i]*y8[i] + x37[i]*y7[i] + x38[i]*y6[i] + x39[i]*y5[i]
            + x40[i]*y4[i] + x41[i]*y3[i] + x42[i]*y2[i] + x43[i]*y1[i]
            + x44[i]*y0[i];
      *s45 += x0[i]*y45[i] + x1[i]*y44[i] + x2[i]*y43[i] + x3[i]*y42[i]
            + x4[i]*y41[i] + x5[i]*y40[i] + x6[i]*y39[i] + x7[i]*y38[i]
            + x8[i]*y37[i] + x9[i]*y36[i] + x10[i]*y35[i] + x11[i]*y34[i]
            + x12[i]*y33[i] + x13[i]*y32[i] + x14[i]*y31[i] + x15[i]*y30[i]
            + x16[i]*y29[i] + x17[i]*y28[i] + x18[i]*y27[i] + x19[i]*y26[i]
            + x20[i]*y25[i] + x21[i]*y24[i] + x22[i]*y23[i] + x23[i]*y22[i]
            + x24[i]*y21[i] + x25[i]*y20[i] + x26[i]*y19[i] + x27[i]*y18[i]
            + x28[i]*y17[i] + x29[i]*y16[i] + x30[i]*y15[i] + x31[i]*y14[i]
            + x32[i]*y13[i] + x33[i]*y12[i] + x34[i]*y11[i] + x35[i]*y10[i]
            + x36[i]*y9[i] + x37[i]*y8[i] + x38[i]*y7[i] + x39[i]*y6[i]
            + x40[i]*y5[i] + x41[i]*y4[i] + x42[i]*y3[i] + x43[i]*y2[i]
            + x44[i]*y1[i] + x45[i]*y0[i];
      *s46 += x0[i]*y46[i] + x1[i]*y45[i] + x2[i]*y44[i] + x3[i]*y43[i]
            + x4[i]*y42[i] + x5[i]*y41[i] + x6[i]*y40[i] + x7[i]*y39[i]
            + x8[i]*y38[i] + x9[i]*y37[i] + x10[i]*y36[i] + x11[i]*y35[i]
            + x12[i]*y34[i] + x13[i]*y33[i] + x14[i]*y32[i] + x15[i]*y31[i]
            + x16[i]*y30[i] + x17[i]*y29[i] + x18[i]*y28[i] + x19[i]*y27[i]
            + x20[i]*y26[i] + x21[i]*y25[i] + x22[i]*y24[i] + x23[i]*y23[i]
            + x24[i]*y22[i] + x25[i]*y21[i] + x26[i]*y20[i] + x27[i]*y19[i]
            + x28[i]*y18[i] + x29[i]*y17[i] + x30[i]*y16[i] + x31[i]*y15[i]
            + x32[i]*y14[i] + x33[i]*y13[i] + x34[i]*y12[i] + x35[i]*y11[i]
            + x36[i]*y10[i] + x37[i]*y9[i] + x38[i]*y8[i] + x39[i]*y7[i]
            + x40[i]*y6[i] + x41[i]*y5[i] + x42[i]*y4[i] + x43[i]*y3[i]
            + x44[i]*y2[i] + x45[i]*y1[i] + x46[i]*y0[i];
      *s47 += x0[i]*y47[i] + x1[i]*y46[i] + x2[i]*y45[i] + x3[i]*y44[i]
            + x4[i]*y43[i] + x5[i]*y42[i] + x6[i]*y41[i] + x7[i]*y40[i]
            + x8[i]*y39[i] + x9[i]*y38[i] + x10[i]*y37[i] + x11[i]*y36[i]
            + x12[i]*y35[i] + x13[i]*y34[i] + x14[i]*y33[i] + x15[i]*y32[i]
            + x16[i]*y31[i] + x17[i]*y30[i] + x18[i]*y29[i] + x19[i]*y28[i]
            + x20[i]*y27[i] + x21[i]*y26[i] + x22[i]*y25[i] + x23[i]*y24[i]
            + x24[i]*y23[i] + x25[i]*y22[i] + x26[i]*y21[i] + x27[i]*y20[i]
            + x28[i]*y19[i] + x29[i]*y18[i] + x30[i]*y17[i] + x31[i]*y16[i]
            + x32[i]*y15[i] + x33[i]*y14[i] + x34[i]*y13[i] + x35[i]*y12[i]
            + x36[i]*y11[i] + x37[i]*y10[i] + x38[i]*y9[i] + x39[i]*y8[i]
            + x40[i]*y7[i] + x41[i]*y6[i] + x42[i]*y5[i] + x43[i]*y4[i]
            + x44[i]*y3[i] + x45[i]*y2[i] + x46[i]*y1[i] + x47[i]*y0[i];
      *s48 += x0[i]*y48[i] + x1[i]*y47[i] + x2[i]*y46[i] + x3[i]*y45[i]
            + x4[i]*y44[i] + x5[i]*y43[i] + x6[i]*y42[i] + x7[i]*y41[i]
            + x8[i]*y40[i] + x9[i]*y39[i] + x10[i]*y38[i] + x11[i]*y37[i]
            + x12[i]*y36[i] + x13[i]*y35[i] + x14[i]*y34[i] + x15[i]*y33[i]
            + x16[i]*y32[i] + x17[i]*y31[i] + x18[i]*y30[i] + x19[i]*y29[i]
            + x20[i]*y28[i] + x21[i]*y27[i] + x22[i]*y26[i] + x23[i]*y25[i]
            + x24[i]*y24[i] + x25[i]*y23[i] + x26[i]*y22[i] + x27[i]*y21[i]
            + x28[i]*y20[i] + x29[i]*y19[i] + x30[i]*y18[i] + x31[i]*y17[i]
            + x32[i]*y16[i] + x33[i]*y15[i] + x34[i]*y14[i] + x35[i]*y13[i]
            + x36[i]*y12[i] + x37[i]*y11[i] + x38[i]*y10[i] + x39[i]*y9[i]
            + x40[i]*y8[i] + x41[i]*y7[i] + x42[i]*y6[i] + x43[i]*y5[i]
            + x44[i]*y4[i] + x45[i]*y3[i] + x46[i]*y2[i] + x47[i]*y1[i]
            + x48[i]*y0[i];
      *s49 += x0[i]*y49[i] + x1[i]*y48[i] + x2[i]*y47[i] + x3[i]*y46[i]
            + x4[i]*y45[i] + x5[i]*y44[i] + x6[i]*y43[i] + x7[i]*y42[i]
            + x8[i]*y41[i] + x9[i]*y40[i] + x10[i]*y39[i] + x11[i]*y38[i]
            + x12[i]*y37[i] + x13[i]*y36[i] + x14[i]*y35[i] + x15[i]*y34[i]
            + x16[i]*y33[i] + x17[i]*y32[i] + x18[i]*y31[i] + x19[i]*y30[i]
            + x20[i]*y29[i] + x21[i]*y28[i] + x22[i]*y27[i] + x23[i]*y26[i]
            + x24[i]*y25[i] + x25[i]*y24[i] + x26[i]*y23[i] + x27[i]*y22[i]
            + x28[i]*y21[i] + x29[i]*y20[i] + x30[i]*y19[i] + x31[i]*y18[i]
            + x32[i]*y17[i] + x33[i]*y16[i] + x34[i]*y15[i] + x35[i]*y14[i]
            + x36[i]*y13[i] + x37[i]*y12[i] + x38[i]*y11[i] + x39[i]*y10[i]
            + x40[i]*y9[i] + x41[i]*y8[i] + x42[i]*y7[i] + x43[i]*y6[i]
            + x44[i]*y5[i] + x45[i]*y4[i] + x46[i]*y3[i] + x47[i]*y2[i]
            + x48[i]*y1[i] + x49[i]*y0[i];
      *s50 += x0[i]*y50[i] + x1[i]*y49[i] + x2[i]*y48[i] + x3[i]*y47[i]
            + x4[i]*y46[i] + x5[i]*y45[i] + x6[i]*y44[i] + x7[i]*y43[i]
            + x8[i]*y42[i] + x9[i]*y41[i] + x10[i]*y40[i] + x11[i]*y39[i]
            + x12[i]*y38[i] + x13[i]*y37[i] + x14[i]*y36[i] + x15[i]*y35[i]
            + x16[i]*y34[i] + x17[i]*y33[i] + x18[i]*y32[i] + x19[i]*y31[i]
            + x20[i]*y30[i] + x21[i]*y29[i] + x22[i]*y28[i] + x23[i]*y27[i]
            + x24[i]*y26[i] + x25[i]*y25[i] + x26[i]*y24[i] + x27[i]*y23[i]
            + x28[i]*y22[i] + x29[i]*y21[i] + x30[i]*y20[i] + x31[i]*y19[i]
            + x32[i]*y18[i] + x33[i]*y17[i] + x34[i]*y16[i] + x35[i]*y15[i]
            + x36[i]*y14[i] + x37[i]*y13[i] + x38[i]*y12[i] + x39[i]*y11[i]
            + x40[i]*y10[i] + x41[i]*y9[i] + x42[i]*y8[i] + x43[i]*y7[i]
            + x44[i]*y6[i] + x45[i]*y5[i] + x46[i]*y4[i] + x47[i]*y3[i]
            + x48[i]*y2[i] + x49[i]*y1[i] + x50[i]*y0[i];
      *s51 += x0[i]*y51[i] + x1[i]*y50[i] + x2[i]*y49[i] + x3[i]*y48[i]
            + x4[i]*y47[i] + x5[i]*y46[i] + x6[i]*y45[i] + x7[i]*y44[i]
            + x8[i]*y43[i] + x9[i]*y42[i] + x10[i]*y41[i] + x11[i]*y40[i]
            + x12[i]*y39[i] + x13[i]*y38[i] + x14[i]*y37[i] + x15[i]*y36[i]
            + x16[i]*y35[i] + x17[i]*y34[i] + x18[i]*y33[i] + x19[i]*y32[i]
            + x20[i]*y31[i] + x21[i]*y30[i] + x22[i]*y29[i] + x23[i]*y28[i]
            + x24[i]*y27[i] + x25[i]*y26[i] + x26[i]*y25[i] + x27[i]*y24[i]
            + x28[i]*y23[i] + x29[i]*y22[i] + x30[i]*y21[i] + x31[i]*y20[i]
            + x32[i]*y19[i] + x33[i]*y18[i] + x34[i]*y17[i] + x35[i]*y16[i]
            + x36[i]*y15[i] + x37[i]*y14[i] + x38[i]*y13[i] + x39[i]*y12[i]
            + x40[i]*y11[i] + x41[i]*y10[i] + x42[i]*y9[i] + x43[i]*y8[i]
            + x44[i]*y7[i] + x45[i]*y6[i] + x46[i]*y5[i] + x47[i]*y4[i]
            + x48[i]*y3[i] + x49[i]*y2[i] + x50[i]*y1[i] + x51[i]*y0[i];
      *s52 += x0[i]*y52[i] + x1[i]*y51[i] + x2[i]*y50[i] + x3[i]*y49[i]
            + x4[i]*y48[i] + x5[i]*y47[i] + x6[i]*y46[i] + x7[i]*y45[i]
            + x8[i]*y44[i] + x9[i]*y43[i] + x10[i]*y42[i] + x11[i]*y41[i]
            + x12[i]*y40[i] + x13[i]*y39[i] + x14[i]*y38[i] + x15[i]*y37[i]
            + x16[i]*y36[i] + x17[i]*y35[i] + x18[i]*y34[i] + x19[i]*y33[i]
            + x20[i]*y32[i] + x21[i]*y31[i] + x22[i]*y30[i] + x23[i]*y29[i]
            + x24[i]*y28[i] + x25[i]*y27[i] + x26[i]*y26[i] + x27[i]*y25[i]
            + x28[i]*y24[i] + x29[i]*y23[i] + x30[i]*y22[i] + x31[i]*y21[i]
            + x32[i]*y20[i] + x33[i]*y19[i] + x34[i]*y18[i] + x35[i]*y17[i]
            + x36[i]*y16[i] + x37[i]*y15[i] + x38[i]*y14[i] + x39[i]*y13[i]
            + x40[i]*y12[i] + x41[i]*y11[i] + x42[i]*y10[i] + x43[i]*y9[i]
            + x44[i]*y8[i] + x45[i]*y7[i] + x46[i]*y6[i] + x47[i]*y5[i]
            + x48[i]*y4[i] + x49[i]*y3[i] + x50[i]*y2[i] + x51[i]*y1[i]
            + x52[i]*y0[i];
      *s53 += x0[i]*y53[i] + x1[i]*y52[i] + x2[i]*y51[i] + x3[i]*y50[i]
            + x4[i]*y49[i] + x5[i]*y48[i] + x6[i]*y47[i] + x7[i]*y46[i]
            + x8[i]*y45[i] + x9[i]*y44[i] + x10[i]*y43[i] + x11[i]*y42[i]
            + x12[i]*y41[i] + x13[i]*y40[i] + x14[i]*y39[i] + x15[i]*y38[i]
            + x16[i]*y37[i] + x17[i]*y36[i] + x18[i]*y35[i] + x19[i]*y34[i]
            + x20[i]*y33[i] + x21[i]*y32[i] + x22[i]*y31[i] + x23[i]*y30[i]
            + x24[i]*y29[i] + x25[i]*y28[i] + x26[i]*y27[i] + x27[i]*y26[i]
            + x28[i]*y25[i] + x29[i]*y24[i] + x30[i]*y23[i] + x31[i]*y22[i]
            + x32[i]*y21[i] + x33[i]*y20[i] + x34[i]*y19[i] + x35[i]*y18[i]
            + x36[i]*y17[i] + x37[i]*y16[i] + x38[i]*y15[i] + x39[i]*y14[i]
            + x40[i]*y13[i] + x41[i]*y12[i] + x42[i]*y11[i] + x43[i]*y10[i]
            + x44[i]*y9[i] + x45[i]*y8[i] + x46[i]*y7[i] + x47[i]*y6[i]
            + x48[i]*y5[i] + x49[i]*y4[i] + x50[i]*y3[i] + x51[i]*y2[i]
            + x52[i]*y1[i] + x53[i]*y0[i];
      *s54 += x0[i]*y54[i] + x1[i]*y53[i] + x2[i]*y52[i] + x3[i]*y51[i]
            + x4[i]*y50[i] + x5[i]*y49[i] + x6[i]*y48[i] + x7[i]*y47[i]
            + x8[i]*y46[i] + x9[i]*y45[i] + x10[i]*y44[i] + x11[i]*y43[i]
            + x12[i]*y42[i] + x13[i]*y41[i] + x14[i]*y40[i] + x15[i]*y39[i]
            + x16[i]*y38[i] + x17[i]*y37[i] + x18[i]*y36[i] + x19[i]*y35[i]
            + x20[i]*y34[i] + x21[i]*y33[i] + x22[i]*y32[i] + x23[i]*y31[i]
            + x24[i]*y30[i] + x25[i]*y29[i] + x26[i]*y28[i] + x27[i]*y27[i]
            + x28[i]*y26[i] + x29[i]*y25[i] + x30[i]*y24[i] + x31[i]*y23[i]
            + x32[i]*y22[i] + x33[i]*y21[i] + x34[i]*y20[i] + x35[i]*y19[i]
            + x36[i]*y18[i] + x37[i]*y17[i] + x38[i]*y16[i] + x39[i]*y15[i]
            + x40[i]*y14[i] + x41[i]*y13[i] + x42[i]*y12[i] + x43[i]*y11[i]
            + x44[i]*y10[i] + x45[i]*y9[i] + x46[i]*y8[i] + x47[i]*y7[i]
            + x48[i]*y6[i] + x49[i]*y5[i] + x50[i]*y4[i] + x51[i]*y3[i]
            + x52[i]*y2[i] + x53[i]*y1[i] + x54[i]*y0[i];
      *s55 += x0[i]*y55[i] + x1[i]*y54[i] + x2[i]*y53[i] + x3[i]*y52[i]
            + x4[i]*y51[i] + x5[i]*y50[i] + x6[i]*y49[i] + x7[i]*y48[i]
            + x8[i]*y47[i] + x9[i]*y46[i] + x10[i]*y45[i] + x11[i]*y44[i]
            + x12[i]*y43[i] + x13[i]*y42[i] + x14[i]*y41[i] + x15[i]*y40[i]
            + x16[i]*y39[i] + x17[i]*y38[i] + x18[i]*y37[i] + x19[i]*y36[i]
            + x20[i]*y35[i] + x21[i]*y34[i] + x22[i]*y33[i] + x23[i]*y32[i]
            + x24[i]*y31[i] + x25[i]*y30[i] + x26[i]*y29[i] + x27[i]*y28[i]
            + x28[i]*y27[i] + x29[i]*y26[i] + x30[i]*y25[i] + x31[i]*y24[i]
            + x32[i]*y23[i] + x33[i]*y22[i] + x34[i]*y21[i] + x35[i]*y20[i]
            + x36[i]*y19[i] + x37[i]*y18[i] + x38[i]*y17[i] + x39[i]*y16[i]
            + x40[i]*y15[i] + x41[i]*y14[i] + x42[i]*y13[i] + x43[i]*y12[i]
            + x44[i]*y11[i] + x45[i]*y10[i] + x46[i]*y9[i] + x47[i]*y8[i]
            + x48[i]*y7[i] + x49[i]*y6[i] + x50[i]*y5[i] + x51[i]*y4[i]
            + x52[i]*y3[i] + x53[i]*y2[i] + x54[i]*y1[i] + x55[i]*y0[i];
      *s56 += x0[i]*y56[i] + x1[i]*y55[i] + x2[i]*y54[i] + x3[i]*y53[i]
            + x4[i]*y52[i] + x5[i]*y51[i] + x6[i]*y50[i] + x7[i]*y49[i]
            + x8[i]*y48[i] + x9[i]*y47[i] + x10[i]*y46[i] + x11[i]*y45[i]
            + x12[i]*y44[i] + x13[i]*y43[i] + x14[i]*y42[i] + x15[i]*y41[i]
            + x16[i]*y40[i] + x17[i]*y39[i] + x18[i]*y38[i] + x19[i]*y37[i]
            + x20[i]*y36[i] + x21[i]*y35[i] + x22[i]*y34[i] + x23[i]*y33[i]
            + x24[i]*y32[i] + x25[i]*y31[i] + x26[i]*y30[i] + x27[i]*y29[i]
            + x28[i]*y28[i] + x29[i]*y27[i] + x30[i]*y26[i] + x31[i]*y25[i]
            + x32[i]*y24[i] + x33[i]*y23[i] + x34[i]*y22[i] + x35[i]*y21[i]
            + x36[i]*y20[i] + x37[i]*y19[i] + x38[i]*y18[i] + x39[i]*y17[i]
            + x40[i]*y16[i] + x41[i]*y15[i] + x42[i]*y14[i] + x43[i]*y13[i]
            + x44[i]*y12[i] + x45[i]*y11[i] + x46[i]*y10[i] + x47[i]*y9[i]
            + x48[i]*y8[i] + x49[i]*y7[i] + x50[i]*y6[i] + x51[i]*y5[i]
            + x52[i]*y4[i] + x53[i]*y3[i] + x54[i]*y2[i] + x55[i]*y1[i]
            + x56[i]*y0[i];
      *s57 += x0[i]*y57[i] + x1[i]*y56[i] + x2[i]*y55[i] + x3[i]*y54[i]
            + x4[i]*y53[i] + x5[i]*y52[i] + x6[i]*y51[i] + x7[i]*y50[i]
            + x8[i]*y49[i] + x9[i]*y48[i] + x10[i]*y47[i] + x11[i]*y46[i]
            + x12[i]*y45[i] + x13[i]*y44[i] + x14[i]*y43[i] + x15[i]*y42[i]
            + x16[i]*y41[i] + x17[i]*y40[i] + x18[i]*y39[i] + x19[i]*y38[i]
            + x20[i]*y37[i] + x21[i]*y36[i] + x22[i]*y35[i] + x23[i]*y34[i]
            + x24[i]*y33[i] + x25[i]*y32[i] + x26[i]*y31[i] + x27[i]*y30[i]
            + x28[i]*y29[i] + x29[i]*y28[i] + x30[i]*y27[i] + x31[i]*y26[i]
            + x32[i]*y25[i] + x33[i]*y24[i] + x34[i]*y23[i] + x35[i]*y22[i]
            + x36[i]*y21[i] + x37[i]*y20[i] + x38[i]*y19[i] + x39[i]*y18[i]
            + x40[i]*y17[i] + x41[i]*y16[i] + x42[i]*y15[i] + x43[i]*y14[i]
            + x44[i]*y13[i] + x45[i]*y12[i] + x46[i]*y11[i] + x47[i]*y10[i]
            + x48[i]*y9[i] + x49[i]*y8[i] + x50[i]*y7[i] + x51[i]*y6[i]
            + x52[i]*y5[i] + x53[i]*y4[i] + x54[i]*y3[i] + x55[i]*y2[i]
            + x56[i]*y1[i] + x57[i]*y0[i];
      *s58 += x0[i]*y58[i] + x1[i]*y57[i] + x2[i]*y56[i] + x3[i]*y55[i]
            + x4[i]*y54[i] + x5[i]*y53[i] + x6[i]*y52[i] + x7[i]*y51[i]
            + x8[i]*y50[i] + x9[i]*y49[i] + x10[i]*y48[i] + x11[i]*y47[i]
            + x12[i]*y46[i] + x13[i]*y45[i] + x14[i]*y44[i] + x15[i]*y43[i]
            + x16[i]*y42[i] + x17[i]*y41[i] + x18[i]*y40[i] + x19[i]*y39[i]
            + x20[i]*y38[i] + x21[i]*y37[i] + x22[i]*y36[i] + x23[i]*y35[i]
            + x24[i]*y34[i] + x25[i]*y33[i] + x26[i]*y32[i] + x27[i]*y31[i]
            + x28[i]*y30[i] + x29[i]*y29[i] + x30[i]*y28[i] + x31[i]*y27[i]
            + x32[i]*y26[i] + x33[i]*y25[i] + x34[i]*y24[i] + x35[i]*y23[i]
            + x36[i]*y22[i] + x37[i]*y21[i] + x38[i]*y20[i] + x39[i]*y19[i]
            + x40[i]*y18[i] + x41[i]*y17[i] + x42[i]*y16[i] + x43[i]*y15[i]
            + x44[i]*y14[i] + x45[i]*y13[i] + x46[i]*y12[i] + x47[i]*y11[i]
            + x48[i]*y10[i] + x49[i]*y9[i] + x50[i]*y8[i] + x51[i]*y7[i]
            + x52[i]*y6[i] + x53[i]*y5[i] + x54[i]*y4[i] + x55[i]*y3[i]
            + x56[i]*y2[i] + x57[i]*y1[i] + x58[i]*y0[i];
      *s59 += x0[i]*y59[i] + x1[i]*y58[i] + x2[i]*y57[i] + x3[i]*y56[i]
            + x4[i]*y55[i] + x5[i]*y54[i] + x6[i]*y53[i] + x7[i]*y52[i]
            + x8[i]*y51[i] + x9[i]*y50[i] + x10[i]*y49[i] + x11[i]*y48[i]
            + x12[i]*y47[i] + x13[i]*y46[i] + x14[i]*y45[i] + x15[i]*y44[i]
            + x16[i]*y43[i] + x17[i]*y42[i] + x18[i]*y41[i] + x19[i]*y40[i]
            + x20[i]*y39[i] + x21[i]*y38[i] + x22[i]*y37[i] + x23[i]*y36[i]
            + x24[i]*y35[i] + x25[i]*y34[i] + x26[i]*y33[i] + x27[i]*y32[i]
            + x28[i]*y31[i] + x29[i]*y30[i] + x30[i]*y29[i] + x31[i]*y28[i]
            + x32[i]*y27[i] + x33[i]*y26[i] + x34[i]*y25[i] + x35[i]*y24[i]
            + x36[i]*y23[i] + x37[i]*y22[i] + x38[i]*y21[i] + x39[i]*y20[i]
            + x40[i]*y19[i] + x41[i]*y18[i] + x42[i]*y17[i] + x43[i]*y16[i]
            + x44[i]*y15[i] + x45[i]*y14[i] + x46[i]*y13[i] + x47[i]*y12[i]
            + x48[i]*y11[i] + x49[i]*y10[i] + x50[i]*y9[i] + x51[i]*y8[i]
            + x52[i]*y7[i] + x53[i]*y6[i] + x54[i]*y5[i] + x55[i]*y4[i]
            + x56[i]*y3[i] + x57[i]*y2[i] + x58[i]*y1[i] + x59[i]*y0[i];
      *s60 += x0[i]*y60[i] + x1[i]*y59[i] + x2[i]*y58[i] + x3[i]*y57[i]
            + x4[i]*y56[i] + x5[i]*y55[i] + x6[i]*y54[i] + x7[i]*y53[i]
            + x8[i]*y52[i] + x9[i]*y51[i] + x10[i]*y50[i] + x11[i]*y49[i]
            + x12[i]*y48[i] + x13[i]*y47[i] + x14[i]*y46[i] + x15[i]*y45[i]
            + x16[i]*y44[i] + x17[i]*y43[i] + x18[i]*y42[i] + x19[i]*y41[i]
            + x20[i]*y40[i] + x21[i]*y39[i] + x22[i]*y38[i] + x23[i]*y37[i]
            + x24[i]*y36[i] + x25[i]*y35[i] + x26[i]*y34[i] + x27[i]*y33[i]
            + x28[i]*y32[i] + x29[i]*y31[i] + x30[i]*y30[i] + x31[i]*y29[i]
            + x32[i]*y28[i] + x33[i]*y27[i] + x34[i]*y26[i] + x35[i]*y25[i]
            + x36[i]*y24[i] + x37[i]*y23[i] + x38[i]*y22[i] + x39[i]*y21[i]
            + x40[i]*y20[i] + x41[i]*y19[i] + x42[i]*y18[i] + x43[i]*y17[i]
            + x44[i]*y16[i] + x45[i]*y15[i] + x46[i]*y14[i] + x47[i]*y13[i]
            + x48[i]*y12[i] + x49[i]*y11[i] + x50[i]*y10[i] + x51[i]*y9[i]
            + x52[i]*y8[i] + x53[i]*y7[i] + x54[i]*y6[i] + x55[i]*y5[i]
            + x56[i]*y4[i] + x57[i]*y3[i] + x58[i]*y2[i] + x59[i]*y1[i]
            + x60[i]*y0[i];
      *s61 += x0[i]*y61[i] + x1[i]*y60[i] + x2[i]*y59[i] + x3[i]*y58[i]
            + x4[i]*y57[i] + x5[i]*y56[i] + x6[i]*y55[i] + x7[i]*y54[i]
            + x8[i]*y53[i] + x9[i]*y52[i] + x10[i]*y51[i] + x11[i]*y50[i]
            + x12[i]*y49[i] + x13[i]*y48[i] + x14[i]*y47[i] + x15[i]*y46[i]
            + x16[i]*y45[i] + x17[i]*y44[i] + x18[i]*y43[i] + x19[i]*y42[i]
            + x20[i]*y41[i] + x21[i]*y40[i] + x22[i]*y39[i] + x23[i]*y38[i]
            + x24[i]*y37[i] + x25[i]*y36[i] + x26[i]*y35[i] + x27[i]*y34[i]
            + x28[i]*y33[i] + x29[i]*y32[i] + x30[i]*y31[i] + x31[i]*y30[i]
            + x32[i]*y29[i] + x33[i]*y28[i] + x34[i]*y27[i] + x35[i]*y26[i]
            + x36[i]*y25[i] + x37[i]*y24[i] + x38[i]*y23[i] + x39[i]*y22[i]
            + x40[i]*y21[i] + x41[i]*y20[i] + x42[i]*y19[i] + x43[i]*y18[i]
            + x44[i]*y17[i] + x45[i]*y16[i] + x46[i]*y15[i] + x47[i]*y14[i]
            + x48[i]*y13[i] + x49[i]*y12[i] + x50[i]*y11[i] + x51[i]*y10[i]
            + x52[i]*y9[i] + x53[i]*y8[i] + x54[i]*y7[i] + x55[i]*y6[i]
            + x56[i]*y5[i] + x57[i]*y4[i] + x58[i]*y3[i] + x59[i]*y2[i]
            + x60[i]*y1[i] + x61[i]*y0[i];
      *s62 += x0[i]*y62[i] + x1[i]*y61[i] + x2[i]*y60[i] + x3[i]*y59[i]
            + x4[i]*y58[i] + x5[i]*y57[i] + x6[i]*y56[i] + x7[i]*y55[i]
            + x8[i]*y54[i] + x9[i]*y53[i] + x10[i]*y52[i] + x11[i]*y51[i]
            + x12[i]*y50[i] + x13[i]*y49[i] + x14[i]*y48[i] + x15[i]*y47[i]
            + x16[i]*y46[i] + x17[i]*y45[i] + x18[i]*y44[i] + x19[i]*y43[i]
            + x20[i]*y42[i] + x21[i]*y41[i] + x22[i]*y40[i] + x23[i]*y39[i]
            + x24[i]*y38[i] + x25[i]*y37[i] + x26[i]*y36[i] + x27[i]*y35[i]
            + x28[i]*y34[i] + x29[i]*y33[i] + x30[i]*y32[i] + x31[i]*y31[i]
            + x32[i]*y30[i] + x33[i]*y29[i] + x34[i]*y28[i] + x35[i]*y27[i]
            + x36[i]*y26[i] + x37[i]*y25[i] + x38[i]*y24[i] + x39[i]*y23[i]
            + x40[i]*y22[i] + x41[i]*y21[i] + x42[i]*y20[i] + x43[i]*y19[i]
            + x44[i]*y18[i] + x45[i]*y17[i] + x46[i]*y16[i] + x47[i]*y15[i]
            + x48[i]*y14[i] + x49[i]*y13[i] + x50[i]*y12[i] + x51[i]*y11[i]
            + x52[i]*y10[i] + x53[i]*y9[i] + x54[i]*y8[i] + x55[i]*y7[i]
            + x56[i]*y6[i] + x57[i]*y5[i] + x58[i]*y4[i] + x59[i]*y3[i]
            + x60[i]*y2[i] + x61[i]*y1[i] + x62[i]*y0[i];
      *s63 += x0[i]*y63[i] + x1[i]*y62[i] + x2[i]*y61[i] + x3[i]*y60[i]
            + x4[i]*y59[i] + x5[i]*y58[i] + x6[i]*y57[i] + x7[i]*y56[i]
            + x8[i]*y55[i] + x9[i]*y54[i] + x10[i]*y53[i] + x11[i]*y52[i]
            + x12[i]*y51[i] + x13[i]*y50[i] + x14[i]*y49[i] + x15[i]*y48[i]
            + x16[i]*y47[i] + x17[i]*y46[i] + x18[i]*y45[i] + x19[i]*y44[i]
            + x20[i]*y43[i] + x21[i]*y42[i] + x22[i]*y41[i] + x23[i]*y40[i]
            + x24[i]*y39[i] + x25[i]*y38[i] + x26[i]*y37[i] + x27[i]*y36[i]
            + x28[i]*y35[i] + x29[i]*y34[i] + x30[i]*y33[i] + x31[i]*y32[i]
            + x32[i]*y31[i] + x33[i]*y30[i] + x34[i]*y29[i] + x35[i]*y28[i]
            + x36[i]*y27[i] + x37[i]*y26[i] + x38[i]*y25[i] + x39[i]*y24[i]
            + x40[i]*y23[i] + x41[i]*y22[i] + x42[i]*y21[i] + x43[i]*y20[i]
            + x44[i]*y19[i] + x45[i]*y18[i] + x46[i]*y17[i] + x47[i]*y16[i]
            + x48[i]*y15[i] + x49[i]*y14[i] + x50[i]*y13[i] + x51[i]*y12[i]
            + x52[i]*y11[i] + x53[i]*y10[i] + x54[i]*y9[i] + x55[i]*y8[i]
            + x56[i]*y7[i] + x57[i]*y6[i] + x58[i]*y5[i] + x59[i]*y4[i]
            + x60[i]*y3[i] + x61[i]*y2[i] + x62[i]*y1[i] + x63[i]*y0[i];
   }
}
