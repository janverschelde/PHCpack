/* Defines the functions with prototypes in vectored_hexa_doubles.h. */

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
