// The file random16_vectors.cpp defines the code for the functions
// specified in random16_vectors.h.

#include "random_numbers.h"
#include "random16_vectors.h"
#include "hexa_double_functions.h"

void random_hexa_double
 ( double *x_hihihihi, double *x_lohihihi,
   double *x_hilohihi, double *x_lolohihi,
   double *x_hihilohi, double *x_lohilohi,
   double *x_hilolohi, double *x_lololohi,
   double *x_hihihilo, double *x_lohihilo,
   double *x_hilohilo, double *x_lolohilo,
   double *x_hihilolo, double *x_lohilolo,
   double *x_hilololo, double *x_lolololo )
{
   const double eps = 2.220446049250313e-16; // 2^(-52)
   const double eps2 = eps*eps;
   const double eps3 = eps*eps2;
   const double eps4 = eps*eps3;
   const double eps5 = eps*eps4;
   const double eps6 = eps*eps5;
   const double eps7 = eps*eps6;
   const double eps8 = eps*eps7;
   const double eps9 = eps*eps8;
   const double eps10 = eps*eps9;
   const double eps11 = eps*eps10;
   const double eps12 = eps*eps11;
   const double eps13 = eps*eps12;
   const double eps14 = eps*eps13;
   const double eps15 = eps*eps14;
   const double r0 = random_double();
   const double r1 = random_double();
   const double r2 = random_double();
   const double r3 = random_double();
   const double r4 = random_double();
   const double r5 = random_double();
   const double r6 = random_double();
   const double r7 = random_double();
   const double r8 = random_double();
   const double r9 = random_double();
   const double r10 = random_double();
   const double r11 = random_double();
   const double r12 = random_double();
   const double r13 = random_double();
   const double r14 = random_double();
   const double r15 = random_double();

   *x_hihihihi = r0;  *x_lohihihi = 0.0;
   *x_hilohihi = 0.0; *x_lolohihi = 0.0;
   *x_hihilohi = 0.0; *x_lohilohi = 0.0;
   *x_hilolohi = 0.0; *x_lololohi = 0.0;
   *x_hihihilo = 0.0; *x_lohihilo = 0.0;
   *x_hilohilo = 0.0; *x_lolohilo = 0.0;
   *x_hihilolo = 0.0; *x_lohilolo = 0.0;
   *x_hilololo = 0.0; *x_lolololo = 0.0;

   hdf_inc_d(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
             x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
             x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
             x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,r1*eps);
   hdf_inc_d(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
             x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
             x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
             x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,r2*eps2);
   hdf_inc_d(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
             x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
             x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
             x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,r3*eps3);
   hdf_inc_d(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
             x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
             x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
             x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,r4*eps4);
   hdf_inc_d(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
             x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
             x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
             x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,r5*eps5);
   hdf_inc_d(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
             x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
             x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
             x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,r6*eps6);
   hdf_inc_d(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
             x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
             x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
             x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,r7*eps7);
   hdf_inc_d(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
             x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
             x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
             x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,r8*eps8);
   hdf_inc_d(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
             x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
             x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
             x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,r9*eps9);
   hdf_inc_d(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
             x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
             x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
             x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,r10*eps10);
   hdf_inc_d(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
             x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
             x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
             x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,r11*eps11);
   hdf_inc_d(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
             x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
             x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
             x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,r12*eps12);
   hdf_inc_d(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
             x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
             x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
             x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,r13*eps13);
   hdf_inc_d(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
             x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
             x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
             x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,r14*eps14);
   hdf_inc_d(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
             x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
             x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
             x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,r15*eps15);
}
