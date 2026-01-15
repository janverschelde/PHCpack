// The file random16_vectors.h specifies functions
// to generate random vectors in hexa double precision.

#ifndef __random16_vectors_h__
#define __random16_vectors_h__

void random_hexa_double
 ( double *x_hihihihi, double *x_lohihihi,
   double *x_hilohihi, double *x_lolohihi,
   double *x_hihilohi, double *x_lohilohi,
   double *x_hilolohi, double *x_lololohi,
   double *x_hihihilo, double *x_lohihilo,
   double *x_hilohilo, double *x_lolohilo,
   double *x_hihilolo, double *x_lohilolo,
   double *x_hilololo, double *x_lolololo );
/*
 * DESCRIPTION :
 *   Returns a random hexa double x in [-1, +1],
 *   with the generation of 16 random doubles
 *   so all 16 parts of the random hexa double are filled.
 *
 * ON RETURN :
 *   x_hihihihi   highest part of the random hexa double x;
 *   x_lohihihi   second highest part of the random hexa double x;
 *   x_hilohihi   third highest part of the random hexa double x;
 *   x_lolohihi   fourth highest part of the random hexa double x;
 *   x_hihilohi   fifth highest part of the random hexa double x;
 *   x_lohilohi   sixth highest part of the random hexa double x;
 *   x_hilolohi   seventh highest part of the random hexa double x;
 *   x_lololohi   eighth highest part of the random hexa double x;
 *   x_hihihilo   eighth lowest part of the random hexa double x;
 *   x_lohihilo   seventh lowest part of the random hexa double x;
 *   x_hilohilo   sixth lowest part of the random hexa double x;
 *   x_lolohilo   fifth lowest part of the random hexa double x;
 *   x_hihilolo   fourth lowest part of the random hexa double x;
 *   x_lohilolo   third lowest part of the random hexa double x;
 *   x_hilololo   second lowest part of the random hexa double x;
 *   x_lolololo   lowest part of the random hexa double x. */

#endif
