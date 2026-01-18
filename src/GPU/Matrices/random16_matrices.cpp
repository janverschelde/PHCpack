// The file random16_matrices.cpp defines the functions specified in
// the file random16_matrices.h.

#include "random16_vectors.h"
#include "random16_matrices.h"

void random_dbl16_matrix
 ( int rows, int cols,
   double **Ahihihihi, double **Alohihihi,
   double **Ahilohihi, double **Alolohihi,
   double **Ahihihilo, double **Alohihilo,
   double **Ahilohilo, double **Alolohilo,
   double **Ahihilohi, double **Alohilohi,
   double **Ahilolohi, double **Alololohi,
   double **Ahihilolo, double **Alohilolo,
   double **Ahilololo, double **Alolololo )
{
   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<cols; j++)
         random_hexa_double
            (&Ahihihihi[i][j],&Alohihihi[i][j],
             &Ahilohihi[i][j],&Alolohihi[i][j],
             &Ahihilohi[i][j],&Alohilohi[i][j],
             &Ahilolohi[i][j],&Alololohi[i][j],
             &Ahihihilo[i][j],&Alohihilo[i][j],
             &Ahilohilo[i][j],&Alolohilo[i][j],
             &Ahihilolo[i][j],&Alohilolo[i][j],
             &Ahilololo[i][j],&Alolololo[i][j]);
   }
}
