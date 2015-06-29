// defines code for execution on the host

#include "mgs_host.h"

void CPU_normalize_and_reduce
 ( CT** v, int rows, int cols, int pivot )
{
   T1 sum(0.0);

   for(int k=0; k<rows; k++)                            // compute 2-norm
      sum = sum + v[pivot][k].real*v[pivot][k].real     // of pivot column
                + v[pivot][k].imag*v[pivot][k].imag;
   sum = sqrt(sum);
   for(int k=0; k<rows; k++)                            // normalize the
      v[pivot][k] = v[pivot][k]/sum;                    // pivot column

   for(int k=pivot+1; k<cols; k++) // reduce k-th column w.r.t. pivot column
   {
      CT inprod(0.0,0.0); // inner product of column k with pivot
      for(int i=0; i<rows; i++)
         inprod = inprod + v[pivot][i].adj()*v[k][i];
      for(int i=0; i<rows; i++)
         v[k][i] = v[k][i] - inprod*v[pivot][i];        // reduction
   }
}

void CPU_QR_normalize_and_reduce
 ( CT** v, CT** R, int rows, int cols, int pivot )
{
   //std::cout << "--------------------------" << std::endl;
   //std::cout << "rows = " << rows << ", cols = " << cols << ", pivot = " << pivot << std::endl;
   T1 sum(0);
   T1 zero(0);

   for(int k=0; k<rows; k++)                            // compute 2-norm
      sum = sum + v[pivot][k].real*v[pivot][k].real     // of pivot column
                + v[pivot][k].imag*v[pivot][k].imag;

   sum = sqrt(sum);
   R[pivot][pivot].real = sum;
   R[pivot][pivot].imag = zero;
   for(int k=0; k<rows; k++)                            // normalize the
      v[pivot][k] = v[pivot][k]/sum;                    // pivot column

   for(int k=pivot+1; k<cols; k++) // reduce k-th column w.r.t. pivot column
   {
      CT inprod(0.0,0.0); // inner product of column k with pivot
      for(int i=0; i<rows; i++)
         inprod = inprod + v[pivot][i].adj()*v[k][i];
      R[pivot][k] = inprod;
      for(int i=0; i<rows; i++)
         v[k][i] = v[k][i] - inprod*v[pivot][i];        // reduction
   }
}

void CPU_backsubstitution
 ( CT** U, CT* rhs, CT* x, int dim )
{
   CT temp;

   for(int i=dim-1; i>=0; i--)
   {
      temp = rhs[i];
      for(int j=i+1; j<dim; j++) temp = temp - U[i][j]*x[j];
      x[i] = temp/U[i][i];
   }
}

void CPU_mgs2 ( CT** v, int rows, int cols )
{
   for(int piv=0; piv<cols; piv++)
      CPU_normalize_and_reduce(v,rows,cols,piv);
}

void CPU_mgs2qr ( CT** v, CT** R, int rows, int cols )
{
   for(int piv=0; piv<cols; piv++)
      CPU_QR_normalize_and_reduce(v,R,rows,cols,piv);
}

void CPU_mgs2qrls
 ( CT** v, CT** R, CT* x, int rows, int cols )
{
	//std::cout << "rows = " << rows << " cols = " << cols << std::endl;
   for(int piv=0; piv<cols-1; piv++) // one extra column with rhs vector
      CPU_QR_normalize_and_reduce(v,R,rows,cols,piv);

   //std::cout << "abc" << std::endl;

   //for(int i=0; i<cols; i++){
   //   for(int j=0; j<cols; j++){
	//	   std::cout << i << " "  << j << " "<< R[i][j];
	//   }
   //}

   CT* rhs = new CT[cols];
   for(int i=0; i<cols; i++) rhs[i] = R[i][cols-1];

   //std::cout << "abcd" << std::endl;

   //for(int i=0; i<cols; i++)
   //   std::cout << i << " " << rhs[i];

   // cout << "the triangular matrix R :" << endl;
   // print_matrix(R,cols,cols);
   CPU_backsubstitution(R,rhs,x,cols-1);
}
