// defines code for execution on the host

#include <cmath>
#include "mgs2_host.h"
#include "gqd_qd_utilT.h"

void print_vector ( complex<T> *v, int dim )
{
   complexH<T1> temp;

   for(int i=0; i<dim; i++)
   {  
      comp1_gqd2qd(&v[i],&temp);
      cout << "v[" << i << "] = " << temp;
   }
}

void print_vector ( complexH<T1> *v, int dim )
{
   for(int i=0; i<dim; i++)
      cout << "v[" << i << "] = " << v[i];
}

void print_matrix ( complexH<T1> **v, int rows, int cols )
{
   for(int i=0; i<cols; i++)
      for(int j=0; j<rows; j++)
         cout << "v[" << i << "," << j << "] = " << v[i][j];
}

void print_matrices ( complexH<T1> **v, complex<T> *v_h, int rows, int cols )
{
   for(int i=0; i<cols; i++)
      for(int j=0; j<rows; j++)
      {
         complexH<T1> temp;
         comp1_gqd2qd(&v_h[i*rows+j],&temp);
         cout << "GPU v[" << i << "," << j << "] = " << temp;
         cout << "CPU v[" << i << "," << j << "] = " << v[i][j];
      }
}

void print_difference ( complexH<T1> **v, complex<T> *v_h, int rows, int cols)
{
   complexH<T1> sum(0.0,0.0);

   for(int i=0; i<cols; i++)
      for(int j=0; j<rows; j++)
      {
         complexH<T1> temp;
         comp1_gqd2qd(&v_h[i*rows+j],&temp);
         temp = temp - v[i][j];
         sum = sum + temp.adj()*temp;
      }
   cout << "2-norm of componentwise difference : " << sqrt(sum.real) << endl;
}

void copy_matrices
 ( complexH<T1>** vfrom, complexH<T1>** vto, int rows, int cols )
{
   for(int i=0; i<cols; i++)
      for(int j=0; j<rows; j++) vto[i][j] = vfrom[i][j];
}

void checkGPUnormal ( complex<T>* v_h, int rows, int cols )
{
   T1 error(0.0);
   T1 one(1.0);
   complexH<T1> temp;

   for(int i=0; i<cols; i++)
   {
      T1 sum(0.0);
      for(int j=0; j<rows; j++)
      {
         comp1_gqd2qd(&v_h[i*rows+j],&temp);
         sum = sum + temp.real*temp.real + temp.imag*temp.imag;
      }
      // cout << "column " << i << " has norm " << sqrt(sum.real) << endl;
      error = error + abs(sum - one);
   }
   cout << "GPU normality error : " << error << endl;
}

void checkCPUnormal ( complexH<T1>** v, int rows, int cols )
{
   T1 error(0.0);
   T1 one(1.0);
   for(int i=0; i<cols; i++)
   {
      T1 sum(0.0);
      for(int j=0; j<rows; j++)
         sum = sum + v[i][j].real*v[i][j].real + v[i][j].imag*v[i][j].imag;
      // cout << "column " << i << " has norm " << sqrt(sum.real) << endl;
      error = error + abs(sum - one);
   }
   cout << "CPU normality error : " << error << endl;
}

void checkGPUorthogonal ( complex<T>* v_h, int rows, int cols )
{
   complexH<T1> error(0.0,0.0);
   complexH<T1> xtemp, ytemp;
   for(int i=0; i<cols; i++)
   {
      for(int j=i+1; j<cols; j++)
      {
         complexH<T1> sum(0.0,0.0);
         for(int k=0; k<rows; k++)
         {
            comp1_gqd2qd(&v_h[i*rows+k],&xtemp);
            comp1_gqd2qd(&v_h[j*rows+k],&ytemp);
            sum = sum + xtemp.adj()*ytemp;
         }
         if(abs(sum.real) > 1.0e-12)
            cout << "< " << i << " , " << j << " > = " << sum;
         sum.real = abs(sum.real);
         sum.imag = abs(sum.imag);
         error = error + sum;
      }
   }
   cout << "GPU orthogonality error : " << error;
}

void checkCPUorthogonal ( complexH<T1>** v, int rows, int cols )
{
   complexH<T1> error(0.0,0.0);
   for(int i=0; i<cols; i++)
   {
      for(int j=i+1; j<cols; j++)
      {
         complexH<T1> sum(0.0,0.0);
         for(int k=0; k<rows; k++)
            sum = sum + v[i][k].adj()*v[j][k];
         // cout << "< " << i << " , " << j << " > = " << sum;
         sum.real = abs(sum.real);
         sum.imag = abs(sum.imag);
         error = error + sum;
      }
   }
   cout << "CPU orthogonality error : " << error;
}

void checkGPUdecomposition
 ( complexH<T1>** A, complex<T>* Q, complex<T>* R,
   int dimR, int rows, int cols )
{
   // first we copy the coalesced stored form of R into RR
   complexH<T1>** RR = new complexH<T1>*[cols];
   for(int i=0; i<cols; i++) RR[i] = new complexH<T1>[cols];
   for(int pivot=0; pivot<cols; pivot++)
   {
      for(int block=0; block<pivot; block++)
         RR[pivot][block].init(0.0,0.0);
      for(int block=0; block<cols-pivot; block++)
      {
         int indR = (dimR-1) - (pivot*(pivot+1))/2
                             - (block*(block+1))/2
                             - block*(pivot+1);
         // cout << "pivot = " << pivot << "  block = " << block;
         // cout << "  indR = " << indR << endl;
         comp1_gqd2qd(&R[indR],&RR[pivot][block+pivot]);
      }
   }
   // cout << "the matrix RR :" << endl; print_matrix(RR,cols,cols);

   T1 error(0.0);
   for(int i=0; i<rows; i++)     // multiply i-th row of Q
   {
      for(int j=0; j<cols; j++)  // with j-th column of R
      {
         complexH<T1> sum(0.0,0.0);
         complexH<T1> temp;

         for(int k=0; k<cols; k++)
         {
            comp1_gqd2qd(&Q[k*rows+i],&temp);
            sum = sum + temp*RR[k][j];
         }
         temp = A[j][i] - sum;   // A is stored column wise
         // cout << i << "," << j << " : " << temp;
         error = error + abs(temp.real) + abs(temp.imag);
      }
   }
   cout << "GPU decomposition error : " << error << endl;
}

void checkCPUdecomposition
 ( complexH<T1>** A, complexH<T1>** Q, complexH<T1>** R, int rows, int cols )
{
   // cout << "the matrix R :" << endl; print_matrix(R,cols,cols);

   T1 error(0.0);

   for(int i=0; i<rows; i++)     // multiply i-th row of Q
   {
      for(int j=0; j<cols; j++)  // with j-th column of R
      {
         complexH<T1> sum(0.0,0.0);
         complexH<T1> temp;

         for(int k=0; k<cols; k++)
            sum = sum + Q[k][i]*R[k][j];
         temp = A[j][i] - sum;   // A is stored column wise
         // cout << i << "," << j << " : " << temp;
         error = error + abs(temp.real) + abs(temp.imag);
      }
   }
   cout << "CPU decomposition error : " << error << endl;
}

void checkGPUsolution
 ( complexH<T1>** A, complex<T>* x, int rows, int cols )
{
   T1 error(0.0);
   complexH<T1> sum;

   for(int i=0; i<rows; i++)
   {
      sum = A[cols-1][i];
      for(int j=0; j<cols-1; j++)
      {
         complexH<T1> temp;

         comp1_gqd2qd(&x[j],&temp);
         sum = sum - A[j][i]*temp;
      }
      error = error + abs(sum.real) + abs(sum.imag);
      // cout << "error at row " << i << " : " << error << endl;
   }
   cout << "GPU solution error : " << error << endl;
}

void checkCPUsolution
 ( complexH<T1>** A, complexH<T1>* x, int rows, int cols )
{
   T1 error(0.0);
   complexH<T1> sum;

   for(int i=0; i<rows; i++)
   {
      sum = A[cols-1][i];
      for(int j=0; j<cols-1; j++)
         sum = sum - A[j][i]*x[j];
      error = error + abs(sum.real) + abs(sum.imag);
      // cout << "error at row " << i << " : " << error << endl;
   }
   cout << "CPU solution error : " << error << endl;
}

void CPU_normalize_and_reduce
 ( complexH<T1>** v, int rows, int cols, int pivot )
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
      complexH<T1> inprod(0.0,0.0); // inner product of column k with pivot
      for(int i=0; i<rows; i++)
         inprod = inprod + v[pivot][i].adj()*v[k][i];
      for(int i=0; i<rows; i++)
         v[k][i] = v[k][i] - inprod*v[pivot][i];        // reduction
   }
}

void CPU_QR_normalize_and_reduce
 ( complexH<T1>** v, complexH<T1>** R, int rows, int cols, int pivot )
{
   T1 sum(0.0);
   T1 zero(0.0);

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
      complexH<T1> inprod(0.0,0.0); // inner product of column k with pivot
      for(int i=0; i<rows; i++)
         inprod = inprod + v[pivot][i].adj()*v[k][i];
      R[pivot][k] = inprod;
      for(int i=0; i<rows; i++)
         v[k][i] = v[k][i] - inprod*v[pivot][i];        // reduction
   }
}

void CPU_backsubstitution
 ( complexH<T1>** U, complexH<T1>* rhs, complexH<T1>* x, int dim )
{
   complexH<T1> temp;

   for(int i=dim-1; i>=0; i--)
   {
      temp = rhs[i];
      for(int j=i+1; j<dim; j++) temp = temp - U[i][j]*x[j];
      x[i] = temp/U[i][i];
   }
}

void CPU_mgs2 ( complexH<T1>** v, int rows, int cols )
{
   for(int piv=0; piv<cols; piv++)
      CPU_normalize_and_reduce(v,rows,cols,piv);
}

void CPU_mgs2qr ( complexH<T1>** v, complexH<T1>** R, int rows, int cols )
{
   for(int piv=0; piv<cols; piv++)
      CPU_QR_normalize_and_reduce(v,R,rows,cols,piv);
}

void CPU_mgs2qrls 
 ( complexH<T1>** v, complexH<T1>** R, complexH<T1>* x, int rows, int cols )
{
   for(int piv=0; piv<cols-1; piv++) // one extra column with rhs vector
      CPU_QR_normalize_and_reduce(v,R,rows,cols,piv);

   complexH<T1>* rhs = new complexH<T1>[rows];
   for(int i=0; i<rows; i++) rhs[i] = R[i][cols-1];

   // cout << "the triangular matrix R :" << endl;
   // print_matrix(R,cols,cols);
   CPU_backsubstitution(R,rhs,x,rows);
}
