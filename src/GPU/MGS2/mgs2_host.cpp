// defines code for execution on the host

#include <cmath>
#include <iomanip>
#include "mgs2_host.h"
#include "gqd_qd_utilT.h"

void print_complex_vector ( complex<T> *v, int dim )
{
   complexH<T1> temp;

   for(int i=0; i<dim; i++)
   {  
      comp1_gqd2qd(&v[i],&temp);
      cout << "v[" << i << "] = " << temp;
   }
}

void print_real_vector ( T *v, int dim )
{
   T1 temp;

   for(int i=0; i<dim; i++)
   {  
      gqd2qd(&v[i],&temp);
      cout << "v[" << i << "] = " << temp << endl;
   }
}

void print_complex_vector ( complexH<T1> *v, int dim )
{
   for(int i=0; i<dim; i++)
      cout << "v[" << i << "] = " << v[i];
}

void print_real_vector2 ( T1 *v, int dim )
{
   for(int i=0; i<dim; i++)
      cout << "v[" << i << "] = " << v[i] << endl;
}

void print_complex_matrix ( complexH<T1> **v, int rows, int cols )
{
   for(int i=0; i<cols; i++)
      for(int j=0; j<rows; j++)
         cout << "v[" << i << "," << j << "] = " << v[i][j];
}

void print_real_matrix ( T1 **v, int rows, int cols )
{
   for(int i=0; i<cols; i++)
      for(int j=0; j<rows; j++)
         cout << scientific << setprecision(16)
              << "v[" << i << "," << j << "] = " << v[i][j] << endl;
}

void print_complex_matrices
 ( complexH<T1> **v, complex<T> *v_h, int rows, int cols )
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

void print_real_matrices ( T1 **v, T *v_h, int rows, int cols )
{
   for(int i=0; i<cols; i++)
      for(int j=0; j<rows; j++)
      {
         T1 temp;
         gqd2qd(&v_h[i*rows+j],&temp);
         cout << "GPU v[" << i << "," << j << "] = " << temp << endl;
         cout << "CPU v[" << i << "," << j << "] = " << v[i][j] << endl;
      }
}

void complex_print_difference
 ( complexH<T1> **v, complex<T> *v_h, int rows, int cols)
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

void real_print_difference ( T1 **v, T *v_h, int rows, int cols)
{
   T1 sum(0.0);

   for(int i=0; i<cols; i++)
      for(int j=0; j<rows; j++)
      {
         T1 temp;
         gqd2qd(&v_h[i*rows+j],&temp);
         temp = temp - v[i][j];
         sum = sum + temp*temp;
      }
   cout << "2-norm of componentwise difference : " << sqrt(sum) << endl;
}

void copy_complex_matrices
 ( complexH<T1>** vfrom, complexH<T1>** vto, int rows, int cols )
{
   for(int i=0; i<cols; i++)
      for(int j=0; j<rows; j++) vto[i][j] = vfrom[i][j];
}

void copy_real_matrices ( T1** vfrom, T1** vto, int rows, int cols )
{
   for(int i=0; i<cols; i++)
      for(int j=0; j<rows; j++) vto[i][j] = vfrom[i][j];
}

void complex_checkGPUnormal ( complex<T>* v_h, int rows, int cols )
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

void real_checkGPUnormal ( T* v_h, int rows, int cols )
{
   T1 error(0.0);
   T1 one(1.0);
   T1 temp;

   for(int i=0; i<cols; i++)
   {
      T1 sum(0.0);
      for(int j=0; j<rows; j++)
      {
         gqd2qd(&v_h[i*rows+j],&temp);
         sum = sum + temp*temp;
      }
      // cout << "column " << i << " has norm " << sqrt(sum.real) << endl;
      error = error + abs(sum - one);
   }
   cout << "GPU normality error : " << error << endl;
}

void complex_checkCPUnormal ( complexH<T1>** v, int rows, int cols )
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

void real_checkCPUnormal ( T1** v, int rows, int cols )
{
   T1 error(0.0);
   T1 one(1.0);
   for(int i=0; i<cols; i++)
   {
      T1 sum(0.0);
      for(int j=0; j<rows; j++)
         sum = sum + v[i][j]*v[i][j];
      // cout << "column " << i << " has norm " << sqrt(sum) << endl;
      error = error + abs(sum - one);
   }
   cout << "CPU normality error : " << error << endl;
}

void complex_checkGPUorthogonal ( complex<T>* v_h, int rows, int cols )
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

void real_checkGPUorthogonal ( T* v_h, int rows, int cols )
{
   T1 error(0.0);
   T1 xtemp, ytemp;
   for(int i=0; i<cols; i++)
   {
      for(int j=i+1; j<cols; j++)
      {
         T1 sum(0.0);
         for(int k=0; k<rows; k++)
         {
            gqd2qd(&v_h[i*rows+k],&xtemp);
            gqd2qd(&v_h[j*rows+k],&ytemp);
            sum = sum + xtemp*ytemp;
         }
         if(abs(sum) > 1.0e-12)
            cout << "< " << i << " , " << j << " > = " << sum << endl;
         error = error + abs(sum);
      }
   }
   cout << "GPU orthogonality error : " << error << endl;
}

void complex_checkCPUorthogonal ( complexH<T1>** v, int rows, int cols )
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

void real_checkCPUorthogonal ( T1** v, int rows, int cols )
{
   T1 error(0.0);
   for(int i=0; i<cols; i++)
   {
      for(int j=i+1; j<cols; j++)
      {
         T1 sum(0.0);
         for(int k=0; k<rows; k++)
            sum = sum + v[i][k]*v[j][k];
         // cout << "< " << i << " , " << j << " > = " << sum;
         error = error + abs(sum);
      }
   }
   cout << "CPU orthogonality error : " << error << endl;
}

void complex_checkGPUdecomposition
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
   // cout << "the matrix RR :" << endl; print_complex_matrix(RR,cols,cols);

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

void real_checkGPUdecomposition
 ( T1** A, T* Q, T* R, int dimR, int rows, int cols )
{
   // first we copy the coalesced stored form of R into RR
   T1** RR = new T1*[cols];
   for(int i=0; i<cols; i++) RR[i] = new T1[cols];
   for(int pivot=0; pivot<cols; pivot++)
   {
      for(int block=0; block<pivot; block++)
         RR[pivot][block] = 0.0;
      for(int block=0; block<cols-pivot; block++)
      {
         int indR = (dimR-1) - (pivot*(pivot+1))/2
                             - (block*(block+1))/2
                             - block*(pivot+1);
         // cout << "pivot = " << pivot << "  block = " << block;
         // cout << "  indR = " << indR << endl;
         gqd2qd(&R[indR],&RR[pivot][block+pivot]);
      }
   }
   // cout << "the matrix RR :" << endl; print_complex_matrix(RR,cols,cols);

   T1 error(0.0);
   for(int i=0; i<rows; i++)     // multiply i-th row of Q
   {
      for(int j=0; j<cols; j++)  // with j-th column of R
      {
         T1 sum(0.0);
         T1 temp;

         for(int k=0; k<cols; k++)
         {
            gqd2qd(&Q[k*rows+i],&temp);
            sum = sum + temp*RR[k][j];
         }
         temp = A[j][i] - sum;   // A is stored column wise
         // cout << i << "," << j << " : " << temp;
         error = error + abs(temp);
      }
   }
   cout << "GPU decomposition error : " << error << endl;
}

void complex_checkCPUdecomposition
 ( complexH<T1>** A, complexH<T1>** Q, complexH<T1>** R, int rows, int cols )
{
   // cout << "the matrix R :" << endl; print_complex_matrix(R,cols,cols);

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

void real_checkCPUdecomposition
 ( T1** A, T1** Q, T1** R, int rows, int cols )
{
   // cout << "the matrix R :" << endl; print_complex_matrix(R,cols,cols);

   T1 error(0.0);

   for(int i=0; i<rows; i++)     // multiply i-th row of Q
   {
      for(int j=0; j<cols; j++)  // with j-th column of R
      {
         T1 sum(0.0);
         T1 temp;

         for(int k=0; k<cols; k++)
            sum = sum + Q[k][i]*R[k][j];
         temp = A[j][i] - sum;   // A is stored column wise
         // cout << i << "," << j << " : " << temp;
         error = error + abs(temp);
      }
   }
   cout << "CPU decomposition error : " << error << endl;
}

void complex_checkGPUsolution
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

void real_checkGPUsolution ( T1** A, T* x, int rows, int cols )
{
   T1 error(0.0);
   T1 sum;

   for(int i=0; i<rows; i++)
   {
      sum = A[cols-1][i];
      for(int j=0; j<cols-1; j++)
      {
         T1 temp;

         gqd2qd(&x[j],&temp);
         sum = sum - A[j][i]*temp;
      }
      error = error + abs(sum);
      // cout << "error at row " << i << " : " << error << endl;
   }
   cout << "GPU solution error : " << error << endl;
}

void complex_checkCPUsolution
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

void real_checkCPUsolution ( T1** A, T1* x, int rows, int cols )
{
   T1 error(0.0);
   T1 sum;

   for(int i=0; i<rows; i++)
   {
      sum = A[cols-1][i];
      for(int j=0; j<cols-1; j++)
         sum = sum - A[j][i]*x[j];
      error = error + abs(sum);
      // cout << "error at row " << i << " : " << error << endl;
   }
   cout << "CPU solution error : " << error << endl;
}

void CPU_complex_normalize_and_reduce
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

void CPU_real_normalize_and_reduce
 ( T1** v, int rows, int cols, int pivot )
{
   T1 sum(0.0);

   for(int k=0; k<rows; k++)                   // compute 2-norm
      sum = sum + v[pivot][k]*v[pivot][k];     // of pivot column
   sum = sqrt(sum);
   for(int k=0; k<rows; k++)                   // normalize the pivot column
      v[pivot][k] = v[pivot][k]/sum;

   for(int k=pivot+1; k<cols; k++)  // reduce k-th column w.r.t. pivot column
   {
      T1 inprod(0.0);               // inner product of column k with pivot
      for(int i=0; i<rows; i++)
         inprod = inprod + v[pivot][i]*v[k][i];
      for(int i=0; i<rows; i++)
         v[k][i] = v[k][i] - inprod*v[pivot][i];        // reduction
   }
}

void CPU_complex_QR_normalize_and_reduce
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

void CPU_real_QR_normalize_and_reduce
 ( T1** v, T1** R, int rows, int cols, int pivot )
{
   T1 sum(0.0);

   for(int k=0; k<rows; k++)                   // compute 2-norm
      sum = sum + v[pivot][k]*v[pivot][k];     // of pivot column
   sum = sqrt(sum);
   R[pivot][pivot] = sum;
   for(int k=0; k<rows; k++)                   // normalize the
      v[pivot][k] = v[pivot][k]/sum;           // pivot column

   for(int k=pivot+1; k<cols; k++) // reduce k-th column w.r.t. pivot column
   {
      T1 inprod(0.0);              // inner product of column k with pivot
      for(int i=0; i<rows; i++)
         inprod = inprod + v[pivot][i]*v[k][i];
      R[pivot][k] = inprod;
      for(int i=0; i<rows; i++)
         v[k][i] = v[k][i] - inprod*v[pivot][i];        // reduction
   }
}

void CPU_complex_backsubstitution
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

void CPU_real_backsubstitution ( T1** U, T1* rhs, T1* x, int dim )
{
   T1 temp;

   for(int i=dim-1; i>=0; i--)
   {
      temp = rhs[i];
      for(int j=i+1; j<dim; j++) temp = temp - U[i][j]*x[j];
      x[i] = temp/U[i][i];
   }
}

void CPU_complex_mgs2 ( complexH<T1>** v, int rows, int cols )
{
   for(int piv=0; piv<cols; piv++)
      CPU_complex_normalize_and_reduce(v,rows,cols,piv);
}

void CPU_real_mgs2 ( T1** v, int rows, int cols )
{
   for(int piv=0; piv<cols; piv++)
      CPU_real_normalize_and_reduce(v,rows,cols,piv);
}

void CPU_complex_mgs2qr
 ( complexH<T1>** v, complexH<T1>** R, int rows, int cols )
{
   for(int piv=0; piv<cols; piv++)
      CPU_complex_QR_normalize_and_reduce(v,R,rows,cols,piv);
}

void CPU_real_mgs2qr ( T1** v, T1** R, int rows, int cols )
{
   for(int piv=0; piv<cols; piv++)
      CPU_real_QR_normalize_and_reduce(v,R,rows,cols,piv);
}

void CPU_complex_mgs2qrls 
 ( complexH<T1>** v, complexH<T1>** R, complexH<T1>* x, int rows, int cols )
{
   for(int piv=0; piv<cols-1; piv++) // one extra column with rhs vector
      CPU_complex_QR_normalize_and_reduce(v,R,rows,cols,piv);

   complexH<T1>* rhs = new complexH<T1>[rows];
   for(int i=0; i<rows; i++) rhs[i] = R[i][cols-1];

   // cout << "the triangular matrix R :" << endl;
   // print_complex_matrix(R,cols,cols);
   CPU_complex_backsubstitution(R,rhs,x,rows);
}

void CPU_real_mgs2qrls ( T1** v, T1** R, T1* x, int rows, int cols )
{
   for(int piv=0; piv<cols-1; piv++) // one extra column with rhs vector
      CPU_real_QR_normalize_and_reduce(v,R,rows,cols,piv);

   T1* rhs = new T1[rows];
   for(int i=0; i<rows; i++) rhs[i] = R[i][cols-1];

   // cout << "the triangular matrix R :" << endl;
   // print_complex_matrix(R,cols,cols);
   CPU_real_backsubstitution(R,rhs,x,rows);
}
