// prototypes of code for execution of modified Gram-Schmidt on the host

#include "mgs_host.h"
#include "DefineType.h"
#include "gqd_qd_utilT.h"

// print functions 

void print ( complex<T>* a, int dim, string var )
{
   complexH<T1> temp;

   for(int i=0; i<dim; i++)
   {
      comp1_gqd2qd(&a[i],&temp);
      cout << "GPU: " << var << "["<< i << "] = " << temp;
   }
}

void print ( complex<T>* a, complexH<T1>* b, int dim, string var )
{
   complexH<T1> temp;

   for(int i=0; i<dim; i++)
   {
      comp1_gqd2qd(&a[i],&temp);
      cout << "GPU: " << var << "["<< i << "] = " << temp;
      cout << "CPU: " << var << "["<< i << "] = " << b[i];
   }
}

void print
 ( complex<T>* a, complexH<T1>** b,
   int dimX, int dimY, int stride, string var )
{
   complexH<T1> temp;

   for(int i=0; i<dimX; i++)
      for(int j=0; j<dimY; j++)
      {
         comp1_gqd2qd(&a[stride+dimY*i+j],&temp);
         cout << "GPU: "
              << var << "["<< i << "]" << "["<< j << "] = " << temp;
         cout << "CPU: "
              << var << "["<< i << "]" << "["<< j << "] = " << b[i][j];
      }
}

void print ( complexH<T1>** b, int dimX, int dimY, string var )
{
   complexH<T1> temp;

   for(int i=0; i<dimX; i++)
      for(int j=0; j<dimY; j++)
         cout << "CPU: "
              << var << "["<< i << "]" << "["<< j << "] = " << b[i][j];
   cout << endl;
}

void print ( complexH<T1>* b, int dim, string var )
{
   complexH<T1> temp;

   for(int i=0; i<dim; i++)
      cout << "CPU: " << var << "["<< i << "] = "  << b[i];
   cout << endl;
}

// library for the CPU version of the solver

inline complexH<T1> inProd ( complexH<T1>* u, complexH<T1>* v, int dim )
{
    complexH<T1> sum(0.0,0.0);

    for(int i=0; i<dim; i++) sum = sum + u[i].adj()*v[i];

    return sum;
}

inline complexH<T1> inProd_nconj(complexH<T1>* u, complexH<T1>* v, int dim )
{
   complexH<T1> sum(0.0,0.0);

   for(int i=0;i<dim;i++) sum=sum+u[i]*v[i];

   return sum;
}

void MMprod
 ( complexH<T1>** M1, complexH<T1>** M2, complexH<T1>** PR,
   int NrM1, int NcM1, int NcM2)
{
   for(int i=0; i<NrM1; i++)
      for(int j=0; j<NcM2; j++)
         PR[i][j] = inProd_nconj(M1[i],M2[j],NcM1);
}

void MVprod
 ( complexH<T1>** M, complexH<T1>* V, complexH<T1>* PR, int NrM, int NcM )
{
   for(int i=0; i<NrM; i++) PR[i] = inProd(M[i],V,NcM);
}

void MVprod_nconj
 ( complexH<T1>** M, complexH<T1>* V, complexH<T1>* PR, int NrM, int NcM )
{
   for(int i=0; i<NrM; i++)
      PR[i] = inProd_nconj(M[i],V,NcM);
}

void matrixTr ( complexH<T1>** M, complexH<T1>** Mtr, int dimX, int dimY )
{
   for(int i=0; i<dimY; i++)
      for(int j=0; j<dimX; j++) Mtr[i][j] = M[j][i];
}

complexH<T1> MatrixError
 ( complexH<T1>** a, complexH<T1>** b, int dimX, int dimY )
{
   complexH<T1> error;

   for(int i=0; i<dimX; i++)
      for(int j=0; j<dimY; j++)
         error = error + (a[i][j]-b[i][j])*(a[i][j]-b[i][j]).adj();

   return error;
}

complexH<T1> VectorError ( complexH<T1>* a, complexH<T1>* b, int dim )
{
   complexH<T1> error;

   for(int i=0; i<dim; i++)
      error = error + (a[i]-b[i])*(a[i]-b[i]).adj();

   return error;
}

complexH<T1> GPUvsCPUVectorError ( complex<T>* a, complexH<T1>* b, int dim )
{
   complexH<T1> temp;
   complexH<T1> error(0.0,0.0);

   for(int i=0; i<dim; i++)
   {
      comp1_gqd2qd(&a[i],&temp);
      // cout << temp;
      // cout << error;
      error=error+(b[i]-temp)*(b[i]-temp).adj();
   }
   return error;
}

void checkGPUnormal ( complex<T>* v_h, int dim )
{
   complexH<T1> error(0.0,0.0);
   complexH<T1> one(1.0,0.0);
   complexH<T1> temp;

   for(int i=0; i<dim; i++)
   {
      complexH<T1> sum(0.0,0.0);
      for(int j=0; j<dim; j++)
      {
         comp1_gqd2qd(&v_h[i*dim+j],&temp);
         sum = sum + temp*temp.adj();
      }
      // cout << "column " << i << " has norm " << sqrt(sum.real) << endl;
      error = error + (sum - one);
   }
   cout << "GPU normality error : " << error.real << endl;
}

void checkCPUnormal ( complexH<T1>** v, int dim )
{
   complexH<T1> error(0.0,0.0);
   complexH<T1> one(1.0,0.0);
   for(int i=0; i<dim; i++)
   {
      complexH<T1> sum(0.0,0.0);
      for(int j=0; j<dim; j++)
         sum = sum + v[j][i]*v[j][i].adj();
      // cout << "column " << i << " has norm " << sqrt(sum.real) << endl;
      error = error + (sum - one);
   }
   cout << "CPU normality error : " << error.real << endl;
}

void checkGPUorthogonal ( complex<T>* v_h, int dim )
{
   complexH<T1> error(0.0,0.0);
   complexH<T1> xtemp, ytemp;
   for(int i=0; i<dim; i++)
   {
      for(int j=i+1; j<dim; j++)
      {
         complexH<T1> sum(0.0,0.0);
         for(int k=0; k<dim; k++)
         {
            comp1_gqd2qd(&v_h[i*dim+k],&xtemp);
            comp1_gqd2qd(&v_h[j*dim+k],&ytemp);
            sum = sum + xtemp*ytemp.adj();
         }
         if(abs(sum.real) > 1.0e-12)
            cout << "< " << i << " , " << j << " > = " << sum;
         error = error + sum;
      }
   }
   cout << "GPU orthogonality error : " << error;
}

void checkCPUorthogonal ( complexH<T1>** v, int dim )
{
   complexH<T1> error(0.0,0.0);
   for(int i=0; i<dim; i++)
   {
      for(int j=i+1; j<dim; j++)
      {
         complexH<T1> sum(0.0,0.0);
         for(int k=0; k<dim; k++)
            sum = sum + v[k][i]*v[k][j].adj();
         // cout << "< " << i << " , " << j << " > = " << sum;
         error = error + sum;
      }
   }
   cout << "CPU orthogonality error : " << error;
}

void checkGPUvsCPU ( complex<T>* v_h, complexH<T1>** v, int dim )
{
   complexH<T1> error(0.0,0.0);
   complexH<T1> temp;
   for(int i=0; i<dim; i++)
   {
      for(int j=0; j<dim; j++)
      {
         comp1_gqd2qd(&v_h[i*dim+j],&temp);
         temp = temp - v[i][j];
         if(abs(temp.real) > 1.0e-12)
            cout << "(" << i << "," << j << ") = " << temp;
         error = error + temp;
      }
   }
   cout << "difference between GPU and CPU matrix :\n" << error;
}

// Serial Modified Gram-Schmidt

void CPU_GS_wcopy
 ( complexH<T1>** v, complexH<T1>** v_copy, complexH<T1>** R, int dim, int k )
{
   for(int i=0; i<k; i++)
   {
      R[i][i].real = sqrt(inProd(v[i],v[i],dim).real);
      R[i][i].imag = 0.0;

      for(int j=0; j<dim; j++)
      {
         v[i][j].real = v[i][j].real/R[i][i].real;
         v[i][j].imag = v[i][j].imag/R[i][i].real;
      }

      complexH<T1>* projection = new complexH<T1>[dim];
      for(int j=i+1; j<k;j++)
      {
         R[i][j]=inProd(v[i],v_copy[j],dim);
         for(int l=0; l<dim; l++)
         {
            projection[l]=R[i][j]*v[i][l];
            v[j][l]=v[j][l]-projection[l];
         }
      }
   }
}

// in place Serial Modified Gram-Schmidt 

void CPU_GS
 ( complexH<T1>** v, complexH<T1>** R, int dim, int k )
{
   for(int i=0; i<k; i++)
   {
      R[i][i].real = sqrt(inProd(v[i],v[i],dim).real);
      R[i][i].imag = 0.0;
      for(int j=0; j<dim; j++)
      {
         v[i][j].real = v[i][j].real/R[i][i].real;
         v[i][j].imag = v[i][j].imag/R[i][i].real;
      }

      complexH<T1>* projection=new complexH<T1>[dim];
      for(int j=i+1; j<k; j++)
      {
         R[i][j]=inProd(v[i],v[j],dim);
         for(int l=0; l<dim; l++)
         {
            projection[l]=R[i][j]*v[i][l];
            v[j][l]=v[j][l]-projection[l];
         }
      }
   }
}


void matrMult ( complexH<T1>** A1, complexH<T1>** A2, int k, int dim )
{
   for(int i=0; i<k; i++)
      for(int j=0; j<k; j++)
         A2[i][j] = inProd(A1[i],A1[j],dim);
}

// Serial Back Substitution

void BackSubsSec
 ( int dim, complexH<T1>** U, complexH<T1>* y, complexH<T1>* x )
{
   complexH<T1> temp;

   for(int i=dim-1; i>-1; i--)
   {
      temp = y[i];
      for(int j=i+1; j<dim; j++) temp = temp - U[i][j]*x[j];
      x[i] = temp/U[i][i];
   }
}
