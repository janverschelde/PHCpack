// defines code for execution on the host

template <class ComplexType, class RealType>
void CPU_normalize_and_reduce
 ( ComplexType** v, int rows, int cols, int pivot )
{
   RealType sum(0.0);

   for(int k=0; k<rows; k++)                            // compute 2-norm
      sum = sum + v[pivot][k].real*v[pivot][k].real     // of pivot column
                + v[pivot][k].imag*v[pivot][k].imag;
   sum = sqrt(sum);
   for(int k=0; k<rows; k++)                            // normalize the
      v[pivot][k] = v[pivot][k]/sum;                    // pivot column

   for(int k=pivot+1; k<cols; k++) // reduce k-th column w.r.t. pivot column
   {
      ComplexType inprod(0.0,0.0); // inner product of column k with pivot
      for(int i=0; i<rows; i++)
         inprod = inprod + v[pivot][i].adj()*v[k][i];
      for(int i=0; i<rows; i++)
         v[k][i] = v[k][i] - inprod*v[pivot][i];        // reduction
   }
}

template <class ComplexType, class RealType>
void CPU_QR_normalize_and_reduce
 ( ComplexType** v, ComplexType** R, int rows, int cols, int pivot )
{
   // std::cout << "--------------------------" << std::endl;
   // std::cout << "rows = " << rows << ", cols = " << cols
   //           << ", pivot = " << pivot << std::endl;

   RealType sum(0);
   RealType zero(0);

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
      ComplexType inprod(0.0,0.0); // inner product of column k with pivot
      for(int i=0; i<rows; i++)
         inprod = inprod + v[pivot][i].adj()*v[k][i];
      R[pivot][k] = inprod;
      for(int i=0; i<rows; i++)
         v[k][i] = v[k][i] - inprod*v[pivot][i];        // reduction
   }
}

template <class ComplexType>
void CPU_backsubstitution
 ( ComplexType** U, ComplexType* rhs, ComplexType* x, int dim )
{
   ComplexType temp;

   for(int i=dim-1; i>=0; i--)
   {
      temp = rhs[i];
      for(int j=i+1; j<dim; j++) temp = temp - U[i][j]*x[j];
      x[i] = temp/U[i][i];
   }
}

template <class ComplexType, class RealType>
void CPU_mgs2 ( ComplexType** v, int rows, int cols )
{
   for(int piv=0; piv<cols; piv++)
      CPU_normalize_and_reduce<ComplexType,RealType>(v,rows,cols,piv);
}

template <class ComplexType, class RealType>
void CPU_mgs2qr ( ComplexType** v, ComplexType** R, int rows, int cols )
{
   for(int piv=0; piv<cols; piv++)
      CPU_QR_normalize_and_reduce<ComplexType,RealType>(v,R,rows,cols,piv);
}

template <class ComplexType, class RealType>
void CPU_mgs2qrls
 ( ComplexType** v, ComplexType** R, ComplexType* x, int rows, int cols,
   ComplexType* rhs )
{
   for(int piv=0; piv<cols-1; piv++) // one extra column with rhs vector
      CPU_QR_normalize_and_reduce<ComplexType,RealType>(v,R,rows,cols,piv);

   for(int i=0; i<cols; i++) rhs[i] = R[i][cols-1];

   CPU_backsubstitution<ComplexType>(R,rhs,x,cols-1);
}

template <class ComplexType, class RealType>
void CPU_mgs2qrls
 ( ComplexType** v, ComplexType** R, ComplexType* x, int rows, int cols )
{
   ComplexType* rhs = new ComplexType[cols];

   CPU_mgs2qrls<ComplexType,RealType>(v,R,x,rows,cols,rhs);

   delete[] rhs;
}
