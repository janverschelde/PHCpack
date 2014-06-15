with Standard_Integer_Numbers;           use Standard_Integer_NUmbers;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;

package Flag_Transformations is

-- DESCRIPTION :
--   This package exports linear algebra operations to transform
--   pairs of flags, in standard complex arithmetic.

  function Flag_Transformation_Coefficient_Matrix
             ( n : integer32;
               f1,f2,g1,g2 : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the coefficient matrix to compute the matrix A
  --   and the upper triangular matrices T1 and T2 in the transformation
  --   defined by A*f1 = g1*T1, A*f2 = g2*T2.
  --   The square matrix on return has dimension 2*n^2.

  function Flag_Transformation_Right_Hand_Side
             ( n : integer32; g : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the right hand side vector in the linear system to
  --   compute a transformation between two flags.

  procedure Extract_Matrices
              ( n : in integer32; sol : in Standard_Complex_Vectors.Vector;
                A,T1,T2 : out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Given a vector sol of size 2*n*n, extracts three n-dimensional 
  --   matrices: A, T1, and T2, where T1 and T2 are upper triangular
  --   and T1 has ones on its diagonal.

end Flag_Transformations;
