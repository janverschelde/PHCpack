with Standard_Integer_Numbers;           use Standard_Integer_NUmbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_VecMats;

package Flag_Transformations is

-- DESCRIPTION :
--   This package exports linear algebra operations to transform
--   pairs of flags, in standard complex arithmetic.
--   A flag in n-space is represented by an n-by-n matrix:
--   consecutive columns span linear spaces of increasing dimensions.
--   The statement of the general problem goes as follows:
--   Given two pairs of generic flags (F1, F2) and (G1, G2) in n-space,
--   compute the matrix A, upper triangular matrices T1 and T2
--   so that A*F1 = G1*T1 and A*F2 = G2*T2.
--   The specific application has the following problem statement:
--   For two pairs of flags (M, I) and (I, F) in n-space
--   where M is the moved flags, I the identity matrix,
--   and F some random n-dimensional matrix,
--   compute the matrix A, upper triangular matrices T1 and T2
--   so that A*F1 = G1*T1 and A*F2 = G2*T2.

-- DATA STRUCTURES :

  type Standard_Stack_of_Flags is 
    array ( integer32 range <> ) of Standard_Complex_VecMats.Link_to_VecMat;
  type DoblDobl_Stack_of_Flags is 
    array ( integer32 range <> ) of DoblDobl_Complex_VecMats.Link_to_VecMat;
  type QuadDobl_Stack_of_Flags is 
    array ( integer32 range <> ) of QuadDobl_Complex_VecMats.Link_to_VecMat;

-- DEFINING A LINEAR SYSTEM :

  function Coefficient_Matrix
             ( n : integer32;
               f1,f2,g1,g2 : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix;
  function Coefficient_Matrix
             ( n : integer32;
               f1,f2,g1,g2 : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Matrices.Matrix;
  function Coefficient_Matrix
             ( n : integer32;
               f1,f2,g1,g2 : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the coefficient matrix to compute the matrix A
  --   and the upper triangular matrices T1 and T2 in the transformation
  --   defined by A*f1 = g1*T1, A*f2 = g2*T2.
  --   The square matrix on return has dimension 2*n^2.

  -- ON ENTRY :
  --   n       ambient dimension;
  --   f1      the first flag of the first pair of flags;
  --   f2      the second flag of the first pair of flags;
  --   g1      the first flag of the second pair of flags;
  --   g2      the second flag of the second pair of flags.

  -- ON RETURN :
  --   The coefficient matrix of the linear system that defines the
  --   matrix A and the transformation matrices T1 and T2.

  function Right_Hand_Side
             ( n : integer32; g : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Vectors.Vector;
  function Right_Hand_Side
             ( n : integer32; g : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Vectors.Vector;
  function Right_Hand_Side
             ( n : integer32; g : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the right hand side vector in the linear system to
  --   compute a transformation between two flags.

-- PROCESSING THE SOLUTION :

  procedure Extract_Matrices
              ( n : in integer32; sol : in Standard_Complex_Vectors.Vector;
                A,T1,T2 : out Standard_Complex_Matrices.Matrix );
  procedure Extract_Matrices
              ( n : in integer32; sol : in DoblDobl_Complex_Vectors.Vector;
                A,T1,T2 : out DoblDobl_Complex_Matrices.Matrix );
  procedure Extract_Matrices
              ( n : in integer32; sol : in QuadDobl_Complex_Vectors.Vector;
                A,T1,T2 : out QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Given a vector sol of size 2*n*n, extracts three n-dimensional 
  --   matrices: A, T1, and T2, where T1 and T2 are upper triangular
  --   and T1 has ones on its diagonal.

-- THE MAIN TRANSFORMATION :

  procedure Transform
              ( n : in integer32; 
                f1,f2,g1,g2 : in Standard_Complex_Matrices.Matrix;
                A,T1,T2 : out Standard_Complex_Matrices.Matrix );
  procedure Transform
              ( n : in integer32; 
                f1,f2,g1,g2 : in DoblDobl_Complex_Matrices.Matrix;
                A,T1,T2 : out DoblDobl_Complex_Matrices.Matrix );
  procedure Transform
              ( n : in integer32; 
                f1,f2,g1,g2 : in QuadDobl_Complex_Matrices.Matrix;
                A,T1,T2 : out QuadDobl_Complex_Matrices.Matrix );

  -- DECRIPTION :
  --   Transforms two pairs of flags in n-space.

  -- ON ENTRY :
  --   n        ambient dimension;
  --   f1       the first flag of the first pair of flags;
  --   f2       the second flag of the first pair of flags;
  --   g1       the first flag of the second pair of flags;
  --   g2       the second flag of the second pair of flags.

  -- ON RETURN :
  --   A        a regular n-by-n matrix: A*(f1, f2) = (g1*T1, g2*T2);
  --   T1       upper triangular matrix to postmultiply g1 with;
  --   T2       upper triangular matrix to postmultiply g2 with.

  function Residual 
              ( f1,f2,g1,g2,A,T1,T2 : Standard_Complex_Matrices.Matrix ) 
              return double_float;
  function Residual 
              ( f1,f2,g1,g2,A,T1,T2 : DoblDobl_Complex_Matrices.Matrix ) 
              return double_double;
  function Residual 
              ( f1,f2,g1,g2,A,T1,T2 : QuadDobl_Complex_Matrices.Matrix ) 
              return quad_double;

  -- DESCRIPTION :
  --   Returns the 1-norm of the differences A*f1 - g1*T1 and A*f2 - g2*T2,
  --   where the input arguments are the output arguments of Transform.

-- GENERAL WRAPPERS :

  function Move_to_Generic_Flag
              ( n : integer32; G : Standard_Complex_Matrices.Matrix )
              return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Does the transformation from the moving flag into a generic flag
  --   based on the given random n-by-n matrix G.
  --   Without providing G as below, the caller has no control over the
  --   choice for random numbers used to make the generic flag.

  function Move_to_Generic_Flag
              ( n : integer32 ) return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   To resolve a triple Schubert condition in n-space with the
  --   Littlewood-Richardson homotopies, three flags are involved:
  --   the identity, the fixed flag, and the moving flag. 
  --   This function turns the moving flag into a generic flag.

  -- ON ENTRY :
  --   n       the ambient dimension.

  -- ON RETURN :
  --   A generic flag computed so the solution to the triple intersection
  --   problem involving the moving flag will satisfy the generic flag.

  procedure Move_to_Generic_Flag
              ( n : in integer32; F : out Standard_Complex_Matrices.Matrix;
                rsd : out double_float );

  -- DESCRIPTION :
  --   In addition to returning the generic flag, returns the residual of
  --   the linear system solved when computing the transformation.
  --   This residual allows to check the accuracy of the result.

  -- ON ENTRY :
  --   n        the ambient dimension.

  -- ON RETURN :
  --   F        plays the same role as Move_to_Generic_Flag(n),
  --            although different random numbers are generated;
  --   rsd      is the outcome of Residual() function from above.

-- FOR APPLICATION TO RESOLVE SCHUBERT PROBLEMS :

  procedure Transform_Sequence_with_Flag
              ( n,i : in integer32;
                flags : in out Standard_Complex_VecMats.VecMat;
                A,invA,sT : out Standard_Complex_Matrices.Matrix );
  procedure Transform_Sequence_with_Flag
              ( n,i : in integer32;
                flags : in out DoblDobl_Complex_VecMats.VecMat;
                A,invA,sT : out DoblDobl_Complex_Matrices.Matrix );
  procedure Transform_Sequence_with_Flag
              ( n,i : in integer32;
                flags : in out QuadDobl_Complex_VecMats.VecMat;
                A,invA,sT : out QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Takes the i-th flag in flags, that is flags(i), to transform
  --   the pair (M, I) to (I, flags(i)).  Multiplies all flags that
  --   appear with higher index than i with the inverse of A, in invA,
  --   in standard double, double double, or quad double precision.

  -- REQUIRED : i in f'flags'range.

  -- ON ENTRY :
  --   n        the ambient dimensions;
  --   flags    a sequence of flags in n-space.

  -- ON RETURN :
  --   flags    flags(j) for j in i+1..flags'last are multiplied
  --            with the inverse of A, in invA;
  --   A        invertible transformation on flags;
  --   invA     the inverse of A;
  --   sT       equals A times the moved flag, for use to transform
  --            the solutions to the transformed problem.

  procedure Create ( n : in integer32;
                     flags : in Standard_Complex_VecMats.VecMat;
                     stack : out Standard_Stack_of_Flags;
                     A,invA,sT : out Standard_Complex_VecMats.VecMat );
  procedure Create ( n : in integer32;
                     flags : in DoblDobl_Complex_VecMats.VecMat;
                     stack : out DoblDobl_Stack_of_Flags;
                     A,invA,sT : out DoblDobl_Complex_VecMats.VecMat );
  procedure Create ( n : in integer32;
                     flags : in QuadDobl_Complex_VecMats.VecMat;
                     stack : out QuadDobl_Stack_of_Flags;
                     A,invA,sT : out QuadDobl_Complex_VecMats.VecMat );

  -- DESCRIPTION :
  --   Makes a stack of flags, successively transformed,
  --   in standard double, double double, or quad double precision.

  -- REQUIRED : stack'range = flags'first..flags'last-1.

  -- ON ENTRY :
  --   n        the ambient dimension;
  --   flags    a sequence of flags.

  -- ON RETURN :
  --   stack    a stack of flags, in stack(i), flags(i) was used;
  --   A        sequence of invertible transformations on flags;
  --   invA     sequence of inverses of the matrices in A;
  --   sT       sequence of A times the moved flag.

  procedure Clear ( s : in out Standard_Stack_of_Flags );
  procedure Clear ( s : in out DoblDobl_Stack_of_Flags );
  procedure Clear ( s : in out QuadDobl_Stack_of_Flags );

  -- DESCRIPTION :
  --   Deallocation of the stack of flags in s.

end Flag_Transformations;
