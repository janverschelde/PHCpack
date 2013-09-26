with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Matrices;
 
package Standard_Nullity_Matrices is

-- DESCRIPTION :
--   This package offers implementation of nullity matrices,
--   used to calculate the dual space at an isolated zero.
--   The dimension of the dual space is the multiplicity.

  procedure Dimensions_of_Nullity_Matrix
              ( nq,nv,k : in natural32; nr,nc : out natural32 );

  -- DESCRIPTION :
  --   Returns the dimensions of the matrix S(k) to determine the
  --   multiplicity structure at stage k.

  -- ON ENTRY :
  --   nq       number of equations in the system;
  --   nv       number of variables in the system;
  --   k        if nullity(S(k)) = nullity(S(k-1)), then multiplicity = k. 

  -- ON RETURN :
  --   nr       number of rows in the matrix S(k);
  --   nc       number of columns in the matrix S(k).

  function Create_Nullity_Matrix
              ( nq,nv,nr,nc,k : natural32; f : Poly_Sys )
              return Standard_Complex_Poly_Matrices.Matrix;
  function Create_Nullity_Matrix
              ( file : file_type;
                nq,nv,nr,nc,k : natural32; f : Poly_Sys )
              return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the symbolic nullity matrix S(k) for the system f.

  -- ON ENTRY :
  --   file     for writing diagnostics;
  --   nq       number of equations in the system f;
  --   nv       number of variables in the system f;
  --   nr       number of rows in the matrix S(k);
  --   nc       number of columns in the matrix S(k);
  --   k        if nullity(S(k)) = nullity(S(k-1)), then multiplicity = k;
  --   f        polynomial system in nq equations and nv variables.

  function Eval0 ( nm : Standard_Complex_Poly_Matrices.Matrix;
                   z : Standard_Complex_Vectors.Vector )
                 return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the value of the entire nullity matrix at z.

  function Eval1 ( nm : Standard_Complex_Poly_Matrices.Matrix;
                   z : Standard_Complex_Vectors.Vector )
                 return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the value of the nullity matrix at z, skipping the first
  --   column, so that the matrix on return has one column less than nm.
 
  function Evaluate_Nullity_Matrix
              ( nq,nv,nr,nc,k : natural32; f : Poly_Sys;
                z : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Matrices.Matrix;
  function Evaluate_Nullity_Matrix
              ( nq,nv,nr,nc,k : natural32;
                a1 : Standard_Complex_Matrices.Matrix; f : Poly_Sys;
                z : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Matrices.Matrix;
  function Evaluate_Nullity_Matrix
              ( file : file_type;
                nq,nv,nr,nc,k : natural32; f : Poly_Sys;
                z : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Evaluates the nullity matrix S(k) for f at z.

  -- ON ENTRY :
  --   file     for writing diagnostics;
  --   nq       number of equations in the system f;
  --   nv       number of variables in the system f;
  --   nr       number of rows in the matrix S(k);
  --   nc       number of columns in the matrix S(k);
  --   k        if nullity(S(k)) = nullity(S(k-1)), then multiplicity = k;
  --   a1       previous nullity matrix S(k-1), evaluated at k-1,
  --            if k = 1, then a1 is the Jacobian matrix of f;
  --   f        polynomial system in nq equations and nv variables;
  --   z        vector of range 1..nv.

end Standard_Nullity_Matrices;
