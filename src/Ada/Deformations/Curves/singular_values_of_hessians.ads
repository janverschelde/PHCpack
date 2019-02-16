with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Double_Double_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with Quad_Double_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Hessians;
with DoblDobl_Complex_Hessians;
with QuadDobl_Complex_Hessians;

package Singular_Values_of_Hessians is

-- DESCRIPTION :
--   Given symbolic definitions of Hessian matrices,
--   evaluates the Hessians at numerical vectors and
--   returns the singular values.

  procedure Singular_Values
             ( A : in out Standard_Complex_Matrices.Matrix;
               s : out Standard_Complex_Vectors.Vector );
  procedure Singular_Values
             ( A : in out DoblDobl_Complex_Matrices.Matrix;
               s : out DoblDobl_Complex_Vectors.Vector );
  procedure Singular_Values
             ( A : in out QuadDobl_Complex_Matrices.Matrix;
               s : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Computes the singular values of A and returns the result in s,
  --   computed in double, double double, or quad double precision.

  function Standard_Singular_Values
             ( h : Standard_Complex_Hessians.Link_to_Hessian;
               x : Standard_Complex_Vectors.Vector )
             return  Standard_Floating_Vectors.Vector;
  function DoblDobl_Singular_Values
             ( h : DoblDobl_Complex_Hessians.Link_to_Hessian;
               x : DoblDobl_Complex_Vectors.Vector )
             return  Double_Double_Vectors.Vector;
  function QuadDobl_Singular_Values
             ( h : QuadDobl_Complex_Hessians.Link_to_Hessian;
               x : QuadDobl_Complex_Vectors.Vector )
             return  Quad_Double_Vectors.Vector;

  -- DESCRIPTION :
  --   Given a Hessian h and a numerical vector x,
  --   returns the singular values computed in double precision,
  --   or double double precision, or quad double precision.

  function Standard_Singular_Values
             ( h : Standard_Complex_Hessians.Array_of_Hessians;
               x : Standard_Complex_Vectors.Vector )
             return  Standard_Floating_Vectors.Vector;
  function DoblDobl_Singular_Values
             ( h : DoblDobl_Complex_Hessians.Array_of_Hessians;
               x : DoblDobl_Complex_Vectors.Vector )
             return  Double_Double_Vectors.Vector;
  function QuadDobl_Singular_Values
             ( h : QuadDobl_Complex_Hessians.Array_of_Hessians;
               x : QuadDobl_Complex_Vectors.Vector )
             return  Quad_Double_Vectors.Vector;

  -- DESCRIPTION :
  --   Given an array of Hessians h and a numerical vector x,
  --   computes the singular values in double, double double,
  --   or quad double precision.  On return is a vector of h'range,
  --   with in the i-th entry the first, largest singular value
  --   of the i-th Hessian in h evaluated at x.

end Singular_Values_of_Hessians;
