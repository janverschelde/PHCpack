with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Double_Double_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with Quad_Double_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Complex_Hessians;
with DoblDobl_Complex_Hessians;
with QuadDobl_Complex_Hessians;

package Singular_Values_of_Hessians is

-- DESCRIPTION :
--   Given symbolic definitions of the Jacobian matrix and Hessian matrices,
--   evaluates the Jacobian and Hessians at numerical vectors and
--   computes an estimate of the distance to the nearest solution.

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

  function Standard_Distance
             ( jm : in Standard_Complex_Jaco_Matrices.Jaco_Mat;
               hs : in Standard_Complex_Hessians.Array_of_Hessians;
               xt : in Standard_Complex_Vectors.Vector )
             return double_float;
  function DoblDobl_Distance
             ( jm : in DoblDobl_Complex_Jaco_Matrices.Jaco_Mat;
               hs : in DoblDobl_Complex_Hessians.Array_of_Hessians;
               xt : in DoblDobl_Complex_Vectors.Vector )
             return double_double;
  function QuadDobl_Distance
             ( jm : in QuadDobl_Complex_Jaco_Matrices.Jaco_Mat;
               hs : in QuadDobl_Complex_Hessians.Array_of_Hessians;
               xt : in QuadDobl_Complex_Vectors.Vector )
             return quad_double;

  -- DESCRIPTION :
  --   Returns an estimate to the distance to the nearest solution,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   jm      Jacobian matrix of a homotopy;
  --   hs      array of Hessians of the polynomials in the homotopy;
  --   xt      solution vector and corresponding t value.

  -- REQUIRED :
  --   The number of variables of all polynomials in jm and hs
  --   correspond to the end of the range of xt.
  --   Moreover, the value for the parameter in xt is at the proper place,
  --   corresponding to the homotopy.

  function Standard_Distance
             ( jm : in Standard_Complex_Jaco_Matrices.Jaco_Mat;
               hs : in Standard_Complex_Hessians.Array_of_Hessians;
               sol : in Standard_Complex_Solutions.Solution )
             return double_float;
  function DoblDobl_Distance
             ( jm : in DoblDobl_Complex_Jaco_Matrices.Jaco_Mat;
               hs : in DoblDobl_Complex_Hessians.Array_of_Hessians;
               sol : in DoblDobl_Complex_Solutions.Solution )
             return double_double;
  function QuadDobl_Distance
             ( jm : in QuadDobl_Complex_Jaco_Matrices.Jaco_Mat;
               hs : in QuadDobl_Complex_Hessians.Array_of_Hessians;
               sol : in QuadDobl_Complex_Solutions.Solution )
             return quad_double;

  -- DESCRIPTION :
  --   Returns an estimate to the distance to the nearest solution,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   jm      Jacobian matrix of a homotopy;
  --   hs      array of Hessians of the polynomials in the homotopy;
  --   sol     solution with the right number of coordinates.

  procedure Standard_Jacobian_Hessians_of_Homotopy
              ( jm : out Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                h : out Standard_Complex_Hessians.Link_to_Array_of_Hessians );
  procedure DoblDobl_Jacobian_Hessians_of_Homotopy
              ( jm : out DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                h : out DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians );
  procedure QuadDobl_Jacobian_Hessians_of_Homotopy
              ( jm : out QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                h : out QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians );

  -- DESCRIPTION :
  --   Given in Standard_Homotopy, DoblDobl_Homotopy, or QuadDobl_Homotopy,
  --   the homotopy defined in double, double double, or quad double precision,
  --   returns the Jacobian matrix and the Hessians for the homotopy.
  --   These procedures are for artificial parameter homotopies.

  -- ON RETURN :
  --   jm       the Jacobian matrix for polynomials with a homotopy
  --            continuation parameter, which is assumed as the last one;
  --   h        Hessians for the polynomials in the homotopy,
  --            with the continuation parameter as the last variable.

  procedure Standard_Jacobian_Hessians_of_Homotopy
              ( k : in integer32;
                jm : out Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                h : out Standard_Complex_Hessians.Link_to_Array_of_Hessians );
  procedure DoblDobl_Jacobian_Hessians_of_Homotopy
              ( k : in integer32;
                jm : out DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                h : out DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians );
  procedure QuadDobl_Jacobian_Hessians_of_Homotopy
              ( k : in integer32;
                jm : out QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                h : out QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians );

  -- DESCRIPTION :
  --   Given in Standard_Homotopy, DoblDobl_Homotopy, or QuadDobl_Homotopy,
  --   the homotopy defined in double, double double, or quad double precision,
  --   returns the Jacobian matrix and the Hessians for the homotopy.
  --   These procedures are for natural parameter homotopies.

  -- ON ENTRY :
  --   k        index of the continuation parameter, as the index of
  --            the variable in the polynomials of the homotopy.

  -- ON RETURN :
  --   jm       the Jacobian matrix for polynomials with a homotopy
  --            continuation parameter, which is assumed at position k;
  --   h        Hessians for the polynomials in the homotopy,
  --            with the continuation parameter as the k-th variable.

end Singular_Values_of_Hessians;
