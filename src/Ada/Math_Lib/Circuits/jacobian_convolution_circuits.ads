with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Jacobian_Convolution_Circuits is

-- DESCRIPTION :
--   Jacobian matrices are computed of convolution circuits, for the constant,
--   leading term in double, double double, or quad double precision.

  function Jacobian ( c : Standard_Speelpenning_Convolutions.Circuits;
                      x : Standard_Complex_Vectors.Vector ) 
                    return Standard_Complex_Matrices.Matrix;
  function Jacobian ( c : DoblDobl_Speelpenning_Convolutions.Circuits;
                      x : DoblDobl_Complex_Vectors.Vector ) 
                    return DoblDobl_Complex_Matrices.Matrix;
  function Jacobian ( c : QuadDobl_Speelpenning_Convolutions.Circuits;
                      x : QuadDobl_Complex_Vectors.Vector ) 
                    return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the Jacobian matrix for the polynomials in the circuit c,
  --   evaluated at the vector x.

  -- REQUIRED : c.dim = x'last.

   procedure Singular_Values
               ( c : in Standard_Speelpenning_Convolutions.Circuits;
                 x : in Standard_Complex_Vectors.Vector;
                 A : out Standard_Complex_Matrices.Matrix;
                 U : out Standard_Complex_Matrices.Matrix;
                 V : out Standard_Complex_Matrices.Matrix;
                 e : out Standard_Complex_Vectors.Vector;
                 s : out Standard_Complex_Vectors.Vector );
   procedure Singular_Values
               ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                 x : in DoblDobl_Complex_Vectors.Vector;
                 A : out DoblDobl_Complex_Matrices.Matrix;
                 U : out DoblDobl_Complex_Matrices.Matrix;
                 V : out DoblDobl_Complex_Matrices.Matrix;
                 e : out DoblDobl_Complex_Vectors.Vector;
                 s : out DoblDobl_Complex_Vectors.Vector );
   procedure Singular_Values
               ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                 x : in QuadDobl_Complex_Vectors.Vector;
                 A : out QuadDobl_Complex_Matrices.Matrix;
                 U : out QuadDobl_Complex_Matrices.Matrix;
                 V : out QuadDobl_Complex_Matrices.Matrix;
                 e : out QuadDobl_Complex_Vectors.Vector;
                 s : out QuadDobl_Complex_Vectors.Vector );

   -- DESCRIPTION :
   --   Evaluates the circuit c at x, at the Jacobian matrix A,
   --   and computes the singular value decomposition of A,
   --   in double, double double, or quad double precision;

end Jacobian_Convolution_Circuits;
