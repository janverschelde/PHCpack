with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_VecMats;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Hessian_Convolution_Circuits is

-- DESCRIPTION :
--   Hessian matrices are computed of convolution circuits, for the constant,
--   leading term in double, double double, or quad double precision.

  function Hessian ( c : Standard_Speelpenning_Convolutions.Circuit;
                     x : Standard_Complex_Vectors.Vector ) 
                   return Standard_Complex_Matrices.Matrix;
  function Hessian ( c : DoblDobl_Speelpenning_Convolutions.Circuit;
                     x : DoblDobl_Complex_Vectors.Vector ) 
                   return DoblDobl_Complex_Matrices.Matrix;
  function Hessian ( c : QuadDobl_Speelpenning_Convolutions.Circuit;
                     x : QuadDobl_Complex_Vectors.Vector ) 
                   return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the Hessian for the polynomial in the circuit c,
  --   evaluated at the vector x.

  -- REQUIRED : c.dim = x'last.

  function Hessian ( c : Standard_Speelpenning_Convolutions.Link_to_Circuit;
                     x : Standard_Complex_Vectors.Vector ) 
                   return Standard_Complex_Matrices.Matrix;
  function Hessian ( c : DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                     x : DoblDobl_Complex_Vectors.Vector ) 
                   return DoblDobl_Complex_Matrices.Matrix;
  function Hessian ( c : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                     x : QuadDobl_Complex_Vectors.Vector ) 
                   return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the Hessian for the polynomial in the circuit c,
  --   evaluated at the vector x.
  --   If c = null, then the zero matrix is returned.

  -- REQUIRED : if c /= null, then c.dim = x'last.

  function Hessians ( c : Standard_Speelpenning_Convolutions.Circuits;
                      x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_VecMats.VecMat;
  function Hessians ( c : DoblDobl_Speelpenning_Convolutions.Circuits;
                      x : DoblDobl_Complex_Vectors.Vector )
                    return DoblDobl_Complex_VecMats.VecMat;
  function Hessians ( c : QuadDobl_Speelpenning_Convolutions.Circuits;
                      x : QuadDobl_Complex_Vectors.Vector )
                    return QuadDobl_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Returns all Hessians of all circuits in c as a vector of matrices
  --   of c'range, in double, double double, or quad double precision.

  -- REQUIRED : for all i in c'range, either c(i) is null,
  --  or c(i).dim = x'last.

   procedure Singular_Values
               ( c : in Standard_Speelpenning_Convolutions.Circuit;
                 x : in Standard_Complex_Vectors.Vector;
                 A : out Standard_Complex_Matrices.Matrix;
                 U : out Standard_Complex_Matrices.Matrix;
                 V : out Standard_Complex_Matrices.Matrix;
                 e : out Standard_Complex_Vectors.Vector;
                 s : out Standard_Complex_Vectors.Vector );
   procedure Singular_Values
               ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                 x : in Standard_Complex_Vectors.Vector;
                 A : out Standard_Complex_Matrices.Matrix;
                 U : out Standard_Complex_Matrices.Matrix;
                 V : out Standard_Complex_Matrices.Matrix;
                 e : out Standard_Complex_Vectors.Vector;
                 s : out Standard_Complex_Vectors.Vector );
   procedure Singular_Values
               ( c : in DoblDobl_Speelpenning_Convolutions.Circuit;
                 x : in DoblDobl_Complex_Vectors.Vector;
                 A : out DoblDobl_Complex_Matrices.Matrix;
                 U : out DoblDobl_Complex_Matrices.Matrix;
                 V : out DoblDobl_Complex_Matrices.Matrix;
                 e : out DoblDobl_Complex_Vectors.Vector;
                 s : out DoblDobl_Complex_Vectors.Vector );
   procedure Singular_Values
               ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                 x : in DoblDobl_Complex_Vectors.Vector;
                 A : out DoblDobl_Complex_Matrices.Matrix;
                 U : out DoblDobl_Complex_Matrices.Matrix;
                 V : out DoblDobl_Complex_Matrices.Matrix;
                 e : out DoblDobl_Complex_Vectors.Vector;
                 s : out DoblDobl_Complex_Vectors.Vector );
   procedure Singular_Values
               ( c : in QuadDobl_Speelpenning_Convolutions.Circuit;
                 x : in QuadDobl_Complex_Vectors.Vector;
                 A : out QuadDobl_Complex_Matrices.Matrix;
                 U : out QuadDobl_Complex_Matrices.Matrix;
                 V : out QuadDobl_Complex_Matrices.Matrix;
                 e : out QuadDobl_Complex_Vectors.Vector;
                 s : out QuadDobl_Complex_Vectors.Vector );
   procedure Singular_Values
               ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                 x : in QuadDobl_Complex_Vectors.Vector;
                 A : out QuadDobl_Complex_Matrices.Matrix;
                 U : out QuadDobl_Complex_Matrices.Matrix;
                 V : out QuadDobl_Complex_Matrices.Matrix;
                 e : out QuadDobl_Complex_Vectors.Vector;
                 s : out QuadDobl_Complex_Vectors.Vector );

   -- DESCRIPTION :
   --   Evaluates the circuit c at x, at the Hessian matrix A,
   --   and computes the singular value decomposition of A,
   --   in double, double double, or quad double precision;

end Hessian_Convolution_Circuits;
