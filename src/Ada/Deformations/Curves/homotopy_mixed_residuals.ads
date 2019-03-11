with text_io;                           use text_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Poly_SysFun;
with DoblDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Poly_SysFun;

package Homotopy_Mixed_Residuals is

-- DESCRIPTION :
--   Computes the mixed residual for the homotopy system in storage,
--   in double, double double, and quad double precision.

  function Standard_AbsVal_Homotopy
             return Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
  function DoblDobl_AbsVal_Homotopy
             return DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
  function QuadDobl_AbsVal_Homotopy
             return QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;

  -- DESCRIPTION :
  --   For the stored homotopy, returns the evaluable form of the
  --   polynomials in homotopy, with absolute values of the coefficients,
  --   in double, double double, or quad double precision.

  function Residual ( abh : Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : Standard_Complex_Vectors.Vector;
                      t : Standard_Complex_Numbers.Complex_Number )
                    return double_float;
  function Residual ( abh : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : DoblDobl_Complex_Vectors.Vector;
                      t : DoblDobl_Complex_Numbers.Complex_Number )
                    return double_double;
  function Residual ( abh : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : QuadDobl_Complex_Vectors.Vector;
                      t : QuadDobl_Complex_Numbers.Complex_Number )
                    return quad_double;

  -- DESCRIPTION :
  --   Returns the mixed residual for the homotopy stored in double,
  --   double double, or quad double precision, for a point on a path.

  -- ON ENTRY :
  --   abh      homotopy system with absolute values as coefficients,
  --            for the stored homotopy stored;
  --   z        a solution point;
  --   t        the corresponding value of the continuation parameter;

  function Residual ( file : file_type;
                      abh : Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : Standard_Complex_Vectors.Vector;
                      t : Standard_Complex_Numbers.Complex_Number )
                    return double_float;
  function Residual ( file : file_type;
                      abh : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : DoblDobl_Complex_Vectors.Vector;
                      t : DoblDobl_Complex_Numbers.Complex_Number )
                    return double_double;
  function Residual ( file : file_type;
                      abh : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : QuadDobl_Complex_Vectors.Vector;
                      t : QuadDobl_Complex_Numbers.Complex_Number )
                    return quad_double;

  -- DESCRIPTION :
  --   Returns the mixed residual for the homotopy stored in double,
  --   double double, or quad double precision, for a point on a path.
  --   For every polynomial, one line is written to the output file.

  -- ON ENTRY :
  --   file     output file to write norms and residuals;
  --   abh      homotopy system with absolute values as coefficients,
  --            for the stored homotopy stored;
  --   z        a solution point;
  --   t        the corresponding value of the continuation parameter;

end Homotopy_Mixed_Residuals;
