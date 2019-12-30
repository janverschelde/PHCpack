with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecMats;
with Standard_Complex_Series;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_Matrices;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Matrices;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Matrices;

package Evaluation_Differentiation_Errors is

-- DESCRIPTION :
--   The functions in this package compute the errors between the
--   results of the evaluation and differentiation methods with
--   series polynomials and convolution circuits.

  function Difference ( s : Standard_Complex_Series.Link_to_Series;
                        c : Standard_Complex_Vectors.Link_to_Vector )
                      return double_float;
  function Difference ( s : DoblDobl_Complex_Series.Link_to_Series;
                        c : DoblDobl_Complex_Vectors.Link_to_Vector )
                      return double_double;
  function Difference ( s : QuadDobl_Complex_Series.Link_to_Series;
                        c : QuadDobl_Complex_Vectors.Link_to_Vector )
                      return quad_double;

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of the differences
  --   between the coefficients of s and the values in c,
  --   in double, double double, or quad double precision.

  -- REQUIRED : s.cff'range = c'range.

  function Difference ( x : Standard_Complex_Vectors.Link_to_Vector;
                        y : Standard_Complex_Vectors.Link_to_Vector )
                      return double_float;
  function Difference ( x : DoblDobl_Complex_Vectors.Link_to_Vector;
                        y : DoblDobl_Complex_Vectors.Link_to_Vector )
                      return double_double;
  function Difference ( x : QuadDobl_Complex_Vectors.Link_to_Vector;
                        y : QuadDobl_Complex_Vectors.Link_to_Vector )
                      return quad_double;

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of the differences
  --   between the coefficients of x and those values of y,
  --   in double, double double, or quad double precision.

  -- REQUIRED : x'range = y'range.

  function Difference ( s : Standard_Complex_Series_Vectors.Vector;
                        c : Standard_Complex_VecVecs.VecVec )
                      return double_float;
  function Difference ( s : DoblDobl_Complex_Series_Vectors.Vector;
                        c : DoblDobl_Complex_VecVecs.VecVec )
                      return double_double;
  function Difference ( s : QuadDobl_Complex_Series_Vectors.Vector;
                        c : QuadDobl_Complex_VecVecs.VecVec )
                      return quad_double;

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of the differences
  --   between the coefficients of s and the values in c,
  --   in double, double double, or quad double precision.

  -- REQUIRED : s'range = c'range.

  function Difference ( x : Standard_Complex_VecVecs.VecVec;
                        y : Standard_Complex_VecVecs.VecVec )
                      return double_float;
  function Difference ( x : DoblDobl_Complex_VecVecs.VecVec;
                        y : DoblDobl_Complex_VecVecs.VecVec )
                      return double_double;
  function Difference ( x : QuadDobl_Complex_VecVecs.VecVec;
                        y : QuadDobl_Complex_VecVecs.VecVec )
                      return quad_double;

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of the differences
  --   between the coefficients of x and those in y,
  --   in double, double double, or quad double precision.

  -- REQUIRED : x'range = y'range.

  function Difference ( jm : Standard_Complex_Series_Matrices.Matrix;
                        vm : Standard_Complex_VecMats.VecMat )
                      return double_float;
  function Difference ( jm : DoblDobl_Complex_Series_Matrices.Matrix;
                        vm : DoblDobl_Complex_VecMats.VecMat )
                      return double_double;
  function Difference ( jm : QuadDobl_Complex_Series_Matrices.Matrix;
                        vm : QuadDobl_Complex_VecMats.VecMat )
                      return quad_double;

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of the differences
  --   between the coefficients of jm and the values in vm.

  -- REQUIRED : vm'range = 0..degree, where degree is the degree
  --   of all series in the matrix jm, for all k in vm'range:
  --   vm(k)'range(1) = jm'range(1) and vm(k)'range(2) = jm'range(2).

  function Difference ( vm1 : Standard_Complex_VecMats.VecMat;
                        vm2 : Standard_Complex_VecMats.VecMat )
                      return double_float;
  function Difference ( vm1 : DoblDobl_Complex_VecMats.VecMat;
                        vm2 : DoblDobl_Complex_VecMats.VecMat )
                      return double_double;
  function Difference ( vm1 : QuadDobl_Complex_VecMats.VecMat;
                        vm2 : QuadDobl_Complex_VecMats.VecMat )
                      return quad_double;

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of the differences
  --   between the coefficients of vm1 and the values in vm2.

  -- REQUIRED : vm1'range = vm2'range and for k in vm1'range:
  --   vm1(k)'range(1) = vm2'range(1) and vm1(k)'range(2) = vm2(k)'range(2).

end Evaluation_Differentiation_Errors;
