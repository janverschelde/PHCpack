with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_VecMats;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_VecVecs;
with TripDobl_Complex_Matrices;
with TripDobl_Complex_VecMats;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_VecMats;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_VecVecs;
with PentDobl_Complex_Matrices;
with PentDobl_Complex_VecMats;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_VecVecs;
with OctoDobl_Complex_Matrices;
with OctoDobl_Complex_VecMats;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_VecVecs;
with DecaDobl_Complex_Matrices;
with DecaDobl_Complex_VecMats;
with Standard_Complex_Series;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_Matrices;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Matrices;
with TripDobl_Complex_Series;
with TripDobl_Complex_Series_Vectors;
with TripDobl_Complex_Series_Matrices;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Matrices;
with PentDobl_Complex_Series;
with PentDobl_Complex_Series_Vectors;
with PentDobl_Complex_Series_Matrices;
with OctoDobl_Complex_Series;
with OctoDobl_Complex_Series_Vectors;
with OctoDobl_Complex_Series_Matrices;
with DecaDobl_Complex_Series;
with DecaDobl_Complex_Series_Vectors;
with DecaDobl_Complex_Series_Matrices;

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
  function Difference ( s : TripDobl_Complex_Series.Link_to_Series;
                        c : TripDobl_Complex_Vectors.Link_to_Vector )
                      return triple_double;
  function Difference ( s : QuadDobl_Complex_Series.Link_to_Series;
                        c : QuadDobl_Complex_Vectors.Link_to_Vector )
                      return quad_double;
  function Difference ( s : PentDobl_Complex_Series.Link_to_Series;
                        c : PentDobl_Complex_Vectors.Link_to_Vector )
                      return penta_double;
  function Difference ( s : OctoDobl_Complex_Series.Link_to_Series;
                        c : OctoDobl_Complex_Vectors.Link_to_Vector )
                      return octo_double;
  function Difference ( s : DecaDobl_Complex_Series.Link_to_Series;
                        c : DecaDobl_Complex_Vectors.Link_to_Vector )
                      return deca_double;

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of the differences
  --   between the coefficients of s and the values in c,
  --   in double, double double, triple double, quad double,
  --   penta double, octo double, or deca double precision.

  -- REQUIRED : s.cff'range = c'range.

  function Difference ( x : Standard_Complex_Vectors.Link_to_Vector;
                        y : Standard_Complex_Vectors.Link_to_Vector )
                      return double_float;
  function Difference ( x : DoblDobl_Complex_Vectors.Link_to_Vector;
                        y : DoblDobl_Complex_Vectors.Link_to_Vector )
                      return double_double;
  function Difference ( x : TripDobl_Complex_Vectors.Link_to_Vector;
                        y : TripDobl_Complex_Vectors.Link_to_Vector )
                      return triple_double;
  function Difference ( x : QuadDobl_Complex_Vectors.Link_to_Vector;
                        y : QuadDobl_Complex_Vectors.Link_to_Vector )
                      return quad_double;
  function Difference ( x : PentDobl_Complex_Vectors.Link_to_Vector;
                        y : PentDobl_Complex_Vectors.Link_to_Vector )
                      return penta_double;
  function Difference ( x : OctoDobl_Complex_Vectors.Link_to_Vector;
                        y : OctoDobl_Complex_Vectors.Link_to_Vector )
                      return octo_double;
  function Difference ( x : DecaDobl_Complex_Vectors.Link_to_Vector;
                        y : DecaDobl_Complex_Vectors.Link_to_Vector )
                      return deca_double;

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of the differences
  --   between the coefficients of x and those values of y,
  --   in double, double double, triple double, quad double,
  --   penta double, octo double, or deca double precision.

  -- REQUIRED : x'range = y'range.

  function Difference ( s : Standard_Complex_Series_Vectors.Vector;
                        c : Standard_Complex_VecVecs.VecVec )
                      return double_float;
  function Difference ( s : DoblDobl_Complex_Series_Vectors.Vector;
                        c : DoblDobl_Complex_VecVecs.VecVec )
                      return double_double;
  function Difference ( s : TripDobl_Complex_Series_Vectors.Vector;
                        c : TripDobl_Complex_VecVecs.VecVec )
                      return triple_double;
  function Difference ( s : QuadDobl_Complex_Series_Vectors.Vector;
                        c : QuadDobl_Complex_VecVecs.VecVec )
                      return quad_double;
  function Difference ( s : PentDobl_Complex_Series_Vectors.Vector;
                        c : PentDobl_Complex_VecVecs.VecVec )
                      return penta_double;
  function Difference ( s : OctoDobl_Complex_Series_Vectors.Vector;
                        c : OctoDobl_Complex_VecVecs.VecVec )
                      return octo_double;
  function Difference ( s : DecaDobl_Complex_Series_Vectors.Vector;
                        c : DecaDobl_Complex_VecVecs.VecVec )
                      return deca_double;

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of the differences
  --   between the coefficients of s and the values in c,
  --   in double, double double, triple double, quad double,
  --   penta double, octo double, or deca double precision.

  -- REQUIRED : s'range = c'range.

  function Difference ( x : Standard_Complex_VecVecs.VecVec;
                        y : Standard_Complex_VecVecs.VecVec )
                      return double_float;
  function Difference ( x : DoblDobl_Complex_VecVecs.VecVec;
                        y : DoblDobl_Complex_VecVecs.VecVec )
                      return double_double;
  function Difference ( x : TripDobl_Complex_VecVecs.VecVec;
                        y : TripDobl_Complex_VecVecs.VecVec )
                      return triple_double;
  function Difference ( x : QuadDobl_Complex_VecVecs.VecVec;
                        y : QuadDobl_Complex_VecVecs.VecVec )
                      return quad_double;
  function Difference ( x : PentDobl_Complex_VecVecs.VecVec;
                        y : PentDobl_Complex_VecVecs.VecVec )
                      return penta_double;
  function Difference ( x : OctoDobl_Complex_VecVecs.VecVec;
                        y : OctoDobl_Complex_VecVecs.VecVec )
                      return octo_double;
  function Difference ( x : DecaDobl_Complex_VecVecs.VecVec;
                        y : DecaDobl_Complex_VecVecs.VecVec )
                      return deca_double;

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of the differences
  --   between the coefficients of x and those in y,
  --   in double, double double, triple double, quad double,
  --   penta double, or octo double precision.

  -- REQUIRED : x'range = y'range.

  function Difference ( jm : Standard_Complex_Series_Matrices.Matrix;
                        vm : Standard_Complex_VecMats.VecMat )
                      return double_float;
  function Difference ( jm : DoblDobl_Complex_Series_Matrices.Matrix;
                        vm : DoblDobl_Complex_VecMats.VecMat )
                      return double_double;
  function Difference ( jm : TripDobl_Complex_Series_Matrices.Matrix;
                        vm : TripDobl_Complex_VecMats.VecMat )
                      return triple_double;
  function Difference ( jm : QuadDobl_Complex_Series_Matrices.Matrix;
                        vm : QuadDobl_Complex_VecMats.VecMat )
                      return quad_double;
  function Difference ( jm : PentDobl_Complex_Series_Matrices.Matrix;
                        vm : PentDobl_Complex_VecMats.VecMat )
                      return penta_double;
  function Difference ( jm : OctoDobl_Complex_Series_Matrices.Matrix;
                        vm : OctoDobl_Complex_VecMats.VecMat )
                      return octo_double;
  function Difference ( jm : DecaDobl_Complex_Series_Matrices.Matrix;
                        vm : DecaDobl_Complex_VecMats.VecMat )
                      return deca_double;

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
  function Difference ( vm1 : TripDobl_Complex_VecMats.VecMat;
                        vm2 : TripDobl_Complex_VecMats.VecMat )
                      return triple_double;
  function Difference ( vm1 : QuadDobl_Complex_VecMats.VecMat;
                        vm2 : QuadDobl_Complex_VecMats.VecMat )
                      return quad_double;
  function Difference ( vm1 : PentDobl_Complex_VecMats.VecMat;
                        vm2 : PentDobl_Complex_VecMats.VecMat )
                      return penta_double;
  function Difference ( vm1 : OctoDobl_Complex_VecMats.VecMat;
                        vm2 : OctoDobl_Complex_VecMats.VecMat )
                      return octo_double;
  function Difference ( vm1 : DecaDobl_Complex_VecMats.VecMat;
                        vm2 : DecaDobl_Complex_VecMats.VecMat )
                      return deca_double;

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of the differences
  --   between the coefficients of vm1 and the values in vm2.

  -- REQUIRED : vm1'range = vm2'range and for k in vm1'range:
  --   vm1(k)'range(1) = vm2'range(1) and vm1(k)'range(2) = vm2(k)'range(2).

  function Sum_of_Errors
             ( x,y : in Standard_Complex_Vectors.Vector )
             return double_float;
  function Sum_of_Errors
             ( x,y : in DoblDobl_Complex_Vectors.Vector )
             return double_double;
  function Sum_of_Errors
             ( x,y : in TripDobl_Complex_Vectors.Vector )
             return triple_double;
  function Sum_of_Errors
             ( x,y : in QuadDobl_Complex_Vectors.Vector )
             return quad_double;
  function Sum_of_Errors
             ( x,y : in PentDobl_Complex_Vectors.Vector )
             return penta_double;
  function Sum_of_Errors
             ( x,y : in OctoDobl_Complex_Vectors.Vector )
             return octo_double;
  function Sum_of_Errors
             ( x,y : in DecaDobl_Complex_Vectors.Vector )
             return deca_double;

  -- DESCRIPTION :
  --   Returns the sum of the component-wise differences between
  --   the vectors x and y.

  -- REQUIRED : x'range = y'range.

  function Sum_of_Errors
             ( A,B : in Standard_Complex_Matrices.Matrix )
             return double_float;
  function Sum_of_Errors
             ( A,B : in DoblDobl_Complex_Matrices.Matrix )
             return double_double;
  function Sum_of_Errors
             ( A,B : in TripDobl_Complex_Matrices.Matrix )
             return triple_double;
  function Sum_of_Errors
             ( A,B : in QuadDobl_Complex_Matrices.Matrix )
             return quad_double;
  function Sum_of_Errors
             ( A,B : in PentDobl_Complex_Matrices.Matrix )
             return penta_double;
  function Sum_of_Errors
             ( A,B : in OctoDobl_Complex_Matrices.Matrix )
             return octo_double;
  function Sum_of_Errors
             ( A,B : in DecaDobl_Complex_Matrices.Matrix )
             return deca_double;

  -- DESCRIPTION :
  --   Returns the sum of the component-wise differences between
  --   the matrices A and B.

  -- REQUIRED : A'range(1) = B'range(1) and A'range(2) = B'range(2).

end Evaluation_Differentiation_Errors;
