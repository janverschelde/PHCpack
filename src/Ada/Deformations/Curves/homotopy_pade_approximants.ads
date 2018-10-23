with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Dense_Series_Vectors;
with Standard_Pade_Approximants;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Pade_Approximants;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Pade_Approximants;

package Homotopy_Pade_Approximants is

-- DESCRIPTION :
--   The procedures in this package encapsulate the creators of Pade
--   approximants for solutions of polynomial systems defined by 
--   homotopies in double, double double, and quad double precision.

  procedure Standard_Pade_Approximant
              ( sol : in Standard_Complex_Vectors.Vector;
                idx,nbequ,numdeg,dendeg : in integer32;
                nbiters : in natural32;
                srv,eva : out Standard_Dense_Series_Vectors.Vector;
                pv : out Standard_Pade_Approximants.Pade_Vector;
                verbose : in boolean := false );
  procedure DoblDobl_Pade_Approximant
              ( sol : in DoblDobl_Complex_Vectors.Vector;
                idx,nbequ,numdeg,dendeg : in integer32;
                nbiters : in natural32;
                srv,eva : out DoblDobl_Dense_Series_Vectors.Vector;
                pv : out DoblDobl_Pade_Approximants.Pade_Vector;
                verbose : in boolean := false );
  procedure QuadDobl_Pade_Approximant
              ( sol : in QuadDobl_Complex_Vectors.Vector;
                idx,nbequ,numdeg,dendeg : in integer32;
                nbiters : in natural32;
                srv,eva : out QuadDobl_Dense_Series_Vectors.Vector;
                pv : out QuadDobl_Pade_Approximants.Pade_Vector;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Given a start solution of a polynomial homotopy,
  --   applies Newton's method to compute a solution series
  --   and then a vector of Pade approximants,
  --   in standard double, double double, and quad double precision.

  -- REQUIRED : the corresponding Homotopy packages Standard_Homotopy,
  --   DoblDobl_Homotopy, and QuadDobl_Homotopy are initialized with
  --   a target and start system.

  -- ON ENTRY :
  --   sol      solution of a start system in a polynomial homotopy;
  --   idx      index of the parameter in the series in the homotopy,
  --            which equals nbequ+1 for an artificial-parameter homotopy;
  --   nbequ    number of equations in the polynomial homotopy;
  --   numdeg   degree of the numerator of the Pade approximant;
  --   dendeg   degree of the denominator of the Pade approximant;
  --   nbiters  upper bound on the number of Newton iterations
  --            to compute a series solution for the solution vector.

  -- ON RETURN :
  --   srv      series solution, with initial coefficients in sol;
  --   eva      evaluated solution series;
  --   pv       vector of Pade approximants.

  function Numerical_Degree
              ( p : Standard_Complex_Vectors.Vector;
                tol : double_float ) return integer32;
  function Numerical_Degree
              ( p : DoblDobl_Complex_Vectors.Vector;
                tol : double_float ) return integer32;
  function Numerical_Degree
              ( p : QuadDobl_Complex_Vectors.Vector;
                tol : double_float ) return integer32;

  -- DESCRIPTION :
  --   The numerical degree of a polynomial with coefficients in p
  --   is the highest index in p for which the coefficient in absolute
  --   value is higher than the given tolerance tol.
  --   If all coefficients are less than tol, then -1 is returned.

  function Standard_Poles
              ( p : Standard_Pade_Approximants.Pade )
              return Standard_Complex_Vectors.Vector;
  function DoblDobl_Poles
              ( p : DoblDobl_Pade_Approximants.Pade )
              return DoblDobl_Complex_Vectors.Vector;
  function QuadDobl_Poles
              ( p : QuadDobl_Pade_Approximants.Pade )
              return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the poles of the Pade approximant,
  --   in standard double, double double, or quad double precision.
  --   The vector or return is of size 1..p'last,
  --   but only the first Numerical_Degree entries of p matter.
  --   The default poles are -1.

  function Standard_Poles
              ( pv : Standard_Pade_Approximants.Pade_Vector )
              return Standard_Complex_VecVecs.VecVec;
  function DoblDobl_Poles
              ( pv : DoblDobl_Pade_Approximants.Pade_Vector )
              return DoblDobl_Complex_VecVecs.VecVec;
  function QuadDobl_Poles
              ( pv : QuadDobl_Pade_Approximants.Pade_Vector )
              return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the poles of the vector of Pade approximants,
  --   in standard double, double double, or quad double precision.

  procedure Smallest_Forward_Pole
              ( v : in Standard_Complex_Vectors.Vector;
                idx : out integer32; minval : out double_float );
  procedure Smallest_Forward_Pole
              ( v : in DoblDobl_Complex_Vectors.Vector;
                idx : out integer32; minval : out double_double );
  procedure Smallest_Forward_Pole
              ( v : in QuadDobl_Complex_Vectors.Vector;
                idx : out integer32; minval : out quad_double );

  -- DESCRIPTION :
  --   Returns the smallest number in the vector v,
  --   in double, double double, and quad double precision,
  --   with a strictly positive real part.

  -- ON ENTRY :
  --   v        a vector of complex numbers.

  -- ON RETURN :
  --   idx      index in v'range to the smallest number in v;
  --   minval   the smallest radius, equals radius(v(idx)),
  --            which is negative if all poles have negative real part.

  procedure Smallest_Forward_Pole
              ( v : in Standard_Complex_VecVecs.VecVec;
                leadidx,idx : out integer32; minval : out double_float );
  procedure Smallest_Forward_Pole
              ( v : in DoblDobl_Complex_VecVecs.VecVec;
                leadidx,idx : out integer32; minval : out double_double );
  procedure Smallest_Forward_Pole
              ( v : in QuadDobl_Complex_VecVecs.VecVec;
                leadidx,idx : out integer32; minval : out quad_double );

  -- DESCRIPTION :
  --   Returns the smallest number in the vector of vectors v,
  --   in double, double double, and quad double precision,
  --   with a strictly positive real part.

  -- ON ENTRY :
  --   v        a vector of vector of complex numbers.

  -- ON RETURN :
  --   leadidx  index in v'range to the smallest number in v;
  --   idx      index in v(leadidx)'range to the smallest number in v;
  --   minval   the smallest radius, equals radius(v(leadidx)(idx)),
  --            which is negative if all poles have negative real part.

end Homotopy_Pade_Approximants;
