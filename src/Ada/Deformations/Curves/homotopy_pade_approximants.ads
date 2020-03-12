with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Series;
with Standard_Complex_Series_Vectors;
with Standard_Pade_Approximants;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Pade_Approximants;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_Vectors;
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
                srv,eva : out Standard_Complex_Series_Vectors.Vector;
                pv : out Standard_Pade_Approximants.Pade_Vector;
                verbose : in boolean := false );
  procedure DoblDobl_Pade_Approximant
              ( sol : in DoblDobl_Complex_Vectors.Vector;
                idx,nbequ,numdeg,dendeg : in integer32;
                nbiters : in natural32;
                srv,eva : out DoblDobl_Complex_Series_Vectors.Vector;
                pv : out DoblDobl_Pade_Approximants.Pade_Vector;
                verbose : in boolean := false );
  procedure QuadDobl_Pade_Approximant
              ( sol : in QuadDobl_Complex_Vectors.Vector;
                idx,nbequ,numdeg,dendeg : in integer32;
                nbiters : in natural32;
                srv,eva : out QuadDobl_Complex_Series_Vectors.Vector;
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

  procedure Standard_Poles
              ( p : in Standard_Pade_Approximants.Pade;
                poles : out Standard_Complex_Vectors.Vector );
  procedure DoblDobl_Poles
              ( p : in DoblDobl_Pade_Approximants.Pade;
                poles : out DoblDobl_Complex_Vectors.Vector );
  procedure QuadDobl_Poles
              ( p : in QuadDobl_Pade_Approximants.Pade;
                poles : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Computes the poles of the Pade approximant,
  --   in standard double, double double, or quad double precision.
  --   The default poles are -1.
  --   Only the first Numerical_Degree entries of poles matter.

  -- REQUIRED :
  --   The vector poles starts at 1 and end at the degree of denominator of p.

  function Allocate_Standard_Poles
             ( dim,deg : in integer32 )
             return Standard_Complex_VecVecs.VecVec;
  function Allocate_DoblDobl_Poles
             ( dim,deg : in integer32 )
             return DoblDobl_Complex_VecVecs.VecVec;
  function Allocate_QuadDobl_Poles
             ( dim,deg : in integer32 )
             return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Allocates space to hold all poles for a Pade vector of
  --   dimension dim and degree of denominator equal to deg,
  --   for standard double, double double, or quad double precision.

  procedure Standard_Poles
              ( pv : in Standard_Pade_Approximants.Pade_Vector;
                poles : in out Standard_Complex_VecVecs.VecVec );
  procedure DoblDobl_Poles
              ( pv : in DoblDobl_Pade_Approximants.Pade_Vector;
                poles : in out DoblDobl_Complex_VecVecs.VecVec );
  procedure QuadDobl_Poles
              ( pv : in QuadDobl_Pade_Approximants.Pade_Vector;
                poles : in out QuadDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Computes the poles of all approximants in the Pade vector pv,
  --   in standard double, double double, or quad double precision.

  -- REQUIRED :
  --   All vectors in poles range from 1 till the degree of the
  --   denominator of the approximant in the Pade vector.

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
  --   The vector on return ranges between 1 and the degree of the denominator
  --   of p, although only the first Numerical_Degree entries in the returned
  --   vector matter.  The default poles are -1.

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

  procedure Closest_Pole
              ( v : in Standard_Complex_Vectors.Vector;
                idx : out integer32; minval : out double_float );
  procedure Closest_Pole
              ( v : in DoblDobl_Complex_Vectors.Vector;
                idx : out integer32; minval : out double_double );
  procedure Closest_Pole
              ( v : in QuadDobl_Complex_Vectors.Vector;
                idx : out integer32; minval : out quad_double );

  -- DESCRIPTION :
  --   Returns the smallest radius in the vector v,
  --   in double, double double, and quad double precision.

  -- ON ENTRY :
  --   v        a vector of complex numbers.

  -- ON RETURN :
  --   idx      index in v'range to the smallest number in v;
  --   minval   the smallest radius, equals radius(v(idx)),
  --            which is negative if all poles have negative real part.

  procedure Closest_Pole
              ( v : in Standard_Complex_VecVecs.VecVec;
                leadidx,idx : out integer32; minval : out double_float );
  procedure Closest_Pole
              ( v : in DoblDobl_Complex_VecVecs.VecVec;
                leadidx,idx : out integer32; minval : out double_double );
  procedure Closest_Pole
              ( v : in QuadDobl_Complex_VecVecs.VecVec;
                leadidx,idx : out integer32; minval : out quad_double );

  -- DESCRIPTION :
  --   Returns the smallest radius in the vector of vectors v,
  --   in double, double double, and quad double precision.

  -- ON ENTRY :
  --   v        a vector of vector of complex numbers.

  -- ON RETURN :
  --   leadidx  index in v'range to the smallest number in v;
  --   idx      index in v(leadidx)'range to the smallest number in v;
  --   minval   the smallest radius, equals radius(v(leadidx)(idx)).

  function Closest_Pole
             ( v : Standard_Complex_VecVecs.VecVec ) return double_float;
  function Closest_Pole
             ( v : DoblDobl_Complex_VecVecs.VecVec ) return double_double;
  function Closest_Pole
             ( v : QuadDobl_Complex_VecVecs.VecVec ) return quad_double;

  -- DESCRIPTION :
  --   Returns the smallest radius in v.

  function Threshold_Index
             ( c : Standard_Complex_Vectors.Vector;
               endidx : integer32; tol : double_float ) return integer32;
  function Threshold_Index
             ( c : DoblDobl_Complex_Vectors.Vector;
               endidx : integer32; tol : double_float ) return integer32;
  function Threshold_Index
             ( c : QuadDobl_Complex_Vectors.Vector;
               endidx : integer32; tol : double_float ) return integer32;

  -- DESCRIPTION :
  --   Returns the index of the first coefficients in c,
  --   starting at the end index endidx and decrementing the index,
  --   that is larger than tol in magnitude.
  --   If all coefficients in c are less than tol in magnitude,
  --   then -1 is returned.

  function Solution_Error_Estimate
             ( sercff : Standard_Complex_Vectors.Link_to_Vector;
               numcff,dencff : Standard_Complex_Vectors.Link_to_Vector )
             return Standard_Complex_Numbers.Complex_Number;
  function Solution_Error_Estimate
             ( sercff : DoblDobl_Complex_Vectors.Link_to_Vector;
               numcff,dencff : DoblDobl_Complex_Vectors.Link_to_Vector )
             return DoblDobl_Complex_Numbers.Complex_Number;
  function Solution_Error_Estimate
             ( sercff : QuadDobl_Complex_Vectors.Link_to_Vector;
               numcff,dencff : QuadDobl_Complex_Vectors.Link_to_Vector )
             return QuadDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Given in sercff are the coefficients of a series,
  --   the numerator and denominator coefficients in numcff and dencff,
  --   of the Pade approximant constructed on the series,
  --   returns an estimate for a component of the solution vector.

  -- REQUIRED : sercff'last = numcff'last + dencff'last + 2.

  procedure Solution_Error
             ( servec : in Standard_Complex_VecVecs.VecVec;
               numcff,dencff : in Standard_Complex_VecVecs.VecVec;
               err : out Standard_Complex_Vectors.Vector );
  procedure Solution_Error
             ( servec : in DoblDobl_Complex_VecVecs.VecVec;
               numcff,dencff : in DoblDobl_Complex_VecVecs.VecVec;
               err : out DoblDobl_Complex_Vectors.Vector );
  procedure Solution_Error
             ( servec : in QuadDobl_Complex_VecVecs.VecVec;
               numcff,dencff : in QuadDobl_Complex_VecVecs.VecVec;
               err : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns in err the estimate for the error on the solution
  --   based on the series with coefficients in servec and on
  --   the coefficients of the Pade approximants in numcff and dencff,
  --   where the degree of each series is numdeg + dendeg + 2,
  --   numdeg = numcff'last and dendeg = dencff'last.

  function Solution_Error_Estimate
             ( s : Standard_Complex_Series.Link_to_Series;
               p : Standard_Pade_Approximants.Pade )
             return Standard_Complex_Numbers.Complex_Number;
  function Solution_Error_Estimate
             ( s : DoblDobl_Complex_Series.Link_to_Series;
               p : DoblDobl_Pade_Approximants.Pade )
             return DoblDobl_Complex_Numbers.Complex_Number;
  function Solution_Error_Estimate
             ( s : QuadDobl_Complex_Series.Link_to_Series;
               p : QuadDobl_Pade_Approximants.Pade )
             return QuadDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns an estimate for a component of the solution vector,
  --   based on the series s and the Pade approximant p.

  function Solution_Error
             ( srv : Standard_Complex_Series_Vectors.Vector;
               pv : Standard_Pade_Approximants.Pade_Vector )
             return Standard_Complex_Vectors.Vector;
  function Solution_Error
             ( srv : DoblDobl_Complex_Series_Vectors.Vector;
               pv : DoblDobl_Pade_Approximants.Pade_Vector )
             return DoblDobl_Complex_Vectors.Vector;
  function Solution_Error
             ( srv : QuadDobl_Complex_Series_Vectors.Vector;
               pv : QuadDobl_Pade_Approximants.Pade_Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Estimates the error on the solution of the predicted Pade 
  --   approximants pv and the series in srv.

  function Solution_Error_Norm
             ( srv : Standard_Complex_Series_Vectors.Vector;
               pv : Standard_Pade_Approximants.Pade_Vector )
             return double_float;
  function Solution_Error_Norm
             ( srv : DoblDobl_Complex_Series_Vectors.Vector;
               pv : DoblDobl_Pade_Approximants.Pade_Vector )
             return double_double;
  function Solution_Error_Norm
             ( srv : QuadDobl_Complex_Series_Vectors.Vector;
               pv : QuadDobl_Pade_Approximants.Pade_Vector )
             return quad_double;

  -- DESCRIPTION :
  --   Returns the 2-norm of the estimate error of the solution vector
  --   of the predicted Pade approximations pv and the series in srv.

end Homotopy_Pade_Approximants;
