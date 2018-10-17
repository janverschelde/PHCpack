with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
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

end Homotopy_Pade_Approximants;
