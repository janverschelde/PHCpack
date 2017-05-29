with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Dense_Series_Vectors;
with Standard_Pade_Approximants;
with DoblDobl_Complex_Vectors;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Pade_Approximants;
with QuadDobl_Complex_Vectors;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Pade_Approximants;

package Homotopy_Pade_Approximants is

-- DESCRIPTION :
--   The procedures in this package encapsulate the creators of Pade
--   approximants for solutions of polynomial systems defined by 
--   homotopies in double, double double, and quad double precision.

  procedure Standard_Pade_Approximant
              ( sol : in Standard_Complex_Vectors.Vector;
                nbequ,numdeg,dendeg : in integer32; nbiters : in natural32;
                srv,eva : out Standard_Dense_Series_Vectors.Vector;
                pv : out Standard_Pade_Approximants.Pade_Vector );
  procedure DoblDobl_Pade_Approximant
              ( sol : in DoblDobl_Complex_Vectors.Vector;
                nbequ,numdeg,dendeg : in integer32; nbiters : in natural32;
                srv,eva : out DoblDobl_Dense_Series_Vectors.Vector;
                pv : out DoblDobl_Pade_Approximants.Pade_Vector );
  procedure QuadDobl_Pade_Approximant
              ( sol : in QuadDobl_Complex_Vectors.Vector;
                nbequ,numdeg,dendeg : in integer32; nbiters : in natural32;
                srv,eva : out QuadDobl_Dense_Series_Vectors.Vector;
                pv : out QuadDobl_Pade_Approximants.Pade_Vector );

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
  --   nbequ    number of equations in the polynomial homotopy;
  --   numdeg   degree of the numerator of the Pade approximant;
  --   dendeg   degree of the denominator of the Pade approximant;
  --   nbiters  upper bound on the number of Newton iterations
  --            to compute a series solution for the solution vector.

  -- ON RETURN :
  --   srv      series solution, with initial coefficients in sol;
  --   eva      evaluated solution series;
  --   pv       vector of Pade approximants.

end Homotopy_Pade_Approximants;
