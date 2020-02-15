with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;

package Multitasked_Pade_Approximations is

-- DESCRIPTION :
--   Provides multitasked algorithms for rational approximations.

  procedure Standard_Multitasking
              ( nbtasks,numdeg,dendeg : in integer32;
                cff : in Standard_Complex_VecVecs.VecVec;
                numcff,dencff : in Standard_Complex_VecVecs.VecVec;
                t : in double_float;
                eva : out Standard_Complex_Vectors.Vector;
                output : in boolean := false );
  procedure DoblDobl_Multitasking
              ( nbtasks,numdeg,dendeg : in integer32;
                cff : in DoblDobl_Complex_VecVecs.VecVec;
                numcff,dencff : in DoblDobl_Complex_VecVecs.VecVec;
                t : in double_double;
                eva : out DoblDobl_Complex_Vectors.Vector;
                output : in boolean := false );
  procedure QuadDobl_Multitasking
              ( nbtasks,numdeg,dendeg : in integer32;
                cff : in QuadDobl_Complex_VecVecs.VecVec;
                numcff,dencff : in QuadDobl_Complex_VecVecs.VecVec;
                t : in quad_double;
                eva : out QuadDobl_Complex_Vectors.Vector;
                output : in boolean := false );

  -- DESCRIPTION :
  --   Constructs rational approximations with multitasking,
  --   in double, double double, and quad double precision.

  -- ON ENTRY :
  --   nbtasks   the number of tasks;
  --   numdeg    degree of the numerators;
  --   dendeg    degree of the denominators;
  --   cff       coefficients of the power series;
  --   numcff    allocated space for the numerator coefficients;
  --   dencff    allocated space for the numerator coefficients;
  --   t         value to evaluate the rational approximations;
  --   output    true if the tasks are verbose,
  --             false if no output during the multitasking.

  -- ON RETURN :
  --   numcff    coefficients of the numerators;
  --   dencff    coefficients of the denominators;
  --   eva       evaluated rational approximants at t.

end Multitasked_Pade_Approximations;
