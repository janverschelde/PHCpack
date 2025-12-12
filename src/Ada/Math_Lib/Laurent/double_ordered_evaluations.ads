with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;

package Double_Ordered_Evaluations is

-- DESCRIPTION :
--   Evaluates power series of first and higher order terms 
--   in Laurent polynomials, with coefficients in double precision.

  procedure Normalize
              ( cf : in out Standard_Complex_Vectors.Vector;
                dg : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Given a sorted sequence of powers in dg,
  --   adds the coefficients in cf which correspond to equal powers. 

  procedure First_Order_Evaluation
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the first order evaluation of a series in a polynomial.

  -- REQUIRED : ycf'range = pcf'first..2*pcf'lst = ydg'range.

  -- ON ENTRY :
  --   pcf      coefficients of a Laurent polynomial;
  --   pct      powers of t of the coefficients;
  --   pdg      supports of a Laurent polynomial;
  --   cff      coefficients of a power series;
  --   pwr      exponents in the power series.
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ycf      coefficients of the evaluated series, up to first order;
  --   ydg      corresponding exponents of the evaluated series.

end Double_Ordered_Evaluations;
