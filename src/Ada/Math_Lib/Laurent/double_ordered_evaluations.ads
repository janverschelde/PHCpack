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

  function Size_Evaluation
             ( dim,ndf,ord,nbr : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size of the evaluation for dim variables,
  --   for all derivatives up to ndf, for ord as the number of terms
  --   (or order) of the series, and for nbr monomials.

-- ON ONE LAURENT POLYNOMIAL :

  procedure First_Derivative_First_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Evaluates a polynomial at a series truncated at the first order,
  --   using constants and the first derivatives.

  -- REQUIRED : ycf'range = 1..(dim+1)*nbr = ydg'range,
  --   where nbr is the number of terms in the polynomial, and
  --   where dim is the number of variables.

  -- ON ENTRY :
  --   pcf      coefficients of a Laurent polynomial;
  --   pct      powers of t of the coefficients;
  --   pdg      supports of a Laurent polynomial;
  --   cf0      constant coefficients of a power series;
  --   cf1      coefficients corresponding to pw1;
  --   pw1      leading exponents in the power series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ycf      coefficients of the evaluated series, up to first order;
  --   ydg      corresponding exponents of the evaluated series.

  procedure First_Derivative_First_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Evaluates a polynomial at a series truncated at the first order,
  --   using constants and the first derivatives.

  -- REQUIRED : ycf'range = 1..(dim+1)*nbr = ydg'range,
  --   where nbr is the number of terms in the polynomial, and
  --   where dim is the number of variables.

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

  procedure First_Derivative_Second_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                cf2 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                pw2 : in Standard_Floating_Vectors.Vector;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Evaluates a polynomial at a series truncated at the second order,
  --   using constants and the first derivatives.

  -- REQUIRED : ycf'range = 1..(2*dim+1)*nbr = ydg'range,
  --   where nbr is the number of terms in the polynomial, and
  --   where dim is the number of variables.

  -- ON ENTRY :
  --   pcf      coefficients of a Laurent polynomial;
  --   pct      powers of t of the coefficients;
  --   pdg      supports of a Laurent polynomial;
  --   cf0      constant coefficients of a power series;
  --   cf1      coefficients corresponding to pw1;
  --   cf2      coefficients corresponding to pw2;
  --   pw1      leading exponents in the power series;
  --   pw2      second exponents in the power series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ycf      coefficients of the evaluated series, up to second order;
  --   ydg      corresponding exponents of the evaluated series.

  procedure First_Derivative_Second_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Evaluates a polynomial at a series truncated at the second order,
  --   using constants and the first derivatives.

  -- REQUIRED : ycf'range = 1..(2*dim+1)*nbr = ydg'range,
  --   where nbr is the number of terms in the polynomial, and
  --   where dim is the number of variables.

  -- ON ENTRY :
  --   pcf      coefficients of a Laurent polynomial;
  --   pct      powers of t of the coefficients;
  --   pdg      supports of a Laurent polynomial;
  --   cff      coefficients of a power series;
  --   pwr      exponents in the power series.
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ycf      coefficients of the evaluated series, up to second order;
  --   ydg      corresponding exponents of the evaluated series.

  procedure Second_Derivative_First_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Evaluates a polynomial at a series truncated at the first order,
  --   using constants, first and second derivatives.

  -- REQUIRED : ycf'range = 1..size*nbr = ydg'range,
  --   where size = 1 + dim + dim*(dim+1)/2,
  --   where nbr is the number of terms in the polynomial, and
  --   where dim is the number of variables.

  -- ON ENTRY :
  --   pcf      coefficients of a Laurent polynomial;
  --   pct      powers of t of the coefficients;
  --   pdg      supports of a Laurent polynomial;
  --   cf0      constant coefficients of a power series;
  --   cf1      coefficients corresponding to pw1;
  --   pw1      leading exponents in the power series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ycf      coefficients of the evaluated series, up to first order;
  --   ydg      corresponding exponents of the evaluated series.

  procedure Second_Derivative_First_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Evaluates a polynomial at a series truncated at the first order,
  --   using constants, first and second derivatives.

  -- REQUIRED : ycf'range = 1..size*nbr = ydg'range,
  --   where size = 1 + dim + dim*(dim+1)/2,
  --   where nbr is the number of terms in the polynomial, and
  --   where dim is the number of variables.

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

  procedure Second_Derivative_Second_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Evaluates a polynomial at a series truncated at the second order,
  --   using constants, first and second derivatives.

  -- REQUIRED : ycf'range = 1..size*nbr = ydg'range,
  --   where size = 1 + 3*dim + 2*dim*(dim+1)/2,
  --   where nbr is the number of terms in the polynomial, and
  --   where dim is the number of variables.

  -- ON ENTRY :
  --   pcf      coefficients of a Laurent polynomial;
  --   pct      powers of t of the coefficients;
  --   pdg      supports of a Laurent polynomial;
  --   cff      coefficients of a power series;
  --   pwr      exponents in the power series.
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ycf      coefficients of the evaluated series, up to second order;
  --   ydg      corresponding exponents of the evaluated series.

  procedure Third_Derivative_First_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Evaluates a polynomial at a series truncated at the first order,
  --   using constants, first, second, and third derivatives.

  -- REQUIRED : ycf'range = 1..size*nbr = ydg'range,
  --   where size = 1 + dim + dim*(dim+1)/2
  --                  + dim + dim*(dim-1) + dim*(dim-1)*(dim-2)/3,
  --   where nbr is the number of terms in the polynomial, and
  --   where dim is the number of variables.

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

  procedure Third_Derivative_Second_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Evaluates a polynomial at a series truncated at the second order,
  --   using constants, first, second, and third derivatives.

  -- REQUIRED : ycf'range = 1..size*nbr = ydg'range,
  --   where size = 1 + 4*dim + dim*(dim+1) + 2*dim + 2*dim*(dim+1)
  --     + dim*(dim-1)*(dim-2)/6,
  --   where nbr is the number of terms in the polynomial, and
  --   where dim is the number of variables.

  -- ON ENTRY :
  --   pcf      coefficients of a Laurent polynomial;
  --   pct      powers of t of the coefficients;
  --   pdg      supports of a Laurent polynomial;
  --   cff      coefficients of a power series;
  --   pwr      exponents in the power series.
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ycf      coefficients of the evaluated series, up to second order;
  --   ydg      corresponding exponents of the evaluated series.

-- INDEXED DERIVATIVES ON ONE LAURENT POLYNOMIAL :

  procedure Fixed_Derivative_First_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                difsum : in integer32; idxnxt : in out integer32;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Updates (ycf, ydg) with all first order evaluations
  --   using indexed derivative with fixed sum difsum.

  -- ON ENTRY :
  --   pcf      coefficients of a Laurent polynomial;
  --   pct      powers of t of the coefficients;
  --   pdg      supports of a Laurent polynomial;
  --   cf0      constant coefficients of a power series;
  --   cf1      coefficients corresponding to pw1;
  --   pw1      leading exponents in the power series;
  --   difsum   sum of all derivative indices;
  --   idxnxt   next index in (ycf, ydg) to be updated;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   idxnxt   next index in (ycf, ydg) to be updated;
  --   ycf      coefficients of the evaluated series, up to first order;
  --   ydg      corresponding exponents of the evaluated series.

  procedure First_Order_Evaluation
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                difmax : in integer32;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Evaluates a Laurent polynomial at a series truncated at first order,
  --   using indexed derivatives sum up to difmax.

  -- REQUIRED :
  --   ycf and ydg have wide enough ranges.

  -- ON ENTRY :
  --   pcf      coefficients of a Laurent polynomial;
  --   pct      powers of t of the coefficients;
  --   pdg      supports of a Laurent polynomial;
  --   cf0      constant coefficients of a power series;
  --   cf1      coefficients corresponding to pw1;
  --   pw1      leading exponents in the power series;
  --   difmax   maximum sum of all derivative indices;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ycf      coefficients of the evaluated series, up to first order;
  --   ydg      corresponding exponents of the evaluated series.

-- ON A LAURENT HOMOTOPY :

  procedure First_Derivative_First_Order
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                cf2 : out Standard_Complex_Vectors.Vector;
                pw2 : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes a Taylor series expansion of the Laurent homotopy
  --   using the constant coefficients of the power series solution,
  --   truncated after the first order, using first derivatives.
  --   Computes the smallest positive exponents of this evaluation.

  -- ON ENTRY :
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hdg      supports of the Laurent homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cf0      constant coefficients of a power series;
  --   cf1      coefficients corresponding to pw1;
  --   pw1      leading exponents in the power series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   cf2      coefficients corresponding to the smallest positive powers;
  --   pw2      smallest positive powers in the evaluated series.

  procedure First_Derivative_First_Order
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                psm : out Standard_Floating_Vectors.Vector;
                csm : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes a Taylor series expansion of the Laurent homotopy
  --   using the constant coefficients of the power series solution,
  --   truncated after the first order, using first derivatives.
  --   Computes the smallest positive exponents of this evaluation.

  -- ON ENTRY :
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hdg      supports of the Laurent homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cff      coefficients of the power series solution;
  --   pwr      exponents of the power series solution;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   psm      smallest positive powers in the evaluated series;
  --   csm      coefficients corresponding to the smallest positive powers.

  procedure First_Derivative_Second_Order
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                cf2 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                pw2 : in Standard_Floating_Vectors.Vector;
                cf3 : out Standard_Complex_Vectors.Vector;
                pw3 : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes a Taylor series expansion of the Laurent homotopy
  --   using the constant coefficients of the power series solution,
  --   truncated after the second order, using first derivatives.
  --   Computes the smallest positive exponents of this evaluation.

  -- ON ENTRY :
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hdg      supports of the Laurent homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cf0      constant coefficients of a power series;
  --   cf1      coefficients corresponding to pw1;
  --   cf2      coefficients corresponding to pw2;
  --   pw1      leading exponents in the power series;
  --   pw2      second exponents in the power series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   cf3      coefficients corresponding to the pw3;
  --   pw3      smallest positive powers in the evaluated series.

  procedure First_Derivative_Second_Order
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                psm : out Standard_Floating_Vectors.Vector;
                csm : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes a Taylor series expansion of the Laurent homotopy
  --   using the constant coefficients of the power series solution,
  --   truncated after the second order, using first derivatives.
  --   Computes the smallest positive exponents of this evaluation.

  -- ON ENTRY :
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hdg      supports of the Laurent homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cff      coefficients of the power series solution;
  --   pwr      exponents of the power series solution;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   psm      smallest positive powers in the evaluated series;
  --   csm      coefficients corresponding to the smallest positive powers.

  procedure Second_Derivative_First_Order
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                cf2 : out Standard_Complex_Vectors.Vector;
                pw2 : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes a Taylor series expansion of the Laurent homotopy
  --   using the constant coefficients of the power series solution,
  --   truncated after the first order, using up to 2nd derivatives.
  --   Computes the smallest positive exponents of this evaluation.

  -- ON ENTRY :
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hdg      supports of the Laurent homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cf0      constant coefficients of a power series;
  --   cf1      coefficients corresponding to pw1;
  --   pw1      leading exponents in the power series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   cf2      coefficients corresponding to the smallest positive powers;
  --   pw2      smallest positive powers in the evaluated series.

  procedure Second_Derivative_First_Order
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                psm : out Standard_Floating_Vectors.Vector;
                csm : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes a Taylor series expansion of the Laurent homotopy
  --   using the constant coefficients of the power series solution,
  --   truncated after the first order, using up to 2nd derivatives.
  --   Computes the smallest positive exponents of this evaluation.

  -- ON ENTRY :
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hdg      supports of the Laurent homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cff      coefficients of the power series solution;
  --   pwr      exponents of the power series solution;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   psm      smallest positive powers in the evaluated series;
  --   csm      coefficients corresponding to the smallest positive powers.

  procedure Second_Derivative_Second_Order
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                psm : out Standard_Floating_Vectors.Vector;
                csm : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes a Taylor series expansion of the Laurent homotopy
  --   using the constant coefficients of the power series solution,
  --   truncated after the second order, using up to 2nd derivatives.
  --   Computes the smallest positive exponents of this evaluation.

  -- ON ENTRY :
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hdg      supports of the Laurent homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cff      coefficients of the power series solution;
  --   pwr      exponents of the power series solution;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   psm      smallest positive powers in the evaluated series;
  --   csm      coefficients corresponding to the smallest positive powers.

  procedure Third_Derivative_First_Order
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                cf2 : out Standard_Complex_Vectors.Vector;
                pw2 : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes a Taylor series expansion of the Laurent homotopy
  --   using the constant coefficients of the power series solution,
  --   truncated after the first order, using up to 3rd derivatives.
  --   Computes the smallest positive exponents of this evaluation.

  -- ON ENTRY :
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hdg      supports of the Laurent homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cf0      constant coefficients of a power series;
  --   cf1      coefficients corresponding to pw1;
  --   pw1      leading exponents in the power series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   cf2      coefficients corresponding to the smallest positive powers;
  --   pw2      smallest positive powers in the evaluated series.

  procedure Third_Derivative_First_Order
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                psm : out Standard_Floating_Vectors.Vector;
                csm : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes a Taylor series expansion of the Laurent homotopy
  --   using the constant coefficients of the power series solution,
  --   truncated after the first order, using up to 3rd derivatives.
  --   Computes the smallest positive exponents of this evaluation.

  -- ON ENTRY :
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hdg      supports of the Laurent homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cff      coefficients of the power series solution;
  --   pwr      exponents of the power series solution;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   psm      smallest positive powers in the evaluated series;
  --   csm      coefficients corresponding to the smallest positive powers.

-- INDEXED DERIVATIVES ON A LAURENT HOMOTOPY :

  procedure First_Order_Evaluation
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                difmax : in integer32;
                cf2 : out Standard_Complex_Vectors.Vector;
                pw2 : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes a Taylor series expansion of the Laurent homotopy
  --   using the constant coefficients of the power series solution,
  --   truncated after the first order, using up to some derivatives.
  --   Computes the smallest positive exponents of this evaluation.

  -- ON ENTRY :
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hdg      supports of the Laurent homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cf0      constant coefficients of a power series;
  --   cf1      coefficients corresponding to pw1;
  --   pw1      leading exponents in the power series;
  --   difmax   maximum sum of all derivative indices;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   psm      smallest positive powers in the evaluated series;
  --   csm      coefficients corresponding to the smallest positive powers.

end Double_Ordered_Evaluations;
