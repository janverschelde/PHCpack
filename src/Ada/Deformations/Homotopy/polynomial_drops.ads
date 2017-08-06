with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;

package Polynomial_Drops is

-- DESCRIPTION :
--   This package offers operations to drop a variable from a polynomial.
--   Dropping variable k simply means that the kth exponent is removed
--   from every term in the polynomial.

  function Drop ( t : Standard_Complex_Polynomials.Term; k : integer32 )
                return Standard_Complex_Polynomials.Term;
  function Drop ( t : Standard_Complex_Laurentials.Term; k : integer32 )
                return Standard_Complex_Laurentials.Term;
  function Drop ( t : DoblDobl_Complex_Polynomials.Term; k : integer32 )
                return DoblDobl_Complex_Polynomials.Term;
  function Drop ( t : DoblDobl_Complex_Laurentials.Term; k : integer32 )
                return DoblDobl_Complex_Laurentials.Term;
  function Drop ( t : QuadDobl_Complex_Polynomials.Term; k : integer32 )
                return QuadDobl_Complex_Polynomials.Term;
  function Drop ( t : QuadDobl_Complex_Laurentials.Term; k : integer32 )
                return QuadDobl_Complex_Laurentials.Term;

  -- DESCRIPTION :
  --   Removes the k-th variable from the term t.

  -- REQUIRED : k is in t.dg'range and t.dg(k) = 0.

  function Drop ( p : Standard_Complex_Polynomials.Poly; k : integer32 )
                return Standard_Complex_Polynomials.Poly;
  function Drop ( p : Standard_Complex_Laurentials.Poly; k : integer32 )
                return Standard_Complex_Laurentials.Poly;
  function Drop ( p : DoblDobl_Complex_Polynomials.Poly; k : integer32 )
                return DoblDobl_Complex_Polynomials.Poly;
  function Drop ( p : DoblDobl_Complex_Laurentials.Poly; k : integer32 )
                return DoblDobl_Complex_Laurentials.Poly;
  function Drop ( p : QuadDobl_Complex_Polynomials.Poly; k : integer32 )
                return QuadDobl_Complex_Polynomials.Poly;
  function Drop ( p : QuadDobl_Complex_Laurentials.Poly; k : integer32 )
                return QuadDobl_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Removes the k-th variable from the polynomial p.

  -- REQUIRED : k is in the range of variables.

  function Drop ( p : Standard_Complex_Poly_Systems.Poly_Sys; k : integer32 )
                return Standard_Complex_Poly_Systems.Poly_Sys;
  function Drop ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys; k : integer32 )
                return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Drop ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys; k : integer32 )
                return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Removes the k-th variable from the polynomial system p.
  --   Only those terms for which t.dg(k) = 0 in p remain!

  -- REQUIRED : k is in the range of variables.

  function Drop ( p : Standard_Complex_Laur_Systems.Laur_Sys; k : integer32 )
                return Standard_Complex_Laur_Systems.Laur_Sys;
  function Drop ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys; k : integer32 )
                return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Drop ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys; k : integer32 )
                return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Removes the k-th variable from the Laurent polynomial system p.
  --   Only those terms for which t.dg(k) = 0 in p remain!

  -- REQUIRED : k is in the range of variables.

  function Remove_Variable
             ( p : Standard_Complex_Laurentials.Poly; k : integer32 )
             return Standard_Complex_Laurentials.Poly;
  function Remove_Variable
             ( p : Standard_Complex_Laur_Systems.Laur_Sys; k : integer32 )
             return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Removes the kth variable from p, regardless what the exponent is.

end Polynomial_Drops;
