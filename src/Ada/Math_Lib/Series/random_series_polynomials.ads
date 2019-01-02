with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Cseries_Polynomials;
with Standard_Cseries_Poly_Systems;
with DoblDobl_Cseries_Polynomials;
with DoblDobl_Cseries_Poly_Systems;
with QuadDobl_Cseries_Polynomials;
with QuadDobl_Cseries_Poly_Systems;

package Random_Series_Polynomials is

-- DESCRIPTION :
--   A random series polynomial is a sum of random monomials, 
--   with random coefficients series up to a fixed degree.

  function Standard_Random_Term
             ( nbvar,degtrm,degcff : natural32 )
             return Standard_CSeries_Polynomials.Term;
  function DoblDobl_Random_Term
             ( nbvar,degtrm,degcff : natural32 )
             return DoblDobl_CSeries_Polynomials.Term;
  function QuadDobl_Random_Term
             ( nbvar,degtrm,degcff : natural32 )
             return QuadDobl_CSeries_Polynomials.Term;

  -- DESCRIPTION :
  --   Returns a random monomial in several variables with 
  --   complex numbers in the series coefficients,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   nbvar     number of variables in the term;
  --   degtrm    upper bound on the degree of the term;
  --   degcff    degree of the coefficient series.

  function Standard_Random_Polynomial
             ( nbvar,nbterms,degpol,degcff : natural32 )
             return Standard_CSeries_Polynomials.Poly;
  function DoblDobl_Random_Polynomial
             ( nbvar,nbterms,degpol,degcff : natural32 )
             return DoblDobl_CSeries_Polynomials.Poly;
  function QuadDobl_Random_Polynomial
             ( nbvar,nbterms,degpol,degcff : natural32 )
             return QuadDobl_CSeries_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns a random polynomial in several variables with
  --   complex numbers in the series coefficients,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   nbvar     number of variables in each term;
  --   nbterms   number of terms in the polynomial on return;
  --   degpol    upper bound on the degree of each term;
  --   degcff    degree of each coefficient series.

  function Standard_Random_System
             ( nbequ,nbvar,nbterms,degpol,degcff : natural32 )
             return Standard_CSeries_Poly_Systems.Poly_Sys;
  function DoblDobl_Random_System
             ( nbequ,nbvar,nbterms,degpol,degcff : natural32 )
             return DoblDobl_CSeries_Poly_Systems.Poly_Sys;
  function QuadDobl_Random_System
             ( nbequ,nbvar,nbterms,degpol,degcff : natural32 )
             return QuadDobl_CSeries_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns a random polynomial system in several variables with
  --   complex numbers in the series coefficients,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   nbequ     number of polynomials in the system;
  --   nbvar     number of variables in each term;
  --   nbterms   number of terms in the polynomial on return;
  --   degpol    upper bound on the degree of each term;
  --   degcff    degree of each coefficient series.

end Random_Series_Polynomials;
