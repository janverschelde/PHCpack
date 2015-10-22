with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Brackets,Bracket_Monomials;         use Brackets,Bracket_Monomials;
with Standard_Bracket_Polynomials;       use Standard_Bracket_Polynomials;

package Straightening_Syzygies is

-- DESCRIPTION :
--   This package generates and uses the so-called van der Waerden syzygies
--   to straighten bracket polynomials w.r.t. the tableau order.

  function Laplace_Expansion ( n,d : natural32 ) return Bracket_Polynomial;

  -- DESCRIPTION :
  --   Returns the Laplace expansion of the determinant of an (n*n)-matrix,
  --   in terms of two blocks of size d and n-d respectively.
  --   The output is in terms of a quadratic bracket polynomial.
  --   The first entry of every bracket for the first block equals 0,
  --   to avoid messing up the structure of the Laplace expansion.

  function Straightening_Syzygy ( b1,b2 : Bracket ) return Bracket_Polynomial;

  -- DESCRIPTION :
  --   Returns the van der Waerden syzygy that can be used to straighten
  --   the nonstandard monomial b1*b2.  If b1*b2 is already standard, then
  --   the polynomial with only the monomial b1*b2 is returned.

  -- REQUIRED : b1 < b2.

  function Straightening_Syzygy ( b : Bracket_Monomial )
                                return Bracket_Polynomial;

  -- DESCRIPTION :
  --   Returns Straightening_Syzygy(b1,b2), with b = b1*b2.

  -- REQUIRED : b = b1*b2.

  function nonStandard_Monomials ( n,d : natural32 ) return Bracket_Polynomial;

  -- DESCRIPTION :
  --   Returns the polynomial of all quadratic nonStandard monomials,
  --   where d is the dimension of the brackets and n the number of 
  --   elements to choose from.

  generic
    with procedure Process ( s : in Bracket_Polynomial;
                             continue : out boolean );
  procedure Enumerate_Syzygies ( p : in Bracket_Polynomial );

  -- DESCRIPTION :
  --   Constructs the straightening syzygy s for every monomial in p
  --   and invokes the procedure Process with s as its argument.
  --   Enumeration stops when continue is set to false.

  function Straighten ( b1,b2 : Bracket ) return Bracket_Polynomial;

  -- DESCRIPTION :
  --   Returns a bracket polynomial, equivalent to b1*b2, that contains
  --   only standard monomials.  This is done by repeatively replacing the
  --   nonstandard monomials b1*b2 by the generated van der Waerden syzygies.

  function Straighten ( b : Bracket_Monomial ) return Bracket_Polynomial;
  function Straighten ( b : Bracket_Term ) return Bracket_Polynomial;
  function Straighten ( b : Bracket_Polynomial ) return Bracket_Polynomial;

  -- DESCRIPTION :
  --   Returns an equivalent bracket polynomial that contains only
  --   standard monomials.

end Straightening_Syzygies;
