with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Matrices;
with Standard_Complex_Poly_Matrices;
with DoblDobl_Complex_Poly_Matrices;
with QuadDobl_Complex_Poly_Matrices;
with Brackets;                          use Brackets;
with Standard_Bracket_Polynomials;      use Standard_Bracket_Polynomials;
with Standard_Bracket_Systems;          use Standard_Bracket_Systems;

package Symbolic_Schubert_Conditions is

-- DESCRIPTION :
--   This package produces the minor equations which arise from
--   Schubert intersection conditions.

  function General_Localization_Map
             ( n,k : in integer32 ) return Standard_Natural_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a general localization map for a k-plane in n-space.
  --   The entries in a localization map are either 0, 1, or 2,
  --   indicating whether the corresponding elements are respectively
  --   zero, one, or any variable.

  function Symbolic_Form_of_Plane
             ( n,k : integer32; locmap : Standard_Natural_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix;
  function Symbolic_Form_of_Plane
             ( n,k : integer32; locmap : Standard_Natural_Matrices.Matrix )
             return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Symbolic_Form_of_Plane
             ( n,k : integer32; locmap : Standard_Natural_Matrices.Matrix )
             return QuadDobl_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a symbolic form of a k-plane in n-space,
  --   following the localization pattern in locmap.

  function Number_of_Equations ( n,k,f,i : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of equations imposed on a k-plane in n-space
  --   that has to intersect an f-plane in an i-space.

  function Number_of_Equations ( n : natural32; b : Bracket ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of minor equations imposed on a k-plane
  --   by a general flag in n-space.  The k-bracket lists the
  --   dimensions of the spaces defined by the flag.

  generic
    with procedure Process ( c : in Bracket; continue : out boolean );
  procedure Enumerate_NotAbove ( n : in natural32; b : in Bracket );

  -- DESCRIPTION :
  --   Enumerates all brackets c that are not (c <= b).
  --   With each new bracket, the procedure Process is executed.
  --   The enumeration stops when the execution of Process results
  --   in the output parameter continue set to false.

  function Number_of_NotAbove ( n : natural32; b : Bracket ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of brackets that are not above b.

  procedure Explain_Equations
              ( n : in natural32; b : in Bracket; nq : out natural32 );

  -- DESCRIPTION :
  --   On return, nq = Number_of_Equations(n,b),
  --   messages are printed to screen, explaining the number of equations.

  function Flag_Minors ( n,k,f,i : natural32 ) return Bracket_Polynomial;

  -- DESCRIPTION :
  --   Returns a symbolic representation of the minor equations on a
  --   k-plane in n-space, imposed by f = b(i), for some bracket.
  --   The equations express Rank([ X | F(b(i)) ] = n + f - i = r.
  --   The minors are represented as products of brackets of size r+1,
  --   indicating a row selection (if first bracket in monomial) or
  --   the corresponding column selection (bracket following first).

  function Flag_Minors ( n : natural32; b : bracket ) return Bracket_System;

  -- DESCRIPTION :
  --   Returns a system of k bracket polynomials (some may be null),
  --   where k = b'last, expressing the intersection conditions on
  --   a k-plane in n-space imposed by the bracket b.

  procedure Flag_Minor_Polynomials
              ( p : in Bracket_Polynomial; bs : in out Bracket_System;
                ind : in out integer32 );

  -- DESCRIPTION :
  --   Stores the minors in the bracket polynomial p into the system bs.

  -- ON ENTRY :
  --   p        the first monomial in each term are row indices,
  --            and the next monomials are column indices;
  --   ind      current index of last defined entry in bs;
  --   bs       bs is filled up to bs(ind), initialize at first to 0.

  -- ON RETURN :
  --   ind      updated current index of last defined entry in bs;
  --   bs       contains all expanded equations for all terms in p.

  function Flag_Minor_System
              ( m : natural32; fm : Bracket_Polynomial ) return Bracket_System;

  -- DESCRIPTION :
  --   Returns a system of m equations, derived from the minors in fm,
  --   calling the Flag_Minor_Polynomials procedure.

  function Flag_Minor_System
              ( m : natural32; bs : Bracket_System ) return Bracket_System;

  -- DESCRIPTION :
  --   Returns a system of m equations, derived from the minors in bs,
  --   where m equals Number_of_Equations(n,b), for some bracket b
  --   representing a k-plane in n-space.
  --   Every equation in the system on return consists of one quadratic
  --   bracket monomial: the rows and columns that determine the minor.

end Symbolic_Schubert_Conditions;
