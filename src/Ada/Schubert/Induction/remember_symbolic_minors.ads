with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Matrices;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Matrices;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Matrices;
with Brackets;                          use Brackets;

package Remember_Symbolic_Minors is

-- DESCRIPTION :
--   This package manages a remember table for all maximal minors
--   of a symbolic n-by-k matrix, where n >= k.
--   The minors are represented as three different types of polynomials,
--   depending on the precision of the coefficients, which can be either
--   standard double, double double, or quad double.

  type Standard_Symbolic_Minors ( m : integer32 ) is record
    b : Array_of_Brackets(1..m);
    p : Standard_Complex_Poly_Systems.Poly_Sys(1..m);
  end record;

  type DoblDobl_Symbolic_Minors ( m : integer32 ) is record
    b : Array_of_Brackets(1..m);
    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..m);
  end record;

  type QuadDobl_Symbolic_Minors ( m : integer32 ) is record
    b : Array_of_Brackets(1..m);
    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..m);
  end record;

  function Number_of_Minors ( n,k : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of k-by-k minors of an n-by-k matrix, n >= k.

  function Create ( n,k : natural32;
                    x : Standard_Complex_Poly_Matrices.Matrix )
                  return Standard_Symbolic_Minors;
  function Create ( n,k : natural32;
                    x : DoblDobl_Complex_Poly_Matrices.Matrix )
                  return DoblDobl_Symbolic_Minors;
  function Create ( n,k : natural32;
                    x : QuadDobl_Complex_Poly_Matrices.Matrix )
                  return QuadDobl_Symbolic_Minors;

  -- DESCRIPTION :
  --   Returns a minor table for all k-by-k minors of a general symbolic
  --   n-by-k matrix.  The size m of the table is Number_of_Minors(n,k).

  -- REQUIRED : x'range(1) = 1..n and x'range(2) = 1..k.

  function Search ( t : Standard_Symbolic_Minors; b : Bracket )
                  return Standard_Complex_Polynomials.Poly;
  function Search ( t : DoblDobl_Symbolic_Minors; b : Bracket )
                  return DoblDobl_Complex_Polynomials.Poly;
  function Search ( t : QuadDobl_Symbolic_Minors; b : Bracket )
                  return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns either the null polynomial if there is no t.b(i) = b,
  --   or else returns the corresponding t.p(i), with b = t.b(i).

  procedure Write ( t : in Standard_Symbolic_Minors );
  procedure Write ( t : in DoblDobl_Symbolic_Minors );
  procedure Write ( t : in QuadDobl_Symbolic_Minors );

  -- DESCRIPTION :
  --   Writes the minor table to screen.

  procedure Query ( t : in Standard_Symbolic_Minors; k : in integer32 );
  procedure Query ( t : in DoblDobl_Symbolic_Minors; k : in integer32 );
  procedure Query ( t : in QuadDobl_Symbolic_Minors; k : in integer32 );

  -- DECRIPTION :
  --   Interactive routine, prompts the user for a k-bracket,
  --   and displays the result of a search in t for the given bracket.

  procedure Clear ( t : in out Standard_Symbolic_Minors );
  procedure Clear ( t : in out DoblDobl_Symbolic_Minors );
  procedure Clear ( t : in out QuadDobl_Symbolic_Minors );

  -- DESCRIPTION :
  --   Deallocation of all memory occupied by the remember table.

end Remember_Symbolic_Minors;
