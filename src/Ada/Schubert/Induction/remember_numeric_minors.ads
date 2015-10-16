with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with Brackets;                          use Brackets;

package Remember_Numeric_Minors is

-- DESCRIPTION :
--   This package manages a remember table for all maximal minors
--   of a numeric n-by-k matrix, where n >= k, computed either in
--   standard double, double double, or quad double precision.

  type Standard_Numeric_Minors ( m : integer32 ) is record
    b : Array_of_Brackets(1..m);
    v : Standard_Complex_Vectors.Vector(1..m);
  end record;

  type DoblDobl_Numeric_Minors ( m : integer32 ) is record
    b : Array_of_Brackets(1..m);
    v : DoblDobl_Complex_Vectors.Vector(1..m);
  end record;

  type QuadDobl_Numeric_Minors ( m : integer32 ) is record
    b : Array_of_Brackets(1..m);
    v : QuadDobl_Complex_Vectors.Vector(1..m);
  end record;

  function Number_of_Minors ( n,k : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of k-by-k minors of an n-by-k matrix, n >= k.

  function Create ( n,k : natural32;
                    x : Standard_Complex_Matrices.Matrix )
                  return Standard_Numeric_Minors;
  function Create ( n,k : natural32;
                    x : DoblDobl_Complex_Matrices.Matrix )
                  return DoblDobl_Numeric_Minors;
  function Create ( n,k : natural32;
                    x : QuadDobl_Complex_Matrices.Matrix )
                  return QuadDobl_Numeric_Minors;

  -- DESCRIPTION :
  --   Returns a minor table for all k-by-k minors of an n-by-k matrix.
  --   The size m of the table is Number_of_Minors(n,k).

  -- REQUIRED : x'range(1) = 1..n and x'range(2) = 1..k.

  function Search ( t : Standard_Numeric_Minors; b : Bracket )
                  return Standard_Complex_Numbers.Complex_Number;
  function Search ( t : DoblDobl_Numeric_Minors; b : Bracket )
                  return DoblDobl_Complex_Numbers.Complex_Number;
  function Search ( t : QuadDobl_Numeric_Minors; b : Bracket )
                  return QuadDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns the corresponding t.p(i), with b = t.b(i).
  --   The default value is zero.

  procedure Write ( t : in Standard_Numeric_Minors );
  procedure Write ( t : in DoblDobl_Numeric_Minors );
  procedure Write ( t : in QuadDobl_Numeric_Minors );

  -- DESCRIPTION :
  --   Writes the remember table t to screen.

  procedure Query ( t : in Standard_Numeric_Minors; k : in integer32 );
  procedure Query ( t : in DoblDobl_Numeric_Minors; k : in integer32 );
  procedure Query ( t : in QuadDobl_Numeric_Minors; k : in integer32 );

  -- DECRIPTION :
  --   Interactive routine, prompts the user for a k-bracket,
  --   and displays the result of a search in t for the given bracket.

  procedure Clear ( t : in out Standard_Numeric_Minors );
  procedure Clear ( t : in out DoblDobl_Numeric_Minors );
  procedure Clear ( t : in out QuadDobl_Numeric_Minors );

  -- DESCRIPTION :
  --   Deallocation of all memory occupied by the remember table.

end Remember_Numeric_Minors;
