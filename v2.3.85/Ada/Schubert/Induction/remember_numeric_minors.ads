with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Brackets;                          use Brackets;

package Remember_Numeric_Minors is

-- DESCRIPTION :
--   This package manages a remember table for all maximal minors
--   of a numeric n-by-k matrix, where n >= k.

  type Numeric_Minor_Table ( m : integer32 ) is record
    b : Array_of_Brackets(1..m);
    v : Vector(1..m);
  end record;

  function Number_of_Minors ( n,k : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of k-by-k minors of an n-by-k matrix, n >= k.

  function Create ( n,k : natural32; x : Matrix ) return Numeric_Minor_Table;

  -- DESCRIPTION :
  --   Returns a minor table for all k-by-k minors of an n-by-k matrix.
  --   The size m of the table is Number_of_Minors(n,k).

  -- REQUIRED : x'range(1) = 1..n and x'range(2) = 1..k.

  function Search ( t : Numeric_Minor_Table;
                    b : Bracket ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the corresponding t.p(i), with b = t.b(i).
  --   The default value is zero.

  procedure Write ( t : in Numeric_Minor_Table );

  -- DESCRIPTION :
  --   Writes the remember table t to screen.

  procedure Query ( t : in Numeric_Minor_Table; k : in integer32 );

  -- DECRIPTION :
  --   Interactive routine, prompts the user for a k-bracket,
  --   and displays the result of a search in t for the given bracket.

  procedure Clear ( t : in out Numeric_Minor_Table );

  -- DESCRIPTION :
  --   Deallocation of all memory occupied by the remember table.

end Remember_Numeric_Minors;
