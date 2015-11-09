with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Standard_Diagonal_Solutions is

-- DESCRIPTION :
--   The functions in this package manipulate solutions for use in the
--   diagonal homotopies in standard double precision.

-- OPERATIONS ON SOLUTIONS :

  function Product ( s1,s2 : Solution ) return Solution;
  function Product ( s1,s2 : Solution_List ) return Solution_List;

  -- DESCRIPTION :
  --   Returns the product of the solutions in s1 and s2.

  function Truncate ( s : Solution; n : integer32 ) return Solution;
  function Truncate ( s : Solution_List; n : integer32 ) return Solution_List;

  -- DESCRIPTION :
  --   Truncates the solution vectors, after the first n entries.

end Standard_Diagonal_Solutions;
