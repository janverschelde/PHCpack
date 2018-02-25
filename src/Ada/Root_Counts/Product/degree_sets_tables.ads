with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Sets_of_Unknowns;                   use Sets_of_Unknowns;

package Degree_Sets_Tables is

-- DESCRIPTION :
--   This package provides facilities for computing generalized permanents,
--   based on a set structure.

-- DATASTRUCTURES :

  type Array_of_Sets is array ( integer32 range <> ) of Set;
  
  type Degree_Sets_Table ( n,m : integer32 ) is record
    s : Array_of_Sets(1..m);
    a : Matrix(1..n,1..m);
  end record;

-- CONSTRUCTOR :

  function Create return Degree_Sets_Table;

  -- DESCRIPTION :
  --   Selects the information from the package Set_Structure to create
  --   a degree set structure table.

-- PERMANENT COMPUTATIONS :

  function Permanent ( dst : Degree_Sets_Table ) return integer32;

  -- DESCRIPTION :
  --   Returns the generalized permanent, based on the set structure.

  function Matching_Permanent ( dst : Degree_Sets_Table ) return integer32;

  -- DESCRIPTION :
  --   Returns the generalized permanent, based on the set structure,
  --   applying the connection with maximum bipartite matching problem.

-- DESTRUCTOR :

  procedure Clear ( ase : in out Array_of_Sets );
  procedure Clear ( dst : in out Degree_Sets_Table );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the array of sets.

end Degree_Sets_Tables;
