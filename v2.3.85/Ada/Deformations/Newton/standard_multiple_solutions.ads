with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Deflation_Trees;

package Standard_Multiple_Solutions is

-- DESCRIPTION :
--   This package offers tools to determine the multiplicity 
--   after deflation, simply by comparing clustered solutions.

  function Equal ( n : natural32; tol : double_float;
                   s1,s2 : Vector ) return boolean;

  -- DESCRIPTION :
  --   Compares the first n components of the two vectors.
  --   Returns true if none of the first n components differs
  --   in absolute value more than the given tolerance.

  procedure Set_Multiplicity
              ( sols : in out Solution_List; s : in Solution;
                tol : in double_float; n,m : in natural32 );

  -- DESCRIPTION :
  --   Every solution in sols close to s within the given tolerance tol
  --   will be given the multiplicity m.

  function Number_of_Occurrences 
              ( sols : Solution_List; s : Solution;
                tol : in double_float; n : in natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of times the solution s occurs in the list.

  function Is_In ( sols : Solution_List; v : Vector;
                   tol : double_float; n : natural32 ) return boolean;

  -- DESCRIPTION :
  --   Return true if the first n components of the solution vector v
  --   belong to the list sols within the given tolerance tol.

  procedure Remove_Duplicates
              ( sols : in out Solution_List;
                tol : in double_float; n : in natural32 );

  -- DESCRIPTION :
  --   Removes all duplicates from the list sols.

  procedure Merge_Multiple_Solutions
              ( sols : in out Solution_List; tol : in double_float );

  -- DESCRIPTION :
  --   Solutions whose vectors are within the tolerance are merged.
  --   When two solution merge, their multiplicity is taken as the
  --   maximum of the two solutions.

  procedure Compute_Multiplicities
              ( sols : in out Solution_List;
                tol : in double_float; n : in natural32 );

  -- DESCRIPTION :
  --   Sets the multiplicity for every solution in the list,
  --   grouping the solutions according to their clusters,
  --   using the tolerance tol as cluster radius.

  procedure Compute_Multiplicities
              ( nd : in out Standard_Deflation_Trees.Node;
                tol : in double_float; n : in natural32 );

  -- DESCRIPTION :
  --   Computes the multiplicities of the solution lists in the tree.
  --   Duplicate entries in the solution lists are removed.

end Standard_Multiple_Solutions;
