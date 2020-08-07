with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Wrapped_Solution_Vectors is

-- DESCRIPTION :
--   The functions below wrap the coordinates of solutions with
--   the continuation parameter t, adding t as last coordinate.

  function Create ( xt : Standard_Complex_Vectors.Vector )
                  return Standard_Complex_Solutions.Solution;
  function Create ( xt : DoblDobl_Complex_Vectors.Vector )
                  return DoblDobl_Complex_Solutions.Solution;
  function Create ( xt : QuadDobl_Complex_Vectors.Vector )
                  return QuadDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Returns the solution representation of the vector xt,
  --   in double, double double, or quad double precision.
  --   The last coordinate of xt is the value of the continuation parameter,
  --   which on return is stored in the field t of the solution.

  function Create ( xt : Standard_Complex_Vectors.Vector ) 
                  return Standard_Complex_Solutions.Solution_List;
  function Create ( xt : DoblDobl_Complex_Vectors.Vector ) 
                  return DoblDobl_Complex_Solutions.Solution_List;
  function Create ( xt : QuadDobl_Complex_Vectors.Vector ) 
                  return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns the solution list representation of the vector xt,
  --   in double, double double, or quad double precision,
  --   wrapping the function which turns the vector xt into a solution.

  function Create ( xt : Standard_Complex_Solutions.Link_to_Solution )
                  return Standard_Complex_Solutions.Link_to_Solution;
  function Create ( xt : DoblDobl_Complex_Solutions.Link_to_Solution )
                  return DoblDobl_Complex_Solutions.Link_to_Solution;
  function Create ( xt : QuadDobl_Complex_Solutions.Link_to_Solution )
                  return QuadDobl_Complex_Solutions.Link_to_Solution;
  function Create ( xt : Standard_Complex_Solutions.Solution_List )
                  return Standard_Complex_Solutions.Solution_List;
  function Create ( xt : DoblDobl_Complex_Solutions.Solution_List )
                  return DoblDobl_Complex_Solutions.Solution_List;
  function Create ( xt : QuadDobl_Complex_Solutions.Solution_List )
                  return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns the solution list representation of the solution(s) xt,
  --   in double, double double, or quad double precision.

-- UPDATING SOLUTION LISTS :

  procedure Update ( sols : in out Standard_Complex_Solutions.Solution_List;
                     xtsols : in Standard_Complex_Solutions.Solution_List );
  procedure Update ( sols : in out DoblDobl_Complex_Solutions.Solution_List;
                     xtsols : in DoblDobl_Complex_Solutions.Solution_List );
  procedure Update ( sols : in out QuadDobl_Complex_Solutions.Solution_List;
                     xtsols : in QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Given in the list xtsols the coordinates of the solutions
  --   jointly with the value of continuation parameter,
  --   the list sols is updated, separating the continuation parameter t
  --   from the last coordinate of the solution vectors in xtsols.

end Wrapped_Solution_Vectors;
