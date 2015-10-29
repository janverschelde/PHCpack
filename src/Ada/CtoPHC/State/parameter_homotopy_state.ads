with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;

package Parameter_Homotopy_State is

-- DESCRIPTION :
--   This packages stores the state of the parameter homotopy:
--   1) the number of equations;
--   2) the number of variables, which includes the parameters;
--   3) the number of parameters;
--   4) the index set of which variables are parameters.

-- CONSTRUCTORS :

  procedure Set_Number_of_Equations ( n : in integer32 );

  -- DESCRIPTION :
  --   Sets the number of equations to the value of n.

  procedure Set_Number_of_Variables ( n : in integer32 );

  -- DESCRIPTION :
  --   Sets the number of variables to the value of n.

  procedure Set_Number_of_Parameters ( n : in integer32 );

  -- DESCRIPTION :
  --   Sets the number of parameters to the value of n.

  procedure Set_Indices ( idx : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Defines the indices to the parameters in the sweep homotopy.

-- SELECTORS :

  function Get_Number_of_Equations return integer32;

  -- DESCRIPTION :
  --   Returns the number of equations.

  function Get_Number_of_Variables return integer32;

  -- DESCRIPTION :
  --   Returns the number of variables.

  function Get_Number_of_Parameters return integer32;

  -- DESCRIPTION :
  --   Returns the number of parameters.

  function Get_Indices return Standard_Integer_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the indices to the parameters.

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   Resets all numbers to zero and clears the index vector
  --   to the parameter variables.

end Parameter_Homotopy_State;
