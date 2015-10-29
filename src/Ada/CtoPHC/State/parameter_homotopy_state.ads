with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;

package Parameter_Homotopy_State is

-- DESCRIPTION :
--   This packages stores the state of the parameter homotopy:
--   1) the number of equations;
--   2) the number of variables, which includes the parameters;
--   3) the number of parameters;
--   4) the index set of which variables are parameters.
--   In addition, the start and target values for the parameters are stored,
--   in standard double, double double, or quad double precision.

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

  -- REQUIRED : idx'range goes from 1 till the number of parameters.

  procedure Set_Start ( v : in Standard_Complex_Vectors.Vector ); 
  procedure Set_Start ( v : in DoblDobl_Complex_Vectors.Vector ); 
  procedure Set_Start ( v : in QuadDobl_Complex_Vectors.Vector ); 

  -- DESCRIPTION :
  --   Sets the start values for the parameters,
  --   in standard double, double double, or quad double precision.

  -- REQUIRED : v'range is the range of the indices to the parameter.

  procedure Set_Target ( v : in Standard_Complex_Vectors.Vector ); 
  procedure Set_Target ( v : in DoblDobl_Complex_Vectors.Vector ); 
  procedure Set_Target ( v : in QuadDobl_Complex_Vectors.Vector ); 

  -- DESCRIPTION :
  --   Sets the target values for the parameters,
  --   in standard double, double double, or quad double precision.

  -- REQUIRED : v'range is the range of the indices to the parameter.

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

  function Get_Start return Standard_Complex_Vectors.Link_to_Vector; 
  function Get_Start return DoblDobl_Complex_Vectors.Link_to_Vector; 
  function Get_Start return QuadDobl_Complex_Vectors.Link_to_Vector; 

  -- DESCRIPTION :
  --   Gets the start values for the parameters,
  --   in standard double, double double, or quad double precision.

  function Get_Target return Standard_Complex_Vectors.Link_to_Vector; 
  function Get_Target return DoblDobl_Complex_Vectors.Link_to_Vector; 
  function Get_Target return QuadDobl_Complex_Vectors.Link_to_Vector; 

  -- DESCRIPTION :
  --   Gets the target values for the parameters,
  --   in standard double, double double, or quad double precision.

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   Resets all numbers to zero and clears the index vector
  --   to the parameter variables.

end Parameter_Homotopy_State;
