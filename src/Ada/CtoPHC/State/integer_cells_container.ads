with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;    use Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Integer_Mixed_Subdivisions;        use Integer_Mixed_Subdivisions;

package Integer_Cells_Container is

-- DESCRIPTION :
--   This package provides a container for regular mixed-cell configurations,
--   with integer valued lifting function, designed for the interface with C.

-- CREATORS :

  procedure Initialize_Supports ( nbr : in natural32 );

  -- DESCRIPTION :
  --   Initializes the container to hold as many different supports
  --   as the value of nbr.

  procedure Initialize 
              ( mixture : in Standard_Integer_Vectors.Link_to_Vector;
                lifting : in Link_to_Array_of_Lists;
                mcc : in Mixed_Subdivision );

  -- DESCRIPTION :
  --   Initializes the container with a mixed-cell configuration.
  --   On entry are the type of mixture, lifted supports and the cells.

  procedure Initialize
              ( mixture : in Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Initializes the container with a type of mixture.

  procedure Initialize ( lifting : in Link_to_Array_of_Lists );

  -- DESCRIPTION :
  --   Initializes the container with lifted supports.

  procedure Initialize ( mcc : in Mixed_Subdivision );

  -- DESCRIPTION :
  --   Initializes the container with a mixed-cell configuration.

  procedure Make_Subdivision;

  -- DESCRIPTION :
  --   After the initialization of the lifted supports,
  --   the mixed cell configuration is computed.
  --   If no mixture is defined, then the supports are
  --   assumed to be fully mixed.

  procedure Generate_Random_Standard_Coefficient_System;
  procedure Generate_Random_DoblDobl_Coefficient_System;
  procedure Generate_Random_QuadDobl_Coefficient_System;

  -- DESCRIPTION :
  --   Generates a random coefficient system for polyhedral continuation.

  -- REQUIRED :
  --   The type of mixture and lifted supports are initialized.

  procedure Initialize_Random_Standard_Coefficient_System
              ( q : in Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Initialize_Random_DoblDobl_Coefficient_System
              ( q : in DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Initialize_Random_QuadDobl_Coefficient_System
              ( q : in QuadDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Initializes the random coefficient system with the given system q.

  -- REQUIRED :
  --   The given system must match the type of mixture and supports.

  procedure Standard_Polyhedral_Homotopy;

  -- DESCRIPTION :
  --   Using the lifting and random coefficient system, sets up the
  --   data for polyhedral homotopy continuation in standard precision.

  procedure DoblDobl_Polyhedral_Homotopy;

  -- DESCRIPTION :
  --   Using the lifting and random coefficient system, sets up the
  --   data for polyhedral homotopy continuation in double double precision.

  procedure QuadDobl_Polyhedral_Homotopy;

  -- DESCRIPTION :
  --   Using the lifting and random coefficient system, sets up the
  --   data for polyhedral homotopy continuation in double double precision.

-- SELECTORS :

  function Length return natural32;

  -- DESCRIPTION :
  --   Returns the number of cells in the container.

  function Dimension return natural32;

  -- DESCRIPTION :
  --   Returns 0 if the container is empty, otherwise the
  --   dimension of lifted points in the cells are returned.
  --   This dimension is the length of the lifted points.

  function Dimension_of_Supports return natural32;

  -- DESCRIPTION :
  --   Returns 0 if the lifted supports are not defined.
  --   Otherwise, returns the dimension of the supports
  --   before the lifting was applied.

  function Type_of_Mixture return Standard_Integer_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the type of mixture of the mixed-cell configuration.

  function Lifted_Supports return Link_to_Array_of_Lists;

  -- DESCRIPTION :
  --   Returns the lifted supports for the mixed-cell configuration.

  procedure Retrieve ( k : in natural32; mic : out Mixed_Cell;
                       fail : out boolean );

  -- DESCRIPTION : 
  --   Returns the k-th cell in the subdivision in mic if not fail.

  procedure Retrieve_Mixed_Cell
             ( k : in natural32; fail : out boolean;
               cnt,lab : out Standard_Integer_Vectors.Link_to_Vector;
               normal : out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Returns the k-th cell in the subdivision if not fail,
  --   in a format most suitable for the C gateway interface.
  --   See the Append_Mixed_Cell below.

  -- ON ENTRY :
  --   k       index number to a mixed cell in the subdivision.

  -- ON RETURN :
  --   fail    true if there is no k-th cell;
  --   cnt     cnt(i) contains the number of points in the i-th support;
  --   lab     lab(i) is the label for the i-th points in the cell;
  --   normal  is the inner normal to the mixed cell.

  function Retrieve return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Returns the mixed subdivision stored in the container.

  function Retrieve_Random_Standard_Coefficient_System
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Retrieve_Random_Standard_Coefficient_System
             return Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
  function Retrieve_Random_DoblDobl_Coefficient_System
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Retrieve_Random_DoblDobl_Coefficient_System
             return DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
  function Retrieve_Random_QuadDobl_Coefficient_System
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Retrieve_Random_QuadDobl_Coefficient_System
             return QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Returns the random coefficient system in the container.

  function Retrieve_Standard_Start_Solution
             ( k,i : natural32 )
             return Standard_Complex_Solutions.Link_to_Solution;
  function Retrieve_Standard_Target_Solution
             ( k,i : natural32 )
             return Standard_Complex_Solutions.Link_to_Solution;

  -- DESCRIPTION :
  --   Returns a pointer to the i-th solution of the k-th cell,
  --   as computed in standard precision.

  -- REQUIRED : k <= Length.

  function Retrieve_DoblDobl_Start_Solution
             ( k,i : natural32 )
             return DoblDobl_Complex_Solutions.Link_to_Solution;
  function Retrieve_DoblDobl_Target_Solution
             ( k,i : natural32 )
             return DoblDobl_Complex_Solutions.Link_to_Solution;

  -- DESCRIPTION :
  --   Returns a pointer to the i-th solution of the k-th cell,
  --   as computed in double double precision.

  -- REQUIRED : k <= Length.

  function Retrieve_QuadDobl_Start_Solution
             ( k,i : natural32 )
             return QuadDobl_Complex_Solutions.Link_to_Solution;
  function Retrieve_QuadDobl_Target_Solution
             ( k,i : natural32 )
             return QuadDobl_Complex_Solutions.Link_to_Solution;

  -- DESCRIPTION :
  --   Returns a pointer to the i-th solution of the k-th cell,
  --   as computed in quad double precision.

  -- REQUIRED : k <= Length.

-- CONSTRUCTORS :

  function Append_to_Support
             ( k : in natural32;
               x : in Standard_Integer_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Adds the point x to the k-th lifted support.
  --   The return value indicates failure and is true if k is out of
  --   the range of the type of mixture, otherwise false is returned.
  
  -- REQUIRED : the type of mixture is initialized.

  procedure Append ( mic : in Mixed_Cell );

  -- DESCRIPTION :
  --   Appends the cell to the current mixed-cell configuration.

  procedure Append_Mixed_Cell
             ( cnt,lab : in Standard_Integer_Vectors.Vector;
               normal : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Creates a mixed cell defined by the input parameters
  --   and appends this cell to the current mixed-cell configuration.
  --   See the Retrieve_Mixed_Cell above.

  -- ON ENTRY :
  --   cnt     cnt(i) contains the number of points in the i-th support;
  --   lab     lab(i) is the label for the i-th points in the cell;
  --   normal  is the inner normal to the mixed cell.

  procedure Solve_Standard_Start_System
              ( k : in natural32; mv : out natural32 );

  -- DESCRIPTION :
  --   Solves the start system defined by the k-th mixed cell in
  --   standard precision.  Returns in "mv" the number of solutions,
  --   which should equal the mixed volume of the cell.

  -- REQUIRED :
  --   The standard polyhedral homotopy has been created.

  procedure Solve_DoblDobl_Start_System
              ( k : in natural32; mv : out natural32 );

  -- DESCRIPTION :
  --   Solves the start system defined by the k-th mixed cell in
  --   double double precision.  Returns in "mv" the number of solutions,
  --   which should equal the mixed volume of the cell.

  -- REQUIRED :
  --   The dobldobl polyhedral homotopy has been created.

  procedure Solve_QuadDobl_Start_System
              ( k : in natural32; mv : out natural32 );

  -- DESCRIPTION :
  --   Solves the start system defined by the k-th mixed cell in
  --   quad double precision.  Returns in "mv" the number of solutions,
  --   which should equal the mixed volume of the cell.

  -- REQUIRED :
  --   The quaddobl polyhedral homotopy has been created.

  procedure Track_Standard_Solution_Path ( k,i,otp : in natural32 );

  -- DESCRIPTION :
  --   Tracks the path originating at the i-th solution of the k-th cell,
  --   in standard double precision.

  -- REQUIRED : 
  --   The start system corresponding to the k-th mixed cell is solved.

  -- ON ENTRY :
  --   k       number of a mixed cell in the configuration;
  --   i       index to a solution of a binomial start system;
  --   otp     output level, no output if zero, otherwise the number
  --           confirms to the output code of path trackers.

  -- ON RETURN :
  --   A target solution corresponding the k-th cell is added.

  procedure Track_DoblDobl_Solution_Path ( k,i,otp : in natural32 );

  -- DESCRIPTION :
  --   Tracks the path originating at the i-th solution of the k-th cell,
  --   in double double precision.

  -- REQUIRED : 
  --   The start system corresponding to the k-th mixed cell is solved.

  -- ON ENTRY :
  --   k       number of a mixed cell in the configuration;
  --   i       index to a solution of a binomial start system;
  --   otp     output level, no output if zero, otherwise the number
  --           confirms to the output code of path trackers.

  -- ON RETURN :
  --   A target solution corresponding the k-th cell is added.

  procedure Track_QuadDobl_Solution_Path ( k,i,otp : in natural32 );

  -- DESCRIPTION :
  --   Tracks the path originating at the i-th solution of the k-th cell,
  --   in quad double precision.

  -- REQUIRED : 
  --   The start system corresponding to the k-th mixed cell is solved.

  -- ON ENTRY :
  --   k       number of a mixed cell in the configuration;
  --   i       index to a solution of a binomial start system;
  --   otp     output level, no output if zero, otherwise the number
  --           confirms to the output code of path trackers.

  -- ON RETURN :
  --   A target solution corresponding the k-th cell is added.

  function Mixed_Volume return natural32;

  -- DESCRIPTION :
  --   If the mixture and lifted supports are defined,
  --   then the mixed volume of the supports is returned.

-- DESTRUCTORS :

  procedure Clear_Cell_Data;

  -- DESCRIPTION :
  --   Clears the cells in the container.

  procedure Clear_Standard_Data;

  -- DESCRIPTION :
  --   Clears data for polyhedral homotopies in standard precision.

  procedure Clear_DoblDobl_Data;

  -- DESCRIPTION :
  --   Clears data for polyhedral homotopies in double double precision.

  procedure Clear;

  -- DESCRIPTION :
  --   Clears all data stored in the container.

end Integer_Cells_Container;
