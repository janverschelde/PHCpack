with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Cells_Interface is

-- DESCRIPTION :
--   The functions below provide an interface to the mixed volume
--   computations and the polyhedral homotopies.

  function Cells_Read_Floating_Mixed_Cells
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a regular mixed cell configuration,
  --   induced by a floating-point lifting function,
  --   and stores the mixed cells in the container.

  function Cells_Write_Floating_Mixed_Cells
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the mixed cells stored in the container
  --   for configurations induced by floating-point lifting.

  function Cells_Number_of_Floating_Mixed_Cells
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of mixed cells stored in the container,
  --   induced by a floating-point lifting.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the number of mixed cells in a configuration
  --           induced by a floating-point lifting.

  function Cells_Dimension_of_Floating_Mixed_Cells
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the dimension of the mixed cells stored in the container,
  --   induced by a floating-point lifting.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the dimension of the mixed cells in a configuration
  --           induced by a floating-point lifting.

  function Cells_Get_Floating_Mixture
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of different supports and 
  --   the number of occurrences of each different support,
  --   for the points lifted by floating-point values.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the number of different supports;
  --   b       as many integers as the number of different supports,
  --           with the number of occurrences of each support.

  function Cells_Floating_Supports_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of points in each lifted support set,
  --   lifted with floating-point values.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the number of lifted supports;
  --   b       as many integers as the value of a[0],
  --           with the size of each support set.

  function Cells_Get_Floating_Support_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns point in the support, lifted with a floating-point value.

  -- ON ENTRY :
  --   a       index to the support;
  --   b       index to a point in the support;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       coordinates of the lifted point.

  function Cells_Floating_Normal
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the inner normal to a mixed cell
  --   with floating-point lifting values.

  -- ON ENTRY :
  --   a       the index to a mixed cell;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       coordinates of the inner normal.

  function Cells_Floating_Cell_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of points in a mixed cell
  --   with floating-point lifting values.

  -- ON ENTRY :
  --   a       the index to a mixed cell;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       contains the number of points in each support set
  --           of the mixed cell.

  function Cells_Get_Floating_Cell_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns a lifted point in a mixed cell
  --   with floating-point lifting values.

  -- ON ENTRY :
  --   a       the index to a mixed cell;
  --   b       in b[0] is the support index;
  --           in b[1] is the index of the point in the support;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       coordinates of the lifted cell point.

  function Cells_Floating_Mixed_Volume
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the mixed volume of a cell with floating-point lifting.

  -- ON ENTRY :
  --   a       the index to a mixed cell;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the mixed volume of the cell.

  function Cells_Set_Floating_Mixture
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the number of different supports and 
  --   the occurrences of each different supports,
  --   for mixed cells induced by floating-point lifting.

  -- ON ENTRY :
  --   a       the number of different supports;
  --   b       as many integers as the value of a[0],
  --           with the occurrence of each different support;
  --   vrblvl  is the verbose level.

  function Cells_Add_Floating_Support_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Adds a point lifted with a floating-point value to a support.

  -- ON ENTRY :
  --   a       in a[0] is the index of the support,
  --   b       in b[0] is the size of the point;
  --   c       the coordinates of the point;
  --   vrblvl  is the verbose level.

  function Cells_Add_Floating_Mixed_Cell
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Appends a mixed cell with floating-point lifting values.

  -- ON ENTRY :
  --   a       in a[0] is the number R of different supports,
  --           in a[1] is the dimension of the lifted points,
  --           in a[2] is the size of b;
  --   b       in b[0] is the total number of points in the cell,
  --           in b[k] is the number of points in the k-th support,
  --           in b[1+R+k] is the label for the k-th point in the cell;
  --   c       the coordinates for the inner normal to the cell;
  --   vrblvl  is the verbose level.

  function Cells_Get_Floating_Mixed_Cell
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Retrieves a mixed cell with floating-point lifting values.

  -- ON ENTRY :
  --   a       in a[0] is the index to the mixed cell;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the total number of points in the cell,
  --           in b[k] is the number of points in the k-th support,
  --           in b[1+R+k] is the label for the k-th point in the cell;
  --   c       the coordinates for the inner normal to the cell.

  function Cells_Make_Standard_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Generates a random coefficient polynomial system in double
  --   precision using the mixtures and supports in the cells container.
  --   The system is stored in the cells container.

  function Cells_Make_DoblDobl_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Generates a random coefficient polynomial system in double double
  --   precision using the mixtures and supports in the cells container.
  --   The system is stored in the cells container.

  function Cells_Make_QuadDobl_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Generates a random coefficient polynomial system in quad double
  --   precision using the mixtures and supports in the cells container.
  --   The system is stored in the cells container.

  function Cells_Read_Standard_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a coefficient polynomial system
  --   in double precision and stores it in the container.

  function Cells_Read_DoblDobl_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a coefficient polynomial system
  --   in double double precision and stores it in the container.

  function Cells_Read_QuadDobl_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a coefficient polynomial system
  --   in quad double precision and stores it in the container.

  function Cells_Write_Standard_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the coefficient system in double precision
  --   to standard output or to the defined output file.

  function Cells_Write_DoblDobl_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the coefficient system in double double precision
  --   to standard output or to the defined output file.

  function Cells_Write_QuadDobl_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the coefficient system in quad double precision
  --   to standard output or to the defined output file.

  function Cells_Standard_System_into_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the coefficient system in double precision
  --   into the systems container.

  function Cells_DoblDobl_System_into_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the coefficient system in double double precision
  --   into the systems container.

  function Cells_QuadDobl_System_into_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the coefficient system in quad double precision
  --   into the systems container.

  function Cells_Standard_System_from_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container for double precision
  --   into the cells container.

  function Cells_DoblDobl_System_from_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container for double double precision
  --   into the cells container.

  function Cells_QuadDobl_System_from_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container for quad double precision
  --   into the cells container.

  function Cells_Standard_Polyhedral_Homotopy
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a polyhedral homotopy to solve a random system
  --   in double precision.

  function Cells_DoblDobl_Polyhedral_Homotopy
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a polyhedral homotopy to solve a random system
  --   in double double precision.

  function Cells_QuadDobl_Polyhedral_Homotopy
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a polyhedral homotopy to solve a random system
  --   in quad double precision.

  function Cells_Standard_Start_Solve 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Solves the start system which corresponds to a mixed cell,
  --   in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index to the mixed cells;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the number of solutions.

  function Cells_DoblDobl_Start_Solve
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Solves the start system which corresponds to a mixed cell,
  --   in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index to the mixed cells;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the number of solutions.

  function Cells_QuadDobl_Start_Solve
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Solves the start system which corresponds to a mixed cell,
  --   in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index to the mixed cells;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the number of solutions.

  function Cells_Standard_Track_One_Path
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Tracks one solution path in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index to a mixed cell;
  --   b       in b[0] is the index to the start solution;
  --           in b[1] is the output code for the trackers;
  --   vrblvl  is the verbose level.

  function Cells_DoblDobl_Track_One_Path
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Tracks one solution path in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index to a mixed cell;
  --   b       in b[0] is the index to the start solution;
  --           in b[1] is the output code for the trackers;
  --   vrblvl  is the verbose level.

  function Cells_QuadDobl_Track_One_Path
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Tracks one solution path in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index to a mixed cell;
  --   b       in b[0] is the index to the start solution;
  --           in b[1] is the output code for the trackers;
  --   vrblvl  is the verbose level.

  function Cells_Standard_TarSol_into_Container
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies a target solution in double precision into the container.

  -- ON ENTRY :
  --   a       index to a mixed cell;
  --   b       index to the i-th solution of the initial form system
  --           with that mixed cell;
  --   vrblvl  is the verbose level.

  function Cells_DoblDobl_TarSol_into_Container
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies a target solution in double double precision
  --   into the container.

  -- ON ENTRY :
  --   a       index to a mixed cell;
  --   b       index to the i-th solution of the initial form system
  --           with that mixed cell;
  --   vrblvl  is the verbose level.

  function Cells_QuadDobl_TarSol_into_Container
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies a target solution in quad double precision
  --   into the container.

  -- ON ENTRY :
  --   a       index to a mixed cell;
  --   b       index to the i-th solution of the initial form system
  --           with that mixed cell;
  --   vrblvl  is the verbose level.

  function Cells_Standard_StaSol_into_Container
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies a start solution in double precision
  --   into the container.

  -- ON ENTRY :
  --   a       index to a mixed cell;
  --   b       index to the i-th solution of the initial form system
  --           with that mixed cell;
  --   vrblvl  is the verbose level.

  function Cells_DoblDobl_StaSol_into_Container
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies a start solution in double double precision
  --   into the container.

  -- ON ENTRY :
  --   a       index to a mixed cell;
  --   b       index to the i-th solution of the initial form system
  --           with that mixed cell;
  --   vrblvl  is the verbose level.

  function Cells_QuadDobl_StaSol_into_Container
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies a start solution in quad double precision
  --   into the container.

  -- ON ENTRY :
  --   a       index to a mixed cell;
  --   b       index to the i-th solution of the initial form system
  --           with that mixed cell;
  --   vrblvl  is the verbose level.

  function Cells_Standard_Permute
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Permutes the target system stored in double precision.
  --   The value of the verbose level is given in vrblvl.

  function Cells_DoblDobl_Permute
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Permutes the target system stored in double double precision.
  --   The value of the verbose level is given in vrblvl.

  function Cells_QuadDobl_Permute
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Permutes the target system stored in quad double precision.
  --   The value of the verbose level is given in vrblvl.

  function Cells_Floating_Mixed_Volume
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the mixed volume of the cells
  --   induced by a floating-point lifting.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the sum of all mixed volumes of all cells.

  function Cells_Set_Floating_Number_of_Supports
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the number of different supports,
  --   lifted with floating-point values.

  -- ON ENTRY :
  --   a       numbers of different supports;
  --   vrblvl  the verbose level.

  function Cells_Read_Integer_Mixed_Cells
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a regular mixed cell configuration,
  --   induced by an integer valued lifting function,
  --   and stores the mixed cells in the container.

  function Cells_Write_Integer_Mixed_Cells
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the mixed cells stored in the container
  --   for configurations induced by floating-point lifting.

  function Cells_Number_of_Integer_Mixed_Cells
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of mixed cells stored in the container,
  --   induced by an integer-valued lifting.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the number of mixed cells in a configuration
  --           induced by an integer lifting.

  function Cells_Dimension_of_Integer_Mixed_Cells
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the dimension of the mixed cells stored in the container,
  --   induced by an integer-valued lifting.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the dimension of mixed cells in a configuration
  --           induced by an integer lifting.

  function Cells_Get_Integer_Mixture
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of different supports and 
  --   the number of occurrences of each different support,
  --   for the points lifted by integer values.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the number of different supports;
  --   b       as many integers as the number of different supports,
  --           with the number of occurrences of each support.

  function Cells_Integer_Supports_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of points in each lifted support set,
  --   lifted with integer values.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the number of lifted supports;
  --   b       as many integers as the value of a[0],
  --           with the size of each support set.

  function Cells_Get_Integer_Support_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns point in the support, lifted with an integer value.

  -- ON ENTRY :
  --   a       index to the support;
  --   b       index to a point in the support;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       coordinates of the lifted point.

  function Cells_Integer_Normal
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the inner normal to a mixed cell
  --   with integer lifting values.

  -- ON ENTRY :
  --   a       the index to a mixed cell;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       coordinates of the inner normal.

  function Cells_Integer_Cell_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of points in a mixed cell
  --   with integer lifting values.

  -- ON ENTRY :
  --   a       the index to a mixed cell;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       contains the number of points in each support set
  --           of the mixed cell.

  function Cells_Get_Integer_Cell_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns a lifted point in a mixed cell
  --   with integer lifting values.

  -- ON ENTRY :
  --   a       the index to a mixed cell;
  --   b       in b[0] is the support index;
  --           in b[1] is the index of the point in the support;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       coordinates of the lifted cell point.

  function Cells_Integer_Mixed_Volume
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the mixed volume of a cell with integer lifting.

  -- ON ENTRY :
  --   a       the index to a mixed cells;
  --   vrblv   is the verbose level.

  -- ON RETURN :
  --   b       the mixed volume of the cell.

  function Cells_Set_Integer_Number_of_Supports
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the number of different supports,
  --   lifted with integer values.

  -- ON ENTRY :
  --   a       numbers of different supports;
  --   vrblvl  the verbose level.

  function Cells_Set_Integer_Mixture
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the number of different supports and 
  --   the occurrences of each different supports,
  --   for mixed cells induced by integer-valued lifting.

  -- ON ENTRY :
  --   a       the number of different supports;
  --   b       as many integers as the value of a[0],
  --           with the occurrence of each different support;
  --   vrblvl  is the verbose level.

  function Cells_Add_Integer_Support_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Adds a point lifted with an integer value to a support.

  -- ON ENTRY :
  --   a       in a[0] is the index of the support,
  --   b       in b[0] is the size of the point;
  --   c       the coordinates of the point;
  --   vrblvl  is the verbose level.

  function Cells_Make_Integer_Subdivision
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes the subdivision for the points with integer liftings.
  --   The value of the verbose level is in vrblvl.

  function Cells_Add_Integer_Mixed_Cell
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Appends a mixed cell with integer lifting values.

  -- ON ENTRY :
  --   a       in a[0] is the number R of different supports,
  --           in a[1] is the dimension of the lifted points,
  --           in a[2] is the size of b;
  --   b       in b[0] is the total number of points in the cell,
  --           in b[k] is the number of points in the k-th support,
  --           in b[1+R+k] is the label for the k-th point in the cell;
  --   c       the coordinates for the inner normal to the cell;
  --   vrblvl  is the verbose level.

  function Cells_Get_Integer_Mixed_Cell
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Retrieves a mixed cell with integer lifting values.

  -- ON ENTRY :
  --   a       in a[0] is the index to the mixed cell;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the total number of points in the cell,
  --           in b[k] is the number of points in the k-th support,
  --           in b[1+R+k] is the label for the k-th point in the cell;
  --   c       the coordinates for the inner normal to the cell.

  function Cells_Is_Stable
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns true or false, depending whether the mixed cells
  --   were computed with a stable lifting bound.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       either 0 or 1, false or true, depending on whether
  --           the cells were computed with a stable lifting bound.

  function Cells_Number_of_Original_Cells
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of original cells,
  --   cells without artificial origin.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the number of original cells.

  function Cells_Number_of_Stable_Cells
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of stable cells,
  --   cells which contribute to the count of
  --   solutions with zero coordinates.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the number of stable cells.

  function Cells_Standard_Stable_Solve
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Solve a stable start system in double precision.

  -- ON ENTRY :
  --   a       the index of a mixed cell;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the number of solutions computed.

  function Cells_DoblDobl_Stable_Solve
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Solve a stable start system in double double precision.

  -- ON ENTRY :
  --   a       the index of a mixed cell;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the number of solutions computed.

  function Cells_QuadDobl_Stable_Solve
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Solve a stable start system in quad double precision.

  -- ON ENTRY :
  --   a       the index of a mixed cell;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the number of solutions computed.

  function Cells_Floating_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears all mixed cells induced by floating-point lifting.
  --   In vrblvl is the value of the verbose level.

  function Cells_Integer_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears all mixed cells induced by integer-valued lifting.
  --   In vrblvl is the value of the verbose level.

end Cells_Interface;
