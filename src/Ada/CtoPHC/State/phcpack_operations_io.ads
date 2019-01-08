with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package PHCpack_Operations_io is

-- DESCRIPTION :
--   Offers input/output routines for data stored in PHCpack_Operations.

-- INPUT OPERATIONS :

  procedure Read_Start_System;
  procedure Read_Start_System ( filename : in string );
  procedure Read_Start_Laurent_System;
  procedure Read_Start_Laurent_System ( filename : in string );

  procedure Read_DoblDobl_Start_System;
  procedure Read_DoblDobl_Start_System ( filename : in string );
  procedure Read_DoblDobl_Start_Laurent_System;
  procedure Read_DoblDobl_Start_Laurent_System ( filename : in string );

  procedure Read_QuadDobl_Start_System;
  procedure Read_QuadDobl_Start_System ( filename : in string );
  procedure Read_QuadDobl_Start_Laurent_System;
  procedure Read_QuadDobl_Start_Laurent_System ( filename : in string );

  procedure Read_Multprec_Start_System ( decimals : in natural32 );
  procedure Read_Multprec_Start_System
              ( filename : in string; decimals : in natural32 );

  -- DESCRIPTION :
  --   Reads the target system using as many decimal places as the value
  --   of decimals to parse the multiprecision numbers.

  procedure Read_Start_System_without_Solutions;
  procedure Read_Start_System_without_Solutions ( filename : in string );

  -- DESCRIPTION :
  --   Prompts the user for a file name and reads the start system,
  --   without the list of solutions.  However, as side effect, the
  --   Solutions_Container.Solution_Input_File has a value, which permits
  --   to read the solutions from the same file as the file name provided
  --   by the user, whenever a solution is needed.
  --   If the filename is provided, then the file with the name will be read.

  procedure Read_Linear_Product_Start_System ( fail : out boolean );
  procedure Read_Linear_Product_Start_System
              ( filename : in string; fail : out boolean );

  -- DESCRIPTION :
  --   Either the user will be prompted to provide a file name
  --   or when filename is provided as parameter,
  --   a linear-product start system will be read from file.
  --   After a successful reading, the data in Standard_Complex_Prod_Planes 
  --   contains a linear-product start system.

  procedure Read_Witness_Set_for_Diagonal_Homotopy
              ( k : in natural32; n,dim,deg : out natural32;
                fail : out boolean );

  -- DESCRIPTION :
  --   Reads the polynomials for the k-th witness set from file,
  --   given by the user.

  -- ON ENTRY :
  --   k        indicates first (k=1) or second (k=2) witness set.

  -- ON RETURN :
  --   n        ambient dimension of the k-th witness set;
  --   dim      dimension of the k-th witness set;
  --   deg      degree of the k-th witness set;
  --   fail     false if the reading went well, otherwise true.

  procedure Reset_Witness_Input_File
              ( k : in natural32; deg,dim : out natural32; fail : out boolean );

  -- DESCRIPTION :
  --   Resets the input file to read the k-th witness set.

  -- ON ENTRY :
  --   k        indicates first (k=1) or second (k=2) witness set.

  -- ON RETURN :
  --   deg      degree of the witness set, #solutions on file;
  --   dim      ambient dimension, length of the solution vectors;
  --   fail     true if something went wrong, false otherwise.

  procedure Read_Start_Solutions;
  procedure Read_Start_Solutions ( filename : in string );

  procedure Read_Target_System;
  procedure Read_Target_System ( filename : in string );
  procedure Read_Target_Laurent_System;
  procedure Read_Target_Laurent_System ( filename : in string );

  procedure Read_Target_System_without_Solutions;
  procedure Read_Target_System_without_Solutions ( filename : in string );

  procedure Read_DoblDobl_Target_System_without_Solutions;

  -- DESCRIPTION :
  --   Prompts the user for a file name and reads from file
  --   a system with double double precision coefficients.
  --   The system is stored as a dobldobl target system.

  procedure Read_QuadDobl_Target_System_without_Solutions;

  -- DESCRIPTION :
  --   Prompts the user for a file name and reads from file
  --   a system with quad double precision coefficients.
  --   The system is stored as a quaddobl target system.

  procedure Read_Target_Solutions;
  procedure Read_Target_Solutions ( filename : in string );

  procedure Read_DoblDobl_Target_System;
  procedure Read_DoblDobl_Target_System ( filename : in string );
  procedure Read_QuadDobl_Target_System;
  procedure Read_QuadDobl_Target_System ( filename : in string );
  procedure Read_DoblDobl_Target_Laurent_System;
  procedure Read_DoblDobl_Target_Laurent_System ( filename : in string );
  procedure Read_QuadDobl_Target_Laurent_System;
  procedure Read_QuadDobl_Target_Laurent_System ( filename : in string );

  procedure Read_Multprec_Target_System ( decimals : in natural32 );
  procedure Read_Multprec_Target_System
              ( filename : in string; decimals : in natural32 );

  -- DESCRIPTION :
  --   Reads the target system using as many decimal places as the value
  --   of decimals to parse the multiprecision numbers.

-- OUTPUT OPERATIONS :

  procedure Write_Start_System;
  procedure Write_Start_System ( filename : in string );
  procedure Write_Start_Laurent_System;
  procedure Write_Start_Laurent_System ( filename : in string );
  procedure Write_Start_Solutions;
  procedure Write_Start_Solutions ( filename : in string );

  procedure Write_DoblDobl_Start_System;
  procedure Write_DoblDobl_Start_System ( filename : in string );
  procedure Write_DoblDobl_Start_Laurent_System;
  procedure Write_DoblDobl_Start_Laurent_System ( filename : in string );
  procedure Write_DoblDobl_Start_Solutions;
  procedure Write_DoblDobl_Start_Solutions ( filename : in string );

  procedure Write_QuadDobl_Start_System;
  procedure Write_QuadDobl_Start_System ( filename : in string );
  procedure Write_QuadDobl_Start_Laurent_System;
  procedure Write_QuadDobl_Start_Laurent_System ( filename : in string );
  procedure Write_QuadDobl_Start_Solutions;
  procedure Write_QuadDobl_Start_Solutions ( filename : in string );

  procedure Write_Multprec_Start_System;
  procedure Write_Multprec_Start_System ( filename : in string );
  procedure Write_Multprec_Start_Solutions;
  procedure Write_Multprec_Start_Solutions ( filename : in string );

  -- DESCRIPTION :
  --   Writes the start system (and solutions) to the defined output file,
  --   or to the given filename defined in the string.

  procedure Write_Target_System;
  procedure Write_Target_System ( filename : in string );
  procedure Write_Target_Laurent_System;
  procedure Write_Target_Laurent_System ( filename : in string );
  procedure Write_Target_Solutions;
  procedure Write_Target_Solutions ( filename : in string );

  procedure Write_DoblDobl_Target_System;
  procedure Write_DoblDobl_Target_System ( filename : in string );
  procedure Write_DoblDobl_Target_Laurent_System;
  procedure Write_DoblDobl_Target_Laurent_System ( filename : in string );
  procedure Write_DoblDobl_Target_Solutions;
  procedure Write_DoblDobl_Target_Solutions ( filename : in string );

  procedure Write_QuadDobl_Target_System;
  procedure Write_QuadDobl_Target_System ( filename : in string );
  procedure Write_QuadDobl_Target_Laurent_System;
  procedure Write_QuadDobl_Target_Laurent_System ( filename : in string );
  procedure Write_QuadDobl_Target_Solutions;
  procedure Write_QuadDobl_Target_Solutions ( filename : in string );

  procedure Write_Multprec_Target_System;
  procedure Write_Multprec_Target_System ( filename : in string );
  procedure Write_Multprec_Target_Solutions;
  procedure Write_Multprec_Target_Solutions ( filename : in string );

  -- DESCRIPTION :
  --   Writes the target system (and solutions) to the defined output file,
  --   or to the given filename defined in the string.

end PHCpack_Operations_io;
