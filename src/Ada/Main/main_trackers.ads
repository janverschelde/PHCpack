with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Main_Trackers is

-- DESCRIPTION :
--   The procedures below give access to jumpstarting homotopies.

  function Ask_for_Start_Type return character;

  -- DESCRIPTION :
  --   Displays the menu of the types of start systems, and asks the
  --   user to make a choice from a menu.

  -- ON RETURN :
  --   A character with a number of a type of start system:
  --     '1' : the start system is based on the total degree;
  --     '2' : a linear-product start system will be used;
  --     '3' : the user will provide start system and solutions;
  --     '4' : polyhedral continuation to solve random system;
  --     '5' : diagonal homotopy to intersect algebraic sets;
  --     '6' : one level down continuation in a homotopy cascade;
  --     '7' : remove last slack variable in a witness set.

  procedure Standard_Track
              ( target,start,output : in file_type;
                kind : in character );

  -- DESCRIPTION :
  --   Runs the cheater homotopy (option 3 of the menu above)
  --   with standard double floating point arithmetic.

  -- REQUIRED : kind is either '1' or '3'.

  procedure DoblDobl_Track
              ( target,start,output : in file_type; kind : in character );

  -- DESCRIPTION :
  --   Tracks paths with a homotopy (option 1 or 3 of the menu above)
  --   with double double floating point arithmetic.

  procedure QuadDobl_Track
              ( target,start,output : in file_type; kind : in character  );

  -- DESCRIPTION :
  --   Tracks paths with a homotopy (option 1 or 3 of the menu above)
  --   with quad double floating point arithmetic.

  procedure Track ( target,start,output : in file_type; kd : in character );

  -- DESCRIPTION :
  --   Prompts the user for the working precision 
  --   and tracks the paths for a homotopy defined by kd.

  procedure Main ( targetfilename,startfilename,outfilename : in string;
                   verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   This is the main interface to the path trackers.
  --   The arguments are the respective names of the files
  --   for the target system, start system, output file;
  --   and at last the verbose level.

end Main_Trackers;
