package Machines is

-- DESCRIPTION :
--   This package offers some routines to obtain some information
--   about the machines one wants to use.

  function Process_ID return integer;

  -- DESCRIPTION :
  --   Returns the process ID of the program being executed.

  function Process_ID return string;

  -- DESCRIPTION :
  --   Returns the process ID in the form of a string of characters.
  --   This is useful to make unique tempory file names, as is done
  --   in the following procedures.

  function User_Name ( pid : string ) return string;

  -- DESCRIPTION :
  --   Returns the login-name of the user.
  --   The parameter pid indicates the process ID of the calling routine.

  function Architecture ( pid : string ) return string;
  function Architecture ( pid : string; machine : string ) return string;

  -- DESCRIPTION :
  --   The architecture of the requested machine will be returned;
  --   pid is the process ID of the process calling the routine.
  --   By default, the type of the current machine is returned.

  function Host_Name ( pid : string ) return string;

  -- DESCRIPTION :
  --   The name of the machine one is currently working on is returned;
  --   pid is the ID of the process calling the routine.

  function Date ( pid : string ) return string;

  -- DESCRIPTION :
  --   Returns the current date; pid contains the ID of the caller.

end Machines;
