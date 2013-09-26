package Unix_Command_Line is

-- DESCRIPTION :
--   This package allows to access the arguments from the command line.

  function Number_of_Arguments return natural;

  -- DESCRIPTION :
  --   Returns the number of arguments of the command that started
  --   the current process.
  --   So, if the command has no arguments, then zero is returned.

  function Argument ( i : natural ) return string;

  -- DESCRIPTION :
  --   Returns the ith argument of the command line.

end Unix_Command_Line;
