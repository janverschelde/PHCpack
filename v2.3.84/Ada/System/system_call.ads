package System_Call is

-- DESCRIPTION :
--   This package provides a routine to call UNIX shell commands.

  procedure Call ( Command: in string );

  -- DESCRIPTION :
  --   Command is passed to and executed by the operating system.

  System_Error: exception;   -- failure to execute Command properly

end System_Call;
