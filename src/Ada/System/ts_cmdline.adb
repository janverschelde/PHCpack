with text_io,integer_io;                 use text_io,integer_io;
with Unix_Command_Line;                  use Unix_Command_Line;

procedure ts_cmdline is

-- DESCRIPTION :
--   Prints the arguments of the command line.

  narg : constant natural := Number_of_Arguments;

begin
  put("The number of arguments : "); put(narg,1); new_line;
  put_line("The arguments :");
  for k in 0..narg loop
    put("arg("); put(k,1); put("): ");
    put_line(Argument(k));
  end loop;
end ts_cmdline;
