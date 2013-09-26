with text_io;                        use text_io;
with Machines;                       use Machines;

procedure ts_mach is

-- DESCRIPTION :
--   Prints system information about the current process.

  pid : constant string := Process_ID;
  username : constant string := User_Name(pid);
  hostname : constant string := Host_Name(pid);
  archi : constant string := Architecture(pid);
  today : constant string := Date(pid);

begin
  put_line("Hello " & username & " !");
  put_line("This is " & hostname & ".");
  put_line(archi);
  put_line("Running process with ID : " & pid & ".");
  put_line("Today date is " & today & ".");
end ts_mach;
