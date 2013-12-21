with text_io;                            use text_io;
with System_Call;

procedure ts_syscall is

-- DESCRIPTION :
--   Performs an "ls" to check the system call.

begin
  new_line;
  put_line("Test of system call: execution of ls.");
  new_line;
  System_Call.Call("ls");
end ts_syscall;
