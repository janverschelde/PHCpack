with text_io,integer_io;                use text_io,integer_io;
with System;
with GNAT.Threads;                      use GNAT.Threads;

package body Multithreading is

-- NOTE :
--   The delay in executing the thread code was necessary to avoid
--   that one thread grabs all CPU cycles and blocks the other threads.
--   The delay of 1 second is enough to have all threads be created.
--   There ought to be a better solution, but it works for now...

  procedure Start_Threads ( n : in natural ) is

    a : array(1..n) of System.Address;
    p : array(1..n) of Void_Ptr;

    function code ( k : Void_Ptr; b : Void_Ptr ) return Void_Ptr is

      id : constant integer := integer(b.all);

    begin
      delay 1.0; -- wait for all threads to be created
      Thread_Code(id);
      return null;
    exception 
      when others =>
        put("exception raised in thread "); put(id,1); new_line;
        raise;
    end code;

  begin
    for i in 1..n loop
      p(i) := new integer'(i);
      a(i) := Create_Thread(code'address,p(i),10000,0);
    end loop;
  end Start_Threads;

end Multithreading;
