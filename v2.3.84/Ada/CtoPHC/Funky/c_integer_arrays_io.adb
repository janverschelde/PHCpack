with Interfaces.C;
with C_Integer_io;                          use C_Integer_io;

package body C_Integer_Arrays_io is

  procedure put ( n : in natural; v : in C_Integer_Array ) is
  begin
    put(Standard_Output,n,v);
  end put;

  procedure put ( file : in file_type;
                  n : in natural; v : in C_Integer_Array ) is
  begin
    for i in 0..Interfaces.C.size_T(n-1) loop
      put(file," "); put(file,v(i),1);
    end loop;
  end put;

end C_Integer_Arrays_io;
