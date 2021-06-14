with Interfaces.C;
with C_Double_io;                          use C_Double_io;

package body C_Double_Arrays_io is

  procedure put ( n : in natural; v : in C_Double_Array ) is
  begin
    put(Standard_Output,n,v);
  end put;

  procedure put ( file : in file_type;
                  n : in natural; v : in C_Double_Array ) is
  begin
    for i in 0..Interfaces.C.size_T(n-1) loop
      put(file,v(i)); new_line(file);
    end loop;
  end put;

end C_Double_Arrays_io;
