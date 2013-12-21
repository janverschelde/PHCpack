with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;

package body Generic_Arrays_of_Vector_Lists_io is

  use Lists_io;

  procedure get ( al : in out Link_to_Array_of_Lists ) is
  begin
    get(Standard_Input,al);
  end get;

  procedure get ( n : in natural32; m : in Standard_Natural_Vectors.Vector;
                  al : out Array_of_Lists ) is
  begin
    get(Standard_Input,n,m,al);
  end get;

  procedure get ( file : in file_type; al : in out Link_to_Array_of_Lists ) is

    n : integer32 := 0;

  begin
    get(file,n);
    al := new Array_of_Lists(1..n);
    declare
      m : Standard_Natural_Vectors.Vector(1..n) := (1..n => 0);
    begin
      get(file,m);
      get(file,natural32(n),m,al.all);
    end;
  end get;

  procedure get ( file : in file_type;
                  n : in natural32; m : in Standard_Natural_Vectors.Vector;
		  al : out Array_of_Lists ) is
  begin
    for i in al'range loop
      get(file,n,m(i),al(i));
    end loop;
  end get;

  procedure put ( al : in Array_of_Lists ) is
  begin
    put(Standard_Output,al);
  end put;

  procedure put ( file : in file_type; al : in Array_of_Lists ) is
  begin
    for i in al'range loop
      put(file,al(i)); new_line(file);
    end loop;
  end put;

end Generic_Arrays_of_Vector_Lists_io;
