with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with DecaDobl_Complex_Vectors_io;

package body DecaDobl_Complex_Vector_Series_io is

  procedure put ( v : in Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( file : in file_type; v : in Vector ) is
  begin
    for i in 0..v.deg loop
      if i > 0
       then new_line(file);
      end if;
      DecaDobl_Complex_Vectors_io.put_line(file,v.cff(i));
    end loop;
  end put;

end DecaDobl_Complex_Vector_Series_io;
