with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Vectors_io;

package body Standard_Dense_Vector_Series_io is

-- DESCRPTION :
--   Provides very basic output for series with vector coefficients.

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
      Standard_Complex_Vectors_io.put_line(file,v.cff(i));
    end loop;
  end put;

end Standard_Dense_Vector_Series_io;
