with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with PentDobl_Complex_Matrices_io;

package body PentDobl_Complex_Matrix_Series_io is

  procedure put ( A : in Matrix ) is
  begin
    put(standard_output,A);
  end put;

  procedure put ( file : in file_type; A : in Matrix ) is
  begin
    for k in 0..A.deg loop
      if k > 0
       then new_line(file);
      end if;
      PentDobl_Complex_Matrices_io.put(A.cff(k).all);
    end loop;
  end put;

end PentDobl_Complex_Matrix_Series_io;
