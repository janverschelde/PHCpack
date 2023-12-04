with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with HexaDobl_Complex_Matrices_io;

package body HexaDobl_Complex_Matrix_Series_io is

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
      HexaDobl_Complex_Matrices_io.put(A.cff(k).all);
    end loop;
  end put;

end HexaDobl_Complex_Matrix_Series_io;
