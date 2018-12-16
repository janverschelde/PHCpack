with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with QuadDobl_Complex_Monomials;        use QuadDobl_Complex_Monomials;
with QuadDobl_Complex_Monomials_io;     use QuadDobl_Complex_Monomials_io;

package body QuadDobl_Monomial_Vectors_io is

  procedure put ( v : in Monomial_Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( file : in file_type; v : in Monomial_Vector ) is
  begin
    for i in v'range loop
      if v(i) /= null then
        put(file,"-> monomial "); put(file,i,1); put_line(file," :");
        put(file,v(i).all);
      end if;
    end loop;
  end put;

  procedure put ( v : in Link_to_Monomial_Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( file : in file_type; v : in Link_to_Monomial_Vector ) is
  begin
    if v /= null then
      for i in v'range loop
        if v(i) /= null then
          put(file,"-> monomial "); put(file,i,1); put_line(file," :");
          put(file,v(i).all);
        end if;
      end loop;
    end if;
  end put;

end QuadDobl_Monomial_Vectors_io;
