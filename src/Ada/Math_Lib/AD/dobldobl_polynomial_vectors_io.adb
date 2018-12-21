with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with DoblDobl_Monomial_Vectors;         use DoblDobl_Monomial_Vectors;
with DoblDobl_Monomial_Vectors_io;      use DoblDobl_Monomial_Vectors_io;

package body DoblDobl_Polynomial_Vectors_io is

  procedure put ( v : in Polynomial_Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( file : in file_type; v : in Polynomial_Vector ) is
  begin
    for i in v'range loop
      if v(i) /= null then
        put(file,"-> polynomial "); put(file,i,1); put_line(file," :");
        put(file,v(i).all);
      end if;
    end loop;
  end put;

  procedure put ( v : in Link_to_Polynomial_Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( file : in file_type; v : in Link_to_Polynomial_Vector ) is
  begin
    if v /= null then
      for i in v'range loop
        if v(i) /= null then
          put(file,"-> polynomial "); put(file,i,1); put_line(file," :");
          put(file,v(i).all);
        end if;
      end loop;
    end if;
  end put;

  procedure put ( s : in System ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( file : in file_type; s : in System ) is
  begin
    put(file,"largest exponents : "); put(file,s.maxexp);
    if s.deg1
     then put_line(file,"  no large exponents");
     else put_line(file,"  exponents larger than one");
    end if;
    put(file,s.pols);
  end put;

  procedure put ( s : in Link_to_System ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( file : in file_type; s : in Link_to_System ) is
  begin
    if s /= null
     then put(file,s.all);
    end if;
  end put;

end DoblDobl_Polynomial_Vectors_io;
