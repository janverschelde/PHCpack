with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with DoblDobl_Complex_Monomials;        use DoblDobl_Complex_Monomials;
with DoblDobl_Complex_Monomials_io;     use DoblDobl_Complex_Monomials_io;

package body DoblDobl_Monomial_Vectors_io is

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
    if v /= null
     then put(file,v.all);
    end if;
  end put;

  procedure put ( p : in Polynomial ) is
  begin
    put(standard_output,p);
  end put;

  procedure put ( file : in file_type; p : in Polynomial ) is
  begin
    put(file,"degree : "); put(file,Degree(p),1); new_line(file);
    put(file,"largest exponents : "); put(file,p.maxexp);
    if p.deg1
     then put_line(file,"  no large exponents");
     else put_line(file,"  exponents larger than one");
    end if;
    put_line(file,"-> monomial 0 :");
    put(file,"coefficient : ");
    put(file,p.cff0); new_line(file);
    put(file,p.mons);
  end put;

  procedure put ( p : in Link_to_Polynomial ) is
  begin
    put(standard_output,p);
  end put;

  procedure put ( file : in file_type; p : in Link_to_Polynomial ) is
  begin
    if p /= null
     then put(file,p.all);
    end if;
  end put;

end DoblDobl_Monomial_Vectors_io;
