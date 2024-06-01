with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;

package body Double_Taylor_Homotopies_io is

  procedure put ( tm : in Taylor_Monomial ) is
  begin
    put(standard_output,tm);
  end put;

  procedure put ( file : in file_type; tm : in Taylor_Monomial ) is
  begin
    put(file,tm.dim,1); put(file,"  ");
    put(file,tm.deg,1); new_line(file);
    put(file,tm.pwr); new_line(file);
    put(file,tm.cst); new_line(file);
    put_line(file,tm.cff);
    put(file,tm.exp); new_line(file);
  end put;

  procedure put ( tm : in Link_to_Taylor_Monomial ) is
  begin
    put(standard_output,tm);
  end put;

  procedure put ( file : in file_type; tm : in Link_to_Taylor_Monomial ) is
  begin
    if tm /= null
     then put(file,tm.all);
    end if;
  end put;

  procedure put ( tmv : in Taylor_Monomial_Vector ) is
  begin
    put(standard_output,tmv);
  end put;

  procedure put ( file : in file_type; tmv : in Taylor_Monomial_Vector ) is
  begin
    for i in tmv'range loop
      put(file,tmv(i));
    end loop;
  end put;

  procedure put ( tmv : in Link_to_Taylor_Monomial_Vector ) is
  begin
    put(standard_output,tmv);
  end put;

  procedure put ( file : in file_type;
                  tmv : in Link_to_Taylor_Monomial_Vector ) is
  begin
    if tmv /= null
     then put(file,tmv.all);
    end if;
  end put;

  procedure put ( th : in Taylor_Homotopy ) is
  begin
    put(standard_output,th);
  end put;

  procedure put ( file : in file_type; th : in Taylor_Homotopy ) is
  begin
    for i in th'range loop
      put(file,th(i));
    end loop;
  end put;

  procedure put ( th : in Link_to_Taylor_Homotopy ) is
  begin
    put(standard_output,th);
  end put;

  procedure put ( file : in file_type; th : in Link_to_Taylor_Homotopy ) is
  begin
    if th /= null
     then put(file,th.all);
    end if;
  end put;

end Double_Taylor_Homotopies_io;
