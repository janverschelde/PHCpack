with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with TripDobl_Complex_Polynomials;      use TripDobl_Complex_Polynomials;
with TripDobl_Complex_Polynomials_io;   use TripDobl_Complex_Polynomials_io;
with TripDobl_Polynomial_Convertors;    use TripDobl_Polynomial_Convertors;
with Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;

package body TripDobl_Complex_Poly_Systems_io is

  procedure get ( p : out Link_to_Poly_Sys ) is

    lp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Multprec_Complex_Polynomials_io.Set_Working_Precision(7);
    Multprec_Complex_Poly_Systems_io.get(lp);
    declare
      dp : constant TripDobl_Complex_Poly_Systems.Poly_Sys
         := Multprec_Poly_Sys_to_TripDobl_Complex(lp.all);
    begin
      p := new TripDobl_Complex_Poly_Systems.Poly_Sys'(dp);
    end;
    Multprec_Complex_Poly_Systems.Clear(lp);
  end get;
  
  procedure get ( file : in file_type; p : out Poly_Sys ) is

    lp : Multprec_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    Multprec_Complex_Polynomials_io.Set_Working_Precision(7);
    Multprec_Complex_Poly_Systems_io.get(file,lp);
    p := Multprec_Poly_Sys_to_TripDobl_Complex(lp);
    Multprec_Complex_Poly_Systems.Clear(lp);
  end get;
  
  procedure get ( file : in file_type; p : out Link_to_Poly_Sys ) is

    lp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Multprec_Complex_Polynomials_io.Set_Working_Precision(7);
    Multprec_Complex_Poly_Systems_io.get(file,lp);
    declare
      dp : constant TripDobl_Complex_Poly_Systems.Poly_Sys
         := Multprec_Poly_Sys_to_TripDobl_Complex(lp.all);
    begin
      p := new TripDobl_Complex_Poly_Systems.Poly_Sys'(dp);
    end;
    Multprec_Complex_Poly_Systems.Clear(lp);
  end get;

  procedure put ( p : in Poly_Sys ) is
  begin
    put(standard_output,p);
  end put;

  procedure put ( file : in file_type; p : in Poly_Sys ) is
  begin
    for i in p'range loop
      put(file,p(i));
      new_line(file);
    end loop;
  end put;

  procedure put ( n : in natural32; p : in Poly_Sys ) is
  begin
    put(standard_output,n,p);
  end put;

  procedure put ( file : in file_type; n : in natural32; p : in Poly_Sys ) is
  begin
    put(file,n,2); new_line(file);
    put(file,p);
  end put;

  procedure put ( n,m : in natural32; p : in Poly_Sys ) is
  begin
    put(standard_output,n,m,p);
  end put;

  procedure put ( file : in file_type;
                  n,m : in natural32; p : in Poly_Sys ) is
  begin
    put(file,n,2); put(file," "); put(file,m,2); new_line(file);
    put(file,p);
  end put;

  procedure put_line ( p : in Poly_Sys ) is
  begin
    put_line(standard_output,p);
  end put_line;

  procedure put_line ( file : in file_type; p : in Poly_Sys ) is

    n : constant natural32 := Number_of_Unknowns(p(p'first));
    nq : constant natural32 := natural32(p'length);

  begin
    put(file,nq,2);
    if n /= nq 
     then put(file," "); put(file,n,1);
    end if;
    new_line(file);
    for i in p'range loop
      put_line(file,p(i));
      new_line(file);
    end loop;
  end put_line;

  procedure put ( p : in Poly_Sys; s : in Array_of_Symbols ) is
  begin
    put(standard_output,p,s);
  end put;

  procedure put ( file : in file_type;
                  p : in Poly_Sys; s : in Array_of_Symbols ) is
  begin
    for i in p'range loop
      put(file,p(i),s);
      new_line(file);
    end loop;
  end put;

end TripDobl_Complex_Poly_Systems_io;
