with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with DecaDobl_Complex_Laurentials;      use DecaDobl_Complex_Laurentials;
with DecaDobl_Complex_Laurentials_io;   use DecaDobl_Complex_Laurentials_io;
with DecaDobl_Polynomial_Convertors;    use DecaDobl_Polynomial_Convertors;
with Multprec_Complex_Laurentials_io;
with Multprec_Complex_Laur_Systems;
with Multprec_Complex_Laur_Systems_io;

package body DecaDobl_Complex_Laur_Systems_io is

  procedure get ( p : out Link_to_Laur_Sys ) is

    lp : Multprec_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    Multprec_Complex_Laurentials_io.Set_Working_Precision(5);
    Multprec_Complex_Laur_Systems_io.get(lp);
    declare
      dp : constant DecaDobl_Complex_Laur_Systems.Laur_Sys
         := Multprec_Laur_Sys_to_DecaDobl_Complex(lp.all);
    begin
      p := new DecaDobl_Complex_Laur_Systems.Laur_Sys'(dp);
    end;
    Multprec_Complex_Laur_Systems.Clear(lp);
  end get;

  procedure get ( file : in file_type; p : out Link_to_Laur_Sys ) is

    lp : Multprec_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    Multprec_Complex_Laurentials_io.Set_Working_Precision(5);
    Multprec_Complex_Laur_Systems_io.get(file,lp);
    declare
      dp : constant DecaDobl_Complex_Laur_Systems.Laur_Sys
         := Multprec_Laur_Sys_to_DecaDobl_Complex(lp.all);
    begin
      p := new DecaDobl_Complex_Laur_Systems.Laur_Sys'(dp);
    end;
    Multprec_Complex_Laur_Systems.Clear(lp);
  end get;

  procedure put ( p : in Laur_Sys ) is
  begin
    put(standard_output,p);
  end put;

  procedure put ( n : natural32; p : in Laur_Sys ) is
  begin
    put(standard_output,n,p);
  end put;

  procedure put ( n,m : natural32; p : in Laur_Sys ) is
  begin
    put(standard_output,n,m,p);
  end put;

  procedure put ( file : in file_type; p : in Laur_Sys ) is
  begin
    for i in p'range loop
      put(file,p(i));
      new_line(file);
    end loop;
  end put;

  procedure put ( file : in file_type; n : in natural32; p : in Laur_Sys ) is
  begin
    put(file,n,2); new_line(file);
    for i in p'range loop
      put(file,p(i));
      new_line(file);
    end loop;
  end put;

  procedure put ( file : in file_type;
                  n,m : in natural32; p : in Laur_Sys ) is
  begin
    put(file,n,2); put(file," ");
    put(file,m,2); new_line(file);
    for i in p'range loop
      put(file,p(i));
      new_line(file);
    end loop;
  end put;

  procedure put_line ( p : in Laur_Sys ) is
  begin
    put_line(standard_output,p);
  end put_line;

  procedure put_line ( file : in file_type; p : in Laur_Sys ) is

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

end DecaDobl_Complex_Laur_Systems_io;
