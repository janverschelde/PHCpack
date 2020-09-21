with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with OctoDobl_CSeries_Polynomials;      use OctoDobl_CSeries_Polynomials;
with OctoDobl_CSeries_Polynomials_io;   use OctoDobl_CSeries_Polynomials_io;

package body OctoDobl_CSeries_Poly_Systems_io is

  procedure put ( p : in Poly_Sys ) is
  begin
    put(standard_output,p);
  end put;

  procedure put ( file : in file_type; p : in Poly_Sys ) is

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
    put(file,nq,1);
    if nv /= nq
     then put(file,"  "); put(file,nv,1);
    end if;
    new_line(file);
    for i in p'range loop
      put(file,p(i));
    end loop;
  end put;

end OctoDobl_CSeries_Poly_Systems_io;
