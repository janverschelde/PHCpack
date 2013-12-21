with Multprec_Complex_to_Real_Poly;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;

package body Multprec_Floating_Poly_Systems_io is

  procedure get ( p : out Link_to_Poly_Sys ) is

    cp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    get(cp);
    declare
      s : Multprec_Floating_Poly_Systems.Poly_Sys(cp'range);
    begin
      s := Multprec_Complex_to_Real_Poly.Convert_Complex_to_Real(cp.all);
      p := new Multprec_Floating_Poly_Systems.Poly_Sys'(s);
    end;
    Multprec_Complex_Poly_Systems.Clear(cp);
  end get;

  procedure put ( p : in Poly_Sys ) is
  begin
    put(standard_output,p);
  end put;

  procedure put ( file : in file_type; p : in Poly_Sys ) is

    cp : Multprec_Complex_Poly_Systems.Poly_Sys(p'range)
       := Multprec_Complex_to_Real_Poly.Convert_Real_to_Complex(p);

  begin
    put(file,cp);
    Multprec_Complex_Poly_Systems.Clear(cp);
  end put;

end Multprec_Floating_Poly_Systems_io;
