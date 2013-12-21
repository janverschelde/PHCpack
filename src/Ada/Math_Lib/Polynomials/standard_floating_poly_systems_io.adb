with Standard_Complex_to_Real_Poly;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;

package body Standard_Floating_Poly_Systems_io is

  procedure get ( p : out Link_to_Poly_Sys ) is

    cp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    get(cp);
   -- declare
   --   s : Standard_Floating_Poly_Systems.Poly_Sys(cp'range);
   -- begin
   --   s := Standard_Complex_to_Real_Poly.Convert_Complex_to_Real(cp.all);
   --   p := new Standard_Floating_Poly_Systems.Poly_Sys'(s);
   -- end;
    p := Standard_Complex_to_Real_Poly.Convert_Complex_to_Real(cp);
    Standard_Complex_Poly_Systems.Clear(cp);
  end get;

  procedure get ( file : in file_type; p : out Link_to_Poly_Sys ) is

    cp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    get(file,cp);
    p := Standard_Complex_to_Real_Poly.Convert_Complex_to_Real(cp);
    Standard_Complex_Poly_Systems.Clear(cp);
  end get;

  procedure put ( p : in Poly_Sys ) is
  begin
    put(standard_output,p);
  end put;

  procedure put ( file : in file_type; p : in Poly_Sys ) is

    cp : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
       := Standard_Complex_to_Real_Poly.Convert_Real_to_Complex(p);

  begin
    put(file,cp);
    Standard_Complex_Poly_Systems.Clear(cp);
  end put;

end Standard_Floating_Poly_Systems_io;
