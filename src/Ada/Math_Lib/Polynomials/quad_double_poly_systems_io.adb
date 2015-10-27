with QuadDobl_Complex_to_Real_Poly;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;

package body Quad_Double_Poly_Systems_io is

  procedure get ( p : out Link_to_Poly_Sys ) is

    cp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    get(cp);
   -- declare
   --   s : Quad_Double_Poly_Systems.Poly_Sys(cp'range);
   -- begin
   --   s := QuadDobl_Complex_to_Real_Poly.Convert_Complex_to_Real(cp.all);
   --   p := new Quad_Double_Poly_Systems.Poly_Sys'(s);
   -- end;
    p := QuadDobl_Complex_to_Real_Poly.Convert_Complex_to_Real(cp);
    QuadDobl_Complex_Poly_Systems.Clear(cp);
  end get;

  procedure get ( file : in file_type; p : out Link_to_Poly_Sys ) is

    cp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    get(file,cp);
    p := QuadDobl_Complex_to_Real_Poly.Convert_Complex_to_Real(cp);
    QuadDobl_Complex_Poly_Systems.Clear(cp);
  end get;

  procedure put ( p : in Poly_Sys ) is
  begin
    put(standard_output,p);
  end put;

  procedure put ( file : in file_type; p : in Poly_Sys ) is

    cp : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range)
       := QuadDobl_Complex_to_Real_Poly.Convert_Real_to_Complex(p);

  begin
    put(file,cp);
    QuadDobl_Complex_Poly_Systems.Clear(cp);
  end put;

end Quad_Double_Poly_Systems_io;
