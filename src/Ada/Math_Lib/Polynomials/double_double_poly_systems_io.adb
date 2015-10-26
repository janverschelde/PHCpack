with DoblDobl_Complex_to_Real_Poly;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;

package body Double_Double_Poly_Systems_io is

  procedure get ( p : out Link_to_Poly_Sys ) is

    cp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    get(cp);
   -- declare
   --   s : Double_Double_Poly_Systems.Poly_Sys(cp'range);
   -- begin
   --   s := DoblDobl_Complex_to_Real_Poly.Convert_Complex_to_Real(cp.all);
   --   p := new Double_Double_Poly_Systems.Poly_Sys'(s);
   -- end;
    p := DoblDobl_Complex_to_Real_Poly.Convert_Complex_to_Real(cp);
    DoblDobl_Complex_Poly_Systems.Clear(cp);
  end get;

  procedure get ( file : in file_type; p : out Link_to_Poly_Sys ) is

    cp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    get(file,cp);
    p := DoblDobl_Complex_to_Real_Poly.Convert_Complex_to_Real(cp);
    DoblDobl_Complex_Poly_Systems.Clear(cp);
  end get;

  procedure put ( p : in Poly_Sys ) is
  begin
    put(standard_output,p);
  end put;

  procedure put ( file : in file_type; p : in Poly_Sys ) is

    cp : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range)
       := DoblDobl_Complex_to_Real_Poly.Convert_Real_to_Complex(p);

  begin
    put(file,cp);
    DoblDobl_Complex_Poly_Systems.Clear(cp);
  end put;

end Double_Double_Poly_Systems_io;
