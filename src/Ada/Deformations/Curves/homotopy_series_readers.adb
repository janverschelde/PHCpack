with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Homotopy;
with Standard_System_and_Solutions_io;
with DoblDobl_Homotopy;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Homotopy;
with QuadDobl_System_and_Solutions_io;
with Standard_Series_Poly_Systems;
with DoblDobl_Series_Poly_Systems;
with QuadDobl_Series_Poly_Systems;
with Series_and_Homotopies;
with Series_and_Predictors;

package body Homotopy_Series_Readers is

  procedure Standard_Reader
              ( nbequ : out integer32;
                sols : out Standard_Complex_Solutions.Solution_List;
                tpow : in natural32;
                gamma : in Standard_Complex_Numbers.Complex_Number ) is

    target,start : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Testing the creation of a homotopies as a series system ...");
    new_line;
    put_line("Reading the target system ..."); get(target);
    nbequ := target'last;
    new_line;
    put_line("Reading the start system and its solutions ...");
    Standard_System_and_Solutions_io.get(start,sols);
    Standard_Homotopy.Create(target.all,start.all,tpow,gamma);
  end Standard_Reader;

  procedure DoblDobl_Reader
              ( nbequ : out integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                tpow : in natural32;
                gamma : in DoblDobl_Complex_Numbers.Complex_Number ) is

    target,start : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Testing the creation of a homotopies as a series system ...");
    new_line;
    put_line("Reading the target system ..."); get(target);
    nbequ := target'last;
    new_line;
    put_line("Reading the start system and its solutions ...");
    DoblDobl_System_and_Solutions_io.get(start,sols);
    DoblDobl_Homotopy.Create(target.all,start.all,tpow,gamma);
  end DoblDobl_Reader;

  procedure QuadDobl_Reader
              ( nbequ : out integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                tpow : in natural32;
                gamma : in QuadDobl_Complex_Numbers.Complex_Number ) is

    target,start : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Testing the creation of a homotopies as a series system ...");
    new_line;
    put_line("Reading the target system ..."); get(target);
    nbequ := target'last;
    new_line;
    put_line("Reading the start system ...");
    QuadDobl_System_and_Solutions_io.get(start,sols);
    QuadDobl_Homotopy.Create(target.all,start.all,tpow,gamma);
  end QuadDobl_Reader;

  procedure Standard_Reader
              ( nbequ : out integer32;
                sols : out Standard_Complex_Solutions.Solution_List;
                tpow : in natural32 := 2 ) is

    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;

  begin
    Standard_Reader(nbequ,sols,tpow,gamma);
  end Standard_Reader;

  procedure DoblDobl_Reader
              ( nbequ : out integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                tpow : in natural32 := 2 ) is

    gamma : constant DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Random_Numbers.Random1;

  begin
    DoblDobl_Reader(nbequ,sols,tpow,gamma);
  end DoblDobl_Reader;

  procedure QuadDobl_Reader
              ( nbequ : out integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                tpow : in natural32 := 2 ) is

    gamma : constant QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Random_Numbers.Random1;

  begin
    QuadDobl_Reader(nbequ,sols,tpow,gamma);
  end QuadDobl_Reader;

  procedure Standard_Series_Newton
              ( sol : in Standard_Complex_Vectors.Vector;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out Standard_Dense_Series_Vectors.Vector ) is

    hom : Standard_Complex_Poly_Systems.Poly_Sys(1..nbequ)
        := Standard_Homotopy.Homotopy_System;
    sys : Standard_Series_Poly_Systems.Poly_Sys(1..nbequ)
        := Series_and_Homotopies.Create(hom,idx);
    nit : constant integer32 := integer32(nbiters);

  begin
    Series_and_Predictors.Newton_Prediction(nit,sys,sol,srv,eva);
   -- Standard_Complex_Poly_Systems.Clear(hom); -- sharing, do not clear!
    Standard_Series_Poly_Systems.Clear(sys);
  end Standard_Series_Newton;

  procedure Standard_Series_Newton
              ( sol : in Standard_Complex_Solutions.Solution;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out Standard_Dense_Series_Vectors.Vector ) is
  begin
    Standard_Series_Newton(sol.v,idx,nbequ,nbterms,nbiters,srv,eva);
  end Standard_Series_Newton;

  procedure DoblDobl_Series_Newton
              ( sol : in DoblDobl_Complex_Vectors.Vector;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out DoblDobl_Dense_Series_Vectors.Vector ) is

    hom : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nbequ)
        := DoblDobl_Homotopy.Homotopy_System;
    sys : DoblDobl_Series_Poly_Systems.Poly_Sys(1..nbequ)
        := Series_and_Homotopies.Create(hom,idx);
    nit : constant integer32 := integer32(nbiters);

  begin
    Series_and_Predictors.Newton_Prediction(nit,sys,sol,srv,eva);
   -- DoblDobl_Complex_Poly_Systems.Clear(hom); -- sharing, do not clear!
    DoblDobl_Series_Poly_Systems.Clear(sys);
  end DoblDobl_Series_Newton;

  procedure DoblDobl_Series_Newton
              ( sol : in DoblDobl_Complex_Solutions.Solution;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out DoblDobl_Dense_Series_Vectors.Vector ) is
  begin
    DoblDobl_Series_Newton(sol.v,idx,nbequ,nbterms,nbiters,srv,eva);
  end DoblDobl_Series_Newton;

  procedure QuadDobl_Series_Newton
              ( sol : in QuadDobl_Complex_Vectors.Vector;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out QuadDobl_Dense_Series_Vectors.Vector ) is

    hom : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nbequ)
        := QuadDobl_Homotopy.Homotopy_System;
    sys : QuadDobl_Series_Poly_Systems.Poly_Sys(1..nbequ)
        := Series_and_Homotopies.Create(hom,idx);
    nit : constant integer32 := integer32(nbiters);

  begin
    Series_and_Predictors.Newton_Prediction(nit,sys,sol,srv,eva);
   -- QuadDobl_Complex_Poly_Systems.Clear(hom); -- sharing, do not clear!
    QuadDobl_Series_Poly_Systems.Clear(sys);
  end QuadDobl_Series_Newton;

  procedure QuadDobl_Series_Newton
              ( sol : in QuadDobl_Complex_Solutions.Solution;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out QuadDobl_Dense_Series_Vectors.Vector ) is
  begin
    QuadDobl_Series_Newton(sol.v,idx,nbequ,nbterms,nbiters,srv,eva);
  end QuadDobl_Series_Newton;

end Homotopy_Series_Readers;
