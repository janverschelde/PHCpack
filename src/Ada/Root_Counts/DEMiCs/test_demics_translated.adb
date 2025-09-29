with Ada.Text_IO;                       use Ada.Text_IO;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Cyclic_Roots_System;
with Cyclic_Laurent_System;
with Floating_Mixed_Subdivisions;       use Floating_Mixed_Subdivisions;
with DEMiCs_Translated;

package body Test_DEMiCs_Translated is

  procedure Test_Cyclic ( dim : in integer32; vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Computes the mixed volume of the cyclic n-roots system
  --   where n equals dim.

    p : constant Standard_Complex_Poly_Systems.Poly_Sys(1..dim)
      := Cyclic_Roots_System.Double_Cyclic_System(dim);
    mv : integer32;

  begin
    if vrblvl > 0 then
      new_line;
      put("the cyclic "); put(dim,1); put_line("-roots polynomials : ");
      put(p);
    end if;
    mv := DEMiCs_Translated.Mixed_Volume(p,vrblvl);
    put("mixed volume of cyclic "); put(dim,1);
    put("-roots : "); put(mv,1); new_line;
  end Test_Cyclic;

  procedure Test_Reformulated
              ( dim : in integer32; vrblvl : in integer32 := 0 ) is

    p : constant Standard_Complex_Laur_Systems.Laur_Sys(1..dim-1)
      := Cyclic_Laurent_System.Cyclic_System(natural32(dim));
    mv : integer32;

  begin
    if vrblvl > 0 then
      new_line;
      put("the reformulated cyclic "); put(dim,1);
      put_line("-roots polynomials : "); put(p);
    end if;
    mv := DEMiCs_Translated.Mixed_Volume(p,vrblvl);
    put("mixed volume of reformulated cyclic "); put(dim,1);
    put("-roots : "); put(mv,1); new_line;
  end Test_Reformulated;

  procedure Test_Labels ( dim : in integer32; vrblvl : in integer32 := 0 ) is

    p : constant Standard_Complex_Poly_Systems.Poly_Sys(1..dim)
      := Cyclic_Roots_System.Double_Cyclic_System(dim);
    mv : integer32;

  begin
    if vrblvl > 0 then
      new_line;
      put("the cyclic "); put(dim,1); put_line("-roots polynomials : ");
      put(p);
    end if;
    mv := DEMiCs_Translated.Mixed_Labels(p,true,vrblvl);
    put("mixed volume of cyclic "); put(dim,1);
    put("-roots : "); put(mv,1); new_line;
  end Test_Labels;

  procedure Test_Cells ( dim : in integer32; vrblvl : in integer32 := 0 ) is

    p : constant Standard_Complex_Poly_Systems.Poly_Sys(1..dim)
      := Cyclic_Roots_System.Double_Cyclic_System(dim);
    mv : integer32;
    mcc : Mixed_Subdivision;

  begin
    if vrblvl > 0 then
      new_line;
      put("the cyclic "); put(dim,1); put_line("-roots polynomials : ");
      put(p);
    end if;
    mv := DEMiCs_Translated.Mixed_Labels(p,true,vrblvl);
    put("mixed volume of cyclic "); put(dim,1);
    put("-roots : "); put(mv,1); new_line;
    mcc := DEMiCs_Translated.Mixed_Cells(vrblvl);
    put("number of mixed cells : ");
    put(integer32(Length_Of(mcc)),1); new_line;
  end Test_Cells;

  procedure Test_Cyclic_Roots ( vrblvl : in integer32 := 0 ) is
  begin
    new_line;
    put_line("-> running tests on the cyclic n-roots system ...");
    for dim in 3..11 loop
      Test_Cyclic(integer32(dim),vrblvl);
    end loop;
  end Test_Cyclic_Roots;

  procedure Test_Reformulated_Cyclic ( vrblvl : in integer32 := 0 ) is
  begin
    new_line;
    put_line("-> running tests on the reformulated cyclic n-roots system ...");
    for dim in 3..11 loop
      Test_Reformulated(integer32(dim),vrblvl);
    end loop;
  end Test_Reformulated_Cyclic;

  procedure Main is

    vrblvl : integer32 := 99; -- default verbose level
    ans : character;
  
  begin
    put_line("Testing the DEMiCs algorithm ...");
    new_line;
    put("Intermediate output wanted ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'n'
     then vrblvl := 0;
    end if;
    new_line;
    put_line("MENU for testing the translated DEMiCs :");
    put_line("  1. run sequence of cyclic n-roots problems");
    put_line("  2. compute labels to points in the mixed cells");
    put_line("  3. convert labels into mixed cells");
    put_line("  4. test on the reformulated cyclic n-roots systems");
    put("Type 1, 2, 3, or 4 to select a test : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' => Test_Cyclic_Roots(vrblvl);
      when '2' => Test_Labels(5,vrblvl);
      when '3' => Test_Cells(5,vrblvl);
      when '4' => Test_Reformulated_Cyclic(vrblvl);
      when others => null;
    end case;
  end Main;

end Test_DEMiCs_Translated;
