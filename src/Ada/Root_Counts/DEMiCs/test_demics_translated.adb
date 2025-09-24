with Ada.Text_IO;                       use Ada.Text_IO;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Cyclic_Roots_System;
with DEMiCs_Translated;

package body Test_DEMiCs_Translated is

  procedure Test_Cyclic ( dim : in integer32; vrblvl : in integer32 := 0 ) is

    p : constant Poly_Sys(1..dim)
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

  procedure Test_Cyclic_Roots ( vrblvl : in integer32 := 0 ) is
  begin
    put_line("-> running tests on the cyclic n-roots system ...");
    for dim in 3..11 loop
      Test_Cyclic(integer32(dim),vrblvl);
    end loop;
  end Test_Cyclic_Roots;

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
    Test_Cyclic_Roots(vrblvl);
  end Main;

end Test_DEMiCs_Translated;
