with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Drivers_for_Set_Structures;         use Drivers_for_Set_Structures;
with Random_Product_Start_Systems;
with Standard_Linear_Product_System;
with Set_Structure;
with Set_Structure_io;
with Degree_Sets_Tables;
with Degree_Sets_Tables_io;

procedure ts_drivss is

-- DESCRIPTION :
--   Reads a polynomial system and calls the driver.

  procedure Test_Driver is

  -- DESCRIPTION :
  --   Calls the main driver to construct a set structure
  --   for a polynomial system.

    file : file_type;
    lp : Link_to_Poly_Sys;
    lpos : List;
    b : natural32 := 0;

  begin
    get(lp);
    declare
      q : Poly_Sys(lp'range);
      qsols : Solution_List;
    begin
      put_line("Reading the output file.");
      Read_Name_and_Create_File(file);
      Driver_for_Set_Structure(file,lp.all,b,lpos,q,qsols);
    end;
  end Test_Driver;

  procedure Test_Permanent_Computation ( dim : in natural32 ) is

  -- DESCRIPTION :
  --   Given a set structure stored internally,
  --   for a square polynomial system of the given dimension dim,
  --   this procedure tests the computation of the permanent in three ways:
  --   1) based on taking unions of the sets;
  --   2) running the maximum bipartite matching algorithms;
  --   3) solving a random linear product start system.

    dst : constant Degree_Sets_Tables.Degree_Sets_Table
        := Degree_Sets_Tables.Create;
    cnt : integer32;
    ans : character;
    nbsols : natural64;

  begin
    put_line("The degree sets table :");
    Degree_Sets_Tables_io.put(dst);
    new_line;
    put("Compute the permanent with sets ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      cnt := Degree_Sets_Tables.Permanent(dst);
      put("The formal root count : "); put(cnt,1); new_line;
    end if;
    put("Compute the permanent with maximum bipartite matching ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      cnt := Degree_Sets_Tables.Matching_Permanent(dst);
      put("The formal root count : "); put(cnt,1); new_line;
    end if;
    put("Compute the permanent solving a linear-product system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Standard_Linear_Product_System.Init(dim);
      Random_Product_Start_Systems.Build_Random_Product_System(dim);
      nbsols := Standard_Linear_Product_System.Count_All_Solutions(1.0E-8);
      put("The number of solutions : "); put(nbsols,1); new_line;
    end if;
  end Test_Permanent_Computation;

  procedure Test_Root_Count is

  -- DESCRIPTION :
  --   Prompts for a system and tests the computation of the formal
  --   root count based on a supporting set structure for the system.

    lp : Link_to_Poly_Sys;

  begin
    get(lp);
    Random_Product_Start_Systems.Build_Set_Structure(lp.all);
    put_line("A supporting set structure :");
    Set_Structure_io.put;
    Test_Permanent_Computation(natural32(lp'last));
  end Test_Root_Count;

  procedure Main is

  -- DESCRIPTION :
  --   Displays the menu of the testing operations
  --   and then calls the selected test.

    ans : character;

  begin
    new_line;
    put_line("MENU to test set structures :");
    put_line("  1. test the main driver;");
    put_line("  2. test root count computation.");
    put("Type 1 or 2 to select : ");
    Ask_Alternative(ans,"12");
    new_line;
    case ans is
      when '1' => Test_Driver;
      when '2' => Test_Root_Count;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_drivss;
