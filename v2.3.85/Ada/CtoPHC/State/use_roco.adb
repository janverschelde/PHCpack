with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Random_Product_Start_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Set_Structure,Set_Structure_io;
with Degree_Sets_Tables;
with Standard_Linear_Product_System;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_PolySys_Container;
with Standard_Solutions_Container;

function use_roco ( job : integer32;
                    a : C_intarrs.Pointer;
                    b : C_intarrs.Pointer;
                    c : C_dblarrs.Pointer ) return integer32 is

  function Job0 return integer32 is  -- construct supporting set structure
  
    p : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;

  begin
    Random_Product_Start_Systems.Build_Set_Structure(p.all);
    return 0;
  end Job0;

  function Job1 return integer32 is -- write supporting set structure
  begin
    Set_Structure_io.put;
    return 0;
  end Job1;

  function Job2 return integer32 is -- compute the Bezout bound

   -- r : constant natural := Set_Structure.B;
    r : constant integer32 
      := Degree_Sets_Tables.Permanent(Degree_Sets_Tables.Create);

  begin
    Assign(r,a);
    return 0;
  end Job2;

  function Job3 return integer32 is -- make random linear-product system

    n : constant natural32 := Standard_PolySys_Container.Dimension;
    q : Link_to_Poly_Sys;

  begin
    Standard_Linear_Product_System.Init(n);
    Random_Product_Start_Systems.Build_Random_Product_System(n);
    q := new Poly_Sys'(Standard_Linear_Product_System.Polynomial_System);
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(q.all);
    return 0;
  end Job3;

  function Job4 return integer32 is -- solve random linear-product system

    sols : Solution_List;
    nb : natural32;

  begin
    Standard_Linear_Product_System.Solve(sols,nb);
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(sols);
    return 0;
  end Job4;

  function Job5 return integer32 is -- clears the set structure
  begin
    Set_Structure.Clear;
    Standard_Linear_Product_System.Clear;
    return 0;
  end Job5;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0;   -- supporting set structure
      when 1 => return Job1;   -- write supporting set structure
      when 2 => return Job2;   -- compute the Bezout bound
      when 3 => return Job3;   -- make random linear-product system
      when 4 => return Job4;   -- solve random linear-product system
      when 5 => return Job5;   -- clears the set structure
      when others => put_line("invalid operation"); return 1;
    end case;
  exception
    when others => put("Exception raised in use_celcon handling job ");
                   put(job,1); put_line(".  Will try to ignore...");
                   return 1;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_roco;
