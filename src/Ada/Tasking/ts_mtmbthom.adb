with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Witness_Sets_io;
with Random_Test_Points;                 use Random_Test_Points;
with Multitasking_Membership_Tests;      use Multitasking_Membership_Tests;

procedure ts_mtmbthom is

-- DESCRIPTION :
--   Development of a multithreaded version of the homotopy membership test.

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a witness set, 
  --   prompts for the number of tasks,
  --   generates a random point in standard double precision
  --   and then launches the multitasked homotopy membership test.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols : Standard_Complex_Solutions.Solution_List;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    tol : constant double_float := 1.0E-6;
    dim,deg,nt,idx : natural32 := 0;

  begin
    Witness_Sets_io.Standard_Read_Embedding(lp,genpts,dim);
    deg := Standard_Complex_Solutions.Length_Of(genpts);
    new_line;
    put("Read a set of degree "); put(deg,1);
    put(" and dimension "); put(dim,1); put_line(".");
    new_line;
    put("Give the number of tasks : "); get(nt);
    new_line;
    put_line("Computing a random point on the solution set ...");
    sols := Standard_Random_Point(standard_output,lp.all,genpts,dim);
    ls := Standard_Complex_Solutions.Head_Of(sols);
    declare
      match : Standard_Complex_Vectors.Vector(ls.v'range);
    begin
      Standard_Membership_Test(nt,dim,tol,lp.all,genpts,ls.v,idx,match);
      if idx = 0 then
        put_line("The test point does not belong to the set???");
      else
        put("The test point found at "); put(idx,1); put_line(".");
        put_line("The test point :"); put_line(ls.v);
        put_line("The matching coordinates :"); put_line(match);
      end if;
    end;
  end Standard_Test;

  procedure Standard_Laurent_Test is

  -- DESCRIPTION :
  --   Prompts the user for a witness set of a Laurent system, 
  --   prompts for the number of tasks,
  --   generates a random point in standard double precision
  --   and then launches the multitasked homotopy membership test.

    lp : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    genpts,sols : Standard_Complex_Solutions.Solution_List;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    tol : constant double_float := 1.0E-6;
    dim,deg,nt,idx : natural32 := 0;

  begin
    Witness_Sets_io.Standard_Read_Embedding(lp,genpts,dim);
    deg := Standard_Complex_Solutions.Length_Of(genpts);
    new_line;
    put("Read a set of degree "); put(deg,1);
    put(" and dimension "); put(dim,1); put_line(".");
    new_line;
    put("Give the number of tasks : "); get(nt);
    new_line;
    put_line("Computing a random point on the solution set ...");
    sols := Standard_Random_Point(standard_output,lp.all,genpts,dim);
    ls := Standard_Complex_Solutions.Head_Of(sols);
    declare
      match : Standard_Complex_Vectors.Vector(ls.v'range);
    begin
      Standard_Membership_Test(nt,dim,tol,lp.all,genpts,ls.v,idx,match);
      if idx = 0 then
        put_line("The test point does not belong to the set???");
      else
        put("The test point found at "); put(idx,1); put_line(".");
        put_line("The test point :"); put_line(ls.v);
        put_line("The matching coordinates :"); put_line(match);
      end if;
    end;
  end Standard_Laurent_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a witness set, 
  --   prompts for the number of tasks,
  --   generates a random point in standard double precision
  --   and then launches the multitasked homotopy membership test
  --   in double double precision.

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols : DoblDobl_Complex_Solutions.Solution_List;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    tol : constant double_float := 1.0E-6;
    dim,deg,nt,idx : natural32 := 0;

  begin
    Witness_Sets_io.DoblDobl_Read_Embedding(lp,genpts,dim);
    deg := DoblDobl_Complex_Solutions.Length_Of(genpts);
    new_line;
    put("Read a set of degree "); put(deg,1);
    put(" and dimension "); put(dim,1); put_line(".");
    new_line;
    put("Give the number of tasks : "); get(nt);
    new_line;
    put_line("Computing a random point on the solution set ...");
    sols := DoblDobl_Random_Point(standard_output,lp.all,genpts,dim);
    ls := DoblDobl_Complex_Solutions.Head_Of(sols);
    declare
      match : DoblDobl_Complex_Vectors.Vector(ls.v'range);
    begin
      DoblDobl_Membership_Test(nt,dim,tol,lp.all,genpts,ls.v,idx,match);
      if idx = 0 then
        put_line("The test point does not belong to the set???");
      else
        put("The test point found at "); put(idx,1); put_line(".");
        put_line("The test point :"); put_line(ls.v);
        put_line("The matching coordinates :"); put_line(match);
      end if;
    end;
  end DoblDobl_Test;

  procedure DoblDobl_Laurent_Test is

  -- DESCRIPTION :
  --   Prompts the user for a witness set of a Laurent system, 
  --   prompts for the number of tasks,
  --   generates a random point in standard double precision
  --   and then launches the multitasked homotopy membership test
  --   in double double precision.

    lp : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    genpts,sols : DoblDobl_Complex_Solutions.Solution_List;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    tol : constant double_float := 1.0E-6;
    dim,deg,nt,idx : natural32 := 0;

  begin
    Witness_Sets_io.DoblDobl_Read_Embedding(lp,genpts,dim);
    deg := DoblDobl_Complex_Solutions.Length_Of(genpts);
    new_line;
    put("Read a set of degree "); put(deg,1);
    put(" and dimension "); put(dim,1); put_line(".");
    new_line;
    put("Give the number of tasks : "); get(nt);
    new_line;
    put_line("Computing a random point on the solution set ...");
    sols := DoblDobl_Random_Point(standard_output,lp.all,genpts,dim);
    ls := DoblDobl_Complex_Solutions.Head_Of(sols);
    declare
      match : DoblDobl_Complex_Vectors.Vector(ls.v'range);
    begin
      DoblDobl_Membership_Test(nt,dim,tol,lp.all,genpts,ls.v,idx,match);
      if idx = 0 then
        put_line("The test point does not belong to the set???");
      else
        put("The test point found at "); put(idx,1); put_line(".");
        put_line("The test point :"); put_line(ls.v);
        put_line("The matching coordinates :"); put_line(match);
      end if;
    end;
  end DoblDobl_Laurent_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a witness set, 
  --   prompts for the number of tasks,
  --   generates a random point in standard double precision
  --   and then launches the multitasked homotopy membership test
  --   in quad double precision.

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols : QuadDobl_Complex_Solutions.Solution_List;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    tol : constant double_float := 1.0E-6;
    dim,deg,nt,idx : natural32 := 0;

  begin
    Witness_Sets_io.QuadDobl_Read_Embedding(lp,genpts,dim);
    deg := QuadDobl_Complex_Solutions.Length_Of(genpts);
    new_line;
    put("Read a set of degree "); put(deg,1);
    put(" and dimension "); put(dim,1); put_line(".");
    new_line;
    put("Give the number of tasks : "); get(nt);
    new_line;
    put_line("Computing a random point on the solution set ...");
    sols := QuadDobl_Random_Point(standard_output,lp.all,genpts,dim);
    ls := QuadDobl_Complex_Solutions.Head_Of(sols);
    declare
      match : QuadDobl_Complex_Vectors.Vector(ls.v'range);
    begin
      QuadDobl_Membership_Test(nt,dim,tol,lp.all,genpts,ls.v,idx,match);
      if idx = 0 then
        put_line("The test point does not belong to the set???");
      else
        put("The test point found at "); put(idx,1); put_line(".");
        put_line("The test point :"); put_line(ls.v);
        put_line("The matching coordinates :"); put_line(match);
      end if;
    end;
  end QuadDobl_Test;

  procedure QuadDobl_Laurent_Test is

  -- DESCRIPTION :
  --   Prompts the user for a witness set of a Laurent system,
  --   prompts for the number of tasks,
  --   generates a random point in standard double precision
  --   and then launches the multitasked homotopy membership test
  --   in quad double precision.

    lp : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    genpts,sols : QuadDobl_Complex_Solutions.Solution_List;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    tol : constant double_float := 1.0E-6;
    dim,deg,nt,idx : natural32 := 0;

  begin
    Witness_Sets_io.QuadDobl_Read_Embedding(lp,genpts,dim);
    deg := QuadDobl_Complex_Solutions.Length_Of(genpts);
    new_line;
    put("Read a set of degree "); put(deg,1);
    put(" and dimension "); put(dim,1); put_line(".");
    new_line;
    put("Give the number of tasks : "); get(nt);
    new_line;
    put_line("Computing a random point on the solution set ...");
    sols := QuadDobl_Random_Point(standard_output,lp.all,genpts,dim);
    ls := QuadDobl_Complex_Solutions.Head_Of(sols);
    declare
      match : QuadDobl_Complex_Vectors.Vector(ls.v'range);
    begin
      QuadDobl_Membership_Test(nt,dim,tol,lp.all,genpts,ls.v,idx,match);
      if idx = 0 then
        put_line("The test point does not belong to the set???");
      else
        put("The test point found at "); put(idx,1); put_line(".");
        put_line("The test point :"); put_line(ls.v);
        put_line("The matching coordinates :"); put_line(match);
      end if;
    end;
  end QuadDobl_Laurent_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision
  --   and then calls the proper test.

    ans : character;

  begin
    new_line;
    put_line("MENU for the multitasking homotopy membership test :");
    put_line("  0. standard double precision, for ordinary polynomial system");
    put_line("  1. double double precision, for ordinary polynomial system");
    put_line("  2. quad double precision, for ordinary polynomial system");
    put_line("  3. standard double precision, for Laurent polynomial system");
    put_line("  4. double double precision, for Laurent polynomial system");
    put_line("  5. quad double precision, for Laurent polynomial system");
    put("Type 0, 1, 2, 3, 4, or 5 to select the test : ");
    Ask_Alternative(ans,"012345");
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when '3' => Standard_Laurent_Test;
      when '4' => DoblDobl_Laurent_Test;
      when '5' => QuadDobl_Laurent_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtmbthom;

