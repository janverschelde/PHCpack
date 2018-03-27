with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Witness_Sets_io;
with Homotopy_Membership_Tests;          use Homotopy_Membership_Tests;
with Homotopy_Membership_Filters;
with Random_Test_Points;                 use Random_Test_Points;

procedure ts_mbthom is

-- DESCRIPTION :
--   Test on membership with homotopy.

  procedure Standard_Membership is

  -- DESCRIPTION :
  --   Prompts the user for a witness set, a list of test points,
  --   and then applies the homotopy membership test
  --   using standard double precision arithmetic.

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols : Standard_Complex_Solutions.Solution_List;
    dim : natural32;
    restol : constant double_float := 1.0E-10;
    homtol : constant double_float := 1.0E-6;
    out2file : boolean;
    ans : character;

  begin
    Witness_Sets_io.Standard_Read_Embedding(lp,genpts,dim);
    new_line;
    put("Do you want the output to a file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans /= 'y' then
      out2file := false;
    else
      out2file := true;
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(file);
    end if;
    new_line;
    put_line("MENU for test points : ");
    put_line("  0. generate a random point on the solution set; or");
    put_line("  1. give a file name with a list of test points.");
    put("Type 0 or 1 for test points : "); Ask_Alternative(ans,"01");
    new_line;
    if ans = '1' then
      Read(sols);
    else
      put_line("Computing a random point on the solution set ...");
      if out2file
       then sols := Standard_Random_Point(file,lp.all,genpts,dim);
       else sols := Standard_Random_Point(standard_output,lp.all,genpts,dim);
      end if;
    end if;
    new_line;
    put_line("See the output file for results...");
    new_line;
    if out2file
     then Homotopy_Membership_Test(file,lp.all,dim,genpts,sols,restol,homtol);
     else Homotopy_Membership_Test(true,lp.all,dim,genpts,sols,restol,homtol);
    end if;
  end Standard_Membership;

  procedure Standard_Membership_Filter is

  -- DESCRIPTION :
  --   Prompts the user for a witness set, a list of test points,
  --   and then applies the homotopy membership test
  --   using standard double precision arithmetic.

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols,mempts,outpts : Standard_Complex_Solutions.Solution_List;
    dim : natural32;
    restol : constant double_float := 1.0E-10;
    rcotol : constant double_float := 1.0E-6;
    homtol : constant double_float := 1.0E-6;
    out2file : boolean;
    ans : character;

  begin
    Witness_Sets_io.Standard_Read_Embedding(lp,genpts,dim);
    new_line;
    put("Do you want the output to a file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans /= 'y' then
      out2file := false;
    else
      out2file := true;
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(file);
    end if;
    new_line;
    put_line("MENU for test points : ");
    put_line("  0. generate a random point on the solution set; or");
    put_line("  1. give a file name with a list of test points.");
    put("Type 0 or 1 for test points : "); Ask_Alternative(ans,"01");
    new_line;
    if ans = '1' then
      Read(sols);
    else
      put_line("Computing a random point on the solution set ...");
      if out2file
       then sols := Standard_Random_Point(file,lp.all,genpts,dim);
       else sols := Standard_Random_Point(standard_output,lp.all,genpts,dim);
      end if;
    end if;
    if out2file then
      new_line;
      put_line("See the output file for results ...");
      new_line;
      Homotopy_Membership_Filters.Filter
        (file,lp.all,genpts,dim,rcotol,restol,homtol,sols,mempts,outpts);
      put(file,"Number of points that are member : ");
      put(file,Standard_Complex_Solutions.Length_Of(mempts),1);
      new_line(file);
      put(file,"Number of points that are outside : ");
      put(file,Standard_Complex_Solutions.Length_Of(outpts),1);
      new_line(file);
    else
      new_line;
      put("Verbose ? (y/n) ");
      Ask_Yes_or_No(ans);
      Homotopy_Membership_Filters.Filter
        (ans = 'y',lp.all,genpts,dim,rcotol,restol,homtol,sols,mempts,outpts);
      put("Number of points that are member : ");
      put(Standard_Complex_Solutions.Length_Of(mempts),1); new_line;
      put("Number of points that are outside : ");
      put(Standard_Complex_Solutions.Length_Of(outpts),1); new_line;
    end if;
  end Standard_Membership_Filter;

  procedure DoblDobl_Membership is

  -- DESCRIPTION :
  --   Prompts the user for a witness set, a list of test points,
  --   and then applies the homotopy membership test
  --   using double double precision arithmetic.

    file : file_type;
    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols : DoblDobl_Complex_Solutions.Solution_List;
    dim : natural32;
    restol : constant double_float := 1.0E-10;
    homtol : constant double_float := 1.0E-6;
    out2file : boolean;
    ans : character;

  begin
    Witness_Sets_io.DoblDobl_Read_Embedding(lp,genpts,dim);
    new_line;
    put("Do you want the output to a file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans /= 'y' then
      out2file := false;
    else
      out2file := true;
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(file);
    end if;
    new_line;
    put_line("MENU for test points : ");
    put_line("  0. generate a random point on the solution set; or");
    put_line("  1. give a file name with a list of test points.");
    put("Type 0 or 1 for test points : "); Ask_Alternative(ans,"01");
    new_line;
    if ans = '1' then
      Read(sols);
    else
      put_line("Computing a random point on the solution set ...");
      if out2file
       then sols := DoblDobl_Random_Point(file,lp.all,genpts,dim);
       else sols := DoblDobl_Random_Point(standard_output,lp.all,genpts,dim);
      end if;
    end if;
    new_line;
    put_line("See the output file for results...");
    new_line;
    if out2file
     then Homotopy_Membership_Test(file,lp.all,dim,genpts,sols,restol,homtol);
     else Homotopy_Membership_Test(true,lp.all,dim,genpts,sols,restol,homtol);
    end if;
  end DoblDobl_Membership;

  procedure QuadDobl_Membership is

  -- DESCRIPTION :
  --   Prompts the user for a witness set, a list of test points,
  --   and then applies the homotopy membership test
  --   using quad double precision arithmetic.

    file : file_type;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols : QuadDobl_Complex_Solutions.Solution_List;
    dim : natural32;
    restol : constant double_float := 1.0E-10;
    homtol : constant double_float := 1.0E-6;
    out2file : boolean;
    ans : character;

  begin
    Witness_Sets_io.QuadDobl_Read_Embedding(lp,genpts,dim);
    new_line;
    put("Do you want the output to a file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans /= 'y' then
      out2file := false;
    else
      out2file := true;
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(file);
    end if;
    new_line;
    put_line("MENU for test points : ");
    put_line("  0. generate a random point on the solution set; or");
    put_line("  1. give a file name with a list of test points.");
    put("Type 0 or 1 for test points : "); Ask_Alternative(ans,"01");
    new_line;
    if ans = '1' then
      Read(sols);
    else
      put_line("Computing a random point on the solution set ...");
      if out2file
       then sols := QuadDobl_Random_Point(file,lp.all,genpts,dim);
       else sols := QuadDobl_Random_Point(standard_output,lp.all,genpts,dim);
      end if;
    end if;
    new_line;
    put_line("See the output file for results...");
    new_line;
    if out2file
     then Homotopy_Membership_Test(file,lp.all,dim,genpts,sols,restol,homtol);
     else Homotopy_Membership_Test(true,lp.all,dim,genpts,sols,restol,homtol);
    end if;
  end QuadDobl_Membership;

  procedure Standard_Laurent_Membership is

  -- DESCRIPTION :
  --   Prompts the user for a witness set, a list of test points,
  --   and then applies the homotopy membership test
  --   using standard double precision arithmetic.

    file : file_type;
    lp : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    genpts,sols : Standard_Complex_Solutions.Solution_List;
    dim : natural32;
    restol : constant double_float := 1.0E-10;
    homtol : constant double_float := 1.0E-6;
    out2file : boolean;
    ans : character;

  begin
    Witness_Sets_io.Standard_Read_Embedding(lp,genpts,dim);
    new_line;
    put("Do you want the output to a file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans /= 'y' then
      out2file := false;
    else
      out2file := true;
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(file);
    end if;
    new_line;
    put_line("MENU for test points : ");
    put_line("  0. generate a random point on the solution set; or");
    put_line("  1. give a file name with a list of test points.");
    put("Type 0 or 1 for test points : "); Ask_Alternative(ans,"01");
    new_line;
    if ans = '1' then
      Read(sols);
    else
      put_line("Computing a random point on the solution set ...");
      if out2file
       then sols := Standard_Random_Point(file,lp.all,genpts,dim);
       else sols := Standard_Random_Point(standard_output,lp.all,genpts,dim);
      end if;
    end if;
    new_line;
    put_line("See the output file for results...");
    new_line;
    if out2file
     then Homotopy_Membership_Test(file,lp.all,dim,genpts,sols,restol,homtol);
     else Homotopy_Membership_Test(true,lp.all,dim,genpts,sols,restol,homtol);
    end if;
  end Standard_Laurent_Membership;

  procedure DoblDobl_Laurent_Membership is

  -- DESCRIPTION :
  --   Prompts the user for a witness set, a list of test points,
  --   and then applies the homotopy membership test
  --   using double double precision arithmetic.

    file : file_type;
    lp : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    genpts,sols : DoblDobl_Complex_Solutions.Solution_List;
    dim : natural32;
    restol : constant double_float := 1.0E-10;
    homtol : constant double_float := 1.0E-6;
    out2file : boolean;
    ans : character;

  begin
    Witness_Sets_io.DoblDobl_Read_Embedding(lp,genpts,dim);
    new_line;
    put("Do you want the output to a file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans /= 'y' then
      out2file := false;
    else
      out2file := true;
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(file);
    end if;
    new_line;
    put_line("MENU for test points : ");
    put_line("  0. generate a random point on the solution set; or");
    put_line("  1. give a file name with a list of test points.");
    put("Type 0 or 1 for test points : "); Ask_Alternative(ans,"01");
    new_line;
    if ans = '1' then
      Read(sols);
    else
      put_line("Computing a random point on the solution set ...");
      if out2file
       then sols := DoblDobl_Random_Point(file,lp.all,genpts,dim);
       else sols := DoblDobl_Random_Point(standard_output,lp.all,genpts,dim);
      end if;
    end if;
    new_line;
    put_line("See the output file for results...");
    new_line;
    if out2file
     then Homotopy_Membership_Test(file,lp.all,dim,genpts,sols,restol,homtol);
     else Homotopy_Membership_Test(true,lp.all,dim,genpts,sols,restol,homtol);
    end if;
  end DoblDobl_Laurent_Membership;

  procedure QuadDobl_Laurent_Membership is

  -- DESCRIPTION :
  --   Prompts the user for a witness set, a list of test points,
  --   and then applies the homotopy membership test
  --   using quad double precision arithmetic.

    file : file_type;
    lp : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    genpts,sols : QuadDobl_Complex_Solutions.Solution_List;
    dim : natural32;
    restol : constant double_float := 1.0E-10;
    homtol : constant double_float := 1.0E-6;
    out2file : boolean;
    ans : character;

  begin
    Witness_Sets_io.QuadDobl_Read_Embedding(lp,genpts,dim);
    new_line;
    put("Do you want the output to a file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans /= 'y' then
      out2file := false;
    else
      out2file := true;
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(file);
    end if;
    new_line;
    put_line("MENU for test points : ");
    put_line("  0. generate a random point on the solution set; or");
    put_line("  1. give a file name with a list of test points.");
    put("Type 0 or 1 for test points : "); Ask_Alternative(ans,"01");
    new_line;
    if ans = '1' then
      Read(sols);
    else
      put_line("Computing a random point on the solution set ...");
      if out2file
       then sols := QuadDobl_Random_Point(file,lp.all,genpts,dim);
       else sols := QuadDobl_Random_Point(standard_output,lp.all,genpts,dim);
      end if;
    end if;
    new_line;
    put_line("See the output file for results...");
    new_line;
    if out2file
     then Homotopy_Membership_Test(file,lp.all,dim,genpts,sols,restol,homtol);
     else Homotopy_Membership_Test(true,lp.all,dim,genpts,sols,restol,homtol);
    end if;
  end QuadDobl_Laurent_Membership;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Membership test with homotopy :");
    put_line("  Input : embedded polynomial system with generic points, and");
    put_line("          list of test points.");
    put_line("  Output : decision whether test point lies on component.");
    new_line;
    put_line("MENU to choose the precision : ");
    put_line("  0. standard double precision homotopy continuation;");
    put_line("  1. double double precision homotopy continuation;");
    put_line("  2. quad double precision homotopy continuation;");
    put_line("  3. Laurent homotopy in standard double precision;");
    put_line("  4. Laurent homotopy in double double precision;");
    put_line("  5. Laurent homotopy in quad double precision;");
    put_line("  6. test membership filter in standard double precision.");
    put("Type 0, 1, 2, 3, 4, 5, or 6 to select the test : ");
    Ask_Alternative(ans,"0123456");
    case ans is
      when '0' => Standard_Membership;
      when '1' => DoblDobl_Membership;
      when '2' => QuadDobl_Membership;
      when '3' => Standard_Laurent_Membership;
      when '4' => DoblDobl_Laurent_Membership;
      when '5' => QuadDobl_Laurent_Membership;
      when '6' => Standard_Membership_Filter;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mbthom;
