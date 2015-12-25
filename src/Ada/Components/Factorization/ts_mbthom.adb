with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Witness_Sets;
with Witness_Sets_io;                    use Witness_Sets_io;
with Sampling_Machine;
with DoblDobl_Sampling_Machine;
with QuadDobl_Sampling_Machine;
with Homotopy_Membership_Tests;          use Homotopy_Membership_Tests;

procedure ts_mbthom is

-- DESCRIPTION :
--   Test on membership with homotopy.

  function Standard_Random_Point
             ( file : file_type;
               p : Standard_Complex_Poly_Systems.Poly_Sys;
               s : Standard_Complex_Solutions.Solution_List;
               d : natural32 )
             return Standard_Complex_Solutions.Solution_List is

  -- DESCRIPTION :
  --   For a given witness set, takes the first point and applies
  --   path tracking to another set of slices to compute another point.

  -- ON ENTRY :
  --   file    for output during sampling;
  --   p       embedded polynomial system;
  --   s       generic points on the solution set;
  --   d       dimension of the solution set;

    res : Standard_Complex_Solutions.Solution_List;
    n : constant natural32 := natural32(p'last);
    hyp : Standard_Complex_VecVecs.VecVec(1..integer32(d))
        := Witness_Sets.Random_Hyperplanes(d,n);
    startsol : Standard_Complex_Solutions.Solution(integer32(n));
    newsol : Standard_Complex_Solutions.Solution(integer32(n));

  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner;
    startsol := Standard_Complex_Solutions.Head_Of(s).all;
    startsol.t := Standard_Complex_Numbers.Create(0.0);
    Sampling_Machine.Sample(file,true,startsol,hyp,newsol);
    Standard_Complex_Solutions.Add(res,newsol);
    Sampling_Machine.Clear;
    return res;
  end Standard_Random_Point;

  function DoblDobl_Random_Point
             ( file : file_type;
               p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               s : DoblDobl_Complex_Solutions.Solution_List;
               d : natural32 )
             return DoblDobl_Complex_Solutions.Solution_List is

  -- DESCRIPTION :
  --   For a given witness set, takes the first point and applies
  --   path tracking to another set of slices to compute another point.

  -- ON ENTRY :
  --   file    for output during sampling;
  --   p       embedded polynomial system;
  --   s       generic points on the solution set;
  --   d       dimension of the solution set;

    res : DoblDobl_Complex_Solutions.Solution_List;
    n : constant natural32 := natural32(p'last);
    hyp : DoblDobl_Complex_VecVecs.VecVec(1..integer32(d))
        := Witness_Sets.Random_Hyperplanes(d,n);
    startsol : DoblDobl_Complex_Solutions.Solution(integer32(n));
    newsol : DoblDobl_Complex_Solutions.Solution(integer32(n));
    ddzero : double_double := Double_Double_Numbers.Create(0.0);

  begin
    DoblDobl_Sampling_Machine.Initialize(p);
    DoblDobl_Sampling_Machine.Default_Tune_Sampler(0);
    DoblDobl_Sampling_Machine.Default_Tune_Refiner;
    startsol := DoblDobl_Complex_Solutions.Head_Of(s).all;
    startsol.t := DoblDobl_Complex_Numbers.Create(ddzero);
    DoblDobl_Sampling_Machine.Sample(file,true,startsol,hyp,newsol);
    DoblDobl_Complex_Solutions.Add(res,newsol);
    DoblDobl_Sampling_Machine.Clear;
    return res;
  end DoblDobl_Random_Point;

  function QuadDobl_Random_Point
             ( file : file_type;
               p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               s : QuadDobl_Complex_Solutions.Solution_List;
               d : natural32 )
             return QuadDobl_Complex_Solutions.Solution_List is

  -- DESCRIPTION :
  --   For a given witness set, takes the first point and applies
  --   path tracking to another set of slices to compute another point.

  -- ON ENTRY :
  --   file    for output during sampling;
  --   p       embedded polynomial system;
  --   s       generic points on the solution set;
  --   d       dimension of the solution set;

    res : QuadDobl_Complex_Solutions.Solution_List;
    n : constant natural32 := natural32(p'last);
    hyp : QuadDobl_Complex_VecVecs.VecVec(1..integer32(d))
        := Witness_Sets.Random_Hyperplanes(d,n);
    startsol : QuadDobl_Complex_Solutions.Solution(integer32(n));
    newsol : QuadDobl_Complex_Solutions.Solution(integer32(n));
    qdzero : quad_double := Quad_Double_Numbers.Create(0.0);

  begin
    QuadDobl_Sampling_Machine.Initialize(p);
    QuadDobl_Sampling_Machine.Default_Tune_Sampler(0);
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
    startsol := QuadDobl_Complex_Solutions.Head_Of(s).all;
    startsol.t := QuadDobl_Complex_Numbers.Create(qdzero);
    QuadDobl_Sampling_Machine.Sample(file,true,startsol,hyp,newsol);
    QuadDobl_Complex_Solutions.Add(res,newsol);
    QuadDobl_Sampling_Machine.Clear;
    return res;
  end QuadDobl_Random_Point;

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
    Standard_Read_Embedding(lp,genpts,dim);
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
      sols := Standard_Random_Point(file,lp.all,genpts,dim);
    end if;
    new_line;
    put_line("See the output file for results...");
    new_line;
    if out2file
     then Homotopy_Membership_Test(file,lp.all,dim,genpts,sols,restol,homtol);
     else Homotopy_Membership_Test(true,lp.all,dim,genpts,sols,restol,homtol);
    end if;
  end Standard_Membership;

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
    DoblDobl_Read_Embedding(lp,genpts,dim);
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
      sols := DoblDobl_Random_Point(file,lp.all,genpts,dim);
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
    QuadDobl_Read_Embedding(lp,genpts,dim);
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
      sols := QuadDobl_Random_Point(file,lp.all,genpts,dim);
    end if;
    new_line;
    put_line("See the output file for results...");
    new_line;
    if out2file
     then Homotopy_Membership_Test(file,lp.all,dim,genpts,sols,restol,homtol);
     else Homotopy_Membership_Test(true,lp.all,dim,genpts,sols,restol,homtol);
    end if;
  end QuadDobl_Membership;

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
    put_line("  0. standard double precision homotopy continuation; or");
    put_line("  1. double double precision homotopy continuation; or");
    put_line("  2. quad double precision homotopy continuation.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Membership;
      when '1' => DoblDobl_Membership;
      when '2' => QuadDobl_Membership;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mbthom;
