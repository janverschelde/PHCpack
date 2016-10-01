with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Standard_Dense_Series_VecVecs;
with DoblDobl_Dense_Series_VecVecs;
with QuadDobl_Dense_Series_VecVecs;
with Regular_Solution_Curves_Series;     use Regular_Solution_Curves_Series;

procedure ts_puiseux is

-- DESCRIPTION :
--   Development of the Newton-Puiseux algorithm,
--   for regular solution curves, defined by complete intersections,
--   in Noether position, with sufficiently general coefficients.

  procedure Tropisms_by_Mixed_Cells
              ( file : in file_type; 
                sup : in out Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32 ) is

  -- DESCRIPTION :
  --   Computes a regular mixed cell configuration for
  --   the supports in sup, with some writing to file.
  --   The mixed volume is in mv.

  begin
    Mixed_Cell_Tropisms(file,sup,mcc,mv);
    new_line(file);
    put(file,"The number of pretropisms : ");
    put(file,Length_Of(mcc),1); new_line(file);
    put(file,"The number of series : ");
    put(file,mv,1); new_line(file);
  end Tropisms_by_Mixed_Cells;

  procedure Tropisms_by_Mixed_Cells
              ( sup : in out Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                report : in boolean ) is

  -- DESCRIPTION :
  --   Computes a regular mixed cell configuration for
  --   the supports in sup, with some writing to screen if report.
  --   The mixed volume is in mv.

  begin
    Mixed_Cell_Tropisms(report,sup,mcc,mv);
    if report then
      new_line;
      put("The number of pretropisms : ");
      put(Length_Of(mcc),1); new_line;
      put("The number of series : "); put(mv,1); new_line;
    end if;
  end Tropisms_by_Mixed_Cells;

  procedure Standard_Test
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and the series,
  --   in standard double precision.  The output is written to file.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;

  begin
    Tropisms_by_Mixed_Cells(file,sup,mcc,mv);
    declare
      s : Standard_Dense_Series_VecVecs.VecVec(1..integer32(mv));
    begin
      s := Series(file,p,mcc,mv,nit);
    end;
  end Standard_Test;

  procedure DoblDobl_Test
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and the series,
  --   in double double precision.  The output is written to file.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;

  begin
    Tropisms_by_Mixed_Cells(file,sup,mcc,mv);
    declare
      s : DoblDobl_Dense_Series_VecVecs.VecVec(1..integer32(mv));
    begin
      s := Series(file,p,mcc,mv,nit);
    end;
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and the series,
  --   in quad double precision.  The output is written to file.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;

  begin
    Tropisms_by_Mixed_Cells(file,sup,mcc,mv);
    declare
      s : QuadDobl_Dense_Series_VecVecs.VecVec(1..integer32(mv));
    begin
      s := Series(file,p,mcc,mv,nit);
    end;
  end QuadDobl_Test;

  procedure Standard_Test 
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                report : in boolean ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and the series,
  --   in standard double precision.
  --   The output is written to screen if report.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;

  begin
    Tropisms_by_Mixed_Cells(sup,mcc,mv,report);
    declare
      s : Standard_Dense_Series_VecVecs.VecVec(1..integer32(mv));
    begin
      s := Series(p,mcc,mv,nit,report);
    end;
  end Standard_Test;

  procedure DoblDobl_Test 
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                report : in boolean ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and the series,
  --   in double double precision.
  --   The output is written to screen if report.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;

  begin
    Tropisms_by_Mixed_Cells(sup,mcc,mv,report);
    declare
      s : DoblDobl_Dense_Series_VecVecs.VecVec(1..integer32(mv));
    begin
      s := Series(p,mcc,mv,nit,report);
    end;
  end DoblDobl_Test;

  procedure QuadDobl_Test 
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                report : in boolean ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and the series,
  --   in quad double precision.
  --   The output is written to screen if report.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;

  begin
    Tropisms_by_Mixed_Cells(sup,mcc,mv,report);
    declare
      s : QuadDobl_Dense_Series_VecVecs.VecVec(1..integer32(mv));
    begin
      s := Series(p,mcc,mv,nit,report);
    end;
  end QuadDobl_Test;

  function Prompt_for_Precision return character is

  -- DESCRIPTION :
  --   Displays the menu for the available precision and
  --   returns '0', '1', or '2' for double, double double,
  --   or quad double precision.

    res : character;

  begin
    new_line;
    put_line("MENU to set the precision level :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(res,"012");
    return res;
  end Prompt_for_Precision;

  procedure Standard_Read
              ( lp : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                nq,nv : out integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent system and returns the system,
  --   along with the number of polynomials in nq
  --   and the number of variables in nv.

    use Standard_Complex_Laurentials;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    get(lp);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("Number of polynomials : "); put(nq,1); new_line;
    put("Number of variables : "); put(nv,1); new_line;
  end Standard_Read;

  procedure DoblDobl_Read
              ( lp : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                nq,nv : out integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent system and returns the system,
  --   along with the number of polynomials in nq
  --   and the number of variables in nv.

    use DoblDobl_Complex_Laurentials;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    get(lp);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("Number of polynomials : "); put(nq,1); new_line;
    put("Number of variables : "); put(nv,1); new_line;
  end DoblDobl_Read;

  procedure QuadDobl_Read
              ( lp : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                nq,nv : out integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent system and returns the system,
  --   along with the number of polynomials in nq
  --   and the number of variables in nv.

    use QuadDobl_Complex_Laurentials;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    get(lp);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("Number of polynomials : "); put(nq,1); new_line;
    put("Number of variables : "); put(nv,1); new_line;
  end QuadDobl_Read;
  
  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system
  --   and checks whether the n polynomials have n+1 variables.
  --   Computations are done in standard double precision.

    lp : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    nq,nv : integer32;
    ans : character;
    report : boolean;
    file : file_type;

  begin
    Standard_Read(lp,nq,nv);
    if nv /= nq+1 then
      put(nv,1); put(" /= "); put(nq,1); put(" + 1");
    else
      new_line;
      put("Do you want intermediate output to file ? (y/n) ");
      Ask_Yes_or_No(ans); report := (ans = 'y');
      if report then
        new_line;
        put_line("Reading the name of the output file ...");
        Read_Name_and_Create_File(file);
        put(file,natural32(nq),natural32(nv),lp.all);
        Standard_Test(file,lp.all);
      else
        new_line;
        put("Do you want intermediate output to screen ? (y/n) ");
        Ask_Yes_or_No(ans); report := (ans = 'y');
        Standard_Test(lp.all,report);
      end if;
    end if;
  end Standard_Main;
  
  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system
  --   and checks whether the n polynomials have n+1 variables.
  --   Computations are done in double double precision.

    lp : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    nq,nv : integer32;
    ans : character;
    report : boolean;
    file : file_type;

  begin
    DoblDobl_Read(lp,nq,nv);
    if nv /= nq+1 then
      put(nv,1); put(" /= "); put(nq,1); put(" + 1");
    else
      new_line;
      put("Do you want intermediate output to file ? (y/n) ");
      Ask_Yes_or_No(ans); report := (ans = 'y');
      if report then
        new_line;
        put_line("Reading the name of the output file ...");
        Read_Name_and_Create_File(file);
        put(file,natural32(nq),natural32(nv),lp.all);
        DoblDobl_Test(file,lp.all);
      else
        new_line;
        put("Do you want intermediate output to screen ? (y/n) ");
        Ask_Yes_or_No(ans); report := (ans = 'y');
        DoblDobl_Test(lp.all,report);
      end if;
    end if;
  end DoblDobl_Main;
  
  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system
  --   and checks whether the n polynomials have n+1 variables.
  --   Computations are done in quad double precision.

    lp : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    nq,nv : integer32;
    ans : character;
    report : boolean;
    file : file_type;

  begin
    QuadDobl_Read(lp,nq,nv);
    if nv /= nq+1 then
      put(nv,1); put(" /= "); put(nq,1); put(" + 1");
    else
      new_line;
      put("Do you want intermediate output to file ? (y/n) ");
      Ask_Yes_or_No(ans); report := (ans = 'y');
      if report then
        new_line;
        put_line("Reading the name of the output file ...");
        Read_Name_and_Create_File(file);
        put(file,natural32(nq),natural32(nv),lp.all);
        QuadDobl_Test(file,lp.all);
      else
        new_line;
        put("Do you want intermediate output to screen ? (y/n) ");
        Ask_Yes_or_No(ans); report := (ans = 'y');
        QuadDobl_Test(lp.all,report);
      end if;
    end if;
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision
  --   and then launches the proper main driver.

    prc : constant character := Prompt_for_Precision;

  begin
    case prc is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_puiseux;
