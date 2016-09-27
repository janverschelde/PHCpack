with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Regular_Solution_Curves_Series;     use Regular_Solution_Curves_Series;

procedure ts_puiseux is

-- DESCRIPTION :
--   Development of the Newton-Puiseux algorithm.

  procedure Test ( file : in file_type; p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and initial terms.
  --   The output is written to file.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    nit : constant integer32 := 7;
    mv : natural32;

  begin
    Mixed_Cell_Tropisms(file,sup,mcc,mv);
    new_line(file);
    put(file,"The number of pretropisms : ");
    put(file,Length_Of(mcc),1); new_line(file);
    put(file,"The number of series : ");
    put(file,mv,1); new_line(file);
    Series(file,p,mcc,nit);
  end Test;

  procedure Test ( p : in Laur_Sys; report : in boolean ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and initial terms.
  --   The output is written to screen if report.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    nit : constant integer32 := 7;
    mv : natural32;

  begin
    Mixed_Cell_Tropisms(report,sup,mcc,mv);
    if report then
      new_line;
      put("The number of pretropisms : ");
      put(Length_Of(mcc),1); new_line;
      put("The number of series : "); put(mv,1); new_line;
    end if;
    Series(p,mcc,nit,report);
  end Test;
  
  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system
  --   and checks whether the n polynomials have n+1 variables.

    use Standard_Complex_Laurentials;

    lp : Link_to_Laur_Sys;
    nq,nv : integer32;
    ans : character;
    report : boolean;
    file : file_type;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    get(lp);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("Number of polynomials : "); put(nq,1); new_line;
    put("Number of variables : "); put(nv,1); new_line;
    if nv /= nq+1 then
      put(nv,1); put(" /= "); put(nq,1); put(" + 1");
    else
      new_line;
      put("Do you want intermediate output to file ? (y/n) ");
      Ask_Yes_or_No(ans);
      report := (ans = 'y');
      if report then
        new_line;
        put_line("Reading the name of the output file ...");
        Read_Name_and_Create_File(file);
        put(file,natural32(nq),natural32(nv),lp.all);
        Test(file,lp.all);
      else
        new_line;
        put("Do you want intermediate output to screen ? (y/n) ");
        Ask_Yes_or_No(ans);
        report := (ans = 'y');
        Test(lp.all,report);
      end if;
    end if;
  end Main;

begin
  Main;
end ts_puiseux;
