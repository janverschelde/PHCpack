with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists;    use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io; use Arrays_of_Integer_Vector_Lists_io;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laur_SysFun;      use Standard_Complex_Laur_SysFun;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Supports_of_Polynomial_Systems;    use Supports_of_Polynomial_Systems;
with Integer_Mixed_Subdivisions;        use Integer_Mixed_Subdivisions;
with Drivers_for_Static_Lifting;        use Drivers_for_Static_Lifting;
with Black_Box_Solvers;                 use Black_Box_Solvers;

procedure ts_puiseux is

-- DESCRIPTION :
--   Development of the Newton-Puiseux algorithm.

  procedure Mixed_Cell_Tropisms
              ( report : in boolean;
                sup : in out Array_of_Lists;
                mcc : out Mixed_Subdivision;
                mv : out natural32 ) is

  -- DESCRIPTION :
  --   Given a system of n Laurent polynomials in n+1 variables,
  --   computes the tropisms, where the last variable is the parameter.
  --   The tropisms are computed via a regular mixed cell configuration,
  --   induced by the lifting defined by the last variable.
  --   As the lifting may differ, even if several supports would be the same,
  --   the type of mixed is assumed to be fully mixed.

  -- ON ENTRY :
  --   report   if intermediate output has to be written to screen;
  --   sup      a list of n supports in n+1 variables.

  -- ON RETURN :
  --   mcc      a mixed cell configuration induced by the lifting
  --            defined by the last exponent in each monomial of p;
  --   mv       the mixed volume of the cells in mcc.

    dim : constant integer32 := sup'last;
    mix : Standard_Integer_Vectors.Vector(1..dim);

  begin
    if report
     then put_line("The supports : "); put(sup);
    end if;
    mix := (mix'range => 1);
    if report then
      Integer_Create_Mixed_Cells(standard_output,dim,mix,false,sup,mcc);
      Integer_Volume_Computation(standard_output,dim,mix,true,sup,mcc,mv);
    else
      Integer_Create_Mixed_Cells(dim,mix,sup,mcc);
      Integer_Volume_Computation(dim,mix,true,sup,mcc,mv);
    end if;
  end Mixed_Cell_Tropisms;

  procedure Mixed_Cell_Tropisms
              ( file : in file_type;
                sup : in out Array_of_Lists;
                mcc : out Mixed_Subdivision;
                mv : out natural32 ) is

  -- DESCRIPTION :
  --   Given a system of n Laurent polynomials in n+1 variables,
  --   computes the tropisms, where the last variable is the parameter.
  --   The tropisms are computed via a regular mixed cell configuration,
  --   induced by the lifting defined by the last variable.
  --   As the lifting may differ, even if several supports would be the same,
  --   the type of mixed is assumed to be fully mixed.

  -- ON ENTRY :
  --   file     for intermediate output;
  --   sup      a list of n supports in n+1 variables.

  -- ON RETURN :
  --   mcc      a mixed cell configuration induced by the lifting
  --            defined by the last exponent in each monomial of p;
  --   mv       the mixed volume of the cells in mcc.

    dim : constant integer32 := sup'last;
    mix : Standard_Integer_Vectors.Vector(1..dim);

  begin
    new_line(file);
    put_line(file,"THE SUPPORTS : ");
    put(file,sup);
    mix := (mix'range => 1);
    Integer_Create_Mixed_Cells(file,dim,mix,false,sup,mcc);
    Integer_Volume_Computation(file,dim,mix,true,sup,mcc,mv);
  end Mixed_Cell_Tropisms;

  procedure Initial_Coefficients
              ( p : in Laur_Sys; mic : in Mixed_Cell;
                sols : out Solution_List ) is

  -- DESCRIPTION :
  --   Calls the blackbox solver on the subsystem of p
  --   supported by the points that span the mixed cell mic.
  --   The solutions are in the list sols on return.

    q : Laur_Sys(p'range) := Select_Terms(p,mic.pts.all);
    idx : constant integer32 := p'last+1;
    one : constant Complex_Number := Create(1.0);
    s : Laur_Sys(q'range) := Eval(q,one,idx);
    rc : natural32;

  begin
    put_line("The initial form system with t :"); put(q);
    put_line("The initial form system with t = 1 :"); put(s);
    Solve(s,false,rc,sols);
    put("Computed "); put(Length_Of(sols),1); put_line(" solutions.");
    put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Initial_Coefficients;

  procedure Initials
              ( p : in Laur_Sys;
                mcc : in Mixed_Subdivision;
                mv : in natural32 ) is

  -- DESCRIPTION :
  --   Solves all initial form systems defined by the cells in mcc.

    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;

  begin
    for k in 1..Length_Of(mcc) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
        sols : Solution_List;
      begin
        put("Tropism "); put(k,1); put(" is ");
        put(mic.nor); new_line;
        Initial_Coefficients(p,mic,sols);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Initials;

  procedure Pretropisms
              ( p : Standard_Complex_Laur_Systems.Laur_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user if output to a file is wanted
  --   and then launches the pretropisms computation.

    sup : Array_of_Lists(p'range) := Create(p);
    ans : character;
    report : boolean;
    file : file_type;
    cells : Mixed_Subdivision;
    mixvol : natural32;

  begin
    new_line;
    put("Do you want intermediate output to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    report := (ans = 'y');
    if report then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(file);
      Mixed_Cell_Tropisms(file,sup,cells,mixvol);
    else
      new_line;
      put("Do you want intermediate output to screen ? (y/n) ");
      Ask_Yes_or_No(ans);
      report := (ans = 'y');
      Mixed_Cell_Tropisms(report,sup,cells,mixvol);
    end if;
    put("The number of tropisms : "); put(Length_Of(cells),1); new_line;
    put("The number of series : "); put(mixvol,1); new_line;
    Initials(p,cells,mixvol);
  end Pretropisms;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system
  --   and checks whether the n polynomials have n+1 variables.

    lp : Link_to_Laur_Sys;
    nq,nv : integer32;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    get(lp);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("Number of polynomials : "); put(nq,1); new_line;
    put("Number of variables : "); put(nv,1); new_line;
    if nv /= nq+1
     then put(nv,1); put(" /= "); put(nq,1); put(" + 1");
     else Pretropisms(lp.all);
    end if;
  end Main;

begin
  Main;
end ts_puiseux;
