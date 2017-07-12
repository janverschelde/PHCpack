with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Sample_Point_Lists;                use Sample_Point_Lists;
with DoblDobl_Sample_Lists;             use DoblDobl_Sample_Lists;
with QuadDobl_Sample_Lists;             use QuadDobl_Sample_Lists;
with Combinatorial_Factorization;       use Combinatorial_Factorization;
with Drivers_to_Breakup_Solutions;      use Drivers_to_Breakup_Solutions;

procedure ts_combfac is

-- DESCRIPTION :
--   Test facility to develop a combinatorial factorization
--   using linear traces.

-- AUXILIARIES FOR COMBINATORIAL TESTS :

  procedure Write_Factors
              ( factors : in Standard_Natural_VecVecs.VecVec;
                count,depth : in natural32 ) is
  begin
    put(count,1); put(" : ");
    Write(Standard_Output,factors);
    put(" at depth "); put(depth,1); new_line;
  end Write_Factors;

  procedure Enum is new Enumerate_Factorizations(Write_Factors);

  function Ask_Factor
             ( f : Standard_Natural_Vectors.Vector ) return boolean is

    ans : character;

  begin
    put("Accept"); put(f); put(" as a factor ? (y/n) ");
    Ask_Yes_or_No(ans);
    return (ans = 'y');
  end Ask_Factor;

  function Seek is new Search_Factorization_with_Output(Ask_Factor);

-- AUXILIARIES FOR LINEAR TRACE CERTIFICATES :

  procedure Standard_Factor is

  -- DESCRIPTION :
  --   Performs a combinatorial factorization on sample point grids
  --   computed in standard double precision, for a witness set defined
  --   by an ordinary polynomial system.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim,n : natural32;
    grid : Array_of_Standard_Sample_Lists(0..2);

  begin
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results...");
    grid := Create(file,lp.all,sols,dim);
    n := Length_Of(sols);
    declare
      factors : constant Standard_Natural_VecVecs.VecVec
              := Factor(file,n,grid);
    begin
      null;
    end;
  end Standard_Factor;

  procedure Standard_Laurent_Factor is

  -- DESCRIPTION :
  --   Performs a combinatorial factorization on sample point grids
  --   computed in standard double precision, for a witness set defined
  --   by a Laurent polynomial system.

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    file : file_type;
    lp : Link_to_Laur_Sys;
    sols : Solution_List;
    dim,n : natural32;
    grid : Array_of_Standard_Sample_Lists(0..2);

  begin
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results...");
    grid := Create(file,lp.all,sols,dim);
    n := Length_Of(sols);
    declare
      factors : constant Standard_Natural_VecVecs.VecVec
              := Factor(file,n,grid);
    begin
      null;
    end;
  end Standard_Laurent_Factor;

  procedure DoblDobl_Factor is

  -- DESCRIPTION :
  --   Performs a combinatorial factorization on sample point grids
  --   computed in double double precision, for a witness set defined
  --   by an ordinary polynomial system.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim,n : natural32;
    grid : Array_of_DoblDobl_Sample_Lists(0..2);

  begin
    DoblDobl_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results...");
    grid := Create(file,lp.all,sols,dim);
    n := Length_Of(sols);
    declare
      factors : constant Standard_Natural_VecVecs.VecVec
              := Factor(file,n,grid);
    begin
      null;
    end;
  end DoblDobl_Factor;

  procedure DoblDobl_Laurent_Factor is

  -- DESCRIPTION :
  --   Performs a combinatorial factorization on sample point grids
  --   computed in double double precision, for a witness set defined
  --   by a Laurent polynomial system.

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    file : file_type;
    lp : Link_to_Laur_Sys;
    sols : Solution_List;
    dim,n : natural32;
    grid : Array_of_DoblDobl_Sample_Lists(0..2);

  begin
    DoblDobl_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results...");
    grid := Create(file,lp.all,sols,dim);
    n := Length_Of(sols);
    declare
      factors : constant Standard_Natural_VecVecs.VecVec
              := Factor(file,n,grid);
    begin
      null;
    end;
  end DoblDobl_Laurent_Factor;

  procedure QuadDobl_Factor is

  -- DESCRIPTION :
  --   Performs a combinatorial factorization on sample point grids
  --   computed in quad double precision, for a witness set defined
  --   by an ordinary polynomial system.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim,n : natural32;
    grid : Array_of_QuadDobl_Sample_Lists(0..2);

  begin
    QuadDobl_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results...");
    grid := Create(file,lp.all,sols,dim);
    n := Length_Of(sols);
    declare
      factors : constant Standard_Natural_VecVecs.VecVec
              := Factor(file,n,grid);
    begin
      null;
    end;
  end QuadDobl_Factor;

  procedure QuadDobl_Laurent_Factor is

  -- DESCRIPTION :
  --   Performs a combinatorial factorization on sample point grids
  --   computed in quad double precision, for a witness set defined
  --   by a Laurent polynomial system.

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    file : file_type;
    lp : Link_to_Laur_Sys;
    sols : Solution_List;
    dim,n : natural32;
    grid : Array_of_QuadDobl_Sample_Lists(0..2);

  begin
    QuadDobl_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results...");
    grid := Create(file,lp.all,sols,dim);
    n := Length_Of(sols);
    declare
      factors : constant Standard_Natural_VecVecs.VecVec
              := Factor(file,n,grid);
    begin
      null;
    end;
  end QuadDobl_Laurent_Factor;

  procedure Call_Factor is

  -- DESCRIPTION :
  --   Prompts the user for a precision level,
  --   calls the corresponding sampler and combinatorial factorization.

    prc,ans : character;

  begin
    new_line;
    put_line("MENU to select the precision of the samples :");
    put_line("  0. default standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(prc,"012");
    new_line;
    put("Laurent polynomial system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      case prc is
        when '0' => Standard_Laurent_Factor;
        when '1' => DoblDobl_Laurent_Factor;
        when '2' => QuadDobl_Laurent_Factor;
        when others => null;
      end case;
    else
      case prc is
        when '0' => Standard_Factor;
        when '1' => DoblDobl_Factor;
        when '2' => QuadDobl_Factor;
        when others => null;
      end case;
    end if;
  end Call_Factor;

  procedure Main is

  -- DESCRIPTION :
  --   Shows the main testing loop.

    n : natural32 := 0;
    ans : character;

  begin
    loop
      new_line;
      put_line("MENU to test combinatorial enumerations :");
      put_line("  0. exit this program;");
      put_line("  1. enumerate all possible factors;");
      put_line("  2. enumerate all possible factorizations;");
      put_line("  3. search for one factorization;");
      put_line("  4. perform combinatorial factorization.");
      put("Type 0, 1, 2, 3, or 4 to select : ");
      Ask_Alternative(ans,"01234");
      if ans /= '4' then
        new_line;
        exit when ans = '0';
        put("Give the number of witness points : "); get(n);
      end if;
      case ans is
        when '1' => Enumerate_Factors(n);
        when '2' => Enum(n);
        when '3' =>
          declare
            factors : constant Standard_Natural_VecVecs.VecVec
                    := Seek(Standard_Output,n);
          begin
            put("Factorization found : ");
            Write(Standard_Output,factors); new_line;
          end;
        when '4' => Call_Factor;
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_combfac;
