with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Natural_VecVecs;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Sampling_Machine; 
with DoblDobl_Sampling_Machine;
with QuadDobl_Sampling_Machine;
with Sample_Point_Lists;                use Sample_Point_Lists;
with DoblDobl_Sample_Lists;             use DoblDobl_Sample_Lists;
with QuadDobl_Sample_Lists;             use QuadDobl_Sample_Lists;
with Monodromy_Partitions;              use Monodromy_Partitions;
with Monodromy_Component_Breakup;       use Monodromy_Component_Breakup;

procedure ts_newfac is

-- DESCRIPTION :
--   Test facility for monodromy breakup interlaced with linear traces.

  procedure Standard_Factor is

  -- DESCRIPTION :
  --   Performs the factorization in standard double precision.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    ans : character;
    grid : Array_of_Standard_Sample_Lists(0..2);
    f : Standard_Natural_VecVecs.Link_to_VecVec;

  begin
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put("Do you want intermediate output to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the output file");
      Read_Name_and_Create_File(file);
      new_line;
      put_line("See the output file for results...");
      new_line;
      Sampling_Machine.Initialize(lp.all);
      Sampling_Machine.Default_Tune_Sampler(0);
      Sampling_Machine.Default_Tune_Refiner;
      grid := Create(file,lp.all,sols,dim);
      Factor(file,lp.all,dim,grid,f);
      Sampling_Machine.Clear;
    else
      new_line;
      Sampling_Machine.Initialize(lp.all);
      Sampling_Machine.Default_Tune_Sampler(0);
      Sampling_Machine.Default_Tune_Refiner;
      grid := Create(lp.all,sols,dim);
      Factor(lp.all,dim,grid,f);
      Sampling_Machine.Clear;
      put_line("The factorization : ");
      Write_Factors(Standard_Output,f.all);
    end if;
  end Standard_Factor;

  procedure DoblDobl_Factor is

  -- DESCRIPTION :
  --   Performs the factorization in standard double precision.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    ans : character;
    grid : Array_of_DoblDobl_Sample_Lists(0..2);
    f : Standard_Natural_VecVecs.Link_to_VecVec;

  begin
    DoblDobl_Read_Embedding(lp,sols,dim);
    new_line;
    put("Do you want intermediate output to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the output file");
      Read_Name_and_Create_File(file);
      new_line;
      put_line("See the output file for results...");
      new_line;
      DoblDobl_Sampling_Machine.Initialize(lp.all);
      DoblDobl_Sampling_Machine.Default_Tune_Sampler(0);
      DoblDobl_Sampling_Machine.Default_Tune_Refiner;
      grid := Create(file,lp.all,sols,dim);
      Factor(file,lp.all,dim,grid,f);
      DoblDobl_Sampling_Machine.Clear;
    else
      new_line;
      DoblDobl_Sampling_Machine.Initialize(lp.all);
      DoblDobl_Sampling_Machine.Default_Tune_Sampler(0);
      DoblDobl_Sampling_Machine.Default_Tune_Refiner;
      grid := Create(lp.all,sols,dim);
      Factor(lp.all,dim,grid,f);
      DoblDobl_Sampling_Machine.Clear;
      put_line("The factorization : ");
      Write_Factors(Standard_Output,f.all);
    end if;
  end DoblDobl_Factor;

  procedure QuadDobl_Factor is

  -- DESCRIPTION :
  --   Performs the factorization in standard double precision.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    ans : character;
    grid : Array_of_QuadDobl_Sample_Lists(0..2);
    f : Standard_Natural_VecVecs.Link_to_VecVec;

  begin
    QuadDobl_Read_Embedding(lp,sols,dim);
    new_line;
    put("Do you want intermediate output to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the output file");
      Read_Name_and_Create_File(file);
      new_line;
      put_line("See the output file for results...");
      new_line;
      QuadDobl_Sampling_Machine.Initialize(lp.all);
      QuadDobl_Sampling_Machine.Default_Tune_Sampler(0);
      QuadDobl_Sampling_Machine.Default_Tune_Refiner;
      grid := Create(file,lp.all,sols,dim);
      Factor(file,lp.all,dim,grid,f);
      QuadDobl_Sampling_Machine.Clear;
    else
      new_line;
      QuadDobl_Sampling_Machine.Initialize(lp.all);
      QuadDobl_Sampling_Machine.Default_Tune_Sampler(0);
      QuadDobl_Sampling_Machine.Default_Tune_Refiner;
      grid := Create(lp.all,sols,dim);
      Factor(lp.all,dim,grid,f);
      QuadDobl_Sampling_Machine.Clear;
      put_line("The factorization : ");
      Write_Factors(Standard_Output,f.all);
    end if;
  end QuadDobl_Factor;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision and then calls
  --   the corresponding driver routine.

    ans : character;

  begin
    new_line;
    put_line("Factorization with monodromy and linear traces.");
    new_line;
    put_line("MENU to select the precision level : ");
    put_line("  0. standard double precision; or");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Factor;
      when '1' => DoblDobl_Factor;
      when '2' => QuadDobl_Factor;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_newfac;
