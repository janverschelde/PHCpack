with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Natural_VecVecs;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Sampling_Machine; 
with Sample_Point_Lists;                use Sample_Point_Lists;
with Monodromy_Partitions;              use Monodromy_Partitions;
with Monodromy_Component_Breakup;       use Monodromy_Component_Breakup;

procedure ts_newfac is

-- DESCRIPTION :
--   Test facility for monodromy breakup interlaced with linear traces.

  procedure Main is

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    ans : character;
    grid : Array_of_Standard_Sample_Lists(0..2);
    f : Standard_Natural_VecVecs.Link_to_VecVec;

  begin
    new_line;
    put_line("Factorization with monodromy and linear traces.");
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
  end Main;

begin
  Main;
end ts_newfac;
