with text_io,integer_io;                use text_io,integer_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Witness_Sets_io;                   use Witness_Sets_io;
with Multitasking_Sampling;

procedure ts_mtsam is

-- DESCRIPTION :
--   Interactive development of a sampler using multiple tasks.

  procedure Main is

    n,d : natural;
    p : Link_to_Poly_Sys;
    sols : Solution_List;
    file : file_type;

  begin
    new_line;
    put_line("Multitasked sampling a positive dimensional solution set ...");
    Standard_Read_Embedding(p,sols,d);
    new_line;
    put("Read witness set of dimension "); put(d,1);
    put(" and degree "); put(Length_Of(sols),1); put_line(".");
    put("give number of tasks : "); get(n); skip_line;
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
   -- Multitasking_Sampling.Driver_to_Sampler(file,n,d,p.all,sols);
    Multitasking_Sampling.Driver_to_Monodromy(file,n,d,p.all,sols);
  end Main;

begin
  Main;
end ts_mtsam;
