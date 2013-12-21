with text_io ;                           use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Lists_of_Integer64_Vectors;
with Lists_of_Integer64_Vectors_io;      use Lists_of_Integer64_Vectors_io;
with Standard_Integer32_Triangulations;
with Standard_Integer32_Triangulations_io;
 use Standard_Integer32_Triangulations_io;
with Standard_Integer64_Triangulations;
with Standard_Integer64_Triangulations_io;
 use Standard_Integer64_Triangulations_io;
with Standard_Dynamic32_Triangulations;  use Standard_Dynamic32_Triangulations;
with Standard_Dynamic64_Triangulations;  use Standard_Dynamic64_Triangulations;

procedure ts_dyntri is

-- DESCRIPTION :
--   Test the creation of a regular triangulation by dynamic lifting.

  procedure Standard_Integer32_Computation ( file : in file_type ) is

  -- DESCRIPTION :
  --   Uses standard 32-bit arithmetic to compute a regular triangulation.

    n,d,vol : natural32 := 0;
    L,lifted,lifted_last : Lists_of_Integer_Vectors.List;
    t : Standard_Integer32_Triangulations.Triangulation;

  begin
    new_line;
    put("Give d, the dimension : "); get(d);
    put("Give n, the number of points : "); get(n);
    put("Give "); put(n,1); put(" vectors of length "); put(d,1);
    put_line(" :"); get(d,n,L);
    put(file,"n = "); put(file,n,1);
    put(file,"  d = "); put(file,d,1); new_line(file);
    put_line(file,"The points : "); put(file,L);
    Dynamic_Lifting(L,true,false,0,lifted,lifted_last,t);
    put_line(file,"The lifted points : "); put(file,lifted);
    put_line(file,"The triangulation : "); put(file,d,t,vol);
    put(file,"The volume : "); put(file,vol,1); new_line(file);
  end Standard_Integer32_Computation;

  procedure Standard_Integer64_Computation ( file : in file_type ) is

  -- DESCRIPTION :
  --   Uses standard 64-bit arithmetic to compute a regular triangulation.

    n,d : natural32 := 0;
    L,lifted,lifted_last : Lists_of_Integer64_Vectors.List;
    t : Standard_Integer64_Triangulations.Triangulation;
    vol : natural64;

  begin
    new_line;
    put("Give d, the dimension : "); get(d);
    put("Give n, the number of points : "); get(n);
    put("Give "); put(n,1); put(" vectors of length "); put(d,1);
    put_line(" :"); get(d,n,L);
    put(file,"n = "); put(file,n,1);
    put(file,"  d = "); put(file,d,1); new_line(file);
    put_line(file,"The points : "); put(file,L);
    Dynamic_Lifting(L,true,false,0,lifted,lifted_last,t);
    put_line(file,"The lifted points : "); put(file,lifted);
    put_line(file,"The triangulation : "); put(file,d,t,vol);
    put(file,"The volume : ");
    Standard_Natural_Numbers_io.put(file,vol,1); new_line(file);
  end Standard_Integer64_Computation;

  procedure Main is

    file : file_type;
    ans : character;

  begin
    new_line;
    put_line("Testing the computation of a regular triangulation");
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put("Use 64-bit arithmetic ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Standard_Integer64_Computation(file);
     else Standard_Integer32_Computation(file);
    end if;
  end Main;

begin
  Main;
end ts_dyntri;
