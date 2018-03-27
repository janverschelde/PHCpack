with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Witness_Sets_io;
with Random_Test_Points;                 use Random_Test_Points;

procedure ts_mtmbthom is

-- DESCRIPTION :
--   Development of a multithreaded version of the homotopy membership test.

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a witness set, 
  --   prompts for the number of tasks,
  --   generates a random point in standard double precision
  --   and then launches the multitasked homotopy membership test.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols : Standard_Complex_Solutions.Solution_List;
    dim,deg,nt : natural32 := 0;

  begin
    Witness_Sets_io.Standard_Read_Embedding(lp,genpts,dim);
    deg := Standard_Complex_Solutions.Length_Of(genpts);
    new_line;
    put("Read a set of degree "); put(deg,1);
    put(" and dimension "); put(dim,1); put_line(".");
    new_line;
    put("Give the number of tasks : "); get(nt);
    new_line;
    put_line("Computing a random point on the solution set ...");
    sols := Standard_Random_Point(standard_output,lp.all,genpts,dim);
  end Standard_Test;

begin
  Standard_Test;
end ts_mtmbthom;

