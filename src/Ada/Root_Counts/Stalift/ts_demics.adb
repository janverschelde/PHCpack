with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;
with Drivers_for_DEMiCs_Algorithm;       use Drivers_for_DEMiCs_Algorithm;

procedure ts_demics is

-- DESCRIPTION :
--   Test on blackbox call to the DEMiCs Algorithm.

  procedure Call_DEMiCs ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Calls DEMiCs on the system p.

    dim : constant natural32 := natural32(p'last);
    mix : Standard_Integer_Vectors.Link_to_Vector;
    lif : Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
    mcc : Mixed_Subdivision;
    mv : natural32;

  begin
    BlackBox_DEMiCs_Algorithm(p,mix,lif,mcc,mv);
    put("The mixed volume : "); put(mv,1); new_line;
    put("The type of mixture : "); put(mix); new_line;
    put_line("The lifted supports : ");
    Floating_Mixed_Subdivisions_io.put(lif.all);
    Floating_Mixed_Subdivisions_io.put(dim,mix.all,mcc,mv);
    put("The mixed volume : "); put(mv,1); new_line;
  end Call_DEMiCs;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and then calls DEMiCs.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
  
  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    new_line;
    Call_DEMiCs(lp.all);
  end Main;

begin
  Main;
end ts_demics;
