with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;
with DEMiCs_Algorithm;                   use DEMiCs_Algorithm;

procedure ts_demicsrun is

-- DESCRIPTION :
--   Development of a run with DEMiCs, invoked by an Ada main program.

  procedure Compute_Mixed_Volume ( p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Computes the mixed volume of p.

    dim : constant integer32 := p'last;
    ans : character;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    verbose : boolean;
    mv : natural32;
    mcc : Mixed_Subdivision;

  begin
    new_line;
    put("Verbose ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    Extract_Supports(p,mix,sup,verbose);
    Call_DEMiCs(mix,sup,verbose);
    Show_Output;
    declare
      lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    begin
      Process_Output(dim,mix,sup,lifsup,mcc,verbose);
      put_line("The lifted supports :");
      Floating_Mixed_Subdivisions_io.put(lifsup);
    end;
    put_line("The mixed-cell configuration :");
    Floating_Mixed_Subdivisions_io.put(natural32(dim),mix.all,mcc,mv);
    put("The mixed volume : "); put(mv,1); new_line;
  end Compute_Mixed_Volume;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and then prepares the input for DEMiCs.

    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    Compute_Mixed_Volume(lp.all);
  end Main;

begin
  Main;
end ts_demicsrun;
