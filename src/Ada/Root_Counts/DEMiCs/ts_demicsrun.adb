with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Floating_Lifting_Functions;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;
with DEMiCs_Algorithm;                   use DEMiCs_Algorithm;
with DEMiCs_Output_Data;
with Drivers_for_Static_Lifting;
with use_c2phc; -- to force the compilation of use_c2phc.adb ...

procedure ts_demicsrun is

-- DESCRIPTION :
--   Development of a run with DEMiCs, invoked by an Ada main program.

  procedure Process_Output 
              ( dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stable : in boolean; stlb : in double_float;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Processes the output of DEMiCs.

    lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    mcc,orgmcc,stbmcc : Mixed_Subdivision;
    mv,smv,tmv,orgcnt,stbcnt : natural32;

  begin
    Process_Output(dim,mix,sup,lifsup,mcc,verbose);
    put_line("The lifted supports :");
    Floating_Mixed_Subdivisions_io.put(lifsup);
    if not stable then
      put_line("The mixed-cell configuration :");
      Floating_Mixed_Subdivisions_io.put(natural32(dim),mix.all,mcc,mv);
      put("The mixed volume : "); put(mv,1); new_line;
    else
      Split_Original_Cells(mcc,stlb,orgmcc,stbmcc,orgcnt,stbcnt);
      new_line;
      put("Number of cells without artificial origin : ");
      put(orgcnt,1); new_line;
      put("#extra stable cells with artificial origin : ");
      put(stbcnt,1); new_line;
      Drivers_for_Static_Lifting.Floating_Volume_Computation
        (standard_output,dim,stlb,mix.all,mcc,mv,smv,tmv);
    end if;
  end Process_Output;

  procedure Compute_Mixed_Volume ( p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Computes the mixed volume of p.

    dim : constant integer32 := p'last;
    ans : character;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    verbose,stable : boolean;
    stlb : double_float;

  begin
    new_line;
    put("Verbose ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put("Monitor the adding of cell indices ? (y/n) ");
    Ask_Yes_or_No(ans);
    DEMiCs_Output_Data.monitor := (ans = 'y');
    new_line;
    put("Do you want the stable mixed volume ? (y/n) ");
    Ask_Yes_or_No(ans);
    stable := (ans = 'y');
    DEMiCs_Output_Data.stable := stable;
    if ans = 'y'
     then stlb := Floating_Lifting_Functions.Lifting_Bound(p);
     else stlb := 0.0;
    end if;
    Extract_Supports(p,mix,sup,verbose);
    Call_DEMiCs(mix,sup,stable,stlb,verbose);
    Show_Output;
    Process_Output(dim,mix,sup,stable,stlb,verbose);
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
