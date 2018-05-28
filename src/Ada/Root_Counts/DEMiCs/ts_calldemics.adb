with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;       use Standard_Floating_VecVecs_io;
with Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;
with DEMiCs_Command_Line;                use DEMiCs_Command_Line;
with DEMiCs_Output_Convertors;           use DEMiCs_Output_Convertors;

procedure ts_calldemics is

-- DESCRIPTION :
--   Tests the basic file based command line interface to DEMiCs.

  procedure Call_DEMiCs 
              ( inpfile,outfile : in string; verbose : in boolean;
                p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Given names of the input and output files, the verbose flag,
  --   the mixed volume of the system p will be computed.

  -- ON ENTRY :
  --   inpfile  name of the input file;
  --   outfile  name of the output file;
  --   verbose  flag to indicate if extra output is wanted;
  --   p        a polynomial system to compute the mixed volume for.

    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    mix,perm : Standard_Integer_Vectors.Link_to_Vector;
    cells : Lists_of_Integer_Vectors.List;
    mcc : Mixed_Subdivision;
    mv : natural32;

  begin
    Prepare_Input(p,inpfile,mix,perm,sup,verbose);
    declare
      lif : Standard_Floating_VecVecs.VecVec(mix'range);
      lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    begin -- assume demics is in the execution path
      Call_DEMiCs(inpfile,outfile,execfile=>"demics",verbose=>verbose);
      Process_Output(p'last,outfile,mix.all,mv,lif,cells,verbose);
      new_line;
      put_line("The lifting values for the supports : "); put(lif);
      new_line;
      put_line("The indices to the mixed cells :"); put(cells);
      new_line;
      put("The mixed volume : "); put(mv,1); new_line;
      if mix'last = p'last
       then lifsup := Apply_Lifting(sup,lif);         -- fully mixed
       else lifsup := Apply_Lifting(mix.all,sup,lif); -- semi- or un-mixed
      end if;
      new_line;
      put_line("The lifting applied to the supports :");
      Floating_Mixed_Subdivisions_io.put(lifsup);
      if mix'last = p'last
       then mcc := Make_Mixed_Cells(sup'last,cells,lifsup);
       else mcc := Make_Mixed_Cells(sup'last,mix.all,cells,lifsup);
      end if;
      new_line;
      put_line("The mixed-cell configuration :");
      Floating_Mixed_Subdivisions_io.put(natural32(sup'last),mix.all,mcc,mv);
      put("The mixed volume : "); put(mv,1); new_line;
    end;
  end Call_DEMiCs;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   prepares the input for demics, calls demics,
  --   and then parses the output file.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    randnbr : constant string := Random_Name("");
    iptname : constant string := "/tmp/demics_input" & randnbr;
    optname : constant string := "/tmp/demics_output" & randnbr;
    ans : character;
    verbose : boolean;

  begin
    new_line;
    put_line("Calling DEMiCs for mixed volume computation ...");
    new_line;
    get(lp);
    new_line;
    put_line("Writing input to " & iptname & ".");
    put_line("Writing output to " & optname & ".");
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    Call_DEMiCs(iptname,optname,verbose,lp.all);
  end Main;

begin
  Main;
end ts_calldemics;
