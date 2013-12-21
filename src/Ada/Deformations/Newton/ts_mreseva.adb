with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with Multprec_Complex_Vectors;           use Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Multprec_Complex_Poly_Systems;      use Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_SysFun;       use Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Solutions;         use Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;
with Multprec_Residual_Evaluations;      use Multprec_Residual_Evaluations;

procedure ts_mreseva is

-- DESCRIPTION : test on multi-precision residual computation.

  procedure Test_Solution_Residuals
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys ) is

    outfile,solsfile : file_type;
    p_eval : Multprec_Complex_Poly_SysFun.Eval_Poly_Sys(p'range) := Create(p);
    stsols : Standard_Complex_Solutions.Solution_List;
    mpsols : Multprec_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    new_line;
    put_line("Reading the name of the file with the solutions.");
    Read_Name_and_Open_File(solsfile);
    get(solsfile,stsols);
    mpsols := Create(stsols);
   -- get(solsfile,mpsols);
    put_line(outfile,"The list of solutions : ");
    put(outfile,mpsols);
    Residuals(outfile,p_eval,mpsols);
  end Test_Solution_Residuals;

  procedure Interactive_Test_Residuals
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys ) is

    p_eval : Multprec_Complex_Poly_SysFun.Eval_Poly_Sys(p'range) := Create(p);
    root,eva : Vector(p'range);
    res : Floating_Number;
    ans : character;

  begin
    loop
      new_line;
      put("Give "); put(root'last,1);
      put_line(" complex numbers for the root : ");
      get(root);
      eva := Eval(p_eval,root);
      put_line("The evaluated root : "); put_line(eva);
      res := Residual(p_eval,root);
      put("The residual : "); put(res); new_line;
      put("Do you want more tests ? (y/n) "); Ask_Yes_or_No(ans);
      Clear(root); Clear(eva); Clear(res);
      exit when (ans /= 'y');
    end loop;
  end Interactive_Test_Residuals;

  procedure Main is
 
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    ans : character;

  begin
    put_line("Choose one of the following : ");
    put_line("  1. Evaluate user-given vectors");
    put_line("  2. Residuals for a solution list.");
    put("Type 1 or 2 to select : "); Ask_Alternative(ans,"12");
    new_line;
    get(lp);
    declare
      mp : Multprec_Complex_Poly_Systems.Poly_Sys(lp'range) := Convert(lp.all);
    begin
      if ans = '1'
       then Interactive_Test_Residuals(mp);
       else Test_Solution_Residuals(mp);
      end if;
    end;
  end Main;

begin
  new_line;
  put_line("Test on the multi-precision residual computation.");
  new_line;
  Main;
end ts_mreseva;
