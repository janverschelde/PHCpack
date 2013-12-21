with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Multprec_Floating_Vectors;
with Multprec_Complex_Solutions_io;
with Continuation_Parameters;
with Standard_Continuation_Data_io;
with Multprec_Continuation_Data;         use Multprec_Continuation_Data;
with Multprec_Path_Trackers;             use Multprec_Path_Trackers;

package body Multprec_IncFix_Continuation is

  tol : constant double_float := 10.0E-14;

  procedure Silent_Continue
               ( sols : in out Solution_List; proj : in boolean;
                 target : in Complex_Number := Create(integer(1)) ) is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    patpp : constant Pred_Pars
          := Convert(Continuation_Parameters.Create_for_Path);
    endpp : constant Pred_Pars
          := Convert(Continuation_Parameters.Create_End_Game);
    patcp : constant Corr_Pars
          := Convert(Continuation_Parameters.Create_for_Path);
    endcp : constant Corr_Pars
          := Convert(Continuation_Parameters.Create_End_Game);
    pathdir : Multprec_Floating_Vectors.Link_to_Vector;
    errv : Floating_Number;
    order : constant integer32
          := integer32(Continuation_Parameters.endext_order);

    procedure Continue_along_Path is
      new Linear_Single_Normal_Silent_Continue(Norm,H,dH,dH);

    procedure Continue_in_End_Game is
      new Linear_Single_Conditioned_Silent_Continue(Norm,H,dH,dH);

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
        s : Solu_Info := Deep_Create(ls.all);
      begin
        Continue_along_Path(s,target,Create(tol),false,patpp,patcp);
        Continue_in_End_Game(s,target,Create(tol),false,
                             order,pathdir,errv,endpp,endcp);
        Append(res,res_last,Deep_Create(s));
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Copy(res,sols); Clear(res);
  end Silent_Continue;

  procedure Write_Solution ( file : in file_type; s : in Solu_Info ) is
  begin
    put(file,"t : "); put(file,s.sol.t); new_line(file);
    put(file,"m : "); put(file,s.sol.m,1);
    put(file,"                  Length of path : ");
    put(file,s.length_path);
    new_line(file);
    put_line(file,"the solution for t : ");
    Multprec_Complex_Solutions_io.put_vector(file,s.sol.all);
    put(file,"==");
    put(file," err : "); put(file,s.cora,2,3,3);  put(file," =");
    put(file," rco : "); put(file,s.rcond,2,3,3); put(file," =");
    put(file," res : "); put(file,s.resa,2,3,3);  put_line(file," ==");
  end Write_Solution;

  procedure Reporting_Continue
               ( file : in file_type; sols : in out Solution_List;
                 proj : in boolean;
                 target : in Complex_Number := Create(integer(1)) ) is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    patpp : constant Pred_Pars
          := Convert(Continuation_Parameters.Create_for_Path);
    endpp : constant Pred_Pars
          := Convert(Continuation_Parameters.Create_End_Game);
    patcp : constant Corr_Pars
          := Convert(Continuation_Parameters.Create_for_Path);
    endcp : constant Corr_Pars
          := Convert(Continuation_Parameters.Create_End_Game);
    pathdir : Multprec_Floating_Vectors.Link_to_Vector;
    errv : Floating_Number;
    order : constant integer32
          := integer32(Continuation_Parameters.endext_order);
    len : constant natural32 := Length_Of(sols);

    procedure Continue_along_Path is
      new Linear_Single_Normal_Reporting_Continue(Norm,H,dH,dH);

    procedure Continue_in_End_Game is
      new Linear_Single_Conditioned_Reporting_Continue(Norm,H,dH,dH);

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,len,1); put(file," ");
    put(file,Head_Of(sols).n,1); new_line(file);
    for i in 1..76 loop put(file,"="); end loop; new_line(file);
    for i in 1..len loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
        s : Solu_Info := Deep_Create(ls.all);
      begin
        Continue_along_Path(file,s,target,Create(tol),false,patpp,patcp);
        Continue_in_End_Game(file,s,target,Create(tol),false,
                             order,pathdir,errv,endpp,endcp);
        Standard_Continuation_Data_io.Write_Statistics
          (file,integer32(i),s.nstep,s.nfail,s.niter,s.nsyst);
        new_line(file);
        Write_Solution(file,s); flush(file);
        Append(res,res_last,Deep_Create(s));
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Copy(res,sols); Clear(res);
  end Reporting_Continue;

end Multprec_IncFix_Continuation;
