with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Vectors;
with Continuation_Parameters;            use Continuation_Parameters;
with DoblDobl_Continuation_Data;         use DoblDobl_Continuation_Data;
with DoblDobl_Continuation_Data_io;      use DoblDobl_Continuation_Data_io;
with DoblDobl_Path_Trackers;             use DoblDobl_Path_Trackers;

package body DoblDobl_IncFix_Continuation is

  tol : constant double_float := 10.0E-14;

  procedure Silent_Continue
               ( sols : in out Solution_List;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(integer(1)) ) is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    patpp : constant Pred_Pars := Continuation_Parameters.Create_for_Path;
    endpp : constant Pred_Pars := Continuation_Parameters.Create_End_Game;
    patcp : constant Corr_Pars := Continuation_Parameters.Create_for_Path;
    endcp : constant Corr_Pars := Continuation_Parameters.Create_End_Game;
    w : integer32 := 1;
    pathdir : Double_Double_Vectors.Link_to_Vector;
    errv : double_double;
    order : constant integer32
          := integer32(Continuation_Parameters.endext_order);

    procedure Continue_along_Path is
      new Linear_Single_Normal_Silent_Continue(Norm,H,dH,dH);

    procedure Continue_End_Game is
      new Linear_Single_Conditioned_Silent_Continue(Norm,H,dH,dH);

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
        s : Solu_Info := Deep_Create(ls.all);
      begin
        Continue_along_Path(s,target,tol,false,patpp,patcp,nbq);
        Continue_End_Game
          (s,target,tol,false,order,w,pathdir,errv,endpp,endcp,nbq);
        Append(res,res_last,Deep_Create(s));
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Copy(res,sols); Clear(res);
  end Silent_Continue;

  procedure Silent_Continue_with_Stop
               ( sols : in out Solution_List;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(integer(1)) ) is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    patpp : constant Pred_Pars := Continuation_Parameters.Create_for_Path;
    endpp : constant Pred_Pars := Continuation_Parameters.Create_End_Game;
    patcp : constant Corr_Pars := Continuation_Parameters.Create_for_Path;
    endcp : constant Corr_Pars := Continuation_Parameters.Create_End_Game;
    w : integer32 := 1;
    pathdir : Double_Double_Vectors.Link_to_Vector;
    errv : double_double;
    order : constant integer32
          := integer32(Continuation_Parameters.endext_order);

    procedure Continue_along_Path is
      new Linear_Single_Normal_Silent_Continue(Norm,H,dH,dH);

    procedure Continue_End_Game is
      new Linear_Single_Conditioned_Silent_Continue(Norm,H,dH,dH);

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
        s : Solu_Info := Deep_Create(ls.all);
      begin
        Continue_along_Path(s,target,tol,false,patpp,patcp,nbq);
        Continue_End_Game
          (s,target,tol,false,order,w,pathdir,errv,endpp,endcp,nbq);
        Append(res,res_last,Deep_Create(s));
        exit when Stop_Test(s.sol.all);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Copy(res,sols); Clear(res);
  end Silent_Continue_with_Stop;

  procedure Reporting_Continue
               ( file : in file_type; sols : in out Solution_List;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(integer(1)) ) is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    patpp : constant Pred_Pars := Continuation_Parameters.Create_for_Path;
    endpp : constant Pred_Pars := Continuation_Parameters.Create_End_Game;
    patcp : constant Corr_Pars := Continuation_Parameters.Create_for_Path;
    endcp : constant Corr_Pars := Continuation_Parameters.Create_End_Game;
    w : integer32 := 1;
    pathdir : Double_Double_Vectors.Link_to_Vector;
    tol_zero : constant double_float := endcp.epsaf;
    errv : double_double;
    order : constant integer32
          := integer32(Continuation_Parameters.endext_order);
    len : constant natural32 := Length_Of(sols);
    cnt,nbfail,nbregu,nbsing,kind : natural32 := 0;

    procedure Continue_along_Path is
      new Linear_Single_Normal_Reporting_Continue(Norm,H,dH,dH);

    procedure Continue_End_Game is
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
        Continue_along_Path(file,s,target,tol,false,patpp,patcp,nbq);
        Continue_End_Game(file,s,target,tol,false,
                          order,w,pathdir,errv,endpp,endcp,nbq);
        Write_Next_Solution
          (file,cnt,s,tol_zero,tol_zero,nbfail,nbregu,nbsing,kind);
        text_io.flush(file);
        Append(res,res_last,Deep_Create(s));
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Copy(res,sols); Clear(res);
    put(file,"Number of path failures      : ");
    put(file,nbfail,1); new_line(file);
    put(file,"Number of regular solutions  : ");
    put(file,nbregu,1); new_line(file);
    put(file,"Number of singular solutions : ");
    put(file,nbsing,1); new_line(file);
    new_line(file);
  end Reporting_Continue;

  procedure Reporting_Continue_with_Stop
               ( file : in file_type; sols : in out Solution_List;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(integer(1)) ) is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    patpp : constant Pred_Pars := Continuation_Parameters.Create_for_Path;
    endpp : constant Pred_Pars := Continuation_Parameters.Create_End_Game;
    patcp : constant Corr_Pars := Continuation_Parameters.Create_for_Path;
    endcp : constant Corr_Pars := Continuation_Parameters.Create_End_Game;
    w : integer32 := 1;
    pathdir : Double_Double_Vectors.Link_to_Vector;
    tol_zero : constant double_float := endcp.epsaf;
    errv : double_double;
    order : constant integer32
          := integer32(Continuation_Parameters.endext_order);
    len : constant natural32 := Length_Of(sols);
    cnt,nbfail,nbregu,nbsing,kind : natural32 := 0;

    procedure Continue_along_Path is
      new Linear_Single_Normal_Reporting_Continue(Norm,H,dH,dH);

    procedure Continue_End_Game is
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
        Continue_along_Path(file,s,target,tol,false,patpp,patcp);
        Continue_End_Game(file,s,target,tol,false,
                          order,w,pathdir,errv,endpp,endcp);
        Write_Next_Solution
          (file,cnt,s,tol_zero,tol_zero,nbfail,nbregu,nbsing,kind);
        text_io.flush(file);
        Append(res,res_last,Deep_Create(s));
        exit when Stop_Test(s.sol.all);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Copy(res,sols); Clear(res);
    put(file,"Number of path failures      : ");
    put(file,nbfail,1); new_line(file);
    put(file,"Number of regular solutions  : ");
    put(file,nbregu,1); new_line(file);
    put(file,"Number of singular solutions : ");
    put(file,nbsing,1); new_line(file);
    new_line(file);
  end Reporting_Continue_with_Stop;

  procedure Silent_Toric_Continue
               ( sols : in out Solution_List; proj : in boolean;
                 w : in out Standard_Integer_Vectors.Vector;
                 v : in out Double_Double_VecVecs.VecVec;
                 errv : in out Double_Double_Vectors.Vector;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(integer(1)) ) is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    patpp : constant Pred_Pars := Continuation_Parameters.Create_for_Path;
    endpp : constant Pred_Pars := Continuation_Parameters.Create_End_Game;
    patcp : constant Corr_Pars := Continuation_Parameters.Create_for_Path;
    endcp : constant Corr_Pars := Continuation_Parameters.Create_End_Game;
    order : constant integer32
          := integer32(Continuation_Parameters.endext_order);

    procedure Continue_along_Path is
      new Linear_Single_Normal_Silent_Continue(Norm,H,dH,dH);

    procedure Continue_End_Game is
      new Linear_Single_Conditioned_Silent_Continue(Norm,H,dH,dH);

  begin
    for i in w'range loop
      exit when Is_Null(tmp);
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
        s : Solu_Info := Deep_Create(ls.all);
      begin
        Continue_along_Path(s,target,tol,false,patpp,patcp,nbq);
        Continue_End_Game(s,target,tol,false,order,
                          w(i),v(i),errv(i),endpp,endcp,nbq);
        Append(res,res_last,Deep_Create(s));
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Copy(res,sols); Clear(res);
  end Silent_Toric_Continue;

  procedure Reporting_Toric_Continue
               ( file : in file_type; sols : in out Solution_List;
                 proj : in boolean;
                 w : in out Standard_Integer_Vectors.Vector;
                 v : in out Double_Double_VecVecs.VecVec;
                 errv : in out Double_Double_Vectors.Vector;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(integer(1)) ) is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    patpp : constant Pred_Pars := Continuation_Parameters.Create_for_Path;
    endpp : constant Pred_Pars := Continuation_Parameters.Create_End_Game;
    patcp : constant Corr_Pars := Continuation_Parameters.Create_for_Path;
    endcp : constant Corr_Pars := Continuation_Parameters.Create_End_Game;
    tol_zero : constant double_float := endcp.epsaf;
    order : constant integer32
          := integer32(Continuation_Parameters.endext_order);
    len : constant natural32 := Length_Of(sols);
    cnt,nbfail,nbregu,nbsing,kind : natural32 := 0;

    procedure Continue_along_Path is
      new Linear_Single_Normal_Reporting_Continue(Norm,H,dH,dH);

    procedure Continue_End_Game is
      new Linear_Single_Conditioned_Reporting_Continue(Norm,H,dH,dH);

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,len,1); put(file," ");
    put(file,Head_Of(sols).n,1); new_line(file);
    for i in 1..76 loop put(file,"="); end loop; new_line(file);
    for i in 1..integer32(len) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
        s : Solu_Info := Deep_Create(ls.all);
      begin
        Continue_along_Path(file,s,target,tol,false,patpp,patcp,nbq);
        Continue_End_Game(file,s,target,tol,false,
                          order,w(i),v(i),errv(i),endpp,endcp,nbq);
        Write_Next_Solution
          (file,cnt,s,tol_zero,tol_zero,nbfail,nbregu,nbsing,kind);
        text_io.flush(file);
        Append(res,res_last,Deep_Create(s));
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Copy(res,sols); Clear(res);
    put(file,"Number of path failures      : ");
    put(file,nbfail,1); new_line(file);
    put(file,"Number of regular solutions  : ");
    put(file,nbregu,1); new_line(file);
    put(file,"Number of singular solutions : ");
    put(file,nbsing,1); new_line(file);
    new_line(file);
  end Reporting_Toric_Continue;

end DoblDobl_IncFix_Continuation;
