with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Valid_Vector_Checks;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Path_Trackers;             use Standard_Path_Trackers;
with Continuation_Parameters;            use Continuation_Parameters;
with Standard_Continuation_Data;         use Standard_Continuation_Data;

package body Standard_IncFix_Continuation is

  sia_size : constant natural32 := 20000;

-- AUXILIAIRIES :

  function At_Infinity ( s : Solution; proj : boolean ) return boolean is

  -- DESCRIPTION :
  --   Decides whether a given solution lies at infinity.

  begin
    if proj then
      if AbsVal(s.v(s.v'last)) < 1.0/tol_endg_at_infinity
       then return true;
       else return false;
      end if;
    else
      for i in 1..s.n loop
        if AbsVal(s.v(i)) > tol_endg_at_infinity
         then return true;
        end if;
      end loop;
      return false;
    end if;
  end At_Infinity;

  function Equals ( s : in Solu_Info_Array; x : in Vector; i : in integer32;
                    d : in double_float; proj : in boolean ) return integer32 is

  -- DESCRIPTION :
  --   Returns the index j in the solution array s(s'first..i) of the
  --   solution which equals x.

    eq : boolean := false;
    j : integer32 := s'first;

  begin
    while j < i loop
      if Valid_Vector_Checks.is_valid(s(j).sol.v) then
        if not At_Infinity(s(j).sol.all,proj) then
          eq := true;
          if proj then
            for k in x'range loop
              if AbsVal(s(j).sol.v(k)/s(j).sol.v(x'last) - x(k)/x(x'last)) > d
               then eq := false; exit;
              end if;
            end loop;
          else
            for k in x'range loop
              if AbsVal(s(j).sol.v(k) - x(k)) > d
               then eq := false; exit;
              end if;
            end loop;
          end if;
          exit when eq;
        end if;
      end if;
      j := j + 1;
    end loop;
    return j;
  end Equals;

  procedure Add_Clustered ( i,n : in integer32; sols : in Solution_List;
                            clusols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Adds the ith start solution to the list clusols.
  --   The multiplicity of the solution equals i.

    ls : constant Link_to_Solution := Retrieve(sols,natural32(i));
    tol : constant double_float := 10.0**(-12);

  begin
    if ls /= null and then not Is_In(clusols,ls.all,tol) then
      ls.m := i;
      Add(clusols,ls.all);
    end if;
  end Add_Clustered;

  procedure Write_Bar ( file : in file_type ) is
  begin
    put(file,"========================================");
    put_line(file,"===================================");
  end Write_Bar;

  procedure Write_Statistics ( file : in file_type; i : in integer32;
                               nstep,nfail,niter,nsyst : in natural32 ) is

  -- DESCRIPTION :
  --   Writes the computing statistics of the ith path on file.

  begin
    put(file,"== "); put(file,i,1); put(file," = ");
    put(file," #step : "); put(file,nstep,3);
    put(file," #fail : "); put(file,nfail,2);
    put(file," #iter : "); put(file,niter,3);
    if nsyst /= niter
     then put(file," #syst : "); put(file,nsyst,3);
    end if;
    put(file," = ");
  end Write_Statistics;

  procedure Write_Diagnostics 
               ( file : in file_type; s : in out Solu_Info_Array;
                 c : in Corr_Pars; tol : in double_float; i : in integer32; 
                 proj : in boolean;
                 ninfi,nregu,nsing,nclus,nfail : in out natural32;
                 sols : in Solution_List; clusols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Writes the diagnostics for the ith solution.
  --   If it is a clustered solution, then the corresponding start solution
  --   will be added from the list sols to the list clusols.

    j : integer32;

  begin
    if At_Infinity(s(i).sol.all,proj) then
      ninfi := ninfi + 1;
      put_line(file,"at infinity ==");
    elsif ((s(i).cora <= c.epsax) or else (s(i).resa <= c.epsaf)
        or else (s(i).corr <= c.epsrx) or else (s(i).resr <= c.epsrf)) then
      j := Equals(s,s(i).sol.v,i,tol,proj);
      if j = i then
        if s(i).rcond > tol_endg_inverse_condition then
          nregu := nregu + 1;
          put_line(file,"regular solution ==");
        else
          nsing := nsing + 1;
          put_line(file,"singular solution ==");
        end if;
      elsif s(i).rcond < tol_endg_inverse_condition then
        nsing := nsing + 1;
        s(j).sol.m := s(j).sol.m + 1;
        s(i).sol.m := s(i).sol.m + 1;
        put(file,"multiple, see ");
        put(file,j,1); put_line(file," ==");
      else
        nclus := nclus + 1;
        put(file,"clustered with ");
        put(file,j,1); put_line(file," ==");
        Add_Clustered(i,s(i).sol.n,sols,clusols);
        Add_Clustered(j,s(j).sol.n,sols,clusols);
      end if;
    elsif s(i).rcond < tol_endg_inverse_condition then
      nfail := nfail + 1;
      put_line(file,"failure ==");
    else
      nfail := nfail + 1;
      put_line(file,"failure ==");
    end if;
  end Write_Diagnostics;

  procedure Write_Solution ( file : in file_type; s : in Solu_Info ) is

  -- DESCRIPTION :
  --   Writes the solution and the length of the path on file.

  begin
    put(file,"t : "); put(file,s.sol.t); new_line(file);
    put(file,"m : "); put(file,s.sol.m,1);
    put(file,"                  Length of path : ");
    put(file,s.length_path);
    new_line(file);
    put_line(file,"the solution for t : ");
    put_vector(file,s.sol.all);
    put(file,"==");
    put(file," err : "); put(file,s.cora,2,3,3);  put(file," =");
    put(file," rco : "); put(file,s.rcond,2,3,3); put(file," =");
    put(file," res : "); put(file,s.resa,2,3,3);  put_line(file," ==");
  end Write_Solution;

  procedure Diagnostics 
               ( s : in out Solu_Info_Array; c : in Corr_Pars;
                 tol : in double_float; i : in integer32; proj : in boolean;
                 ninfi,nregu,nsing,nclus,nfail : in out natural32;
                 sols : in Solution_List; clusols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Analyzes the ith solution.  If it is a clustered solution, then the
  --   corresponding start solution from the list sols will be added to the
  --   list clusols.

    j : integer32;

  begin
    if At_Infinity(s(i).sol.all,proj) then
      ninfi := ninfi + 1;
    elsif ((s(i).cora <= c.epsax) or else (s(i).resa <= c.epsaf)
        or else (s(i).corr <= c.epsrx) or else (s(i).resr <= c.epsrf)) then
      j := Equals(s,s(i).sol.v,i,tol,proj);
      if j = i then
        if s(i).rcond > tol_endg_inverse_condition
         then nregu := nregu + 1;
         else nsing := nsing + 1;
        end if;
      elsif s(i).rcond < tol_endg_inverse_condition then
        nsing := nsing + 1;
        s(j).sol.m := s(j).sol.m + 1;
        s(i).sol.m := s(i).sol.m + 1;
      else
        nclus := nclus + 1;
        Add_Clustered(i,s(i).sol.n,sols,clusols);
        Add_Clustered(j,s(j).sol.n,sols,clusols);
      end if;
    elsif s(i).rcond < tol_endg_inverse_condition then
      nfail := nfail + 1;
    else
      nfail := nfail + 1;
    end if;
  end Diagnostics;

  procedure Write_Summary_Diagnostics 
               ( file : in file_type;
                 ninfi,nregu,nsing,nfail,nclus : in natural32 ) is

  -- DESCRIPTION :
  --   Writes a summary after the continuation.

  begin
    put(file,"== ");
    put(file,"#regu : "); put(file,nregu,1); put(file," = " );
    put(file,"#sing : "); put(file,nsing,1); put(file," = " );
    put(file,"#clus : "); put(file,nclus,1); put(file," = " );
    put(file,"#infi : "); put(file,ninfi,1); put(file," = " );
    put(file,"#fail : "); put(file,nfail,1);
    put_line(file," == " );
  end Write_Summary_Diagnostics;

  procedure Merge_Clustered
               ( s : in out Solu_Info_Array; clusols : in Solution_List ) is

  -- DESCRIPTION :
  --   The new solutions, which were clustered before, are merged with
  --   the solution array, by using there multiplicity.

    tmp : Solution_List := clusols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      s(ls.m).sol := new Solution'(ls.all);
      tmp := Tail_Of(tmp);
    end loop;
  end Merge_Clustered;

-- TARGET ROUTINES :

  procedure Silent_Small_Continue
               ( sols : in out Solution_List; proj : in boolean;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) ) is

    sia : Solu_Info_Array(1..integer32(Length_Of(sols))) := Deep_Create(sols);
    ppa,pen : Pred_Pars;
    cpa,cen : Corr_Pars;
    tol : constant double_float := 10.0**(-10);
    w : integer32 := 1;
    dumv : Standard_Floating_Vectors.Link_to_Vector := null;
    err : double_float := 0.0;

    procedure LCont1 is
      new Linear_Single_Normal_Silent_Continue(Norm,H,dH,dH);
    procedure LCont2 is
      new Linear_Single_Conditioned_Silent_Continue(Norm,H,dH,dH);
    procedure LContN1 is
      new Linear_Multiple_Normal_Silent_Continue(Norm,H,dH,dH);

    procedure Rerun_Clustered
                 ( s : in out Solu_Info_Array;
                   clusols : in out Solution_List ) is

      oldmax : constant natural32 := max_reruns;
      oldblk : constant natural32 := block_size;

    begin
      condition := condition + 1;
      Continuation_Parameters.Tune(condition);
      max_reruns := oldmax - 1;
      block_size := Length_Of(clusols);
      Silent_Small_Continue(clusols,proj,nbq,target);
      block_size := oldblk;
      Merge_Clustered(s,clusols);
      Deep_Clear(clusols);
    end Rerun_Clustered;

    procedure Sequential_Continue
                 ( s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p1,p2 : in Pred_Pars; c_path,c_end : in Corr_Pars ) is

      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      for i in s'range loop
        LCont1(s(i),target,tol,proj,p1,c_path,nbq);
        LCont2(s(i),target,tol,proj,0,w,dumv,err,p2,c_end,nbq);
        Diagnostics(s,c_end,tol,i,proj,
                    ninfi,nregu,nsing,nclus,nfail,sols,clusols);
      end loop;
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(s,clusols);
      end if;
    end Sequential_Continue;

    procedure Continue_End_Game
                 ( s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p : in Pred_Pars; c : in Corr_Pars ) is

    -- DESCRIPTION :
    --   End game for the simultaneous path following.

      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      for i in s'range loop
        LCont2(s(i),target,tol,proj,0,w,dumv,err,p,c);
      end loop;
      for i in s'range loop
        Diagnostics(s,c,tol,i,proj,ninfi,nregu,nsing,nclus,nfail,sols,clusols);
      end loop;
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(s,clusols);
      end if;
    end Continue_end_Game;

    procedure Parallel_Continue
                 ( s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p_path,p_end : in Pred_Pars; c_path,c_end : in Corr_Pars ) is

    -- DESCRIPTION :
    --   This procedure implements the simultaneous continuation of
    --   different solution paths.

      nb,index : integer32;
      blck : integer32 := integer32(block_size);

    begin
      nb := 1; index := 0;
      while index < s'last loop
        if blck > s'last - index
         then blck := s'last - index;
        end if;
        declare
          sbk : Solu_Info_Array(1..blck) := s(index+1..index+blck);
        begin
          LContN1(sbk,target,tol,tol_path_distance,proj,p_path,c_path);
          Continue_end_Game(sbk,target,tol,p_end,c_end);
          s(index+1..index+blck) := sbk;
        end;
        nb := nb + 1;
        index := index + blck;
      end loop;
    end Parallel_Continue;

  begin
    ppa := Continuation_Parameters.Create_for_Path;
    pen := Continuation_Parameters.Create_End_Game;
    cpa := Continuation_Parameters.Create_for_Path;
    cen := Continuation_Parameters.Create_End_Game;
    if block_size = 1
     then Sequential_Continue(sia,target,tol,ppa,pen,cpa,cen);
     else Parallel_Continue(sia,target,tol,ppa,pen,cpa,cen);
    end if;
    Deep_Clear(sols);
    sols := Shallow_Create(sia);
  end Silent_Small_Continue;

  procedure Silent_Continue
               ( sols : in out Solution_List; proj : in boolean;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) ) is

    sials : constant Solu_Info_Array_List := Create(sols,sia_size);
    tmp : Solu_Info_Array_List := sials;
    ppa,pen : Pred_Pars;
    cpa,cen : Corr_Pars;
    tol : constant double_float := 1.0E-10;
    w : integer32 := 1;
    dumv : Standard_Floating_Vectors.Link_to_Vector := null;
    err : double_float := 0.0;

    procedure LCont1 is
      new Linear_Single_Normal_Silent_Continue(Norm,H,dH,dH);
    procedure LCont2 is
      new Linear_Single_Conditioned_Silent_Continue(Norm,H,dH,dH);
    procedure LContN1 is
      new Linear_Multiple_Normal_Silent_Continue(Norm,H,dH,dH);

    procedure Rerun_Clustered
                 ( s : in out Solu_Info_Array;
                   clusols : in out Solution_List ) is

      oldmax : constant natural32 := max_reruns;
      oldblk : constant natural32 := block_size;

    begin
      condition := condition + 1;
      Continuation_Parameters.Tune(condition);
      max_reruns := oldmax - 1;
      block_size := Length_Of(clusols);
      Silent_Continue(clusols,proj,nbq,target);
      block_size := oldblk;
      Merge_Clustered(s,clusols);
      Deep_Clear(clusols);
    end Rerun_Clustered;

    procedure Sequential_Continue
                 ( s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p1,p2 : in Pred_Pars; c_path,c_end : in Corr_Pars ) is

      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      for i in s'range loop
        LCont1(s(i),target,tol,proj,p1,c_path,nbq);
        LCont2(s(i),target,tol,proj,0,w,dumv,err,p2,c_end,nbq);
        Diagnostics(s,c_end,tol,i,proj,
                    ninfi,nregu,nsing,nclus,nfail,sols,clusols);
      end loop;
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(s,clusols);
      end if;
    end Sequential_Continue;

    procedure Continue_End_Game
                 ( s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p : in Pred_Pars; c : in Corr_Pars ) is

    -- DESCRIPTION :
    --   End game for the simultaneous path following.

      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      for i in s'range loop
        LCont2(s(i),target,tol,proj,0,w,dumv,err,p,c);
      end loop;
      for i in s'range loop
        Diagnostics(s,c,tol,i,proj,ninfi,nregu,nsing,nclus,nfail,sols,clusols);
      end loop;
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(s,clusols);
      end if;
    end Continue_end_Game;

    procedure Parallel_Continue
                 ( s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p_path,p_end : in Pred_Pars;
                   c_path,c_end : in Corr_Pars ) is

    -- DESCRIPTION :
    --   This procedure implements the simultaneous continuation of
    --   different solution paths.

      nb,index : integer32;
      blck : integer32 := integer32(block_size);

    begin
      nb := 1; index := 0;
      while index < s'last loop
        if blck > s'last - index
         then blck := s'last - index;
        end if;
        declare
          sbk : Solu_Info_Array(1..blck) := s(index+1..index+blck);
        begin
          LContN1(sbk,target,tol,tol_path_distance,proj,p_path,c_path);
          Continue_end_Game(sbk,target,tol,p_end,c_end);
          s(index+1..index+blck) := sbk;
        end;
        nb := nb + 1;
        index := index + blck;
      end loop;
    end Parallel_Continue;

  begin
    ppa := Continuation_Parameters.Create_for_Path;
    pen := Continuation_Parameters.Create_End_Game;
    cpa := Continuation_Parameters.Create_for_Path;
    cen := Continuation_Parameters.Create_End_Game;
    while not Is_Null(tmp) loop
      declare
        sia : constant Link_to_Solu_Info_Array := Head_Of(tmp);
      begin
        if block_size = 1
	 then Sequential_Continue(sia.all,target,tol,ppa,pen,cpa,cen);
         else Parallel_Continue(sia.all,target,tol,ppa,pen,cpa,cen);
        end if;
        Set_Head(tmp,sia);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Deep_Clear(sols);
    sols := Concat(sials);
  --exception
  --  when others => put_line("exception in Silent_Continue"); raise;
  end Silent_Continue;

  procedure Silent_Continue_with_Stop
               ( sols : in out Solution_List; proj : in boolean;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) ) is

    sia : Solu_Info_Array(1..integer32(Length_Of(sols))) := Deep_Create(sols);
    ppa,pen : Pred_Pars;
    cpa,cen : Corr_Pars;
    tol : constant double_float := 1.0E-10;
    w : integer32 := 1;
    dumv : Standard_Floating_Vectors.Link_to_Vector := null;
    err : double_float := 0.0;

    procedure LCont1 is
      new Linear_Single_Normal_Silent_Continue(Norm,H,dH,dH);
    procedure LCont2 is
      new Linear_Single_Conditioned_Silent_Continue(Norm,H,dH,dH);
    procedure LContN1 is
      new Linear_Multiple_Normal_Silent_Continue(Norm,H,dH,dH);

    procedure Rerun_Clustered
                 ( s : in out Solu_Info_Array;
                   clusols : in out Solution_List ) is

      oldmax : constant natural32 := max_reruns;
      oldblk : constant natural32 := block_size;

    begin
      condition := condition + 1;
      Continuation_Parameters.Tune(condition);
      max_reruns := oldmax - 1;
      block_size := Length_Of(clusols);
      Silent_Continue_with_Stop(clusols,proj,nbq,target);
      block_size := oldblk;
      Merge_Clustered(s,clusols);
      Deep_Clear(clusols);
    end Rerun_Clustered;

    procedure Sequential_Continue
                 ( s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p1,p2 : in Pred_Pars; c_path,c_end : in Corr_Pars ) is

      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      for i in s'range loop
        LCont1(s(i),target,tol,proj,p1,c_path,nbq);
        LCont2(s(i),target,tol,proj,0,w,dumv,err,p2,c_end,nbq);
        Diagnostics(s,c_end,tol,i,proj,
                    ninfi,nregu,nsing,nclus,nfail,sols,clusols);
        exit when Stop_Test(s(i).sol.all);
      end loop;
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(s,clusols);
      end if;
    end Sequential_Continue;

    procedure Continue_End_Game
                 ( s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p : in Pred_Pars; c : in Corr_Pars ) is

    -- DESCRIPTION :
    --   End game for the simultaneous path following.

      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      for i in s'range loop
        LCont2(s(i),target,tol,proj,0,w,dumv,err,p,c);
      end loop;
      for i in s'range loop
        Diagnostics(s,c,tol,i,proj,ninfi,nregu,nsing,nclus,nfail,sols,clusols);
      end loop;
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(s,clusols);
      end if;
    end Continue_end_Game;

    procedure Parallel_Continue
                 ( s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p_path,p_end : in Pred_Pars; c_path,c_end : in Corr_Pars ) is

    -- DESCRIPTION :
    --   This procedure implements the simultaneous continuation of
    --   different solution paths.

      nb,index : integer32;
      blck : integer32 := integer32(block_size);

    begin
      nb := 1; index := 0;
      while index < s'last loop
        if blck > s'last - index
         then blck := s'last - index;
        end if;
        declare
          sbk : Solu_Info_Array(1..blck) := s(index+1..index+blck);
        begin
          LContN1(sbk,target,tol,tol_path_distance,proj,p_path,c_path);
          Continue_end_Game(sbk,target,tol,p_end,c_end);
          s(index+1..index+blck) := sbk;
        end;
        nb := nb + 1;
        index := index + blck;
      end loop;
    end Parallel_Continue;

  begin
    ppa := Continuation_Parameters.Create_for_Path;
    pen := Continuation_Parameters.Create_End_Game;
    cpa := Continuation_Parameters.Create_for_Path;
    cen := Continuation_Parameters.Create_End_Game;
    if block_size = 1
     then Sequential_Continue(sia,target,tol,ppa,pen,cpa,cen);
     else Parallel_Continue(sia,target,tol,ppa,pen,cpa,cen);
    end if;
    Deep_Clear(sols);
    sols := Shallow_Create(sia);
  end Silent_Continue_with_Stop;

  procedure Reporting_Small_Continue
               ( file : in file_type; sols : in out Solution_List;
                 proj : in boolean; nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) ) is

    sia : Solu_Info_Array(1..integer32(Length_Of(sols))) := Deep_Create(sols);
    ppa,pen : Pred_Pars;
    cpa,cen : Corr_Pars;
    tol : constant double_float := 1.0E-10;
    w : integer32 := 1;
    dumv : Standard_Floating_Vectors.Link_to_Vector := null;
    err : double_float := 0.0;

    procedure LCont1 is
      new Linear_Single_Normal_Reporting_Continue(Norm,H,dH,dH);
    procedure LCont2 is
      new Linear_Single_Conditioned_Reporting_Continue(Norm,H,dH,dH);
    procedure LContN1 is
      new Linear_Multiple_Normal_Reporting_Continue(Norm,H,dH,dH);
   -- procedure CCont2 is
   --   new Circular_Single_Conditioned_Reporting_Continue(Norm,H,dH,dH);

    procedure Rerun_Clustered
                 ( file : in file_type; s : in out Solu_Info_Array;
                   clusols : in out Solution_List ) is

      oldmax : constant natural32 := max_reruns;
      oldblk : constant natural32 := block_size;

    begin
      condition := condition + 1;
      Continuation_Parameters.Tune(condition);
      max_reruns := oldmax - 1;
      block_size := Length_Of(clusols);
      Reporting_Small_Continue(file,clusols,proj,nbq,target);
      block_size := oldblk;
      Merge_Clustered(s,clusols);
      Deep_Clear(clusols);
    end Rerun_Clustered;

    procedure Sequential_Continue
                 ( file : in file_type; s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p1,p2 : in Pred_Pars; c_path,c_end : in Corr_Pars ) is

      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      Write_Bar(file);
      for i in s'range loop
        LCont1(file,s(i),target,tol,proj,p1,c_path,nbq);
        LCont2(file,s(i),target,tol,proj,0,w,dumv,err,p2,c_end,nbq);
        Write_Statistics(file,i,s(i).nstep,s(i).nfail,s(i).niter,s(i).nsyst);
        Write_Diagnostics(file,s,c_end,tol,i,proj,
                          ninfi,nregu,nsing,nclus,nfail,sols,clusols);
        Write_Solution(file,s(i));
      end loop;
      Write_Summary_Diagnostics(file,ninfi,nregu,nsing,nfail,nclus);
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(file,s,clusols);
      end if;
    end Sequential_Continue;

    procedure Continue_End_Game 
                 ( file : in file_type; s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p : in Pred_Pars; c : in Corr_Pars ) is
  
      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      for i in s'range loop
        LCont2(file,s(i),target,tol,proj,0,w,dumv,err,p,c);
      end loop;
      Write_Bar(file);
      for i in s'range loop
        Write_Statistics(file,i,s(i).nstep,s(i).nfail,s(i).niter,s(i).nsyst);
        Write_Diagnostics(file,s,c,tol,i,proj,
                          ninfi,nregu,nsing,nclus,nfail,sols,clusols);
        Write_Solution(file,s(i));
      end loop;
      put_line(file,"The computed solutions :");
      declare
        solus : Solution_List := Deep_Create(s);
      begin
        put(file,solus); Deep_Clear(solus);
      end;
      Write_Summary_Diagnostics(file,ninfi,nregu,nsing,nfail,nclus);
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(file,s,clusols);
      end if;
    end Continue_end_Game;

    procedure Parallel_Continue 
                 ( file : in file_type; s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p_path,p_end : in Pred_Pars; c_path,c_end : in Corr_Pars ) is

    -- DESCRIPTION :
    --   This procedure implements the simultaneous continuation of
    --   different solution paths.

      nb,index : integer32;
      blck : integer32 := integer32(block_size);

    begin
      nb := 1; index := 0;
      while index < s'last loop
        if blck > s'last - index
         then blck := s'last - index;
        end if;
        declare
          sbk : Solu_Info_Array(1..blck) := s(index+1..index+blck);
        begin
          LContN1(file,sbk,target,tol,tol_path_distance,proj,p_path,c_path);
          Continue_end_Game(file,sbk,target,tol,p_end,c_end);
          s(index+1..index+blck) := sbk;
        end;
        nb := nb + 1;
        index := index + blck;
      end loop;
    end Parallel_Continue;

  begin
    ppa := Continuation_Parameters.Create_for_Path;
    pen := Continuation_Parameters.Create_End_Game;
    cpa := Continuation_Parameters.Create_for_Path;
    cen := Continuation_Parameters.Create_End_Game;
    if block_size = 1
     then Sequential_Continue(file,sia,target,tol,ppa,pen,cpa,cen);
     else Parallel_Continue(file,sia,target,tol,ppa,pen,cpa,cen);
    end if;
    Deep_Clear(sols);
    sols := Shallow_Create(sia);
  end Reporting_Small_Continue;

  procedure Reporting_Continue
               ( file : in file_type; sols : in out Solution_List;
                 proj : in boolean; nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) ) is

    sials : constant Solu_Info_Array_List := Create(sols,sia_size);
    tmp : Solu_Info_Array_List := sials;
    ppa,pen : Pred_Pars;
    cpa,cen : Corr_Pars;
    tol : constant double_float := 1.0E-10;
    w : integer32 := 1;
    dumv : Standard_Floating_Vectors.Link_to_Vector := null;
    err : double_float := 0.0;

    procedure LCont1 is
      new Linear_Single_Normal_Reporting_Continue(Norm,H,dH,dH);
    procedure LCont2 is
      new Linear_Single_Conditioned_Reporting_Continue(Norm,H,dH,dH);
    procedure LContN1 is
      new Linear_Multiple_Normal_Reporting_Continue(Norm,H,dH,dH);
   -- procedure CCont2 is
   --   new Circular_Single_Conditioned_Reporting_Continue(Norm,H,dH,dH);

    procedure Rerun_Clustered
                 ( file : in file_type; s : in out Solu_Info_Array;
                   clusols : in out Solution_List ) is

      oldmax : constant natural32 := max_reruns;
      oldblk : constant natural32 := block_size;

    begin
      condition := condition + 1;
      Continuation_Parameters.Tune(condition);
      max_reruns := oldmax - 1;
      block_size := Length_Of(clusols);
      Reporting_Continue(file,clusols,proj,nbq,target);
      block_size := oldblk;
      Merge_Clustered(s,clusols);
      Deep_Clear(clusols);
    end Rerun_Clustered;

    procedure Sequential_Continue
                 ( file : in file_type; s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p1,p2 : in Pred_Pars; c_path,c_end : in Corr_Pars ) is

      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      Write_Bar(file);
      for i in s'range loop
        LCont1(file,s(i),target,tol,proj,p1,c_path,nbq);
        LCont2(file,s(i),target,tol,proj,0,w,dumv,err,p2,c_end,nbq);
        Write_Statistics(file,i,s(i).nstep,s(i).nfail,s(i).niter,s(i).nsyst);
        Write_Diagnostics(file,s,c_end,tol,i,proj,
                          ninfi,nregu,nsing,nclus,nfail,sols,clusols);
        Write_Solution(file,s(i));
      end loop;
      Write_Summary_Diagnostics(file,ninfi,nregu,nsing,nfail,nclus);
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(file,s,clusols);
      end if;
    --exception
    --  when others => put_line("exception in sequential continue");
    --                 flush(file); raise;
    end Sequential_Continue;

    procedure Continue_End_Game 
                 ( file : in file_type; s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p : in Pred_Pars; c : in Corr_Pars ) is
  
      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      for i in s'range loop
        LCont2(file,s(i),target,tol,proj,0,w,dumv,err,p,c);
      end loop;
      Write_Bar(file);
      for i in s'range loop
        Write_Statistics(file,i,s(i).nstep,s(i).nfail,s(i).niter,s(i).nsyst);
        Write_Diagnostics(file,s,c,tol,i,proj,
                          ninfi,nregu,nsing,nclus,nfail,sols,clusols);
        Write_Solution(file,s(i));
      end loop;
      put_line(file,"The computed solutions :");
      declare
        solus : Solution_List := Deep_Create(s);
      begin
        put(file,solus); Deep_Clear(solus);
      end;
      Write_Summary_Diagnostics(file,ninfi,nregu,nsing,nfail,nclus);
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(file,s,clusols);
      end if;
    --exception 
    --  when others => put_line("exception in continue end Game"); raise;
    end Continue_end_Game;

    procedure Parallel_Continue 
                 ( file : in file_type; s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p_path,p_end : in Pred_Pars; c_path,c_end : in Corr_Pars ) is

    -- DESCRIPTION :
    --   This procedure implements the simultaneous continuation of
    --   different solution paths.

      nb,index : integer32;
      blck : integer32 := integer32(block_size);

    begin
      nb := 1; index := 0;
      while index < s'last loop
        if blck > s'last - index
         then blck := s'last - index;
        end if;
        declare
          sbk : Solu_Info_Array(1..blck) := s(index+1..index+blck);
        begin
          LContN1(file,sbk,target,tol,tol_path_distance,proj,p_path,c_path);
          Continue_end_Game(file,sbk,target,tol,p_end,c_end);
          s(index+1..index+blck) := sbk;
        end;
        nb := nb + 1;
        index := index + blck;
      end loop;
    end Parallel_Continue;

  begin
    ppa := Continuation_Parameters.Create_for_Path;
    pen := Continuation_Parameters.Create_End_Game;
    cpa := Continuation_Parameters.Create_for_Path;
    cen := Continuation_Parameters.Create_End_Game;
    while not Is_Null(tmp) loop
      declare
        sia : constant Link_to_Solu_Info_Array := Head_Of(tmp);
      begin
        if block_size = 1
         then Sequential_Continue(file,sia.all,target,tol,ppa,pen,cpa,cen);
         else Parallel_Continue(file,sia.all,target,tol,ppa,pen,cpa,cen);
        end if;
        Set_Head(tmp,sia);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Deep_Clear(sols);
    sols := Concat(sials);
  --exception
  --  when others => put_line("exception in reporting continue"); raise;
  end Reporting_Continue;

  procedure Reporting_Continue_with_Stop
               ( file : in file_type; sols : in out Solution_List;
                 proj : in boolean; nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) ) is

    sia : Solu_Info_Array(1..integer32(Length_Of(sols))) := Deep_Create(sols);
    ppa,pen : Pred_Pars;
    cpa,cen : Corr_Pars;
    tol : constant double_float := 1.0E-10;
    w : integer32 := 1;
    dumv : Standard_Floating_Vectors.Link_to_Vector := null;
    err : double_float := 0.0;

    procedure LCont1 is
      new Linear_Single_Normal_Reporting_Continue(Norm,H,dH,dH);
    procedure LCont2 is
      new Linear_Single_Conditioned_Reporting_Continue(Norm,H,dH,dH);
    procedure LContN1 is
      new Linear_Multiple_Normal_Reporting_Continue(Norm,H,dH,dH);
   -- procedure CCont2 is
   --   new Circular_Single_Conditioned_Reporting_Continue(Norm,H,dH,dH);

    procedure Rerun_Clustered
                 ( file : in file_type; s : in out Solu_Info_Array;
                   clusols : in out Solution_List ) is

      oldmax : constant natural32 := max_reruns;
      oldblk : constant natural32 := block_size;

    begin
      condition := condition + 1;
      Continuation_Parameters.Tune(condition);
      max_reruns := oldmax - 1;
      block_size := Length_Of(clusols);
      Reporting_Continue_with_Stop(file,clusols,proj,nbq,target);
      block_size := oldblk;
      Merge_Clustered(s,clusols);
      Deep_Clear(clusols);
    end Rerun_Clustered;

    procedure Sequential_Continue
                 ( file : in file_type; s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p1,p2 : in Pred_Pars; c_path,c_end : in Corr_Pars ) is

      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      Write_Bar(file);
      for i in s'range loop
        LCont1(file,s(i),target,tol,proj,p1,c_path);
        LCont2(file,s(i),target,tol,proj,0,w,dumv,err,p2,c_end);
        Write_Statistics(file,i,s(i).nstep,s(i).nfail,s(i).niter,s(i).nsyst);
        Write_Diagnostics(file,s,c_end,tol,i,proj,
                          ninfi,nregu,nsing,nclus,nfail,sols,clusols);
        Write_Solution(file,s(i));
        exit when Stop_Test(s(i).sol.all);
      end loop;
      Write_Summary_Diagnostics(file,ninfi,nregu,nsing,nfail,nclus);
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(file,s,clusols);
      end if;
    end Sequential_Continue;

    procedure Continue_End_Game 
                 ( file : in file_type; s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p : in Pred_Pars; c : in Corr_Pars ) is
  
      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      for i in s'range loop
        LCont2(file,s(i),target,tol,proj,0,w,dumv,err,p,c);
      end loop;
      Write_Bar(file);
      for i in s'range loop
        Write_Statistics(file,i,s(i).nstep,s(i).nfail,s(i).niter,s(i).nsyst);
        Write_Diagnostics(file,s,c,tol,i,proj,
                          ninfi,nregu,nsing,nclus,nfail,sols,clusols);
        Write_Solution(file,s(i));
      end loop;
      put_line(file,"The computed solutions :");
      declare
        solus : Solution_List := Deep_Create(s);
      begin
        put(file,solus); Deep_Clear(solus);
      end;
      Write_Summary_Diagnostics(file,ninfi,nregu,nsing,nfail,nclus);
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(file,s,clusols);
      end if;
    end Continue_end_Game;

    procedure Parallel_Continue 
                 ( file : in file_type; s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p_path,p_end : in Pred_Pars; c_path,c_end : in Corr_Pars ) is

    -- DESCRIPTION :
    --   This procedure implements the simultaneous continuation of
    --   different solution paths.

      nb,index : integer32;
      blck : integer32 := integer32(block_size);

    begin
      nb := 1; index := 0;
      while index < s'last loop
        if blck > s'last - index
         then blck := s'last - index;
        end if;
        declare
          sbk : Solu_Info_Array(1..blck) := s(index+1..index+blck);
        begin
          LContN1(file,sbk,target,tol,tol_path_distance,proj,p_path,c_path);
          Continue_end_Game(file,sbk,target,tol,p_end,c_end);
          s(index+1..index+blck) := sbk;
        end;
        nb := nb + 1;
        index := index + blck;
      end loop;
    end Parallel_Continue;

  begin
    ppa := Continuation_Parameters.Create_for_Path;
    pen := Continuation_Parameters.Create_End_Game;
    cpa := Continuation_Parameters.Create_for_Path;
    cen := Continuation_Parameters.Create_End_Game;
    if block_size = 1
     then Sequential_Continue(file,sia,target,tol,ppa,pen,cpa,cen);
     else Parallel_Continue(file,sia,target,tol,ppa,pen,cpa,cen);
    end if;
    Deep_Clear(sols);
    sols := Shallow_Create(sia);
  end Reporting_Continue_with_Stop;

-- CONTINUATION WITH ESTIMATION OF PATH DIRECTIONS :

  procedure Silent_Toric_Continue
               ( sols : in out Solution_List; proj : in boolean;
                 w : in out Standard_Integer_Vectors.Vector;
                 v : in out Standard_Floating_VecVecs.VecVec;
                 errv : in out Standard_Floating_Vectors.Vector;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) ) is

    rtoric : constant integer32
           := integer32(Continuation_Parameters.endext_order);
    sia : Solu_Info_Array(1..integer32(Length_Of(sols))) := Deep_Create(sols);
    ppa,pen : Pred_Pars;
    cpa,cen : Corr_Pars;
    tol : constant double_float := 10.0**(-10);

    procedure LCont1 is
      new Linear_Single_Normal_Silent_Continue(Norm,H,dH,dH);
    procedure LCont2 is
      new Linear_Single_Conditioned_Silent_Continue(Norm,H,dH,dH);
    procedure LContN1 is
      new Linear_Multiple_Normal_Silent_Continue(Norm,H,dH,dH);

    procedure Rerun_Clustered
                 ( s : in out Solu_Info_Array;
                   clusols : in out Solution_List ) is

      oldmax : constant natural32 := max_reruns;
      oldblk : constant natural32 := block_size;

    begin
      condition := condition + 1;
      Continuation_Parameters.Tune(condition);
      max_reruns := oldmax - 1;
      block_size := Length_Of(clusols);
      Silent_Toric_Continue(clusols,proj,w,v,errv,nbq,target);
      block_size := oldblk;
      Merge_Clustered(s,clusols);
      Deep_Clear(clusols);
    end Rerun_Clustered;

    procedure Sequential_Continue
                 ( s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p1,p2 : in Pred_Pars; c_path,c_end : in Corr_Pars ) is

      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      for i in s'range loop
        LCont1(s(i),target,tol,proj,p1,c_path,nbq);
        LCont2(s(i),target,tol,proj,rtoric,w(i),v(i),errv(i),p2,c_end,nbq);
        Diagnostics(s,c_end,tol,i,proj,
                    ninfi,nregu,nsing,nclus,nfail,sols,clusols);
      end loop;
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(s,clusols);
      end if;
    end Sequential_Continue;

    procedure Continue_End_Game
                 ( s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p : in Pred_Pars; c : in Corr_Pars ) is

      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      for i in s'range loop
        LCont2(s(i),target,tol,proj,rtoric,w(i),v(i),errv(i),p,c);
      end loop;
      for i in s'range loop
        Diagnostics(s,c,tol,i,proj,ninfi,nregu,nsing,nclus,nfail,sols,clusols);
      end loop;
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(s,clusols);
      end if;
    end Continue_end_Game;

    procedure Parallel_Continue
                 ( s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p_path,p_end : in Pred_Pars; c_path,c_end : in Corr_Pars ) is

    -- DESCRIPTION :
    --   This procedure implements the simultaneous continuation of
    --   different solution paths.

      nb,index : integer32;
      blck : integer32 := integer32(block_size);

    begin
      nb := 1; index := 0;
      while index < s'last loop
        if blck > s'last - index
         then blck := s'last - index;
        end if;
        declare
          sbk : Solu_Info_Array(1..blck) := s(index+1..index+blck);
        begin
          LContN1(sbk,target,tol,tol_path_distance,proj,p_path,c_path);
          Continue_end_Game(sbk,target,tol,p_end,c_end);
          s(index+1..index+blck) := sbk;
        end;
        nb := nb + 1;
        index := index + blck;
      end loop;
    end Parallel_Continue;

  begin
    ppa := Continuation_Parameters.Create_for_Path;
    pen := Continuation_Parameters.Create_End_Game;
    cpa := Continuation_Parameters.Create_for_Path;
    cen := Continuation_Parameters.Create_End_Game;
    if block_size = 1
     then Sequential_Continue(sia,target,tol,ppa,pen,cpa,cen);
     else Parallel_Continue(sia,target,tol,ppa,pen,cpa,cen);
    end if;
    Deep_Clear(sols);
    sols := Shallow_Create(sia);
  end Silent_Toric_Continue;

  procedure Reporting_Toric_Continue
               ( file : in file_type; sols : in out Solution_List;
                 proj : in boolean;
                 w : in out Standard_Integer_Vectors.Vector;
                 v : in out Standard_Floating_VecVecs.VecVec;
                 errv : in out Standard_Floating_Vectors.Vector;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) ) is

    rtoric : constant integer32
           := integer32(Continuation_Parameters.endext_order);
    sia : Solu_Info_Array(1..integer32(Length_Of(sols))) := Deep_Create(sols);
    ppa,pen : Pred_Pars;
    cpa,cen : Corr_Pars;
    tol : constant double_float := 10.0**(-10);

    procedure LCont1 is
      new Linear_Single_Normal_Reporting_Continue(Norm,H,dH,dH);
    procedure LCont2 is
      new Linear_Single_Conditioned_Reporting_Continue(Norm,H,dH,dH);
    procedure LContN1 is
      new Linear_Multiple_Normal_Reporting_Continue(Norm,H,dH,dH);
   -- procedure CCont2 is
   --   new Circular_Single_Conditioned_Reporting_Continue(Norm,H,dH,dH);

    procedure Rerun_Clustered
                 ( file : in file_type; s : in out Solu_Info_Array;
                   clusols : in out Solution_List ) is

      oldmax : constant natural32 := max_reruns;
      oldblk : constant natural32 := block_size;

    begin
      condition := condition + 1;
      Continuation_Parameters.Tune(condition);
      max_reruns := oldmax - 1;
      block_size := Length_Of(clusols);
      Reporting_Toric_Continue(file,clusols,proj,w,v,errv,nbq,target);
      block_size := oldblk;
      Merge_Clustered(s,clusols);
      Deep_Clear(clusols);
    end Rerun_Clustered;

    procedure Sequential_Continue
                 ( file : in file_type; s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p1,p2 : in Pred_Pars; c_path,c_end : in Corr_Pars ) is

      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      Write_Bar(file);
      for i in s'range loop
        LCont1(file,s(i),target,tol,proj,p1,c_path);
        new_line(file);
        put(file,"running polyhedral end game on path ");
        put(file,i,1); put_line(file," :");
        LCont2(file,s(i),target,tol,proj,rtoric,w(i),v(i),errv(i),p2,c_end);
        Write_Statistics(file,i,s(i).nstep,s(i).nfail,s(i).niter,s(i).nsyst);
        Write_Diagnostics(file,s,c_end,tol,i,proj,
                          ninfi,nregu,nsing,nclus,nfail,sols,clusols);
        Write_Solution(file,s(i));
      end loop;
      Write_Summary_Diagnostics(file,ninfi,nregu,nsing,nfail,nclus);
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(file,s,clusols);
      end if;
    end Sequential_Continue;

    procedure Continue_End_Game 
                 ( file : in file_type; s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p : in Pred_Pars; c : in Corr_Pars ) is
  
      ninfi,nregu,nsing,nfail,nclus : natural32 := 0;
      clusols : Solution_List;

    begin
      for i in s'range loop
        LCont2(file,s(i),target,tol,proj,rtoric,w(i),v(i),errv(i),p,c);
      end loop;
      Write_Bar(file);
      for i in s'range loop
        Write_Statistics(file,i,s(i).nstep,s(i).nfail,s(i).niter,s(i).nsyst);
        Write_Diagnostics(file,s,c,tol,i,proj,
                          ninfi,nregu,nsing,nclus,nfail,sols,clusols);
        Write_Solution(file,s(i));
      end loop;
      put_line(file,"The computed solutions :");
      declare
        solus : Solution_List := Deep_Create(s);
      begin
        put(file,solus); Deep_Clear(solus);
      end;
      Write_Summary_Diagnostics(file,ninfi,nregu,nsing,nfail,nclus);
      if (nclus > 0) and then (max_reruns > 0)
       then Rerun_Clustered(file,s,clusols);
      end if;
    end Continue_end_Game;

    procedure Parallel_Continue 
                 ( file : in file_type; s : in out Solu_Info_Array;
                   target : in Complex_Number; tol : in double_float;
                   p_path,p_end : in Pred_Pars; c_path,c_end : in Corr_Pars ) is

    -- DESCRIPTION :
    --   This procedure implements the simultaneous continuation of
    --   different solution paths.

      nb,index : integer32;
      blck : integer32 := integer32(block_size);

    begin
      nb := 1; index := 0;
      while index < s'last loop
        if blck > s'last - index
         then blck := s'last - index;
        end if;
        declare
          sbk : Solu_Info_Array(1..blck) := s(index+1..index+blck);
        begin
          LContN1(file,sbk,target,tol,tol_path_distance,proj,p_path,c_path);
          Continue_end_Game(file,sbk,target,tol,p_end,c_end);
          s(index+1..index+blck) := sbk;
        end;
        nb := nb + 1;
        index := index + blck;
      end loop;
    end Parallel_Continue;

  begin
    ppa := Continuation_Parameters.Create_for_Path;
    pen := Continuation_Parameters.Create_End_Game;
    cpa := Continuation_Parameters.Create_for_Path;
    cen := Continuation_Parameters.Create_End_Game;
    if block_size = 1
     then Sequential_Continue(file,sia,target,tol,ppa,pen,cpa,cen);
     else Parallel_Continue(file,sia,target,tol,ppa,pen,cpa,cen);
    end if;
    Deep_Clear(sols);
    sols := Shallow_Create(sia);
  end Reporting_Toric_Continue;

end Standard_IncFix_Continuation;
