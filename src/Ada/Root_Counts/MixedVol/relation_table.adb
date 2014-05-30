with text_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Simplex_Pivoting;                  use Simplex_Pivoting;

package body Relation_Table is

  procedure RelTable
               ( nVar,nSpt : in integer32;
                 Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                 SptIdx : in Standard_Integer_Vectors.Link_to_Vector;
                 lft : in Standard_Floating_Vectors.Link_to_Vector;
                 RelTab : out Boolean_Matrix;
                 L0 : in out Link_to_L0_IML ) is

    nv1 : constant integer32 := nVar+1;
    ynNuPt : boolean;
    jj,info,ell,TstPt,Strt1PtTst,nPtIn : integer32;
    dtmp : double_float;
    labels : Standard_Integer_Vectors.Link_to_Vector
           := new Standard_Integer_Vectors.Vector(0..1);
    Bidx_s : Standard_Integer_Vectors.Link_to_Vector
           := new Standard_Integer_Vectors.Vector(0..nVar-1);
    Bidx : Standard_Integer_Vectors.Link_to_Vector
         := new Standard_Integer_Vectors.Vector(0..nVar+1);
    LPidx : Standard_Integer_Vectors.Link_to_Vector
          := new Standard_Integer_Vectors.Vector(0..SptIdx(nSpt));
    x_s : Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector(0..nVar-1);
    c : Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector(0..nVar+1);
    x : Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector(0..nVar+1);
    Binv_s : Standard_Floating_Matrices.Matrix(0..nVar-1,0..nVar-1);
    a : Standard_Floating_Matrices.Link_to_Matrix
      := new Standard_Floating_Matrices.Matrix
               (0..SptIdx(nSpt),0..nVar+1);
    Binv : Standard_Floating_Matrices.Link_to_Matrix
         := new Standard_Floating_Matrices.Matrix(0..nVar+1,0..nVar+1);

  begin
    for i in 0..(SptIdx(nSpt)-1) loop
      for j in 0..(SptIdx(nSpt)-1) loop
        RelTab(i,j) := false;
      end loop;
    end loop;
    for j in SptIdx(1)..(SptIdx(nSpt)-1) loop                -- [0] changing
      for i in 0..(nVar-1) loop
        a(j,i) := -double_float(Spt(j)(i));
      end loop;
      a(j,nVar) := 1.0;                                        -- WRT alpha0
      a(j,nVar+1) := lft(j);
    end loop;
    for CurSpt in 0..(nSpt-1) loop
      for FixPt in SptIdx(CurSpt)..(SptIdx(CurSpt+1)-1) loop
        nPtIn := -1;                            -- #constraints involved - 1
        for j in SptIdx(CurSpt)..(FixPt-1) loop
          if RelTab(FixPt,j) then
            nPtIn := nPtIn + 1;
            LPidx(nPtIn) := j;
            for i in 0..(nVar-1) loop
              a(j,i) := double_float(Spt(FixPt)(i) - Spt(j)(i));
            end loop;
            a(j,nVar) := 1.0;                                  -- for alpha0
            a(j,nv1) := lft(j) - lft(FixPt);            -- the constant term
          end if;
        end loop;
        Strt1PtTst := nPtIn + 1;      -- starting index of pts for 1-pt test
        for j in (FixPt+1)..(SptIdx(CurSpt+1)-1) loop
          nPtIn := nPtIn + 1;
          LPidx(nPtIn) := j;
          for i in 0..(nVar-1) loop
            a(j,i) := double_float(Spt(FixPt)(i) - Spt(j)(i));
          end loop;
          a(j,nVar) := 1.0;
          a(j,nv1) := lft(j) - lft(FixPt);
        end loop;
       -- To extend the point by using 1-point test
       -- To find a nondegenerate solution for 1-point tests of the points
        RlTbLP2_e(nPtIn+1,nVar,SptIdx(nSpt),a,nv1,LPidx,Bidx,x,Binv,info);
        if FixPt < SptIdx(CurSpt+1)-1 then
                  -- record any possible points in Bidx passed 1-point test
          if CurSpt = 0 then
            ynNuPt := false;
            for i in 0..(nVar-1) loop
              if (Bidx(i) > FixPt) and then (not RelTab(FixPt,Bidx(i))) then  
                ynNuPt := true;
                RelTab(Bidx(i),FixPt) := true;
                RelTab(FixPt,Bidx(i)) := true;
              end if;
            end loop;
            if ynNuPt then
              labels(0) := FixPt;
              for i in 0..(nVar-1) loop
                if Bidx(i) > FixPt then
                  labels(1) := Bidx(i);
                  L0_Add2(L0,labels,nVar,Bidx,x,Binv);  -- add 2 points
                end if;
              end loop;
            end if;
          else
            for i in 0..(nVar-1) loop
              if Bidx(i) > FixPt then
                RelTab(Bidx(i),FixPt) := true;
                RelTab(FixPt,Bidx(i)) := true;
                for j in i+1..(nVar-1) loop
                  if Bidx(j) > FixPt then
                    RelTab(Bidx(i),Bidx(j)) := true;
                    RelTab(Bidx(j),Bidx(i)) := true;
                  end if;
                end loop;
              end if;
            end loop;
          end if;
             -- x, Bidx, Binv are available for 1-point test on other points
          if info >= 0 then
            for j in Strt1PtTst..nPtIn loop
              TstPt := LPidx(j);
              if not RelTab(FixPt,TstPt) then
                for i in 0..(nVar-1) loop
                  c(i) := -a(TstPt,i);
                end loop;
                if CurSpt = 0
                 then dlp1_1pt_s(nPtIn+1,nVar,a,nv1,c,LPidx,FixPt,TstPt,
                                 Bidx,x,Binv,RelTab,L0);
                 else dlp1_1pt_i(nPtIn+1,nVar,a,nv1,c,LPidx,FixPt,TstPt,
                                 Bidx,x,Binv,RelTab);
                end if;
              end if;
            end loop;
          else
            for j in Strt1PtTst..nPtIn loop
              TstPt := LPidx(j);
              if not RelTab(FixPt,TstPt) then
                for i in 0..(nVar-1) loop
                  c(i) := -a(TstPt,i);
                end loop;
                if(CurSpt = 0)
                 then dlp2_1pt_s(nPtIn+1,nVar,a,nv1,c,LPidx,FixPt,TstPt, 
                                 Bidx,x,Binv,RelTab,L0);
                 else dlp2_1pt_i(nPtIn+1,nVar,a,nv1,c,LPidx,FixPt,TstPt,
                                 Bidx,x,Binv,RelTab);  
                end if;
              end if;
            end loop;
          end if;
        end if;
        for i in 0..(nVar-1) loop    -- save the starting point for next LP
          for j in 0..(nVar-1) loop
            Binv_s(i,j) := Binv(i,j);
          end loop;
          Bidx_s(i) := Bidx(i);
          x_s(i) := x(i);
        end loop;
           -- check point FixPt in S_CurSpt and each point in S_i (i>CurSpt)
        nPtIn := -1;       -- To get all constraints of pts related to FixPt
        for i in SptIdx(CurSpt)..(SptIdx(CurSpt+1)-1) loop
          if RelTab(FixPt,i) then
            nPtIn := nPtIn + 1;
            LPidx(nPtIn) := i;
            a(i,nVar) := 0.0;                                   -- no alpha0
          end if;  
        end loop;
        Strt1PtTst := nPtIn + 1;      -- starting index of pts for 1-pt test
        for NxSpt in CurSpt+1..(nSpt-1) loop
          nPtIn := Strt1PtTst - 1;
          for i in SptIdx(NxSpt)..(SptIdx(NxSpt+1)-1) loop
            nPtIn := nPtIn + 1;
            LPidx(nPtIn) := i;
          end loop;
        -- The part of Ax<=B for this support was formed at the beginning
        -- to do the 1-point test for the points in S_(NxSpt)
        -- to form the starting point from saved starting Bidx, x, Binv
          for i in 0..(nVar-1) loop
            Bidx(i) := Bidx_s(i);
            x(i) := x_s(i);
          end loop;
          ell := -1;
          x(nVar) := DBL_MAX;  
          for j in Strt1PtTst..nPtIn loop
            jj := LPidx(j);
            dtmp := a(jj,nv1);
            for i in 0..(nVar-1) loop
              dtmp := dtmp - a(jj,i)*x_s(i);
            end loop;
            if dtmp < x(nVar) then
              x(nVar) := dtmp;  
              ell := jj;
            end if;
          end loop;
          Bidx(nVar) := ell;
          for i in 0..(nVar-1) loop
            for j in 0..(nVar-1) loop
              Binv(i,j) := Binv_s(i,j);
            end loop;
            Binv(i,nVar) := 0.0;
          end loop;
          for i in 0..(nVar-1) loop
            Binv(nVar,i) := 0.0;
          end loop;
          for j in 0..(nVar-1) loop
            Binv(j,nVar) := 0.0;
            for i in 0..(nVar-1) loop
              Binv(j,nVar) := Binv(j,nVar) - Binv(j,i)*a(ell,i);
            end loop;
          end loop;
          Binv(nVar,nVar) := 1.0;
             -- find a nondegenerate solution for 1-point tests of the points
          TstPt := LPidx(Strt1PtTst);
          for i in 0..nVar loop
            c(i) := -a(TstPt,i);
          end loop;
          RlTbLP2_a(nPtIn+1,nVar+1,a,nv1,c,LPidx,Bidx,x,Binv,info);
                    -- record any possible points in Bidx passed 1-point test
          for i in 0..nVar loop
            if Bidx(i) >= SptIdx(NxSpt) then
              if Bidx(i) >= TstPt then
                RelTab(Bidx(i),FixPt) := true;
                RelTab(FixPt,Bidx(i)) := true;
              end if;
              for j in (i+1)..nVar loop
                if Bidx(j) >= SptIdx(NxSpt) then
                  RelTab(Bidx(i),Bidx(j)) := true;
                  RelTab(Bidx(j),Bidx(i)) := true;
                end if;
              end loop;
            end if;
          end loop;
          if info >= 0 then           -- do the 1-point test for other points
            for j in (Strt1PtTst+1)..nPtIn loop
              TstPt := LPidx(j);
              if not RelTab(FixPt,TstPt) then
                for i in 0..nVar loop
                  c(i) := -a(TstPt,i);
                end loop;
                dlp1_1pt_i(nPtIn+1,nVar+1,a,nv1,c,LPidx,FixPt,TstPt,
                           Bidx,x,Binv,RelTab);
              end if;
            end loop;
          else
            for j in (Strt1PtTst+1)..nPtIn loop
              TstPt := LPidx(j);
              if not RelTab(FixPt,TstPt) then
                for i in 0..nVar loop
                  c(i) := -a(TstPt,i);
                end loop;
                dlp2_1pt_i(nPtIn+1,nVar+1,a,nv1,c,LPidx,FixPt,TstPt,
                           Bidx,x,Binv,RelTab);
              end if;
            end loop;
          end if;
        end loop;
      end loop;
    end loop;
    Standard_Integer_Vectors.Clear(labels);
    Standard_Integer_Vectors.Clear(Bidx_s);
    Standard_Integer_Vectors.Clear(Bidx);
    Standard_Integer_Vectors.Clear(LPidx);
    Standard_Floating_Vectors.Clear(c);
    Standard_Floating_Vectors.Clear(x);
    Standard_Floating_Vectors.Clear(x_s);
    Standard_Floating_Matrices.Clear(a);
    Standard_Floating_Matrices.Clear(Binv);
  end RelTable;

  procedure RlTbLP2_a
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nv1 : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 LPidx : in Standard_Integer_Vectors.Link_to_Vector;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 info : out integer32 ) is

    ibrk : boolean;
    ell,k,k1,flag : integer32;
    sigj,vmax,bot,dtmp : double_float;
    v : Standard_Floating_Vectors.Vector(0..na-1);

  begin
    loop
      ibrk := false;              -- step 1 : search for outgoing constraint
      for i in 0..(na-1) loop
        if Bidx(i) = -1
         then ibrk := true; exit;
        end if;
      end loop;
      if ibrk then
        vmax := -1.0;
        for i in 0..(na-1) loop
          if Bidx(i) = -1 then
            dtmp := c(0)*Binv(i,0);
            for j in 1..(na-1) loop
              dtmp := dtmp + c(j)*Binv(i,j);
            end loop;
            if abs(dtmp) > vmax then
              bot := dtmp;
              vmax := abs(bot);
              k := i;
            end if;
          end if;
        end loop;
        if vmax > eps then
          if bot >= 0.0 then
            for i in 0..(na-1) loop
              v(i) := Binv(k,i);
            end loop;
          else
            for i in 0..(na-1) loop
              v(i) := -Binv(k,i);
            end loop;
          end if;
        else
          ibrk := false;
        end if;
      end if;
      if not ibrk then
        vmax := -DBL_MAX;
        for i in 0..(na-1) loop
          if Bidx(i) /= -1 then
            dtmp := c(0)*Binv(i,0);
            for j in 1..(na-1) loop
              dtmp := dtmp + c(j)*Binv(i,j);
            end loop;
            if dtmp > vmax then
              vmax := dtmp;
              k := i;
            end if;
          end if;
        end loop;
        if vmax < eps
	 then exit;              -- found optimal solution, solve smaller LP
        end if;
        for i in 0..(na-1) loop
	  v(i) := Binv(k,i);
        end loop;
      end if;
      Search_Incoming(ma,na,nv1,k,LPidx,Bidx,x,v,a,Binv,sigj,ell);
      if ell < 0
       then raise unbounded_LP;
      end if;
      for i in 0..(na-1) loop               -- step 3 : update x, Bidx, Binv
        x(i) := x(i) - sigj*v(i);
      end loop;
      Update_Base(Bidx,Binv,na,k,ell,a);
    end loop;
    k1 := -1;
    info := 0;
    for i in 0..(na-1) loop
      if Bidx(i) > -1
       then info := info + 1;
      end if;
    end loop;
    if info = na
     then return;                                       -- case of full rank
    end if;
    k := 0;
    while k < na  loop
      if (Bidx(k) > -1) or (k = k1) then
        k := k + 1;                                             -- continue;
      else
        for i in 0..(na-1) loop
          v(i) := Binv(k,i);
        end loop;
        flag := 1;
        ell := -1;                    -- step 2 : search incoming constraint
        while ell = -1 loop
          Search_Incoming(ma,na,nv1,k,LPidx,Bidx,x,v,a,Binv,sigj,ell);
          ibrk := false;
          if ell < 0 then
            if flag = 1 then
              for i in 0..(na-1) loop
                v(i) := -Binv(k,i);
              end loop;
              flag := 0;
            else
              k := k + 1;
              ibrk := true;
              exit;
            end if;
          end if;
        end loop;
        if not ibrk then                               -- if(ibrk) continue;
          for i in 0..(na-1) loop           -- step 3 : update x, Bidx, Binv
            x(i) := x(i) - sigj*v(i);
          end loop;
          Update_Base(Bidx,Binv,na,k,ell,a);
          info := info + 1;
          if info = na
           then return;
           else k1 := k; k := 0;
          end if;
        end if;
      end if;
    end loop;
    info := -info;                              -- the case of not full rank
  end RlTbLP2_a;

  procedure RlTbLP2_e
               ( ma,na,NumCol : in integer32;
                 a : in out Standard_Floating_Matrices.Link_to_Matrix;
                 nv1 : in integer32; 
                 LPidx : in Standard_Integer_Vectors.Link_to_Vector;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 info : out integer32 ) is

    ibrk : boolean;
    ell,ii,k,k1,flag : integer32;
    sigj,vmax,bot,dtmp : double_float;
    nap1 : constant integer32 := na + 1;
    map1 : constant integer32 := ma + 1;
    v : Standard_Floating_Vectors.Vector(0..nap1);

  begin                          -- Expand the system (add variable epsilon)
            --  Assumed a, b, Binv, x have extra 1 dimension to hold epsilon 
    info := 0;
    for i in 0..(na-1) loop
      a(NumCol,i) := 0.0;                                     -- -epsilon<=0
    end loop;
    a(NumCol,na) := -1.0;
    a(NumCol,nv1) := 0.0;
    LPidx(ma) := NumCol;
    dtmp := 0.0;
    ell := NumCol;                     -- the index of the equality equation
    for i in 0..(ma-1) loop
      ii := LPidx(i);
      if a(ii,nv1) < -eps then
        a(ii,nap1-1) := -1.0;
        if -a(ii,nv1) > dtmp then
          dtmp := -a(ii,nv1);
          ell := ii;
        end if;
      else
        a(ii,nap1-1) := 0.0;
      end if;
    end loop;
    for i in 0..(na-1) loop                  -- form the first x, Bidx, Binv
      Bidx(i) := -1;               
    end loop;
    Bidx(nap1-1) := ell;
    for i in 0..(na-1) loop
      x(i) := 0.0;
    end loop;
    x(na) := dtmp;
    for i in 0..(nap1-1) loop
      for j in 0..(nap1-1) loop
        Binv(i,j) := 0.0;
      end loop;
      Binv(i,i) := 1.0;
    end loop;
    for j in 0..(na-1) loop
      Binv(j,nap1-1) := a(ell,j);
    end loop;
    Binv(nap1-1,nap1-1) := -1.0;
    loop                          -- apply the LP algorithm to the larger LP
      ibrk := false;        --  step 1: To find which constraint to kick out
      for i in 0..nap1-1 loop
        if Bidx(i) = -1
         then ibrk := true; exit;
        end if;
      end loop;
      if ibrk then
        vmax := -1.0;
        for i in 0..(nap1-1) loop
          if Bidx(i) = -1 then
            dtmp := Binv(i,nap1-1);
            if abs(dtmp) > vmax then
              bot := dtmp;
              vmax := abs(bot);
              k := i;
            end if;
          end if;
        end loop;
        if vmax > eps then
          if bot >= 0.0 then
            for i in 0..(nap1-1) loop
              v(i) := Binv(k,i);
            end loop;
          else
            for i in 0..(nap1-1) loop
              v(i) := -Binv(k,i);
            end loop;
          end if;
        else
          ibrk := false;
        end if;
      end if;
      if not ibrk then
        vmax := -DBL_MAX;
        for i in 0..(nap1-1) loop
          if Bidx(i) /= -1 then
            dtmp := Binv(i,nap1-1);
            if dtmp > vmax then
              vmax := dtmp;
              k := i;
            end if;
          end if;
        end loop;
        if vmax < eps
	 then exit;              -- found optimal solution, solve smaller LP
        end if;
        for i in 0..(nap1-1) loop
          v(i) := Binv(k,i);
        end loop;
      end if;
      Search_Incoming(map1,nap1,nv1,k,LPidx,Bidx,x,v,a,Binv,sigj,ell);
      if ell < 0
       then raise unbounded_LP;
      end if;
      for i in 0..(nap1-1) loop            -- step 3 : update x, Bidx, Binv
        x(i) := x(i) - sigj*v(i);
      end loop;
      Update_Base(Bidx,Binv,nap1,k,ell,a);
    end loop;
    if Bidx(nap1-1) /= NumCol then     -- change x, Bidx, Binv into the ones
      k := -1;                       -- without the x(na) = epsilon variable
      ii := -1;
      while ii < nap1-1 loop
        ii := ii + 1;
        if Bidx(ii) = NumCol then
          k := ii;
          ii := nap1-1;
          exit;
        end if;
      end loop;
      if k = -1 then                                          -- for testing
        text_io.put_line("RlTbLP2_e : no index ma in Bidx"); return;
      end if;                               -- should not happen in our case
      Bidx(k) := Bidx(nap1-1);
      for i in 0..(na-1) loop
        Binv(k,i) := Binv(nap1-1,i);
      end loop;
    end if;
    for i in 0..(ma-1) loop
      a(LPidx(i),nap1-1) := 0.0;                 -- reset values WRT epsilon
    end loop;
    for i in 0..(na-1) loop
      if Bidx(i) > -1
       then info := info + 1;
      end if;
    end loop;
    if info = na
     then return;                                       -- case of full rank
    end if;
    k1 := -1;
    k := 0;
    while k < na loop
      if (Bidx(k) > -1) or (k = k1) then
        k := k + 1;                                             -- continue;
      else
        for i in 0..(na-1) loop
          v(i) := Binv(k,i);
        end loop;
        flag := 1;
        ell := -1;                    -- step 2 : search incoming constraint
        while ell = -1 loop
          Search_Incoming(ma,na,nv1,k,LPidx,Bidx,x,v,a,Binv,sigj,ell);
          ibrk := false;
          if ell < 0 then
            if flag = 1 then
              for i in 0..(na-1) loop
                v(i) := -Binv(k,i);
              end loop;
              flag := 0;
            else
              k := k + 1;
              ibrk := true;
              exit;
            end if;
          end if;
        end loop;
        if not ibrk then                               -- if(ibrk) continue;
          for i in 0..(na-1) loop           -- step 3 : update x, Bidx, Binv
            x(i) := x(i) - sigj*v(i); 
          end loop;
          Update_Base(Bidx,Binv,na,k,ell,a);
          info := info + 1;
          if info = na
           then return;
           else k1 := k; k := 0;
          end if;
        end if;
      end if;
    end loop;
    info := -info;                              -- the case of not full rank
  end RlTbLP2_e;

  procedure dlp2_1pt_i
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nv1 : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 LPidx : in Standard_Integer_Vectors.Link_to_Vector;
                 FixPt,TstPt : in integer32;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 RelTab : in out Boolean_Matrix ) is

    k,ell : integer32;
    vmax,sigj : double_float;

  begin
    loop
      Search_Outgoing(na,c,Bidx,Binv,vmax,k);
      if vmax < eps
       then return;                           -- leave with optimal solution
      end if;
      Search_Incoming(ma,na,nv1,k,LPidx,Bidx,x,a,Binv,sigj,ell);
      for i in 0..(na-1) loop                        -- update x, Bidx, Binv
        x(i) := x(i) - sigj*Binv(k,i);
      end loop;
      Update_Base(Bidx,Binv,na,k,ell,a);
      if ell >= TstPt then
        for i in 0..(na-1) loop
          if (((i /= k) and (Bidx(i) > -1)) 
	     and then not RelTab(ell,Bidx(i))) then
            RelTab(ell,Bidx(i)) := true;
            RelTab(Bidx(i),ell) := true;
          end if;
        end loop;
        if not RelTab(ell,FixPt) then
          RelTab(ell,FixPt) := true;
          RelTab(FixPt,ell) := true;
        end if;
      end if;
    end loop;
  end dlp2_1pt_i;

  procedure dlp2_1pt_s
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nv1 : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 LPidx : in Standard_Integer_Vectors.Link_to_Vector;
                 FixPt,TstPt : in integer32;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 RelTab : in out Boolean_Matrix;
                 L0 : in out Link_to_L0_IML ) is

    k,ell : integer32;
    vmax,sigj : double_float;
    labels : Standard_Integer_Vectors.Link_to_Vector
           := new Standard_Integer_Vectors.Vector(0..1);

  begin 
    loop
      Search_Outgoing(na,c,Bidx,Binv,vmax,k);
      if vmax < eps then
        Standard_Integer_Vectors.Clear(labels);
        return;                               -- leave with optimal solution
      end if;
      Search_Incoming(ma,na,nv1,k,LPidx,Bidx,x,a,Binv,sigj,ell);
      for i in 0..(na-1) loop                        -- update x, Bidx, Binv
        x(i) := x(i) - sigj*Binv(k,i);
      end loop;
      Update_Base(Bidx,Binv,na,k,ell,a);
      if (ell >= TstPt) and not RelTab(ell,FixPt) then -- save Bidx, x, Binv 
         RelTab(ell,FixPt) := true;
         RelTab(FixPt,ell) := true;
         labels(0) := FixPt;
         labels(1) := ell;               -- labels(0) < labels(1) is assumed
         L0_Add2(L0,labels,na,Bidx,x,Binv);                  -- add 2 points
      end if;
    end loop;
  end dlp2_1pt_s;

  procedure dlp1_1pt_i
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nv1 : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 LPidx : in Standard_Integer_Vectors.Link_to_Vector;
                 FixPt,TstPt : in integer32;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 RelTab : in out Boolean_Matrix ) is

    k,ell : integer32;
    vmax,sigj : double_float;

  begin
    loop
      Search_Outgoing(na,c,Binv,vmax,k);
      if vmax < eps
       then return;                             -- out with optimal solution
      end if;
      Search_Incoming(ma,na,nv1,k,LPidx,Bidx,x,a,Binv,sigj,ell);
      for i in 0..(na-1) loop                        -- update x, Bidx, Binv 
        x(i) := x(i) - sigj*Binv(k,i);
      end loop;
      Update_Base(Bidx,Binv,na,k,ell,a);
      if ell >= TstPt then
        for i in 0..(na-1) loop
          if (i /= k) and then not RelTab(ell,Bidx(i)) then
            RelTab(ell,Bidx(i)) := true;
            RelTab(Bidx(i),ell) := true;
          end if;
        end loop;
        if not RelTab(ell,FixPt) then
          RelTab(ell,FixPt) := true;
          RelTab(FixPt,ell) := true;
        end if;
      end if;
    end loop;
  end dlp1_1pt_i;

  procedure dlp1_1pt_s
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nv1 : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 LPidx : in Standard_Integer_Vectors.Link_to_Vector;
                 FixPt,TstPt : in integer32;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 RelTab : in out Boolean_Matrix;
                 L0 : in out Link_to_L0_IML ) is

    k,ell : integer32;
    vmax,sigj : double_float;
    labels : Standard_Integer_Vectors.Link_to_Vector
           := new Standard_Integer_Vectors.Vector(0..1);

  begin
    loop
      Search_Outgoing(na,c,Binv,vmax,k); 
      if vmax < eps then
        Standard_Integer_Vectors.Clear(labels);
        return;                                 -- out with optimal solution
      end if;
      Search_Incoming(ma,na,nv1,k,LPidx,Bidx,x,a,Binv,sigj,ell);
      for i in 0..(na-1) loop                        -- update x, Bidx, Binv
        x(i) := x(i) - sigj*Binv(k,i);
      end loop;
      Update_Base(Bidx,Binv,na,k,ell,a);
      if (ell >= TstPt) and not RelTab(ell,FixPt) then -- save Bidx, x, Binv 
        RelTab(ell,FixPt) := true;
        RelTab(FixPt,ell) := true;
        labels(0) := FixPt;
        labels(1) := ell;                -- labels[0] < labels[1] is assumed
        L0_Add2(L0,labels,na,Bidx,x,Binv);                   -- add 2 points 
      end if;
    end loop;
  end dlp1_1pt_s;

end Relation_Table;
