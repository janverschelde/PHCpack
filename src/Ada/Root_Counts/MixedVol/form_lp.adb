with text_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;

-- for testing purposes:
with integer_io; use text_io,integer_io;

package body Form_LP is

  procedure Form_LP
                ( nVar,nSpt : in integer32;
                  SptType,SptIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  RelTab : in Boolean_Matrix;
                  ElmMtrx : in Standard_Floating_Matrices.Link_to_Matrix;
                  NonZero : in out Standard_Integer_Vectors.Link_to_Vector;
                  Lvl2CoDim : in Standard_Integer_Vectors.Link_to_Vector;
                  A : in Standard_Floating_VecMats.Link_to_VecMat;
                  Lvl : in integer32; ptr : in Link_to_LPdata;
                  Cell : in Standard_Integer_Vectors.Link_to_Vector;
                  Cell_Orig : in out Standard_Integer_Vectors.Link_to_Vector;
                  Lvl2LPdim : in Standard_Integer_Vectors.Link_to_Vector;
                  Lvl2Spt : in Standard_Integer_Vectors.Link_to_Vector;
                  FixFrstPt : in Standard_Integer_Vectors.Link_to_Vector;
                  FixLstPt : in Standard_Integer_Vectors.Link_to_Vector;
                  MinNumPt : in Standard_Integer_Vectors.Link_to_Vector;
                  PtIn : in out Standard_Integer_VecVecs.Link_to_VecVec;
                  Strt1Pt,End1Pt : out integer32;
                  ToOrig : in out Standard_Integer_VecVecs.Link_to_VecVec;
                  Pre2Cur : in out Standard_Integer_VecVecs.Link_to_VecVec; 
                  LPdim : out integer32;
                  FrstPt : in out Standard_Integer_Vectors.Link_to_Vector;
                  LstPt : in out Standard_Integer_Vectors.Link_to_Vector;
                  Cur2Pre : in out Standard_Integer_Vectors.Link_to_Vector;
                  x : in out Standard_Floating_Vectors.Link_to_Vector;
                  Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                  Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                  info : out integer32 ) is

    ibrk : boolean;
    ell,ii,kout,EndOldPart,FixPt_Orig,FixPt,nPtIn : integer32;
    smallnum : constant double_float := 1.0E-12;
    DBL_MAX : constant double_float := 1.0E+20;
    dtmp : double_float;
    Elm1Var : Standard_Floating_Vectors.Vector(0..nVar);
    Lvl1 : constant integer32 := Lvl-1;
    ynFixLstPt : constant integer32 := FixLstPt(Lvl);
    ynFixFrstPt : constant integer32 := FixFrstPt(Lvl);
    CurSpt : constant integer32 := Lvl2Spt(Lvl);   -- index of current support
    CoDim : constant integer32 := Lvl2CoDim(Lvl);
    PreLPdim : constant integer32 := Lvl2LPdim(Lvl-1);

  begin
    LPdim := Lvl2LPdim(Lvl);           -- Step 1 : search for LP constraints
    FixPt := Cell(Lvl);                                  -- last point fixed
    if Lvl = 2 then
      FixPt_Orig := FixPt;
      FixPt := Pre2Cur(Lvl1)(FixPt);
    else
      FixPt_Orig := ToOrig(Lvl1)(FixPt);
    end if;
    nPtIn := -1;                  -- #points involved in the constraints - 1
    for i in 0..(FrstPt(Lvl1)-1) loop
      ii := ToOrig(Lvl1)(i); 
      if RelTab(FixPt_Orig,ii) then                             -- always in
        nPtIn := nPtIn + 1;
        ToOrig(Lvl)(nPtIn) := ii;
        Cur2Pre(nPtIn) := i;
        Pre2Cur(Lvl)(i) := nPtIn;
      else
        Pre2Cur(Lvl)(i) := -1;   
      end if;
    end loop;
    FrstPt(Lvl) := nPtIn + 1;
    for i in FrstPt(Lvl1)..(FixPt-1) loop
      ii := ToOrig(Lvl1)(i);
      if (RelTab(FixPt_Orig,ii) and (PtIn(Lvl1)(i) = 1)) then
        nPtIn := nPtIn + 1;
        ToOrig(Lvl)(nPtIn) := ii;
        Cur2Pre(nPtIn) := i;
        Pre2Cur(Lvl)(i) := nPtIn;
      else
        Pre2Cur(Lvl)(i) := -1;
      end if;
    end loop;
    Strt1Pt := nPtIn + 1;
    for i in (FixPt+1)..LstPt(Lvl1) loop
      ii := ToOrig(Lvl1)(i);
      if (RelTab(FixPt_Orig,ii) and (PtIn(Lvl1)(i) = 1)) then
        nPtIn := nPtIn + 1;
        ToOrig(Lvl)(nPtIn) := ii;
        Cur2Pre(nPtIn) := i;
        Pre2Cur(Lvl)(i) := nPtIn;
      else
        Pre2Cur(Lvl)(i) := -1;
      end if;
    end loop;
    if((ynFixLstPt = 0) and nPtIn-Strt1Pt < MinNumPt(Lvl)-1) then
      info := 1; End1pt := Strt1Pt; return;     -- not enough pts to extend
    end if;
    End1Pt := nPtIn;
    EndOldPart := nPtIn;
    LstPt(Lvl) := nPtIn;
    if Lvl = 2
     then Cell_Orig(2) := Cell(2);
     else Cell_Orig(Lvl) := ToOrig(Lvl-1)(Cell(Lvl));
    end if;
    if ynFixLstPt = 1 then            -- the last pt fixed is last pt needed
      Strt1Pt := nPtIn + 1;                            -- extend to next spt
      FrstPt(Lvl) := nPtIn + 1;
      for i in SptIdx(CurSpt+1)..(SptIdx(CurSpt+2)-1) loop
        ibrk := false;
        for j in 1..Lvl loop
          if not RelTab(i,Cell_Orig(j))
           then ibrk := true; exit;
          end if;
        end loop;
        if not ibrk then
          nPtIn := nPtIn + 1;
          ToOrig(Lvl)(nPtIn) := i;
          Cur2Pre(nPtIn) := i;
          Pre2Cur(Lvl)(i) := nPtIn;
        else
          Pre2Cur(Lvl)(i) := -1;
        end if;
      end loop;
      End1Pt := nPtIn;
      LstPt(Lvl) := nPtIn;
      if End1Pt - Strt1Pt < SptType(CurSpt+1)
       then info := 1; return;         -- not enough pts; need irep(j)+1 pts
      end if;
    end if;
   -- Step 2 : eliminate a variable by the eq. of the last fixed point
   --          form the equation used for elimination
    for i in 0..(PreLPdim-1) loop
      Elm1Var(i) := A(Lvl1)(FixPt,i);
    end loop;
    Elm1Var(nVar) := A(Lvl1)(FixPt,nVar);
            -- eliminate the kout-th variable and shift columns of a to left
    if ynFixFrstPt = 0 then                         -- still in same support
      ibrk := false;
      kout := PreLPdim-1;
      for i in reverse 0..PreLPdim-1 loop
        if abs(Elm1Var(i)) > smallnum
         then ibrk := true; kout := i; exit;
        end if;
      end loop;
      if not ibrk then
        text_io.put_line("form_LP : cannot find a variable to eliminate.");
        return;
      end if;
      for i in 0..(kout-1) loop
        Elm1Var(i) := Elm1Var(i)/Elm1Var(kout);
      end loop;
      Elm1Var(nVar) := Elm1Var(nVar)/Elm1Var(kout);
      NonZero(CoDim) := kout;
      for i in 0..(kout-1) loop
        ElmMtrx(CoDim,i) := Elm1Var(i);
      end loop;
      ElmMtrx(CoDim,nVar) := Elm1Var(nVar);
      for i in 0..EndOldPart loop
        ii := Cur2Pre(i);
        if abs(A(Lvl1)(ii,kout)) > smallnum then             -- to eliminate
          for j in 0..(kout-1) loop
            A(Lvl)(i,j) := A(Lvl1)(ii,j) - A(Lvl1)(ii,kout)*Elm1Var(j);
          end loop;
          A(Lvl)(i,nVar) := A(Lvl1)(ii,nVar) - A(Lvl1)(ii,kout)*Elm1Var(nVar);
        else
          for j in 0..(kout-1) loop
            A(Lvl)(i,j) := A(Lvl1)(ii,j);
          end loop;
          A(Lvl)(i,nVar) := A(Lvl1)(ii,nVar);
        end if;
        for j in kout..PreLPdim-2 loop                           -- to shift
          A(Lvl)(i,j) := A(Lvl1)(ii,j+1);
        end loop;
      end loop;
      if ynFixLstPt = 1 then       -- set the values for the variable alpha0
        for i in 0..EndOldPart loop            -- constraints without alpha0
          A(Lvl)(i,PreLPdim-1) := 0.0;
        end loop;
        for i in FrstPt(Lvl)..nPtIn loop          -- constraints with alpha0
          ii := ToOrig(Lvl)(i);
          for j in 0..nVar-1 loop
            A(Lvl)(i,j) := A(0)(ii,j);
          end loop;
          A(Lvl)(i,nVar) := A(0)(ii,nVar);
        end loop;
        for k in 1..CoDim loop
          kout := NonZero(k);
          for j in 0..nVar loop
            Elm1Var(j) := ElmMtrx(k,j);
          end loop;
          for i in FrstPt(Lvl)..nPtIn loop
            if abs(A(Lvl)(i,kout)) > smallnum then            -- to eliminate
              for j in 0..(kout-1) loop
                A(Lvl)(i,j) := A(Lvl)(i,j) - A(Lvl)(i,kout)*Elm1Var(j);
              end loop;
              A(Lvl)(i,nVar) := A(Lvl)(i,nVar) - A(Lvl)(i,kout)*Elm1Var(nVar);
            end if;
            for j in kout..(nVar-k-1) loop                       -- to shift
              A(Lvl)(i,j) := A(Lvl)(i,j+1);
            end loop;
          end loop;
        end loop;               -- kout = NonZero[CoDim] will be used below!
        for i in FrstPt(Lvl)..End1Pt loop         -- constraints with alpha0
          A(Lvl)(i,PreLPdim-1) := 1.0;          -- for added variable alpha0
        end loop;
      end if;
    else          -- fixing the first point from a support, eliminate alpha0
      kout := PreLPdim - 1;
      if abs(Elm1Var(kout)) > smallnum then              -- eliminate alpha0
        for i in 0..(kout-1) loop
          Elm1Var(i) := Elm1Var(i)/Elm1Var(kout);
        end loop;
        Elm1Var(nVar) := Elm1Var(nVar)/Elm1Var(kout);
        for i in 0..(FrstPt(Lvl)-1) loop     -- copy the previous a, b saved
          ii := Cur2Pre(i);
          for j in 0..(PreLPdim-2) loop
            A(Lvl)(i,j) := A(Lvl1)(ii,j);
          end loop;
          A(Lvl)(i,nVar) := A(Lvl1)(ii,nVar);
        end loop;
        for i in FrstPt(Lvl)..End1Pt loop                -- eliminate alpha0
          ii := Cur2Pre(i);
          for j in 0..(PreLPdim-2) loop
            A(Lvl)(i,j) := A(Lvl1)(ii,j) - Elm1Var(j);
          end loop;
          A(Lvl)(i,nVar) := A(Lvl1)(ii,nVar) - Elm1Var(nVar);
        end loop;
      end if;
    end if;
    ibrk := false;
    if Lvl = 2 then
      for i in 0..(PreLPdim-1) loop
        Bidx(i) := ptr.JJJ(i);
        if Bidx(i) = FixPt_Orig
         then ii := i; ibrk := true; exit;
        end if;
      end loop;
    else
      for i in 0..(PreLPdim-1) loop
        Bidx(i) := ptr.JJJ(i);
        if Bidx(i) = FixPt
         then ii := i; ibrk := true; exit;
        end if;
      end loop;
    end if;
    if not ibrk then
      text_io.put_line("form_LP : no index match for reused info"); return;
    end if;
  --  put("PreLPdim = "); put(PreLPdim,1); new_line;
  --  put("  Bidx'first = "); put(Bidx'first,1);
  --  put("  Bidx'last = "); put(Bidx'last,1); new_line;
  --  put("  ptr.JJJ'first = "); put(ptr.JJJ'first,1);
  --  put("  ptr.JJJ'last = "); put(ptr.JJJ'last,1); new_line;
    for i in (ii+1)..(PreLPdim-1) loop
     -- put("i = "); put(i,1); new_line;
      Bidx(i-1) := ptr.JJJ(i);
    end loop;
    if Lvl > 2 then
      for i in 0..(PreLPdim-2) loop
        if Bidx(i) >-1 
         then Bidx(i) := Pre2Cur(Lvl)(Bidx(i));    -- may = -1 if eliminated
        end if;
      end loop;
    else 
      for i in 0..(PreLPdim-2) loop
        if Bidx(i) > -1 then
          Bidx(i) := Pre2Cur(1)(Bidx(i));         -- from level 0 to level 1
          if Bidx(i) > -1 
           then Bidx(i) := Pre2Cur(Lvl)(Bidx(i));      -- level 1 to level 2
          end if;
        end if;
      end loop;
    end if;
    for i in 0..(kout-1) loop
      x(i) := ptr.xxx(i);
    end loop;
    for i in (kout+1)..(PreLPdim-1) loop
      x(i-1) := ptr.xxx(i);
    end loop;
    for j in 0..(ii-1) loop
      for i in 0..(kout-1) loop
        Binv(j,i) := ptr.INV(j,i);
      end loop;
      for i in (kout+1)..(PreLPdim-1) loop
        Binv(j,i-1) := ptr.INV(j,i);
      end loop;
    end loop;
    for j in (ii+1)..(PreLPdim-1) loop
      for i in 0..(kout-1) loop
        Binv(j-1,i) := ptr.INV(j,i); 
      end loop;
      for i in (kout+1)..(PreLPdim-1) loop
        Binv(j-1,i-1) := ptr.INV(j,i);
      end loop;
    end loop;
    if ynFixLstPt = 1 then       -- the 1st round 1-pt test for next support
      x(LPdim-1) := DBL_MAX;                         --  the variable alpha0
      ell := -1;
      for i in Strt1Pt..End1Pt loop
        dtmp := A(Lvl)(i,nVar);
        for j in 0..(LPdim-2) loop
          dtmp := dtmp - A(Lvl)(i,j)*x(j);
        end loop;
        if x(LPdim-1) > dtmp then
          ell := i;
          x(LPdim-1) := dtmp;
        end if;
      end loop;
      if ell < 0
       then text_io.put_line("FormLP : no ell"); return;
      end if;
      Bidx(LPdim-1) := ell;       -- constraint becomes last element of base
      for i in 0..(LPdim-1) loop                -- update the inverse matrix
        Binv(LPdim-1,i) := 0.0;
      end loop;
      for j in 0..(LPdim-2) loop
        dtmp := A(Lvl)(ell,0)*Binv(j,0);
        for i in 1..LPdim-2 loop
          dtmp := dtmp + A(Lvl)(ell,i)*Binv(j,i);
        end loop;
        Binv(j,LPdim-1) := -dtmp;
      end loop;
      Binv(LPdim-1,LPdim-1) := 1.0;
    end if;
    for i in FrstPt(Lvl)..(Strt1Pt-1) loop
      PtIn(Lvl)(i) := 1;
    end loop;
    for i in Strt1Pt..End1Pt loop
      PtIn(Lvl)(i) := 0;
    end loop;
    info := 0;
  end Form_LP;

  procedure Form_LP1
                ( nVar,nSpt : in integer32;
                  SptType,SptIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  RelTab : in Boolean_Matrix;
                  A : in Standard_Floating_VecMats.Link_to_VecMat;
                  Lvl : in integer32;
                  Cell : in Standard_Integer_Vectors.Link_to_Vector;
                  Lvl2LPdim : in Standard_Integer_Vectors.Link_to_Vector;
                  FixLstPt : in Standard_Integer_Vectors.Link_to_Vector;
                  MinNumPt : in Standard_Integer_Vectors.Link_to_Vector;
                  PtIn : in out Standard_Integer_VecVecs.Link_to_VecVec;
                  ToOrig : in out Standard_Integer_VecVecs.Link_to_VecVec;
                  Pre2Cur : in out Standard_Integer_VecVecs.Link_to_VecVec;
                  FrstPt : in out Standard_Integer_Vectors.Link_to_Vector;
                  LstPt : in out Standard_Integer_Vectors.Link_to_Vector;
                  info : out integer32 ) is

    ind,FixPt,nPtIn,Strt1Pt : integer32;
    Lvl1 : constant integer32 := Lvl - 1;
    ynFixLstPt : constant integer32 := FixLstPt(Lvl);
    PreLPdim : constant integer32 := Lvl2LPdim(Lvl-1);

  begin
    FixPt := Cell(Lvl);                                  -- last point fixed
    nPtIn := -1;             -- #points involved to form the constraints - 1
    FrstPt(Lvl) := 0;
    for i in SptIdx(0)..(FixPt-1) loop
      if RelTab(FixPt,i) then
        nPtIn := nPtIn + 1;
        ToOrig(Lvl)(nPtIn) := i;
        Pre2Cur(Lvl)(i) := nPtIn;        -- need for saved Bidx for next lvl
      else
        Pre2Cur(Lvl)(i) := -1;
      end if;
    end loop;
    Pre2Cur(Lvl)(FixPt) := -1;
    Strt1Pt := nPtIn + 1;
    for i in (FixPt+1)..(SptIdx(1)-1) loop
      if RelTab(FixPt,i) then
        nPtIn := nPtIn + 1;
        ToOrig(Lvl)(nPtIn) := i;
        Pre2Cur(Lvl)(i) := nPtIn;        -- need for saved Bidx for next lvl
      else
        Pre2Cur(Lvl)(i) := -1;
      end if;
    end loop;
    if((ynFixLstPt = 0) and (nPtIn-Strt1Pt < MinNumPt(Lvl)-1)) then
      info := 1; return;                         -- not enough pts to extend
    end if;
    LstPt(Lvl) := nPtIn;                               -- to form the matrix
    for i in 0..LstPt(Lvl) loop              -- copy the previous a, b saved
      ind := ToOrig(Lvl)(i);
      for j in 0..(PreLPdim-1) loop
        A(Lvl)(i,j) := A(Lvl1)(ind,j) - A(Lvl1)(FixPt,j);
      end loop;
      A(Lvl)(i,nVar) := A(Lvl1)(ind,nVar) - A(Lvl1)(FixPt,nVar);
    end loop;
    for i in 0..LstPt(Lvl) loop
      PtIn(Lvl)(i) := 1;
    end loop;
    info := 0;
  end Form_LP1;
  
end Form_LP;
