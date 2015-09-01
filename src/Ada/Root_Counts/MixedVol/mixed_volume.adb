with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Multprec_Integer_Numbers;
with Standard_Common_Divisors;
with Standard_Integer_Matrices;
with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;
with Multprec_Integer_Matrices;
with Multprec_Integer_Linear_Solvers;    use Multprec_Integer_Linear_Solvers;
with Standard_Floating_VecMats;
with Relation_Table;                     use Relation_Table;
with Zero_Index_Tree;                    use Zero_Index_Tree;
with Index_Tree_LP;                      use Index_Tree_LP;
with One_Level_LP,Form_LP;

-- for testing :
--with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Multprec_Integer_Matrices_io;       use Multprec_Integer_Matrices_io;

package body Mixed_Volume is

--  procedure Write_Relation_Table ( t : in Boolean_Matrix ) is
--
--  -- DESCRIPTION :
--  --   Writes the matrix as 1 for true and 0 for false.
--
--  begin
--    for i in t'range(1) loop
--      for j in t'range(2) loop
--        if t(i,j)
--         then put(" 1");
--         else put(" 0");
--        end if;
--      end loop;
--      new_line;
--    end loop;
--  end Write_Relation_Table;

  procedure MixedVol 
               ( nVar,nSpt,CellSize : in integer32;
                 SptType,SptIdx : in Standard_Integer_Vectors.Link_to_Vector;
                 Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : in Standard_Floating_Vectors.Link_to_Vector;
                 nbCells : out integer32; MCells : out CellStack;
                 MVol : out natural32;
                 multprec_hermite : in boolean := false ) is
  
    k,info,Strt1Pt,End1Pt,LPdim,Lvl : integer32;
    labels : Standard_Integer_Vectors.Link_to_Vector
           := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    NonZero : Standard_Integer_Vectors.Link_to_Vector
            := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    Lvl2CoDim : Standard_Integer_Vectors.Link_to_Vector
              := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    Cell : Standard_Integer_Vectors.Link_to_Vector
         := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    Cell_Orig : Standard_Integer_Vectors.Link_to_Vector
               := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    Lvl2LPdim : Standard_Integer_Vectors.Link_to_Vector
              := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    Lvl2Spt : Standard_Integer_Vectors.Link_to_Vector
            := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    MinNumPt : Standard_Integer_Vectors.Link_to_Vector
             := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    FixFrstPt : Standard_Integer_Vectors.Link_to_Vector
              := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    FixLstPt : Standard_Integer_Vectors.Link_to_Vector
             := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    FrstPt : Standard_Integer_Vectors.Link_to_Vector
           := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    LstPt : Standard_Integer_Vectors.Link_to_Vector
          := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    Cur2Pre : Standard_Integer_Vectors.Link_to_Vector
            := new Standard_Integer_Vectors.Vector(0..SptIdx(nSpt)-1);
    Bidx : Standard_Integer_Vectors.Link_to_Vector
         := new Standard_Integer_Vectors.Vector(0..nVar-1);
    PtIn : Standard_Integer_VecVecs.Link_to_VecVec
         := new Standard_Integer_VecVecs.VecVec(0..CellSize-1);
    ToOrig : Standard_Integer_VecVecs.Link_to_VecVec
           := new Standard_Integer_VecVecs.VecVec(0..CellSize-1);
    Pre2Cur : Standard_Integer_VecVecs.Link_to_VecVec
            := new Standard_Integer_VecVecs.VecVec(0..CellSize-1);
    RelTab : Boolean_Matrix(0..(SptIdx(nSpt)-1),0..(SptIdx(nSpt)-1));
    x : Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector(0..nVar-1);
    ElmMtrx : Standard_Floating_Matrices.Link_to_Matrix
            := new Standard_Floating_Matrices.Matrix(0..CellSize-1,0..nVar);
    Binv : Standard_Floating_Matrices.Link_to_Matrix
         := new Standard_Floating_Matrices.Matrix(0..nVar-1,0..nVar-1);
    A : Standard_Floating_VecMats.Link_to_VecMat
      := new Standard_Floating_VecMats.VecMat(0..CellSize-1);
    ItLp : Link_to_IT_LP := new IT_LP;
    L0 : Link_to_L0_IML := new L0_IML;
    ptr : Link_to_LPdata;
    itcell : Standard_Integer_Vectors.Link_to_Vector;

  begin
    nbCells := 0; MVol := 0;
    IT_LP_Init(ItLp,nSpt,SptType);
    L0_IML_Init(L0);   
    for i in 0..(CellSize-1) loop
      PtIn(i) := new Standard_Integer_Vectors.Vector(0..SptIdx(nSpt)-1);
      ToOrig(i) := new Standard_Integer_Vectors.Vector(0..SptIdx(nSpt)-1);
      Pre2Cur(i) := new Standard_Integer_Vectors.Vector(0..SptIdx(nSpt)-1);
    end loop;
   -- put("creating relation table ... ");
    RelTable(nVar,nSpt,Spt,SptIdx,lft,RelTab,L0);
   -- put_line("done, the relation table :");
   -- Write_Relation_Table(RelTab);
    info := SptIdx(nSpt);                                      -- #rows in A
    for i in 0..(nSpt-2) loop
      for j in 0..(SptType(i)-1) loop
        info := info + SptIdx(i+1);
      end loop;
      info := info + SptIdx(i+2);
    end loop;
    for j in 0..(SptType(nSpt-1)-1) loop
      info := info + SptIdx(nSpt);
    end loop;
    for i in 0..(CellSize-1) loop
      A(i) := new Standard_Floating_Matrices.Matrix(0..info-1,0..nVar);
    end loop;
    info := 0;             -- Note that the level 0 is the artificial head
    Lvl2LPdim(0) := nVar + 1;  -- for conveniece, should = nVar, see below
    Lvl2CoDim(0) := 0;              -- #variables eliminated. [0] not used
    k := 0;
    for i in 0..(nSpt-1) loop
      k := k + 1;
      Lvl2LPdim(k) := Lvl2LPdim(k-1) - 1;
      Lvl2CoDim(k) := Lvl2CoDim(k-1);
      Lvl2Spt(k) := i;
      MinNumPt(k) := SptType(i);
      FixFrstPt(k) := 1;
      FixLstPt(k) := 0;
      for j in 1..SptType(i) loop
        k := k + 1;
        exit when (k = CellSize);
        Lvl2LPdim(k) := Lvl2LPdim(k-1) - 1;
        Lvl2CoDim(k) := Lvl2CoDim(k-1) + 1;
        Lvl2Spt(k) := i;
        MinNumPt(k) := MinNumPt(k-1) - 1;
        FixFrstPt(k) := 0;
        FixLstPt(k) := 0;
      end loop;
      exit when (k = CellSize);
      Lvl2LPdim(k) := Lvl2LPdim(k) + 1;         -- add alpha0 for next spt
      if i < nSpt-1
       then MinNumPt(k) := SptType(i+1) + 1;
      end if;
      FixLstPt(k) := 1;
    end loop;
    Lvl2LPdim(0) := nVar; -- =nVar+1-1, in RelTab, add alpha0, fix 1 point
   -- put("Lvl2LPdim = "); put(Lvl2LPdim); new_line;
    for i in SptIdx(0)..(SptIdx(nSpt)-1) loop               -- define A[0]
      for j in 0..(nVar-1) loop
        A(0)(i,j) := -double_float(Spt(i)(j));
      end loop;
      A(0)(i,nVar) := lft(i);                               -- bIdx = nVar
    end loop;
    for i in SptIdx(0)..(SptIdx(1)-1) loop
      PtIn(0)(i) := 1;                                   -- define PtIn[0]
    end loop;
   -- put_line("entering the main loop ...");
    loop               --  while(L0_Migrate(L0,IT_ResetCurLevelTo1(ItLp)))
      declare
        lnd : Link_to_IndexNode;
        status : integer32;
      begin
        IT_ResetCurLevelTo1(ItLp,lnd);
        L0_Migrate(L0,lnd,status);
       -- put("status = "); put(status,1); new_line;
        exit when (status = 0);
      end;
      IT_RenewNP1(ItLp);
     -- put_line("entering the while loop ...");
      while not IT_IsEmpty(ItLp) loop
        if not IT_IsFull(ItLp) then
          Lvl := IT_Level(ItLp);         -- level-0 is the artificial head
          Cell(Lvl) := IT_FixedIdxNdPtr(ItLp).idx; 
                                    --  The index of the point to be fixed
          ptr := IT_FixedIdxNdPtr(ItLp).info;
                                       --  The ptr to saved x, Binv and jj
          if Lvl <= 1 then        -- eliminate and form new system Ax <= b
            Cell_Orig(Lvl) := Cell(Lvl);
           -- put_line("  executing Form_LP1");
            Form_LP.Form_LP1
              (nVar,nSpt,SptType,SptIdx,RelTab,A,Lvl,Cell,Lvl2LPdim,
               FixLstPt,MinNumPt,PtIn,ToOrig,Pre2Cur,FrstPt,LstPt,info);
           -- put("    Form_LP1 returns info = "); put(info,1); new_line;
          else                -- eliminate, form new Ax <= b and solve the LP
           -- put_line("  executing Form_LP");
            Form_LP.Form_LP
              (nVar,nSpt,SptType,SptIdx,RelTab,ElmMtrx,NonZero,Lvl2CoDim,A,
               Lvl,ptr,Cell,Cell_Orig,Lvl2LPdim,Lvl2Spt,FixFrstPt,FixLstPt,
               MinNumPt,PtIn,Strt1Pt,End1Pt,ToOrig,Pre2Cur,LPdim,FrstPt,LstPt,
               Cur2Pre,x,Binv,Bidx,info);
            if info = 0 then
              One_Level_LP.one_level_LP
                (Strt1Pt,End1Pt,PtIn(Lvl),LPdim,A(Lvl),nVar,x,Binv,Bidx,ItLp);
            end if;
          end if;
        end if;                                      -- This level finished
       -- put_line("entering the nextlevel loop ...");
        loop -- while(not IT_NextLevel(ItLp) and not IT_IsEmpty(ItLp)) loop
         -- put("  level is "); put(IT_Level(ItLp),1); new_line;
          declare
            rcase : integer32;
          begin
            IT_NextLevel(ItLp,rcase);
           -- put("  IT_NextLevel, level is "); put(IT_Level(ItLp),1); new_line;
           -- put("                rcase is "); put(rcase,1); new_line;
           -- if IT_IsEmpty(ItLp)
           --  then put_line("     ItLp is empty");
           --  else put_line("     ItLp is not empty");
           -- end if;
            exit when ((rcase = 1) or IT_IsEmpty(ItLp));
          end;
          if IT_IsFull(ItLp) then
            itcell := IT_Cell(ItLp);
            labels(0) := itcell(1);
            labels(1) := itcell(2);
            for i in 2..(CellSize-1) loop
              labels(i) := ToOrig(i)(itcell(i+1));
            end loop;
           -- put_line("  pushing a mixed cell to the stack");
            Cs_Push(MCells,labels);                       -- store the cell
           -- important: crash on uliftmv, leave out for tropisms
            CellVol(nVar,nSpt,Spt,SptType,labels,MVol,multprec_hermite);
            nbCells := nbCells + 1; 
          end if;
         -- put_line("  executing IT_StepBack");
          IT_StepBack(ItLp);
         -- put_line("  leaving the next level loop");
        end loop; 
      end loop;                       -- End of while ( !IT_IsEmpty(ItLp) )
    end loop;                           -- End of while ( L0.Migrate(inp) )
    Standard_Integer_Vectors.Clear(labels);
    Standard_Integer_Vectors.Clear(NonZero);
    Standard_Integer_Vectors.Clear(Lvl2CoDim);
    Standard_Integer_Vectors.Clear(Cell);
    Standard_Integer_Vectors.Clear(Cell_Orig);
    Standard_Integer_Vectors.Clear(Lvl2LPdim);
    Standard_Integer_Vectors.Clear(Lvl2Spt);
    Standard_Integer_Vectors.Clear(MinNumPt);
    Standard_Integer_Vectors.Clear(FixFrstPt);
    Standard_Integer_Vectors.Clear(FixLstPt);
    Standard_Integer_Vectors.Clear(FrstPt);
    Standard_Integer_Vectors.Clear(LstPt);
    Standard_Integer_Vectors.Clear(Cur2Pre);
    Standard_Integer_Vectors.Clear(Bidx);
    Standard_Integer_VecVecs.Deep_Clear(PtIn);
    Standard_Integer_VecVecs.Deep_Clear(ToOrig);
    Standard_Integer_VecVecs.Deep_Clear(Pre2Cur);
    Standard_Floating_Vectors.Clear(x);
    Standard_Floating_Matrices.Clear(Binv);
    Standard_Floating_Matrices.Clear(ElmMtrx);
    Standard_Floating_VecMats.Deep_Clear(A);
  end MixedVol;

  procedure MixedVol_with_Callback
               ( nVar,nSpt,CellSize : in integer32;
                 SptType,SptIdx : in Standard_Integer_Vectors.Link_to_Vector;
                 Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : in Standard_Floating_Vectors.Link_to_Vector;
                 nbCells : out integer32; MCells : out CellStack;
                 MVol : out natural32;
                 multprec_hermite : in boolean := false;
                 next_cell : access procedure
                   ( idx : Standard_Integer_Vectors.Link_to_Vector )
                   := null ) is
  
    k,info,Strt1Pt,End1Pt,LPdim,Lvl : integer32;
    labels : Standard_Integer_Vectors.Link_to_Vector
           := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    NonZero : Standard_Integer_Vectors.Link_to_Vector
            := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    Lvl2CoDim : Standard_Integer_Vectors.Link_to_Vector
              := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    Cell : Standard_Integer_Vectors.Link_to_Vector
         := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    Cell_Orig : Standard_Integer_Vectors.Link_to_Vector
               := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    Lvl2LPdim : Standard_Integer_Vectors.Link_to_Vector
              := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    Lvl2Spt : Standard_Integer_Vectors.Link_to_Vector
            := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    MinNumPt : Standard_Integer_Vectors.Link_to_Vector
             := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    FixFrstPt : Standard_Integer_Vectors.Link_to_Vector
              := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    FixLstPt : Standard_Integer_Vectors.Link_to_Vector
             := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    FrstPt : Standard_Integer_Vectors.Link_to_Vector
           := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    LstPt : Standard_Integer_Vectors.Link_to_Vector
          := new Standard_Integer_Vectors.Vector(0..CellSize-1);
    Cur2Pre : Standard_Integer_Vectors.Link_to_Vector
            := new Standard_Integer_Vectors.Vector(0..SptIdx(nSpt)-1);
    Bidx : Standard_Integer_Vectors.Link_to_Vector
         := new Standard_Integer_Vectors.Vector(0..nVar-1);
    PtIn : Standard_Integer_VecVecs.Link_to_VecVec
         := new Standard_Integer_VecVecs.VecVec(0..CellSize-1);
    ToOrig : Standard_Integer_VecVecs.Link_to_VecVec
           := new Standard_Integer_VecVecs.VecVec(0..CellSize-1);
    Pre2Cur : Standard_Integer_VecVecs.Link_to_VecVec
            := new Standard_Integer_VecVecs.VecVec(0..CellSize-1);
    RelTab : Boolean_Matrix(0..(SptIdx(nSpt)-1),0..(SptIdx(nSpt)-1));
    x : Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector(0..nVar-1);
    ElmMtrx : Standard_Floating_Matrices.Link_to_Matrix
            := new Standard_Floating_Matrices.Matrix(0..CellSize-1,0..nVar);
    Binv : Standard_Floating_Matrices.Link_to_Matrix
         := new Standard_Floating_Matrices.Matrix(0..nVar-1,0..nVar-1);
    A : Standard_Floating_VecMats.Link_to_VecMat
      := new Standard_Floating_VecMats.VecMat(0..CellSize-1);
    ItLp : Link_to_IT_LP := new IT_LP;
    L0 : Link_to_L0_IML := new L0_IML;
    ptr : Link_to_LPdata;
    itcell : Standard_Integer_Vectors.Link_to_Vector;

  begin
    nbCells := 0; MVol := 0;
    IT_LP_Init(ItLp,nSpt,SptType);
    L0_IML_Init(L0);   
    for i in 0..(CellSize-1) loop
      PtIn(i) := new Standard_Integer_Vectors.Vector(0..SptIdx(nSpt)-1);
      ToOrig(i) := new Standard_Integer_Vectors.Vector(0..SptIdx(nSpt)-1);
      Pre2Cur(i) := new Standard_Integer_Vectors.Vector(0..SptIdx(nSpt)-1);
    end loop;
   -- put("creating relation table ... ");
    RelTable(nVar,nSpt,Spt,SptIdx,lft,RelTab,L0);
   -- put_line("done, the relation table :");
   -- Write_Relation_Table(RelTab);
    info := SptIdx(nSpt);                                      -- #rows in A
    for i in 0..(nSpt-2) loop
      for j in 0..(SptType(i)-1) loop
        info := info + SptIdx(i+1);
      end loop;
      info := info + SptIdx(i+2);
    end loop;
    for j in 0..(SptType(nSpt-1)-1) loop
      info := info + SptIdx(nSpt);
    end loop;
    for i in 0..(CellSize-1) loop
      A(i) := new Standard_Floating_Matrices.Matrix(0..info-1,0..nVar);
    end loop;
    info := 0;             -- Note that the level 0 is the artificial head
    Lvl2LPdim(0) := nVar + 1;  -- for conveniece, should = nVar, see below
    Lvl2CoDim(0) := 0;              -- #variables eliminated. [0] not used
    k := 0;
    for i in 0..(nSpt-1) loop
      k := k + 1;
      Lvl2LPdim(k) := Lvl2LPdim(k-1) - 1;
      Lvl2CoDim(k) := Lvl2CoDim(k-1);
      Lvl2Spt(k) := i;
      MinNumPt(k) := SptType(i);
      FixFrstPt(k) := 1;
      FixLstPt(k) := 0;
      for j in 1..SptType(i) loop
        k := k + 1;
        exit when (k = CellSize);
        Lvl2LPdim(k) := Lvl2LPdim(k-1) - 1;
        Lvl2CoDim(k) := Lvl2CoDim(k-1) + 1;
        Lvl2Spt(k) := i;
        MinNumPt(k) := MinNumPt(k-1) - 1;
        FixFrstPt(k) := 0;
        FixLstPt(k) := 0;
      end loop;
      exit when (k = CellSize);
      Lvl2LPdim(k) := Lvl2LPdim(k) + 1;         -- add alpha0 for next spt
      if i < nSpt-1
       then MinNumPt(k) := SptType(i+1) + 1;
      end if;
      FixLstPt(k) := 1;
    end loop;
    Lvl2LPdim(0) := nVar; -- =nVar+1-1, in RelTab, add alpha0, fix 1 point
   -- put("Lvl2LPdim = "); put(Lvl2LPdim); new_line;
    for i in SptIdx(0)..(SptIdx(nSpt)-1) loop               -- define A[0]
      for j in 0..(nVar-1) loop
        A(0)(i,j) := -double_float(Spt(i)(j));
      end loop;
      A(0)(i,nVar) := lft(i);                               -- bIdx = nVar
    end loop;
    for i in SptIdx(0)..(SptIdx(1)-1) loop
      PtIn(0)(i) := 1;                                   -- define PtIn[0]
    end loop;
   -- put_line("entering the main loop ...");
    loop               --  while(L0_Migrate(L0,IT_ResetCurLevelTo1(ItLp)))
      declare
        lnd : Link_to_IndexNode;
        status : integer32;
      begin
        IT_ResetCurLevelTo1(ItLp,lnd);
        L0_Migrate(L0,lnd,status);
       -- put("status = "); put(status,1); new_line;
        exit when (status = 0);
      end;
      IT_RenewNP1(ItLp);
     -- put_line("entering the while loop ...");
      while not IT_IsEmpty(ItLp) loop
        if not IT_IsFull(ItLp) then
          Lvl := IT_Level(ItLp);         -- level-0 is the artificial head
          Cell(Lvl) := IT_FixedIdxNdPtr(ItLp).idx; 
                                    --  The index of the point to be fixed
          ptr := IT_FixedIdxNdPtr(ItLp).info;
                                       --  The ptr to saved x, Binv and jj
          if Lvl <= 1 then        -- eliminate and form new system Ax <= b
            Cell_Orig(Lvl) := Cell(Lvl);
           -- put_line("  executing Form_LP1");
            Form_LP.Form_LP1
              (nVar,nSpt,SptType,SptIdx,RelTab,A,Lvl,Cell,Lvl2LPdim,
               FixLstPt,MinNumPt,PtIn,ToOrig,Pre2Cur,FrstPt,LstPt,info);
           -- put("    Form_LP1 returns info = "); put(info,1); new_line;
          else                -- eliminate, form new Ax <= b and solve the LP
           -- put_line("  executing Form_LP");
            Form_LP.Form_LP
              (nVar,nSpt,SptType,SptIdx,RelTab,ElmMtrx,NonZero,Lvl2CoDim,A,
               Lvl,ptr,Cell,Cell_Orig,Lvl2LPdim,Lvl2Spt,FixFrstPt,FixLstPt,
               MinNumPt,PtIn,Strt1Pt,End1Pt,ToOrig,Pre2Cur,LPdim,FrstPt,LstPt,
               Cur2Pre,x,Binv,Bidx,info);
            if info = 0 then
              One_Level_LP.one_level_LP
                (Strt1Pt,End1Pt,PtIn(Lvl),LPdim,A(Lvl),nVar,x,Binv,Bidx,ItLp);
            end if;
          end if;
        end if;                                      -- This level finished
       -- put_line("entering the nextlevel loop ...");
        loop -- while(not IT_NextLevel(ItLp) and not IT_IsEmpty(ItLp)) loop
         -- put("  level is "); put(IT_Level(ItLp),1); new_line;
          declare
            rcase : integer32;
          begin
            IT_NextLevel(ItLp,rcase);
           -- put("  IT_NextLevel, level is "); put(IT_Level(ItLp),1); new_line;
           -- put("                rcase is "); put(rcase,1); new_line;
           -- if IT_IsEmpty(ItLp)
           --  then put_line("     ItLp is empty");
           --  else put_line("     ItLp is not empty");
           -- end if;
            exit when ((rcase = 1) or IT_IsEmpty(ItLp));
          end;
          if IT_IsFull(ItLp) then
            itcell := IT_Cell(ItLp);
            labels(0) := itcell(1);
            labels(1) := itcell(2);
            for i in 2..(CellSize-1) loop
              labels(i) := ToOrig(i)(itcell(i+1));
            end loop;
           -- put_line("  pushing a mixed cell to the stack");
            Cs_Push(MCells,labels);                       -- store the cell
            next_cell(labels);
           -- important: crash on uliftmv, leave out for tropisms
            CellVol(nVar,nSpt,Spt,SptType,labels,MVol,multprec_hermite);
            nbCells := nbCells + 1; 
          end if;
         -- put_line("  executing IT_StepBack");
          IT_StepBack(ItLp);
         -- put_line("  leaving the next level loop");
        end loop; 
      end loop;                       -- End of while ( !IT_IsEmpty(ItLp) )
    end loop;                           -- End of while ( L0.Migrate(inp) )
    Standard_Integer_Vectors.Clear(labels);
    Standard_Integer_Vectors.Clear(NonZero);
    Standard_Integer_Vectors.Clear(Lvl2CoDim);
    Standard_Integer_Vectors.Clear(Cell);
    Standard_Integer_Vectors.Clear(Cell_Orig);
    Standard_Integer_Vectors.Clear(Lvl2LPdim);
    Standard_Integer_Vectors.Clear(Lvl2Spt);
    Standard_Integer_Vectors.Clear(MinNumPt);
    Standard_Integer_Vectors.Clear(FixFrstPt);
    Standard_Integer_Vectors.Clear(FixLstPt);
    Standard_Integer_Vectors.Clear(FrstPt);
    Standard_Integer_Vectors.Clear(LstPt);
    Standard_Integer_Vectors.Clear(Cur2Pre);
    Standard_Integer_Vectors.Clear(Bidx);
    Standard_Integer_VecVecs.Deep_Clear(PtIn);
    Standard_Integer_VecVecs.Deep_Clear(ToOrig);
    Standard_Integer_VecVecs.Deep_Clear(Pre2Cur);
    Standard_Floating_Vectors.Clear(x);
    Standard_Floating_Matrices.Clear(Binv);
    Standard_Floating_Matrices.Clear(ElmMtrx);
    Standard_Floating_VecMats.Deep_Clear(A);
  end MixedVol_with_Callback;

  procedure gcd ( r1,r2 : in integer32; k,l,d : out integer32 ) is

    wr1,wr2,tempxx,tempyy,xx,yy,q,r3 : integer32;

  begin
    k := 0; l := 1; xx := 1; yy := 0; 
    wr1 := r1; wr2 := r2;
    q := wr1/wr2; r3 := wr1 mod wr2;
    while r3 /= 0 loop
      tempxx := xx - q*k;
      tempyy := yy - q*l;
      xx := k; yy := l;
      k := tempxx; l := tempyy;
      wr1 := wr2; wr2 := r3;
      q := wr1/wr2; r3 := wr1 mod wr2;
    end loop;
    if wr2 > 0 then
      d := wr2;
    else
      d := -wr2;
      k := -k; l := -l;
    end if;
  end gcd;

  function cell_size
               ( nSpt : integer32;
                 SptType : Standard_Integer_Vectors.Link_to_Vector )
               return integer32 is

    CellSize : integer32 := nSpt;

  begin
    for i in 0..(nSpt-1) loop
      CellSize := CellSize + SptType(i);
    end loop;
    return CellSize;
  end cell_size;

  procedure Multprec_CellVol
               ( nVar,nSpt : in integer32;
                 Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                 SptType,Cell : in Standard_Integer_Vectors.Link_to_Vector;
                 Vol : in out natural32 ) is

  -- NOTE : this version uses a multiprecision Hermite normal form
  --   to compute the determinant.

    ell,m,n,d,tmp : integer32;
    iwork : Multprec_Integer_Matrices.Matrix(1..nVar,1..nVar);
    mp_d : Multprec_Integer_Numbers.Integer_Number;

  begin
   -- put("calling Multprec_CellVol, vol = "); put(Vol,1);
    d := -1;
    n := 0;
    for i in 0..(nSpt-1) loop
      d := d + 1;
      ell := Cell(d);
      for j in 0..(SptType(i)-1) loop
        d := d + 1;
        m := Cell(d);
        n := n + 1;
        for k in 0..(nVar-1) loop
          tmp := Spt(m)(k) - Spt(ell)(k);
          iwork(n,k+1) := Multprec_Integer_Numbers.Create(tmp);
        end loop;
      end loop;
    end loop;
    mp_d := Det(iwork);
    d := Multprec_Integer_Numbers.Create(mp_d);
    if d < 0
     then Vol := Vol + natural32(-d); -- put(" + "); put(-d,1);
     else Vol := Vol + natural32(d); -- put(" + "); put(d,1); 
    end if;
   -- put(" = "); put(Vol,1); new_line;
  exception
    when others => 
       put_line("exception raised in mixed_volume.Multprec_CellVol");
       put_line("iwork is "); put(iwork); 
       Vol := 0; return;
  end Multprec_CellVol;

  procedure Standard_CellVol
               ( nVar,nSpt : in integer32;
                 Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                 SptType,Cell : in Standard_Integer_Vectors.Link_to_Vector;
                 Vol : in out natural32 ) is

  -- NOTE : This version of CellVol does not use multiprecision
  --   arithmetic and may crash or give wrong results...

    jj,i1,kk,ell,m,n,d,tmp : integer32;
    iwork : Standard_Integer_Matrices.Matrix(0..nVar-1,0..nVar-1);

  begin
   -- put_line("calling CellVol, vol = "); put(Vol,1); new_line;
    d := -1;
    n := -1;
    for i in 0..(nSpt-1) loop
      d := d + 1;
      ell := Cell(d);
      for j in 0..(SptType(i)-1) loop
        d := d + 1;
        m := Cell(d);
        n := n + 1;
        for k in 0..(nVar-1) loop
          iwork(n,k) := Spt(m)(k) - Spt(ell)(k);
        end loop;
      end loop;
    end loop;
    for i in 0..(nVar-2) loop
      jj := i;
      while ((jj<nVar) and (iwork(jj,i) = 0)) loop
        jj := jj + 1;
      end loop;
      if jj = nVar then
        text_io.put_line("Volume of a mixed cell is zero!");
      else
        if jj /= i then
          for w in 0..(nVar-1) loop
            tmp := iwork(i,w);
            iwork(i,w) := iwork(jj,w);
            iwork(jj,w) := tmp;
          end loop;
        end if;
        for j in (i+1)..(nVar-1) loop
          if iwork(j,i) /= 0 then
           -- gcd(iwork(i,i),iwork(j,i),kk,ell,d);
            Standard_Common_Divisors.gcd(iwork(i,i),iwork(j,i),kk,ell,d);
           -- IMPORTANT: the gcd in Standard_Common_Divisors causes a crash
            m := -(iwork(j,i)/d);
            n := iwork(i,i)/d;
            for j1 in i..(nVar-1) loop
              i1 := iwork(i,j1);
              iwork(i,j1) := kk*i1 + ell*iwork(j,j1);
              iwork(j,j1) := m*i1 + n*iwork(j,j1);
            end loop;
          end if;
        end loop;
      end if;
    end loop; 
    d := iwork(nVar-1,nVar-1);
    for i in 0..(nVar-2) loop
      d := d*iwork(i,i);
    end loop;
    if d < 0
     then Vol := Vol + natural32(-d);
     else Vol := Vol + natural32(d);  
    end if;
  exception
    when others => 
       put_line("exception raised in mixed_volume.Standard_CellVol");
       put_line("iwork is "); put(iwork); 
       put_line("most likely multiprecision Hermite normal form will help...");
       Vol := 0; return;
  end Standard_CellVol;

  procedure CellVol
               ( nVar,nSpt : in integer32;
                 Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                 SptType,Cell : in Standard_Integer_Vectors.Link_to_Vector;
                 Vol : in out natural32;
                 multprec_hermite : in boolean := false ) is
  begin
    if multprec_hermite
     then Multprec_CellVol(nVar,nSpt,Spt,SptType,Cell,Vol);
     else Standard_CellVol(nVar,nSpt,Spt,SptType,Cell,Vol);
    end if;
  end CellVol;

  procedure solve_linear_system
               ( n : in integer32;
                 A : in Standard_Floating_Matrices.Matrix;
                 b : in out Standard_Floating_Vectors.Vector;
                 okay : out integer32 ) is

    wrk : Standard_Floating_Matrices.Matrix(1..n,1..n);
    rhs : Standard_Floating_Vectors.Vector(1..n);
    piv : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;

  begin
    for i in wrk'range(1) loop
      rhs(i) := b(i-1);
      for j in wrk'range(2) loop
        wrk(i,j) := A(i-1,j-1);
      end loop;
    end loop;
    lufac(wrk,n,piv,info);
    lusolve(wrk,n,piv,rhs);
    for i in rhs'range loop
      b(i-1) := rhs(i);
    end loop;
    if info = 0
     then okay := 1;
     else okay := 0;
    end if;
  end solve_linear_system;

  function Inner_Normal
               ( n,r : integer32;
                 m,c : Standard_Integer_Vectors.Link_to_Vector;
                 Spt : Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : Standard_Floating_Vectors.Link_to_Vector )
               return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the inner normal spanned by those points in the supports
  --   determined by the labels.

  -- ON ENTRY :
  --   n         number of variables;
  --   r         number of different supports;
  --   m         type of mixture;
  --   c         labels to the points in Spt which span the mixed cell;
  --   Spt       coordinates of the supports;
  --   lft       values of the lifting on the supports.

  -- ON RETURN :
  --   vector of range 0..n-1, normal to the points in the mixed cell.

    res : Standard_Floating_Vectors.Vector(0..n-1);
    A : Standard_Floating_Matrices.Matrix(0..n-1,0..n-1);
    okay : integer32;
    offset : integer32 := 0;
    ind : integer32 := -1;

  begin
    for i in 0..(r-1) loop
      for j in 0..(m(i)-1) loop
        ind := ind + 1;
        for k in 0..(n-1) loop
          A(ind,k) := double_float(Spt(c(offset+j+1))(k)
                                 - Spt(c(offset+j))(k));
        end loop;
        res(ind) := lft(c(offset+j)) - lft(c(offset+j+1));
      end loop;
      offset := offset + m(i) + 1;
    end loop;
    solve_linear_system(n,A,res,okay);
    return res;
  end Inner_Normal;

  procedure write_cells
               ( file : in string; nVar,nSpt : in integer32;
                 SptType : in Standard_Integer_Vectors.Link_to_Vector;
                 Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : in Standard_Floating_Vectors.Link_to_Vector;
                 CellSize,nMCells : in integer32;
                 MCells : in out CellStack ) is

    fp : file_type;
    ind : integer32;
    normal : Standard_Floating_Vectors.Vector(0..nVar-1);
    labels : Standard_Integer_Vectors.Link_to_Vector;
    name : Link_to_String;

  begin
    Create_Output_File(fp,file,name);
    put(fp,nVar,1); new_line(fp);
    put(fp,nSpt,1); new_line(fp);
    for i in 0..(nSpt-1) loop
      put(fp," "); put(fp,SptType(i),1);
    end loop;
    new_line(fp);
    put(fp,nMCells,1); new_line(fp);
    for kk in 0..(nMCells-1) loop
      labels := Cs_Cur(MCells);
      normal := Inner_Normal(nVar,nSpt,SptType,labels,Spt,lft);
      for i in 0..(nVar-1) loop
        put(fp,normal(i)); new_line(fp);
      end loop;
      put(fp,"1.0"); new_line(fp);
      ind := -1;
      for i in 0..(nSpt-1) loop
        put(fp,SptType(i)+1,1); new_line(fp);
        for j in 0..SptType(i) loop
          ind := ind + 1;
          for k in 0..(nVar-1) loop
            put(fp," "); put(fp,Spt(labels(ind))(k),1);
          end loop;
          put(fp," "); put(fp,lft(labels(ind))); new_line(fp);
        end loop;
      end loop;
      put_line(fp,"0");
      if kk /= nMCells - 1            
       then Cs_Pop(MCells); 
      end if;
    end loop;
    close(fp);
    Standard_Integer_Vectors.Clear(labels);
  end write_cells;

end Mixed_Volume;
