with text_io;                           use text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Random_Numbers;

-- for testing :
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;

package body Prepare_for_MV is

  procedure ReArrangeOrder
               ( nVar,nSup : in integer32; 
                 SptType : in Standard_Integer_Vectors.Link_to_Vector;
                 perm : in out Standard_Integer_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   In case nSup < nVar, this procedure determines how the original
  --   equations must be permuted according to the type of mixture.

    nxt,tmp : integer32;

  begin
    for i in 0..(nVar-1) loop
      if SptType(i) < 0 then
        nxt := -1;                    -- nxt points to next original support
        for j in (i+1)..(nVar-1) loop
          if nxt < 0 and then SptType(j) < 0
           then nxt := j;
          end if;
          if SptType(j) = i then
            if nxt >= 0 then                      -- permutations are needed
              tmp := perm(j);       
              perm(j) := perm(nxt);
              perm(nxt) := tmp;
              nxt := j;
            end if;
          end if;
        end loop;
      end if;
    end loop;
  end ReArrangeOrder;

  procedure Pre4MV
               ( nVar,nSpt : in integer32; nS : out integer32; 
                 SptType : in out Standard_Integer_Vectors.Link_to_Vector;
                 Spt : in out Standard_Integer_VecVecs.Link_to_VecVec;
                 SptIdx : in out Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : in out Standard_Integer_VecVecs.Link_to_VecVec;
                 VtxIdx : in out Standard_Integer_Vectors.Link_to_Vector;
                 OldIdx : in out Standard_Integer_Vectors.Link_to_Vector;
                 perm : out Standard_Integer_Vectors.Link_to_Vector ) is

   ibreak,ynReOrder : boolean;   
   icnt,tmp,nSup : integer32;
   OrderSpt,iTmp,iTmp1 : Standard_Integer_Vectors.Link_to_Vector;
   NuIdx2OldIdxCopy : Standard_Integer_Vectors.Link_to_Vector;
   SptCopy : Standard_Integer_VecVecs.Link_to_VecVec;

  begin
    SortSpt(nVar,nSpt,Spt,SptIdx);          -- sort points lexicographically
    for i in 0..(nVar-1) loop
      for j in SptIdx(i)..(SptIdx(i+1)-2) loop
        tmp := 0;
        for k in 0..(nVar-1) loop
	  tmp := tmp + abs(Spt(j)(k) - Spt(j+1)(k));
        end loop;
        if tmp = 0 then
          put("points "); put(j+1,1); put(" and "); put(j+2,1);
          put(" of support "); put(i+1,1); put_line(" are the same!");
          nS := 0; return;
	end if;
      end loop;
    end loop;
    NonVertex(nVar,nSpt,SptIdx,Spt,OldIdx);      -- locate nonvertex points
    icnt := -1;
    for i in 0..(nVar-1) loop
      for j in SptIdx(i)..(SptIdx(i+1)-1) loop
        if OldIdx(j) = 1 then
          icnt := icnt + 1;
          for k in 0..(nVar-1) loop
            Vtx(icnt)(k) := Spt(j)(k);
          end loop;
          OldIdx(icnt) := j;   -- OldIdx[i] = j means that the current i-th 
	end if;               -- point is the j-th point in original system
      end loop;
      VtxIdx(i+1) := icnt+1;
    end loop;
    VtxIdx(0) := 0;
    ynReOrder := false;                -- reorder the supports if necessary
    OrderSpt := new Standard_Integer_Vectors.Vector(0..nVar-1);
    iTmp := new Standard_Integer_Vectors.Vector(0..nVar);
    iTmp1 := new Standard_Integer_Vectors.Vector(0..nVar);
    SptCopy := new Standard_Integer_VecVecs.VecVec(0..(SptIdx(nVar)-1));
    for i in 0..(SptIdx(nVar)-1) loop
      SptCopy(i) := new Standard_Integer_Vectors.Vector(0..nVar-1);
    end loop;
    NuIdx2OldIdxCopy
      := new Standard_Integer_Vectors.Vector(0..SptIdx(nVar)-1);
    for i in 0..(nVar-1) loop
      iTmp(i) := 0;
    end loop;
    for i in 0..(nVar-1) loop
      tmp := 1000000;
      for j in 0..(nVar-1) loop
        if(iTmp(j) = 0) and (VtxIdx(j+1)-VtxIdx(j) < tmp) then
          tmp := VtxIdx(j+1) - VtxIdx(j);
          OrderSpt(i) := j;
        end if;
      end loop;
      if OrderSpt(i) /= i
       then ynReOrder := true;
      end if;
      iTmp(OrderSpt(i)) := 1;
    end loop;
    if ynReOrder then    -- order of the supports : (OrderSpt[i],i=0,nVar-1)
      iTmp(0) := 0;                       -- change Spt, Vtx, SptIdx, VtxIdx
      iTmp1(0) := 0;
      icnt := -1;
      for i in 0..(nVar-1) loop
        for j in SptIdx(OrderSpt(i))..(SptIdx(OrderSpt(i)+1)-1) loop
          icnt := icnt + 1;
          for k in 0..(nVar-1) loop
            SptCopy(icnt)(k) := Spt(j)(k);
          end loop;
        end loop;
        iTmp(i+1) := icnt+1;
      end loop;
      for j in 0..icnt loop
        for i in 0..(nVar-1) loop
          Spt(j)(i) := SptCopy(j)(i);
        end loop;
      end loop;
      icnt := -1;
      for i in 0..(nVar-1) loop
        for j in VtxIdx(OrderSpt(i))..(VtxIdx(OrderSpt(i)+1)-1) loop
          icnt := icnt + 1;
          for k in 0..(nVar-1) loop
            SptCopy(icnt)(k) := Vtx(j)(k);
          end loop;
          NuIdx2OldIdxCopy(icnt)
            := OldIdx(j) - (SptIdx(OrderSpt(i))-iTmp(i));
	end loop;
        iTmp1(i+1) := icnt+1;
      end loop;
      for j in 0..icnt loop
        for i in 0..(nVar-1) loop
          Vtx(j)(i) := SptCopy(j)(i);
        end loop;
        OldIdx(j) := NuIdx2OldIdxCopy(j);
      end loop;
      for i in 0..nVar loop
        SptIdx(i) := iTmp(i);
      end loop;
      for i in 0..nVar loop
        VtxIdx(i) := iTmp1(i);
      end loop;
    end if;
    for i in 0..(nVar-1) loop              -- look for equal vertex supports
      SptType(i) := -1;
    end loop;
    nSup := 0;
    for i in 0..(nVar-2) loop
      if SptType(i) < 0  then
        nSup := nSup + 1;
        for j in (i+1)..(nVar-1) loop
          if SptType(j) < 0 then
             -- if(VtxIdx[j+1]-VtxIdx[j] != VtxIdx[i+1]-VtxIdx[i]) continue;
            if VtxIdx(j+1)-VtxIdx(j) = VtxIdx(i+1)-VtxIdx(i) then
              ibreak := false;
              for k in VtxIdx(i)..(VtxIdx(i+1)-1) loop
                for ll in 0..(nVar-1) loop
                  if Vtx(k)(ll) /= Vtx(k+VtxIdx(j)-VtxIdx(i))(ll)
                   then ibreak := true; exit;
                  end if;
                end loop;
                exit when ibreak;
              end loop;
              if not ibreak                          -- if(ibreak) continue;
               then SptType(j) := i;
              end if;
            end if;
          end if;
        end loop;
      end if;
    end loop;
    if SptType(nVar-1) < 0
     then nSup := nSup + 1;
    end if;
  -- SptType[i] = -1 : the i-th support is the base support 
  --            =  j : the i-th support is the same as the j-th base support
    if nSup < nVar then -- rearrange the supports according to the type
     -- put("SptType before ReArrange : "); put(SptType); new_line;
      ReArrangeOrder(nVar,nSup,SptType,OrderSpt);
      ReArrangeSpt(nVar,Spt,SptIdx,Vtx,VtxIdx,SptType,OldIdx);
     -- put("SptType after ReArrange : "); put(SptType); new_line;
    else
      nSup := nVar;
      for i in 0..(nVar-1) loop
        SptType(i) := 1;
      end loop;
    end if;
    nS := nSup;
    perm := OrderSpt; -- only for fully mixed !!!
    Standard_Integer_Vectors.Clear(iTmp);
    Standard_Integer_Vectors.Clear(iTmp1);
    Standard_Integer_Vectors.Clear(NuIdx2OldIdxCopy);
    Standard_Integer_VecVecs.Deep_Clear(SptCopy);
  end Pre4MV;

  procedure SortSpt
               ( nVar,nSpt : in integer32;
                 Spt : in out Standard_Integer_VecVecs.Link_to_VecVec;
                 SptIdx : in Standard_Integer_Vectors.Link_to_Vector ) is

    idx,tmp : integer32;

  begin
    for iSpt in 0..(nSpt-1) loop
      for i in SptIdx(iSpt)..(SptIdx(iSpt+1)-2) loop
        idx := i;
        for j in (i+1)..(SptIdx(iSpt+1)-1) loop
          for k in 0..(nVar-1) loop
            if Spt(j)(k) < Spt(idx)(k) then
              exit;                                            -- LexiGt = 0
            else
              if Spt(j)(k) > Spt(idx)(k)
               then idx := j; exit;                            -- LexiGt = 1
              end if;
            end if;
          end loop;
        end loop;
        if idx > i then
          for j in 0..(nVar-1) loop
            tmp := Spt(i)(j);
            Spt(i)(j) := Spt(idx)(j);
            Spt(idx)(j) := tmp;
          end loop;
        end if;
      end loop;
    end loop;
  end SortSpt;

  procedure LowerTriangular
               ( A : in out Standard_Floating_Matrices.Link_to_Matrix;
                 jStrt,jEnd,na : in integer32; rnk : out integer32;
                 ib : in out Standard_Integer_Vectors.Link_to_Vector ) is
   
    k,itmp : integer32;
    dtmp : double_float;

  begin
    ib(0) := jStrt;                  -- the jStrt-th row of A is (1,0,...,0)
    for i in 1..(na-1) loop
      ib(i) := -1;
    end loop;
    rnk := 1;
    k := jStrt+1;
    while ((rnk < na) and (k <= jEnd)) loop
      dtmp := 1.0E-13;                       -- search for largest component
      itmp := -1;                            -- from A[k][rnk] to A[k][na-1]
      for i in rnk..(na-1) loop
        if abs(A(k,i)) > dtmp then
          dtmp := abs(A(k,i));
          itmp := i;
        end if;
      end loop;
      if(itmp >= 0) then                                      -- nonzero row
         for i in 0..(itmp-1) loop          -- to change aa(k,*) into e_itmp
           A(k,i) := A(k,i)/A(k,itmp);
         end loop;
         for i in (itmp+1)..(na-1) loop
           A(k,i) := A(k,i)/A(k,itmp);
         end loop;
         for j in (k+1)..jEnd loop
           for i in 0..(itmp-1) loop
             A(j,i) := A(j,i) - A(k,i)*A(j,itmp);
           end loop;
           for i in (itmp+1)..(na-1) loop
             A(j,i) := A(j,i) - A(k,i)*A(j,itmp);
           end loop;
           A(j,itmp) := A(j,itmp)/A(k,itmp);
         end loop;
        if itmp /= rnk then         -- interchange columns rnk and itmp of A
          for j in jStrt..jEnd loop
            dtmp := A(j,rnk);
            A(j,rnk) := A(j,itmp);
            A(j,itmp) := dtmp;
          end loop;
        end if;
        for i in 0..rnk loop
          A(k,i) := 0.0;
        end loop;
        A(k,rnk) := 1.0;
        for i in (rnk+1)..(na-1) loop
          A(k,i) := 0.0;
        end loop;
        rnk := rnk + 1;
        ib(rnk-1) := k;
      end if;
      k := k + 1;
    end loop;
  end LowerTriangular;

  procedure RSimplex
               ( A : in Standard_Floating_Matrices.Link_to_Matrix;
                 b : in Standard_Floating_Vectors.Link_to_Vector;
                 jStrt,jEnd,m : in integer32;
                 xb,sb : in out Standard_Floating_Vectors.Link_to_Vector;
                 ib : in out Standard_Integer_Vectors.Link_to_Vector;
                 binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 info : out integer32 ) is

    ibrk : boolean;
    k,ell : integer32;
    sigb,sum,vkmin : double_float;
    eps : constant double_float := 1.0E-6;
    tol : constant double_float := 1.0E-10;
    DBL_MAX : constant double_float := 1.0E+20;
                
  begin
    info := -1;  
    for i in 0..(m-1) loop
      if ib(i) = jStrt
       then info := i; exit;
      end if;
    end loop;
    loop                                             -- compute u and find k
      if info = -1 then                     -- $e_{n+1}$ is not in the basis
        k := jStrt;
        vkmin := -1.0;
      else                                      -- $e_{n+1}$ is in the basis
        k := -1;             -- find vkmin=v_k=min{c_i+a_i^Tu | i\notin I_b}
        vkmin := DBL_MAX;
        for i in jStrt..jEnd loop                     -- check if i is in ib
          ibrk := false;
          for j in 0..(m-1) loop
            if i = ib(j)
             then ibrk := true; exit;
            end if;
          end loop;
          if not ibrk then                             -- if(ibrk) continue;
            sum := 0.0;                                 -- if i is not in ib
            for j in 0..(m-1) loop
              sum := sum + A(i,j)*binv(info,j);
            end loop;
            if sum < vkmin 
             then vkmin := sum; k := i;
            end if;
          end if;
        end loop; 
        if vkmin > -eps
	 then return;                          -- found the optimal solution
        end if;
      end if;           -- k-th column will get into the basis -- to form sb
      for i in 0..(m-1) loop
        sum := 0.0;
        for j in 0..(m-1) loop
          sum := sum + binv(i,j)*A(k,j);
        end loop;
        sb(i) := sum;
      end loop;
      ell := 0;
      sigb := DBL_MAX;
      for i in 0..(m-1) loop
        if (sb(i) > eps) and then ((xb(i)/sb(i)) < sigb)
         then sigb := xb(i)/sb(i);
              ell := i;
        end if;
      end loop;        -- ell-th row gets out of the basis and k-th row gets
      sum := 0.0;                 -- into the basis.  to find the new B^{-1}
      for i in 0..(m-1) loop
        sum := sum + binv(ell,i)*A(k,i);
      end loop;
      if abs(sum) <= eps
       then return;                                              -- abort();
      end if;         
      sum := 1.0/sum;
      for i in 0..(m-1) loop
        binv(ell,i) := binv(ell,i)*sum;
      end loop;
      for j in 0..(m-1) loop
        if j /= ell then
          sum := 0.0;
          for i in 0..(m-1) loop
            sum := sum + binv(j,i)*A(k,i);
          end loop;
          for i in 0..(m-1) loop
            binv(j,i) := binv(j,i) - sum*binv(ell,i);
          end loop;
        end if;
      end loop;
      for i in 0..(m-1) loop         -- form the new basic feasible solution
         sum := 0.0;
         for j in 0..(m-1) loop
           sum := sum + binv(i,j)*b(j);
         end loop;
         xb(i) := sum;
      end loop;
      if ib(ell) = jStrt
       then info := 0;                       -- column jStrt is out of basis
      end if;
      ib(ell) := k;
      if k = jStrt then
        info := ell;
        if xb(info) > tol
         then return;                         -- epsilon is already positive
        end if;
      end if;
    end loop;
  end RSimplex;

  procedure ReArrangeSpt
               ( nVar : in integer32;
                 Spt : in out Standard_Integer_VecVecs.Link_to_VecVec;
                 SptIdx : in out Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : in out Standard_Integer_VecVecs.Link_to_VecVec;
                 VtxIdx : in out Standard_Integer_Vectors.Link_to_Vector;
                 SptType : in out Standard_Integer_Vectors.Link_to_Vector;
                 OldIdx : in out Standard_Integer_Vectors.Link_to_Vector ) is

    icnt,iSpt,ith : integer32;
    iTmp : Standard_Integer_Vectors.Link_to_Vector
         := new Standard_Integer_Vectors.Vector(0..nVar-1);
    SptIdxCopy,Nu2OldCopy : Standard_Integer_Vectors.Link_to_Vector;
    SptCopy : Standard_Integer_VecVecs.Link_to_VecVec;
   
  begin
    iSpt := -1;
    for i in 0..(nVar-1) loop
      if SptType(i) < 0 then                 -- the i-th support is original
        iSpt := iSpt + 1;
        iTmp(iSpt) := i;
      end if;
    end loop;                              -- iSpt counts #original supports
    ith := -1;
    for i in 0..(nVar-1) loop
      if SptType(i) < 0 then            -- count occurrences of i-th support
        ith := ith + 1;
        icnt := 1;
        for j in (i+1)..(nVar-1) loop
          if SptType(j) = i then
            iSpt := iSpt + 1;
            icnt := icnt + 1;
            iTmp(iSpt) := j;
          end if;
        end loop;
        SptType(ith) := icnt;
      end if;        -- SptType[i] = k : #supports which are the same as the 
    end loop;        --             i-th base support; rearrange Spt, SptIdx
   -- put("iTmp = "); put(iTmp); new_line;
    SptCopy := new Standard_Integer_VecVecs.VecVec(0..SptIdx(nVar)-1);
    for i in 0..(SptIdx(nVar)-1) loop
      SptCopy(i) := new Standard_Integer_Vectors.Vector(0..nVar-1);
    end loop;
    SptIdxCopy := new Standard_Integer_Vectors.Vector(0..nVar);
    for j in 0..(SptIdx(nVar)-1) loop
      for i in 0..(nVar-1) loop
        SptCopy(j)(i) := Spt(j)(i);
      end loop;
    end loop;
    icnt := -1;
   -- put_line("Rearranging the supports ...");
    for i in 0..(nVar-1) loop
      for j in SptIdx(iTmp(i))..(SptIdx(iTmp(i)+1)-1) loop
        icnt := icnt + 1;
        for i1 in 0..(nVar-1) loop
          Spt(icnt)(i1) := SptCopy(j)(i1);
        end loop;
      end loop;
      SptIdxCopy(i+1) := icnt+1;
    end loop;
    SptIdxCopy(0) := 0;                                  -- rearrange OldIdx
    Nu2OldCopy := new Standard_Integer_Vectors.Vector(0..SptIdx(nVar)-1);
    for j in 0..(VtxIdx(nVar)-1) loop
      Nu2OldCopy(j) := OldIdx(j);
    end loop;
    icnt := -1;
    for i in 0..(nVar-1) loop
      for j in VtxIdx(iTmp(i))..(VtxIdx(iTmp(i)+1)-1) loop
        icnt := icnt + 1;
        OldIdx(icnt) := Nu2OldCopy(j) - (SptIdx(iTmp(i)) - SptIdxCopy(i));
      end loop;
    end loop;
    for i in 0..nVar loop
      SptIdx(i) := SptIdxCopy(i);
    end loop;
   -- put_line("Rearranging the vertices ...");
    for j in 0..(VtxIdx(nVar)-1) loop               -- rearrange Vtx, VtxIdx
      for i in 0..(nVar-1) loop
        SptCopy(j)(i) := Vtx(j)(i);
      end loop;
    end loop;
    icnt := -1;
    for i in 0..(nVar-1) loop
      for j in VtxIdx(iTmp(i))..VtxIdx(iTmp(i)+1)-1 loop
        icnt := icnt + 1;
        for i1 in 0..(nVar-1) loop
          Vtx(icnt)(i1) := SptCopy(j)(i1);
        end loop;
      end loop;
      SptIdxCopy(i+1) := icnt+1;
    end loop;
    VtxIdx(0) := 0;
    for i in 1..nVar loop
      VtxIdx(i) := SptIdxCopy(i);
    end loop;
    Standard_Integer_Vectors.Clear(iTmp);
    Standard_Integer_Vectors.Clear(Nu2OldCopy);
    Standard_Integer_Vectors.Clear(SptIdxCopy);
    Standard_Integer_VecVecs.Deep_Clear(SptCopy);
  end ReArrangeSpt;

  procedure NonVertex
               ( nVar,nSpt : in integer32; 
                 SptIdx : in Standard_Integer_Vectors.Link_to_Vector;
                 Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                 ynVtx : in out Standard_Integer_Vectors.Link_to_Vector ) is

   info,itmp,jStrt,jEnd,PtIn,rnk : integer32;
   tmp : double_float;
   ib : Standard_Integer_Vectors.Link_to_Vector
      := new Standard_Integer_Vectors.Vector(0..nVar+1);
   ib0 : Standard_Integer_Vectors.Link_to_Vector
       := new Standard_Integer_Vectors.Vector(0..nVar+1);
   sb : Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector(0..nVar+1);
   xb : Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector(0..nVar+1);
   b : Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector(0..SptIdx(nSpt));
   binv : Standard_Floating_Matrices.Link_to_Matrix
        := new Standard_Floating_Matrices.Matrix(0..nVar+1,0..nVar+1);
   A : Standard_Floating_Matrices.Link_to_Matrix
     := new Standard_Floating_Matrices.Matrix(0..SptIdx(nSpt),0..nVar+1);

  begin
    for i in 0..(SptIdx(nSpt)-1) loop
      ynVtx(i) := 1;         -- = 0: nonvertex point; to form the matrix for
    end loop;                --      simplex method to choose a random point
    for i in 0..(nVar-1) loop
      b(i) := Standard_Random_Numbers.Random;
    end loop;
    for j in 0..(SptIdx(nSpt)-1) loop                     -- lift the points 
      tmp := 0.0;
      for k in 0..(nVar-1) loop
        tmp := tmp + (double_float(Spt(j)(k)) - b(k))
                    *(double_float(Spt(j)(k)) - b(k));
      end loop;
      A(j+1,0) := -Sqrt(tmp);  -- j+1: shift a slot for the added 1st column
      for k in 0..(nVar-1) loop
        A(j+1,k+1) := double_float(Spt(j)(k));
      end loop;
      A(j+1,nVar+1) := 1.0;
    end loop;
    for i in 0..(nSpt-1) loop                            -- for each support
      jStrt := SptIdx(i);          -- The first row is the lifting direction
      jEnd := SptIdx(i+1);      -- +1: shift a slot for the added 1st column 
      A(jStrt,0) := 1.0;
      for k in 1..(nVar+1) loop 
        A(jStrt,k) := 0.0;
      end loop;
      rnk := 0;                            -- lower-triangulate the matrix A 
      LowerTriangular(A,jStrt,jEnd,nVar+2,rnk,ib0);
      for j in 0..(rnk-1) loop
         for k in 0..(rnk-1) loop
           binv(j,k) := 0.0;
         end loop;
         binv(j,j) := 1.0;
      end loop;
      for j in (jStrt+1)..jEnd loop                        -- for each point
        PtIn := j;
        for k in 0..(rnk-1) loop
          ib(k) := -1;
        end loop;
        for k in 0..(rnk-1) loop
          if ib0(k) = j
           then ib(k) := PtIn; PtIn := -PtIn; exit;
          end if;  
        end loop;
        if PtIn >= 0 then        -- decide which column will leave the basis
          for k in 0..(rnk-1) loop       -- j-th row will get into the basis
             xb(k) := 0.0;
             for i1 in 0..(rnk-1) loop
               xb(k) := xb(k) + binv(k,i1)*A(j,i1);
             end loop;
          end loop;
          tmp := 0.0; -- find the biggest component from xb[0] to xb[*rnk-1]
          itmp := -1;
          for k in 1..(rnk-1) loop
            if (abs(xb(k)) > tmp) and (ib(k) = -1) then
              tmp := abs(xb(k));
              itmp := k;
            end if;
          end loop;
          if itmp = -1 then
            info := -1; return;       -- abort();  impossible in our problem 
          end if;                   -- itmp-th row will get out of the basis 
          for k in 0..(itmp-1) loop        -- find the new B^{-1} after j in 
            xb(k) := xb(k)/xb(itmp);
          end loop;
          for k in (itmp+1)..(rnk-1) loop
            xb(k) := xb(k)/xb(itmp);
          end loop;
          for i1 in 0..(rnk-1) loop
            tmp := binv(itmp,i1);
            for k in 0..(itmp-1) loop
              binv(k,i1) := binv(k,i1) - xb(k)*tmp;   
            end loop;
            for k in (itmp+1)..(rnk-1) loop
              binv(k,i1) := binv(k,i1) - xb(k)*tmp;
            end loop;
            binv(itmp,i1) := tmp/xb(itmp);
          end loop;
          ib(itmp) := j;
        end if;
        for k in 0..(rnk-1) loop
          if ib(k) /= -1 then
            xb(k) := 1.0;
          else
            xb(k) := 0.0;
            ib(k) := ib0(k);
          end if;
        end loop;
        for k in 0..(rnk-1) loop
          b(k) := A(j,k);
        end loop;
        RSimplex(A,b,jStrt,jEnd,rnk,xb,sb,ib,binv,info);
        if info > -1 then
          if abs(xb(info)) > 1.0e-10                   -- it is a non-vertex
           then ynVtx(j-1) := 0;           -- -1: since the added 1st column
          end if;
        end if;
        for k in 0..(rnk-1) loop
          ib0(k) := ib(k);
        end loop;
      end loop;
    end loop;
    Standard_Integer_Vectors.Clear(ib);
    Standard_Integer_Vectors.Clear(ib0);
    Standard_Floating_Vectors.Clear(sb);
    Standard_Floating_Vectors.Clear(xb);
    Standard_Floating_Vectors.Clear(b);
    Standard_Floating_Matrices.Clear(binv);
    Standard_Floating_Matrices.Clear(A);
  end NonVertex;

end Prepare_for_MV;
