with text_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Simplex_Pivoting;                   use Simplex_Pivoting;

-- for testing :
use text_io;
with Standard_Integer_Vectors_io;  use Standard_Integer_Vectors_io;

package body One_Level_LP is

  procedure one_level_LP 
               ( Strt1Pt,End1Pt : in integer32;
                 PtIn : in out Standard_Integer_Vectors.Link_to_Vector;
                 LPdim : in integer32;
                 A : in Standard_Floating_Matrices.Link_to_Matrix;
                 nVar : in integer32;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 ItLp : in out Link_to_IT_LP ) is

    j,info,TstPt : integer32;
    labels : Standard_Integer_Vectors.Link_to_Vector
           := new Standard_Integer_Vectors.Vector(0..nVar-1);
    c : Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector(0..nVar+1);

  begin                        
   -- put_line("entered one_level_LP");
    TstPt := Strt1Pt;            -- extend the cells by using 1-point test 
    for i in 0..(LPdim-1) loop
      c(i) := -A(TstPt,i);
    end loop;
    dnulp2_a(End1Pt+1,LPdim,A,nVar,c,Bidx,x,Binv,info);
    j := -1;    --  record any possible points in Bidx passed 1-point test
    for i in 0..(LPdim-1) loop
      if Bidx(i) >= TstPt then
        PtIn(Bidx(i)) := 1;
        j := j + 1;
        labels(j) := Bidx(i);
      end if;
    end loop;
    j := j + 1;
    if j > 0 then
     -- put("The labels before sort : "); put(labels(0..j-1)); new_line;
      Sort(j,labels);
     -- put("The labels after sort : "); put(labels(0..j-1)); new_line;
      IT_Add1(ItLp,j,labels,LPdim,Bidx,x,Binv);
    end if;
    if info >= 0 then           -- perform 1-point test for other points
      for TstPt in (Strt1Pt+1)..End1Pt loop
        if PtIn(TstPt) = 0 then
          for i in 0..(LPdim-1) loop
            c(i) := -A(TstPt,i);
          end loop;
          dlp1_1pts(End1Pt+1,LPdim,A,nVar,c,
                    TstPt,Bidx,x,Binv,PtIn,ItLp);
         end if;
      end loop;
    else
      for TstPt in (Strt1Pt+1)..End1Pt loop
        if PtIn(TstPt) = 0 then 
          for i in 0..(LPdim-1) loop
            c(i) := -A(TstPt,i);
          end loop;
          dlp2_1pts(End1Pt+1,LPdim,A,nVar,c,
                    TstPt,Bidx,x,Binv,PtIn,ItLp);
        end if; 
      end loop;
    end if;
    Standard_Integer_Vectors.Clear(labels);
    Standard_Floating_Vectors.Clear(c);
  end one_level_LP;

  procedure dnulp2_a
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nVar : in integer32; 
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 info : out integer32 ) is

    ibrk : boolean;
    ell,k,k1,flag : integer32;
    sigj,vmax,bot,top,dtmp : double_float;
    v : Standard_Floating_Vectors.Vector(0..na-1);

  begin
    info := 0;
    loop
      ibrk := false;             -- step 1 : search outgoing constraint
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
        Search_Outgoing(na,c,Bidx,Binv,vmax,k);
        if vmax < eps
         then exit;          -- found optimal solution, solve smaller LP
        end if;
        for i in 0..(na-1) loop
          v(i) := Binv(k,i);
        end loop;
      end if; 
      sigj := DBL_MAX;            -- step 2 : search incoming constraint
      ell := -1;
      for i in 0..(ma-1) loop
        ibrk := false;
        for j in 0..(na-1) loop
          if Bidx(j) = i
           then ibrk := true; exit;
          end if;
        end loop;
        if not ibrk then                           -- if(ibrk) continue;
          bot := a(i,0)*v(0);
          for j in 1..(na-1) loop
            bot := bot + a(i,j)*v(j);
          end loop;
          if bot < -eps then                -- if(bot >= -eps) continue;
            top := -a(i,nVar);
            for j in 0..(na-1) loop
              top := top + a(i,j)*x(j);
            end loop;
            top := top/bot;
            if top < sigj then              -- if(top >= sigj) continue;
              sigj := top;
              ell := i;
            end if;
          end if;
        end if;
      end loop;
      if ell < 0 then
        text_io.put_line("dnulp2_a : LP unbounded"); return;
      end if;
      for i in 0..(na-1) loop          -- step 3: update x, Bidx, Binv
        x(i) := x(i) - sigj*v(i);
      end loop;
      Update_Base(Bidx,Binv,na,k,ell,a);
    end loop;
    k1 := -1;
    for i in 0..(na-1) loop
      if Bidx(i) > -1
       then info := info + 1;
      end if;
    end loop;
    if info = na
     then return;                               -- case of full rank
    end if;
    k := 0;
    while k < na loop
      if((Bidx(k) > -1) or (k = k1)) then
        k := k + 1;                                     -- continue;
      else
        for i in 0..(na-1) loop
          v(i) := Binv(k,i);
        end loop;
        flag := 1;
        ell := -1;              -- step 2 : search incoming constraint
        while ell = -1 loop
          sigj := DBL_MAX;
          for i in 0..(ma-1) loop
            ibrk := false;
            for j in 0..(na-1) loop
              if Bidx(j) = i
               then ibrk := true; exit;
              end if;
            end loop;
            if not ibrk then                     -- if(ibrk) continue;
              bot := a(i,0)*v(0);
              for j in 1..(na-1) loop
                bot := bot + a(i,j)*v(j);
              end loop;
              if bot < -eps then          -- if(bot >= -eps) continue;
                top := -a(i,nVar);
                for j in 0..(na-1) loop
                  top := top + a(i,j)*x(j);
                end loop;
                top := top/bot;
                if top < sigj then        -- if(top >= sigj) continue;
                  sigj := top;
                  ell := i;
                end if;
              end if;
            end if;
          end loop;
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
        if not ibrk then                           -- if(ibrk) continue;
          for i in 0..(na-1) loop        -- step 3: update x, Bidx, Binv
            x(i) := x(i) - sigj*v(i);
          end loop;
          Update_Base(Bidx,Binv,na,k,ell,a);
          info := info + 1;
          if info = na 
           then return;
           else k1 := k;
                k := 0;
          end if;
        end if;
      end if;
    end loop; 
    info := -info;                        --  the case of not full rank
  end dnulp2_a;

  procedure dlp2_1pts
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nVar : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 TstPt : in integer32;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 PtIn : in out Standard_Integer_Vectors.Link_to_Vector;
                 ItLp  : in out Link_to_IT_LP ) is

    k,ell : integer32;
    sigj,vmax : double_float;

  begin
    loop
      Search_Outgoing(na,c,Bidx,Binv,vmax,k);
      if vmax < eps
       then return;                        -- leave with optimal solution
      end if;
      Search_Incoming(ma,na,nVar,k,Bidx,x,a,Binv,sigj,ell);
      for i in 0..(na-1) loop            -- step 3 : update x, Bidx, Binv
        x(i) := x(i) - sigj*Binv(k,i);
      end loop;
      Update_Base(Bidx,Binv,na,k,ell,a);
      if((ell >= TstPt) and (PtIn(ell) = 0)) then   -- save Bidx, x, Binv
        PtIn(ell) := 1;
        IT_Add2(ItLp,ell,na,Bidx,x,Binv);
      end if;
    end loop;
  end dlp2_1pts;

  procedure dlp1_1pts
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nVar : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector; 
                 TstPt : in integer32;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 PtIn : in out Standard_Integer_Vectors.Link_to_Vector;
                 ItLp : in out Link_to_IT_LP ) is

    k,ell : integer32;
    vmax,sigj : double_float;

  begin
    loop
      Search_Outgoing(na,c,Binv,vmax,k);
      if vmax < eps 
       then return;                        --  leave with optimal solution
      end if;
      Search_Incoming(ma,na,nVar,k,Bidx,x,a,Binv,sigj,ell);
      for i in 0..(na-1) loop              -- step 3: update x, Bidx, Binv
        x(i) := x(i) - sigj*Binv(k,i);
      end loop;
      Update_Base(Bidx,Binv,na,k,ell,a);
      if ((ell >= TstPt) and (PtIn(ell) = 0)) then   -- save Bidx, x, Binv
        PtIn(ell) := 1;
        IT_Add2(ItLp,ell,na,Bidx,x,Binv);
      end if;
    end loop;
  end dlp1_1pts;

  procedure Sort ( n : in integer32; 
                   a : in out Standard_Integer_Vectors.Link_to_Vector ) is

    itmp,j : integer32;

  begin
    for i in 1..(n-1) loop
      itmp := a(i);
      j := i;
      while ((j > 0) and then (itmp < a(j-1))) loop
        a(j) := a(j-1);
        j := j - 1;
      end loop;
      a(j) := itmp;
    end loop;
  end Sort;

end One_Level_LP;
