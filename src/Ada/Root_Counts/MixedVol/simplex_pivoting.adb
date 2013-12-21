package body Simplex_Pivoting is

  procedure Search_Outgoing
               ( na : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in Standard_Floating_Matrices.Link_to_Matrix;
                 vmax : out double_float; k : out integer32 ) is

    dtmp : double_float;

  begin
    vmax := -DBL_MAX;
    for i in 0..(na-1) loop
      dtmp := c(0)*Binv(i,0);
      for j in 1..(na-1) loop
        dtmp := dtmp + c(j)*Binv(i,j);
      end loop;
      if dtmp > vmax then
        vmax := dtmp;
        k := i;
      end if; 
    end loop;
  end Search_Outgoing;

  procedure Search_Outgoing
               ( na : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 Bidx : in Standard_Integer_Vectors.Link_to_Vector;
                 Binv : in Standard_Floating_Matrices.Link_to_Matrix;
                 vmax : out double_float; k : out integer32 ) is

    dtmp : double_float;

  begin
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
  end Search_Outgoing;

  procedure Search_Incoming
               ( ma,na,nv1,k : in integer32;
                 LPidx,Bidx : in Standard_Integer_Vectors.Link_to_Vector;
                 x : in Standard_Floating_Vectors.Link_to_Vector;
                 a,Binv : in Standard_Floating_Matrices.Link_to_Matrix;
                 sigj : out double_float; ell : out integer32 ) is

    ii : integer32;
    ibrk : boolean;
    bot,top : double_float;
  
  begin
    sigj := DBL_MAX;
    ell := -1;
    for i in 0..(ma-1) loop
      ii := LPidx(i);
      ibrk := false;
      for j in 0..(na-1) loop
        if Bidx(j) = ii
         then ibrk := true; exit;
        end if;
      end loop;
      if not ibrk then                               -- if(ibrk) continue;
        bot := a(ii,0)*Binv(k,0);
        for j in 1..(na-1) loop 
          bot := bot + a(ii,j)*Binv(k,j);
        end loop;
        if bot < -eps then                    -- if(bot >= -eps) continue;
          top := -a(ii,nv1);
          for j in 0..(na-1) loop
            top := top + a(ii,j)*x(j);
          end loop;
          top := top/bot;
          if top < sigj then                  -- if(top >= sigj) continue;
            sigj := top;
            ell := ii;
          end if;
        end if;
      end if;
    end loop;
    if ell < 0
     then raise unbounded_LP;
    end if;                    
  end Search_Incoming;

  procedure Search_Incoming
               ( ma,na,nv1,k : in integer32;
                 LPidx,Bidx : in Standard_Integer_Vectors.Link_to_Vector;
                 x : in Standard_Floating_Vectors.Link_to_Vector;
                 v : in Standard_Floating_Vectors.Vector;
                 a,Binv : in Standard_Floating_Matrices.Link_to_Matrix;
                 sigj : out double_float; ell : out integer32 ) is

    ibrk : boolean;
    ii : integer32;
    bot,top : double_float;

  begin
    sigj := DBL_MAX;
    ell := -1;
    for i in 0..(ma-1) loop
      ii := LPidx(i);
      ibrk := false;
      for j in 0..(na-1) loop
        if Bidx(j) = ii
         then ibrk := true; exit;
        end if;
      end loop;
      if not ibrk then                               -- if(ibrk) continue;
        bot := a(ii,0)*v(0);
        for j in 1..(na-1) loop
          bot := bot + a(ii,j)*v(j);
        end loop;
        if bot < -eps then                    -- if(bot >= -eps) continue;
          top := -a(ii,nv1);
          for j in 0..(na-1) loop
            top := top + a(ii,j)*x(j);
          end loop;
          top := top/bot;
          if top < sigj then                  -- if(top >= sigj) continue;
            sigj := top;
            ell := ii;
          end if;
        end if;
      end if;
    end loop;
  end Search_Incoming;

  procedure Search_Incoming
               ( ma,na,nVar,k : in integer32;
                 Bidx : in Standard_Integer_Vectors.Link_to_Vector;
                 x : in Standard_Floating_Vectors.Link_to_Vector;
                 a,Binv : in Standard_Floating_Matrices.Link_to_Matrix;
                 sigj : out double_float; ell : out integer32 ) is

    ibrk : boolean;
    bot,top : double_float;

  begin
    sigj := DBL_MAX;
    ell := -1;
    for i in 0..(ma-1) loop
      ibrk := false;
      for j in 0..(na-1) loop
        if Bidx(j) = i
         then ibrk := true; exit;
        end if; 
      end loop;
      if not ibrk then                             -- if(ibrk) continue;
        bot := a(i,0)*Binv(k,0);
        for j in 1..(na-1) loop
          bot := bot + a(i,j)*Binv(k,j);
        end loop;
        if bot <= -eps then                  -- if(bot > -eps) continue;
          top := -a(i,nVar);
          for j in 0..(na-1) loop
            top := top + a(i,j)*x(j);
          end loop;
          top := top/bot;
          if top <= sigj then                -- if(top > sigj) continue;
            sigj := top;
            ell := i;
          end if;
        end if;
      end if;
    end loop;
    if ell < 0
     then raise unbounded_LP;
    end if;
  end Search_Incoming;

  procedure Update_Base
               ( Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 na,k,ell : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix ) is

    top : double_float;

  begin
    top := a(ell,0)*Binv(k,0);
    for i in 1..(na-1) loop
      top := top + a(ell,i)*Binv(k,i);
    end loop;
    if abs(top) < eps
     then raise singular_base;
    end if;
    top := 1.0/top;
    for i in 0..(na-1) loop
      Binv(k,i) := Binv(k,i)*top;
    end loop;
    for j in 0..(na-1) loop
      if j /= k then
        top := a(ell,0)*Binv(j,0);
        for i in 1..(na-1) loop
          top := top + a(ell,i)*Binv(j,i);
        end loop;
        for i in 0..(na-1) loop
          Binv(j,i) := Binv(j,i) - top*Binv(k,i);
        end loop;
      end if;
    end loop;
    Bidx(k) := ell;
  end Update_Base;

end Simplex_Pivoting;
