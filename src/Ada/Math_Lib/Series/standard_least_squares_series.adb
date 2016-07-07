with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Mathematical_Functions;     use Standard_Mathematical_Functions;
with Standard_Dense_Series;               use Standard_Dense_Series;
with Standard_Algebraic_Series;
with Standard_Dense_Series_Norms;         use Standard_Dense_Series_Norms;

package body Standard_Least_Squares_Series is

-- AUXILIARIES :

  function min0 ( a,b : integer32 ) return integer32 is

  -- DESCRIPTION : returns the minimum of a and b.

  begin
    if a <= b
     then return a;
     else return b;
    end if;
  end min0;

  function dmax1 ( a,b : double_float ) return double_float is

  -- DESCRIPTION : returns the maximum of a and b.

  begin
    if a >= b
     then return a;
     else return b;
    end if;
  end dmax1;

  function Safe_Norm ( s : Series ) return Series is

  -- DESCRIPTION :
  --   Computing the norm with small negative coefficients
  --   about machine precision generates an exception.
  --   For series with too small coefficients,
  --   the zero series is returned.

    res : Series;
    iszero : boolean := true;

  begin
    for i in 0..s.deg loop
      if abs(REAL_PART(s.cff(i))) > 1.0E-13 then
        iszero := false;
      elsif abs(IMAG_PART(s.cff(i))) > 1.0E-13 then
        iszero := false;
      end if;
      exit when not iszero;
    end loop;
    if iszero
     then res := Create(0.0);
     else res := Norm(s);
    end if;
    return res;
  end Safe_Norm;

  function cdabs ( s : Series ) return Series is

  -- DESCRIPTION :
  --   Computes the 2-norm of the series s.

    res : constant Series := Safe_Norm(s);

  begin
    return res;
  end cdabs;

  function csign ( a,b : Series ) return Series is

  -- DESCRIPTION : translated from
  --       csign(zdum1,zdum2) = cdabs(zdum1)*(zdum2/cdabs(zdum2)) 

    res : Series;
    fac : constant Series := cdabs(a)/cdabs(b);

  begin
    res := fac*b;
    return res;
  end csign;

  procedure zswap ( a : in out Standard_Dense_Series_Matrices.Matrix;
                    k1,k2 : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps the columns k1 and k2 in the matrix a.

    tmp : Series;

  begin
    for i in a'range(1) loop
      tmp := a(i,k1); a(i,k1) := a(i,k2); a(i,k2) := tmp;
    end loop;
  end zswap;

  procedure zcopy ( n,start : in integer32;
                    x : in Standard_Dense_Series_Vectors.Vector;
                    y : out Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Copies the n entries of x to y, beginning at start.

  begin
    for i in start..start+n-1 loop
      y(i) := x(i);
    end loop;
  end zcopy;
 
  function znrm2 ( a : Standard_Dense_Series_Matrices.Matrix;
                   row,col : integer32 ) -- return double_float is
                 return Series is

  -- DESCRIPTION :
  --   Computes the 2-norm of the vector in the column col of the matrix,
  --   starting at the given row.

    sum : Series := Create(0.0);

  begin
    for i in row..a'last(1) loop
      sum := sum + Conjugate(a(i,col))*a(i,col);
    end loop;
    return Standard_Algebraic_Series.sqrt(sum,0);
  end znrm2;

  function zdot ( a : Standard_Dense_Series_Matrices.Matrix;
                  row,k1,k2 : integer32 ) return Series is

  -- DESCRIPTION :
  --   Returns the inner product of the vectors in the columns k1 and k2,
  --   starting at the given row.

    res : Series := Create(0.0);

  begin
    for i in row..a'last(1) loop
      res := res + Conjugate(a(i,k1))*a(i,k2);
    end loop;
    return res;
  end zdot;

  function zdotc ( row : integer32;
                   x : Standard_Dense_Series_Matrices.Matrix;
                   y : Standard_Dense_Series_Vectors.Vector )
                 return Series is

  -- DESCRIPTION :
  --   Dot product of two vectors : x(row..x'last(1),row)*y(row..y'last).

    res : Series := Create(0.0);

  begin
    for i in row..y'last loop
      res := res + Conjugate(x(i,row))*y(i);
    end loop;
    return res;
  end zdotc;

  procedure zaxpy ( a : in out Standard_Dense_Series_Matrices.Matrix; 
                    f : in Series; row,k1,k2 : in integer32 ) is

  -- DESCRIPTION :
  --   The column k2 is added with f times the column k1, starting at row.

  begin
    for i in row..a'last(1) loop
      a(i,k2) := a(i,k2) + f*a(i,k1);
    end loop;
  end zaxpy;

  procedure zscal ( a : in out Standard_Dense_Series_Matrices.Matrix;
                    f : in Series; row,col : in integer32 ) is

  -- DESCRIPTION :
  --   Multiplies the column col of the matrix with f, starting at row.

  begin
    for i in row..a'last(1) loop
      a(i,col) := f*a(i,col);
    end loop;
  end zscal;

  procedure zaxpy ( n,row,col : in integer32; f : in Series;
                    x : in Standard_Dense_Series_Matrices.Matrix;
                    y : in out Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   y(i) := y(i) + f*x(i,col) for i in row..row+n-1.

  begin
    for i in row..row+n-1 loop
      y(i) := y(i) + f*x(i,col);
    end loop;
  end zaxpy;

-- TARGET PROCEDURES :

  procedure QRD ( x : in out Standard_Dense_Series_Matrices.Matrix;
                  qraux : in out Standard_Dense_Series_Vectors.Vector;
                  jpvt : in out Standard_Integer_Vectors.Vector;
                  piv : in boolean ) is

    n : constant integer32 := x'length(1);  -- number of rows
    p : constant integer32 := x'length(2);  -- number of columns
    work : Standard_Dense_Series_Vectors.Vector(x'range(2));
    jj,jp,lp1,lup,maxj,pl,pu : integer32;
    maxnrm,tt : double_float;
    nrmxl,t : Series;
    fac : Complex_Number;
    negj,swapj : boolean;
    one : constant Series := Create(1.0);

  begin
    pl := 1; pu := 0;
    if piv then
      for j in x'range(2) loop       -- rearrange columns according to jpvt
        swapj := (jpvt(j) > 0);
        negj := (jpvt(j) < 0);
        jpvt(j) := j;
        if negj
         then jpvt(j) := -j;
        end if;
        if (swapj and then (j /= pl)) then
          zswap(x,pl,j);
          jpvt(j) := jpvt(pl);
          jpvt(pl) := j;
          pl := pl + 1;
        end if;
      end loop;
      pu := p;
      for j in 1..p loop
        jj := p - j + 1;
        if jpvt(jj) < 0 then
          jpvt(jj) := -jpvt(jj);
          if jj /= pu then
            zswap(x,pu,jj);
            jp := jpvt(pu);
            jpvt(pu) := jpvt(jj);
            jpvt(jj) := jp;
          end if;
          pu := pu - 1;
        end if;
      end loop; 
    end if;
    for j in pl..pu loop                  -- compute norms of the free columns
      qraux(j) := znrm2(x,1,j); -- Create(znrm2(x,1,j));
      work(j) := qraux(j);
    end loop;
    lup := min0(n,p);                -- perform the householder reduction of x
    for ell in 1..lup loop
      if (ell >= pl and ell < pu) then
        maxnrm := 0.0;                  -- locate column with largest norm and
        maxj := ell;                         -- bring it in the pivot position
        for j in ell..pu loop
          if REAL_PART(qraux(j).cff(0)) > maxnrm
           then maxnrm := REAL_PART(qraux(j).cff(0)); maxj := j;
          end if;
        end loop;
        if maxj /= ell then
          zswap(x,ell,maxj);
          qraux(maxj) := qraux(ell);
          work(maxj) := work(ell);
          jp := jpvt(maxj);
          jpvt(maxj) := jpvt(ell);
          jpvt(ell) := jp;
        end if;
      end if;
      qraux(ell) := Create(0.0);
      if ell /= n then
        nrmxl := znrm2(x,ell,ell);              -- householder transformation
        if AbsVal(nrmxl.cff(0)) /= 0.0 then                 -- for column ell
          if REAL_PART(cdabs(x(ell,ell)).cff(0)) /= 0.0
           then nrmxl := csign(nrmxl,x(ell,ell));
          end if;
          zscal(x,Inverse(nrmxl),ell,ell);
          x(ell,ell) := one + x(ell,ell);
          lp1 := ell + 1;       --  apply the transformation to the remaining
          for j in lp1..p loop                --  columns, updating the norms
            t := -zdot(x,ell,ell,j)/x(ell,ell);
            zaxpy(x,t,ell,ell,j);
            if (j >= pl) and (j <= pu)
              and (AbsVal(qraux(j).cff(0)) /= 0.0) then
              tt := 1.0 - (REAL_PART(cdabs(x(ell,j)).cff(0))/
                           REAL_PART(qraux(j).cff(0)))**2;
              tt := dmax1(tt,0.0);
              t := Create(tt);
              tt := 1.0 + 0.05*tt*(REAL_PART(qraux(j).cff(0))/
                                   REAL_PART(work(j).cff(0)))**2;
              if tt /= 1.0 then
                fac := Create(SQRT(REAL_PART(t.cff(0))));
                qraux(j) := qraux(j)*fac;
              else 
                qraux(j) := znrm2(x,ell+1,j);
                work(j) := qraux(j);
              end if;
            end if;
          end loop;
          qraux(ell) := x(ell,ell);               -- save the transformation
          x(ell,ell) := -nrmxl;
        end if;
      end if;
    end loop;
  end QRD;

  procedure Basis ( qr : in out Standard_Dense_Series_Matrices.Matrix;
                    x : in Standard_Dense_Series_Matrices.Matrix ) is

    sum : Series;
    wrk : Standard_Dense_Series_Vectors.Vector(qr'range(1));

  begin
    for j in x'range(2) loop               -- compute jth column of q
      for i in qr'range(1) loop
        sum := x(i,j);
        for k in qr'first(2)..(j-1) loop
          sum := sum - qr(i,k)*qr(k,j);
        end loop;
        wrk(i) := sum/qr(j,j);
      end loop;
      for i in qr'range(1) loop
        qr(i,j) := wrk(i);
      end loop;
    end loop;
  end Basis;

  procedure QRLS ( x : in out Standard_Dense_Series_Matrices.Matrix;
                   n,k : in integer32;
                   qraux,y : in Standard_Dense_Series_Vectors.Vector;
                   qy,qty,b,rsd,xb : out Standard_Dense_Series_Vectors.Vector;
                   job : in integer32; info : out integer32 ) is

    cb,cqy,cqty,cr,cxb : boolean;
    jj,ju,kp1 : integer32;
    t,temp : Series;

  begin
    info := 0;                                               -- set info flag
    cqy := (job/10000 /= 0);                     -- determine what to compute
    cqty := (job mod 10000 /= 0);
    cb := ((job mod 1000)/100 /= 0);
    cr := ((job mod 100)/10 /= 0);
    cxb := ((job mod 10) /= 0);
    ju := min0(k,n-1);
    if ju = 0 then                               -- special action when n = 1
      if cqy then qy(1) := y(1); end if;
      if cqty then qty(1) := y(1); end if;
      if cxb then xb(1) := y(1); end if;
      if cb then
        if AbsVal(x(1,1).cff(0)) = 0.0
         then info := 1;
         else b(1) := y(1)/x(1,1);
        end if;
      end if;
      if cr then rsd(1) := Create(0.0); end if;
      return;
    end if;
    if cqy                                     -- set up to compute qy or qty
     then zcopy(n,y'first,y,qy);
    end if;
    if cqty
     then zcopy(n,y'first,y,qty);
    end if;
    if cqy then                                                -- compute qy
      for j in 1..ju loop
        jj := ju - j + 1;
        if AbsVal(qraux(jj).cff(0)) /= 0.0 then
          temp := x(jj,jj);
          x(jj,jj) := qraux(jj);
          t := -zdotc(jj,x,qy)/x(jj,jj);
          zaxpy(n-jj+1,jj,jj,t,x,qy);
          x(jj,jj) := temp;
        end if;
      end loop;
    end if;
    if cqty then                                        -- compute trans(q)*y
      for j in 1..ju loop
        if AbsVal(qraux(j).cff(0)) /= 0.0 then
          temp := x(j,j);
          x(j,j) := qraux(j);
          t := -zdotc(j,x,qty)/x(j,j);
          zaxpy(n-j+1,j,j,t,x,qty);
          x(j,j) := temp;
        end if;
      end loop;
    end if;
    if cb                                   -- set up to compute b,rsd, or xb
     then zcopy(k,qty'first,qty,b);
    end if;
    kp1 := k + 1;
    if cxb then zcopy(k,qty'first,qty,xb); end if;
    if (cr and (k < n))
     then zcopy(n-k,kp1,qty,rsd);
    end if;
    if (cxb and (kp1 <= n)) then
      for i in kp1..n loop
        xb(i) := Create(0.0);
      end loop;
    end if;
    if cr then
      for i in 1..k loop
        rsd(i) := Create(0.0);
      end loop;
    end if;
    if cb then                                                   -- compute b
      for j in 1..k loop
        jj := k - j + 1;
        if AbsVal(x(jj,jj).cff(0)) = 0.0
         then info := jj; exit;
        end if;
        b(jj) := b(jj)/x(jj,jj);
        if jj /= 1
         then t := -b(jj); zaxpy(jj-1,1,jj,t,x,b);
        end if;
      end loop;
    end if;
    if cr or cxb then                        -- compute rsd or xb as requested
      for j in 1..ju loop
        jj := ju - j + 1;
        if AbsVal(qraux(jj).cff(0)) /= 0.0 then
          temp := x(jj,jj);
          x(jj,jj) := qraux(jj);
          if cr then
            t := -zdotc(jj,x,rsd)/x(jj,jj);
            zaxpy(n-jj+1,jj,jj,t,x,rsd);
          end if;
          if cxb then
            t := -zdotc(jj,x,xb)/x(jj,jj);
            zaxpy(n-jj+1,jj,jj,t,x,xb);
          end if;
          x(jj,jj) := temp;
        end if;
      end loop;
    end if;
  end QRLS;

end Standard_Least_Squares_Series;
