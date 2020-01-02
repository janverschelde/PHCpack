with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;            use QuadDobl_Complex_Numbers;
with QuadDobl_Mathematical_Functions;     use QuadDobl_Mathematical_Functions;
with QuadDobl_Complex_Series;             use QuadDobl_Complex_Series;
with QuadDobl_Complex_Algebraic_Series;
with QuadDobl_Complex_Series_Norms;       use QuadDobl_Complex_Series_Norms;

package body QuadDobl_Series_Least_Squares is

-- AUXILIARIES :

  function min0 ( a,b : integer32 ) return integer32 is

  -- DESCRIPTION : returns the minimum of a and b.

  begin
    if a <= b
     then return a;
     else return b;
    end if;
  end min0;

  function dmax1 ( a,b : quad_double ) return quad_double is

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

    res : Series(s.deg);
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
     then res := Create(0,res.deg);
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

    fac : constant Series := cdabs(a)/cdabs(b);
    res : constant Series := fac*b;

  begin
    return res;
  end csign;

  procedure zswap ( a : in out QuadDobl_Complex_Series_Matrices.Matrix;
                    k1,k2 : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps the columns k1 and k2 in the matrix a.

    tmp : Link_to_Series;

  begin
    for i in a'range(1) loop
      tmp := a(i,k1);
      a(i,k1) := a(i,k2);
      a(i,k2) := tmp;
    end loop;
  end zswap;

  procedure zcopy ( n,start : in integer32;
                    x : in QuadDobl_Complex_Series_Vectors.Vector;
                    y : out QuadDobl_Complex_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Copies the n entries of x to y, beginning at start.

  begin
    for i in start..start+n-1 loop
      Copy(x(i),y(i));               -- y(i) := x(i);
    end loop;
  end zcopy;
 
  function znrm2 ( a : QuadDobl_Complex_Series_Matrices.Matrix;
                   row,col : integer32 ) -- return double_float is
                 return Series is

  -- DESCRIPTION :
  --   Computes the 2-norm of the vector in the column col of the matrix,
  --   starting at the given row.

    deg : constant integer32 := a(a'first(1),a'first(2)).deg;
    sum : Series(deg) := Create(0,deg);
    wrk : Series(deg);

  begin
    for i in row..a'last(1) loop -- sum := sum + Conjugate(a(i,col))*a(i,col);
      wrk := Conjugate(a(i,col).all)*a(i,col).all;
      sum := sum + wrk;
    end loop;
    return QuadDobl_Complex_Algebraic_Series.sqrt(sum,0);
  end znrm2;

  function zdot ( a : QuadDobl_Complex_Series_Matrices.Matrix;
                  row,k1,k2 : integer32 ) return Series is

  -- DESCRIPTION :
  --   Returns the inner product of the vectors in the columns k1 and k2,
  --   starting at the given row.

    deg : constant integer32 := a(a'first(1),a'first(2)).deg;
    res : Series(deg) := Create(0,deg);
    wrk : Series(deg);

  begin
    for i in row..a'last(1) loop -- res := res + Conjugate(a(i,k1))*a(i,k2);
      wrk := Conjugate(a(i,k1).all)*a(i,k2).all;
      res := res + wrk;
    end loop;
    return res;
  end zdot;

  function zdotc ( row : integer32;
                   x : QuadDobl_Complex_Series_Matrices.Matrix;
                   y : QuadDobl_Complex_Series_Vectors.Vector )
                 return Series is

  -- DESCRIPTION :
  --   Dot product of two vectors : x(row..x'last(1),row)*y(row..y'last).

    deg : constant integer32 := x(x'first(1),x'first(2)).deg;
    res : Series(deg) := Create(0,deg);
    wrk : Series(deg);

  begin
    for i in row..y'last loop -- res := res + Conjugate(x(i,row))*y(i);
      if y(i) /= null then
        wrk := Conjugate(x(i,row).all)*y(i).all;
        res := res + wrk;
      end if;
    end loop;
    return res;
  end zdotc;

  procedure zaxpy ( a : in out QuadDobl_Complex_Series_Matrices.Matrix; 
                    f : in Series; row,k1,k2 : in integer32 ) is

  -- DESCRIPTION :
  --   The column k2 is added with f times the column k1, starting at row.

    deg : constant integer32 := a(a'first(1),a'first(1)).deg;
    wrk : Series(deg);

  begin
    for i in row..a'last(1) loop
      wrk := f*a(i,k1).all;         -- a(i,k2) := a(i,k2) + f*a(i,k1);
      Add(a(i,k2).all,wrk);
    end loop;
  end zaxpy;

  procedure zscal ( a : in out QuadDobl_Complex_Series_Matrices.Matrix;
                    f : in Series; row,col : in integer32 ) is

  -- DESCRIPTION :
  --   Multiplies the column col of the matrix with f, starting at row.

  begin
    for i in row..a'last(1) loop
      Mul(a(i,col).all,f);          -- a(i,col) := f*a(i,col);
    end loop;
  end zscal;

  procedure zaxpy ( n,row,col : in integer32; f : in Series;
                    x : in QuadDobl_Complex_Series_Matrices.Matrix;
                    y : in out QuadDobl_Complex_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   y(i) := y(i) + f*x(i,col) for i in row..row+n-1.

    deg : constant integer32 := x(x'first(1),x'first(2)).deg;
    wrk : Series(deg);

  begin
    for i in row..row+n-1 loop
      wrk := f*x(i,col).all;      -- y(i) := y(i) + f*x(i,col);
      if y(i) = null
       then y(i) := new Series'(wrk);
       else Add(y(i).all,wrk);
      end if;
    end loop;
  end zaxpy;

  function Create ( x : quad_double; deg : integer32 ) return Series is

    cx : constant Complex_Number := Create(x);

  begin
    return Create(cx,deg);
  end Create;

-- TARGET PROCEDURES :

  procedure QRD ( x : in out QuadDobl_Complex_Series_Matrices.Matrix;
                  qraux : in out QuadDobl_Complex_Series_Vectors.Vector;
                  jpvt : in out Standard_Integer_Vectors.Vector;
                  piv : in boolean ) is

    zero : constant quad_double := create(0.0);
    dd_one : constant quad_double := create(1.0);
    absq : quad_double;
    n : constant integer32 := x'length(1);  -- number of rows
    p : constant integer32 := x'length(2);  -- number of columns
    deg : constant integer32 := x(x'first(1),x'first(2)).deg;
    work : QuadDobl_Complex_Series_Vectors.Vector(x'range(2));
    jj,jp,lp1,lup,maxj,pl,pu : integer32;
    maxnrm,tt,quot : quad_double;
    nrmxl,t : Series(deg);
    fac : Complex_Number;
    negj,swapj : boolean;
    one : constant Series(deg) := Create(1,deg);

  begin
    for i in work'range loop
      work(i) := new Series'(QuadDobl_Complex_Series.Create(0,deg));
    end loop;
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
      qraux(j).all := znrm2(x,1,j); -- Create(znrm2(x,1,j));
      work(j).cff := qraux(j).cff;
    end loop;
    lup := min0(n,p);                -- perform the householder reduction of x
    for ell in 1..lup loop
      if (ell >= pl and ell < pu) then
        maxnrm := zero;                 -- locate column with largest norm and
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
      qraux(ell) := Create(0,deg);
      if ell /= n then
        nrmxl := znrm2(x,ell,ell);              -- householder transformation
        absq := AbsVal(nrmxl.cff(0));
        if absq /= zero then                                -- for column ell
          if REAL_PART(cdabs(x(ell,ell).all).cff(0)) /= zero
           then nrmxl := csign(nrmxl,x(ell,ell).all);
          end if;
          zscal(x,Inverse(nrmxl),ell,ell);
          Add(x(ell,ell).all,one);  -- x(ell,ell) := one + x(ell,ell);
          lp1 := ell + 1;       --  apply the transformation to the remaining
          for j in lp1..p loop                --  columns, updating the norms
            t := -zdot(x,ell,ell,j)/x(ell,ell).all;
            zaxpy(x,t,ell,ell,j);
            if (j >= pl) and (j <= pu)
              and (AbsVal(qraux(j).cff(0)) /= zero) then
              quot := REAL_PART(cdabs(x(ell,j).all).cff(0))/
                      REAL_PART(qraux(j).cff(0));
              quot := quot*quot;
              tt := 1.0 - quot;
              tt := dmax1(tt,zero);
              t := Create(tt,deg);
              quot := REAL_PART(qraux(j).cff(0))/REAL_PART(work(j).cff(0));
              quot := quot*quot;
              tt := 1.0 + 0.05*tt*quot;
              if tt /= dd_one then
                fac := Create(SQRT(REAL_PART(t.cff(0))));
                Mul(qraux(j),fac); -- qraux(j) := qraux(j)*fac;
              else 
                qraux(j).all := znrm2(x,ell+1,j);
                work(j).cff := qraux(j).cff;
              end if;
            end if;
          end loop;
          qraux(ell).cff := x(ell,ell).cff;      -- save the transformation
          x(ell,ell).all := -nrmxl;
        end if;
      end if;
    end loop;
    QuadDobl_Complex_Series_Vectors.Clear(work);
  end QRD;

  procedure Basis ( qr : in out QuadDobl_Complex_Series_Matrices.Matrix;
                    x : in QuadDobl_Complex_Series_Matrices.Matrix ) is

    deg : constant integer32 := qr(qr'first(1),qr'first(2)).deg;
    sum,acc : Series(deg);
    wrk : QuadDobl_Complex_Series_Vectors.Vector(qr'range(1));

  begin
    for i in wrk'range loop
      wrk(i) := new Series'(QuadDobl_Complex_Series.Create(0,deg));
    end loop;
    for j in x'range(2) loop               -- compute jth column of q
      for i in qr'range(1) loop
        sum.cff := x(i,j).cff;
        for k in qr'first(2)..(j-1) loop
          acc := qr(i,k).all*qr(k,j).all;  -- sum := sum - qr(i,k)*qr(k,j);
          sum := sum - acc;
        end loop;
        acc := sum/qr(j,j).all;            -- wrk(i) := sum/qr(j,j);
        wrk(i).cff := acc.cff;
      end loop;
      for i in qr'range(1) loop
        Copy(wrk(i),qr(i,j));              -- qr(i,j) := wrk(i);
      end loop;
    end loop;
    QuadDobl_Complex_Series_Vectors.Clear(wrk);
  end Basis;

  procedure QRLS ( x : in out QuadDobl_Complex_Series_Matrices.Matrix;
                   n,k : in integer32;
                   qraux,y : in QuadDobl_Complex_Series_Vectors.Vector;
                   qy,qty,b,rsd,xb : out QuadDobl_Complex_Series_Vectors.Vector;
                   job : in integer32; info : out integer32 ) is

    zero : constant quad_double := create(1.0);
    absq : quad_double;
    deg : constant integer32 := x(x'first(1),x'first(2)).deg;
    cb,cqy,cqty,cr,cxb : boolean;
    jj,ju,kp1 : integer32;
    t,temp : Series(deg);

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
        absq := AbsVal(x(1,1).cff(0));
        if absq = zero
         then info := 1;
         else b(1) := y(1)/x(1,1);
        end if;
      end if;
      if cr then rsd(1) := Create(0,deg); end if;
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
        absq := AbsVal(qraux(jj).cff(0));
        if absq /= zero then
          temp.cff := x(jj,jj).cff;
          x(jj,jj).cff := qraux(jj).cff;
          t := -zdotc(jj,x,qy)/x(jj,jj).all;
          zaxpy(n-jj+1,jj,jj,t,x,qy);
          x(jj,jj).cff := temp.cff;
        end if;
      end loop;
    end if;
    if cqty then                                        -- compute trans(q)*y
      for j in 1..ju loop
        absq := AbsVal(qraux(j).cff(0));
        if absq /= zero then
          temp.cff := x(j,j).cff;
          x(j,j).cff := qraux(j).cff;
          t := -zdotc(j,x,qty)/x(j,j).all;
          zaxpy(n-j+1,j,j,t,x,qty);
          x(j,j).cff := temp.cff;
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
        xb(i) := Create(0,deg);
      end loop;
    end if;
    if cr then
      for i in 1..k loop
        rsd(i) := Create(0,deg);
      end loop;
    end if;
    if cb then                                                   -- compute b
      for j in 1..k loop
        jj := k - j + 1;
        absq := AbsVal(x(jj,jj).cff(0));
        if absq = zero
         then info := jj; exit;
        end if;
        Div(b(jj),x(jj,jj));     --  b(jj) := b(jj)/x(jj,jj);
        if jj /= 1
         then t := -b(jj).all; zaxpy(jj-1,1,jj,t,x,b);
        end if;
      end loop;
    end if;
    if cr or cxb then                        -- compute rsd or xb as requested
      for j in 1..ju loop
        jj := ju - j + 1;
        absq := AbsVal(qraux(jj).cff(0));
        if absq /= zero then
          temp.cff := x(jj,jj).cff;
          x(jj,jj).cff := qraux(jj).cff;
          if cr then
            t := -zdotc(jj,x,rsd)/x(jj,jj).all;
            zaxpy(n-jj+1,jj,jj,t,x,rsd);
          end if;
          if cxb then
            t := -zdotc(jj,x,xb)/x(jj,jj).all;
            zaxpy(n-jj+1,jj,jj,t,x,xb);
          end if;
          x(jj,jj).cff := temp.cff;
        end if;
      end loop;
    end if;
  end QRLS;

end QuadDobl_Series_Least_Squares;
