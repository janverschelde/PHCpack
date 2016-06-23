with Quad_Double_Numbers;                 use Quad_Double_Numbers;  
with QuadDobl_Complex_Numbers;            use QuadDobl_Complex_Numbers;
with QuadDobl_Mathematical_Functions;     use QuadDobl_Mathematical_Functions;

package body QuadDobl_Complex_QR_Least_Squares is

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

  function cdabs ( a : Complex_Number ) return quad_double is

  -- DESCRIPTION :
  --   Computes the modulus of the complex number, hopefully this
  --   corresponds to the 'cdabs' fortran function.

    res : constant quad_double := SQRT(sqr(REAL_PART(a)) + sqr(IMAG_PART(a)));

  begin
    return res;
  end cdabs;

  function csign ( a,b : Complex_Number ) return Complex_Number is

  -- DESCRIPTION : translated from
  --       csign(zdum1,zdum2) = cdabs(zdum1)*(zdum2/cdabs(zdum2)) 

  begin
    return (Create(cdabs(a)/cdabs(b))*b);
  end csign;

  procedure zswap ( a : in out QuadDobl_Complex_Matrices.Matrix;
                    k1,k2 : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps the columns k1 and k2 in the matrix a.

    tmp : Complex_Number;

  begin
    for i in a'range(1) loop
      tmp := a(i,k1); a(i,k1) := a(i,k2); a(i,k2) := tmp;
    end loop;
  end zswap;

  procedure zcopy ( n,start : in integer32;
                    x : in QuadDobl_Complex_Vectors.Vector;
                    y : out QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Copies the n entries of x to y, beginning at start.

  begin
    for i in start..start+n-1 loop
      y(i) := x(i);
    end loop;
  end zcopy;
 
  function znrm2 ( a : QuadDobl_Complex_Matrices.Matrix;
                   row,col : integer32 ) return quad_double is

  -- DESCRIPTION :
  --   Computes the 2-norm of the vector in the column col of the matrix,
  --   starting at the given row.

    zero : constant quad_double := create(0.0);
    sum : Complex_Number := Create(zero);

  begin
    for i in row..a'last(1) loop
      sum := sum + Conjugate(a(i,col))*a(i,col);
    end loop;
    return SQRT(REAL_PART(sum));
  end znrm2;

  function zdot ( a : QuadDobl_Complex_Matrices.Matrix;
                  row,k1,k2 : integer32 ) 
                return Complex_Number is

  -- DESCRIPTION :
  --   Returns the inner product of the vectors in the columns k1 and k2,
  --   starting at the given row.

    zero : constant quad_double := create(0.0);
    res : Complex_Number := Create(zero);

  begin
    for i in row..a'last(1) loop
      res := res + Conjugate(a(i,k1))*a(i,k2);
    end loop;
    return res;
  end zdot;

  function zdotc ( row : integer32;
                   x : QuadDobl_Complex_Matrices.Matrix;
                   y : QuadDobl_Complex_Vectors.Vector )
                 return Complex_Number is

  -- DESCRIPTION :
  --   Dot product of two vectors : x(row..x'last(1),row)*y(row..y'last).

    zero : constant quad_double := create(0.0);
    res : Complex_Number := Create(zero);

  begin
    for i in row..y'last loop
      res := res + Conjugate(x(i,row))*y(i);
    end loop;
    return res;
  end zdotc;

  procedure zaxpy ( a : in out QuadDobl_Complex_Matrices.Matrix; 
                    f : in Complex_Number; row,k1,k2 : in integer32 ) is

  -- DESCRIPTION :
  --   The column k2 is added with f times the column k1, starting at row.

  begin
    for i in row..a'last(1) loop
      a(i,k2) := a(i,k2) + f*a(i,k1);
    end loop;
  end zaxpy;

  procedure zscal ( a : in out QuadDobl_Complex_Matrices.Matrix;
                    f : in Complex_Number; row,col : in integer32 ) is

  -- DESCRIPTION :
  --   Multiplies the column col of the matrix with f, starting at row.

  begin
    for i in row..a'last(1) loop
      a(i,col) := f*a(i,col);
    end loop;
  end zscal;

  procedure zaxpy ( n,row,col : in integer32; f : in Complex_Number;
                    x : in QuadDobl_Complex_Matrices.Matrix;
                    y : in out QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   y(i) := y(i) + f*x(i,col) for i in row..row+n-1.

  begin
    for i in row..row+n-1 loop
      y(i) := y(i) + f*x(i,col);
    end loop;
  end zaxpy;

  procedure QRD ( x : in out QuadDobl_Complex_Matrices.Matrix;
                  qraux : in out QuadDobl_Complex_Vectors.Vector;
                  jpvt : in out Standard_Integer_Vectors.Vector;
                  piv : in boolean ) is

    n : constant integer32 := x'length(1);  -- number of rows
    p : constant integer32 := x'length(2);  -- number of columns
    work : QuadDobl_Complex_Vectors.Vector(x'range(2));
    jj,jp,lp1,lup,maxj,pl,pu : integer32;
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);
    maxnrm,tt : quad_double;
    nrmxl,t : Complex_Number;
    negj,swapj : boolean;

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
      qraux(j) := Create(znrm2(x,1,j));
      work(j) := qraux(j);
    end loop;
    lup := min0(n,p);                -- perform the householder reduction of x
    for ell in 1..lup loop
      if (ell >= pl and ell < pu) then
        maxnrm := zero;                 -- locate column with largest norm and
        maxj := ell;                         -- bring it in the pivot position
        for j in ell..pu loop
          if REAL_PART(qraux(j)) > maxnrm
           then maxnrm := REAL_PART(qraux(j)); maxj := j;
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
      qraux(ell) := Create(zero);
      if ell /= n then
        nrmxl := Create(znrm2(x,ell,ell));      -- householder transformation
        if AbsVal(nrmxl) /= zero then                       -- for column ell
          if cdabs(x(ell,ell)) /= zero
           then nrmxl := csign(nrmxl,x(ell,ell));
          end if;
          zscal(x,Create(one)/nrmxl,ell,ell);
          x(ell,ell) := Create(one) + x(ell,ell);
          lp1 := ell + 1;       --  apply the transformation to the remaining
          for j in lp1..p loop                --  columns, updating the norms
            t := -zdot(x,ell,ell,j)/x(ell,ell);
            zaxpy(x,t,ell,ell,j);
            if (j >= pl) and (j <= pu) and (AbsVal(qraux(j)) /= zero) then
              tt := 1.0 - sqr(cdabs(x(ell,j))/REAL_PART(qraux(j)));
              tt := dmax1(tt,zero);
              t := Create(tt);
              tt := 1.0 + 0.05*tt*sqr(REAL_PART(qraux(j))/REAL_PART(work(j)));
              if tt /= one then
                qraux(j) := qraux(j)*Create(SQRT(REAL_PART(t)));
              else 
                qraux(j) := Create(znrm2(x,ell+1,j));
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

  procedure Basis ( qr : in out QuadDobl_Complex_Matrices.Matrix;
                    x : in QuadDobl_Complex_Matrices.Matrix ) is

    sum : Complex_Number;
    wrk : QuadDobl_Complex_Vectors.Vector(qr'range(1));

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

  procedure QRLS ( x : in out QuadDobl_Complex_Matrices.Matrix;
                   n,k : in integer32;
                   qraux,y : in QuadDobl_Complex_Vectors.Vector;
                   qy,qty,b,rsd,xb : out QuadDobl_Complex_Vectors.Vector;
                   job : in integer32; info : out integer32 ) is

    cb,cqy,cqty,cr,cxb : boolean;
    jj,ju,kp1 : integer32;
    zero : constant quad_double := create(0.0);
    t,temp : Complex_Number;

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
        if AbsVal(x(1,1)) = zero
         then info := 1;
         else b(1) := y(1)/x(1,1);
        end if;
      end if;
      if cr then rsd(1) := Create(zero); end if;
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
        if AbsVal(qraux(jj)) /= zero then
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
        if AbsVal(qraux(j)) /= zero then
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
        xb(i) := Create(zero);
      end loop;
    end if;
    if cr then
      for i in 1..k loop
        rsd(i) := Create(zero);
      end loop;
    end if;
    if cb then                                                   -- compute b
      for j in 1..k loop
        jj := k - j + 1;
        if AbsVal(x(jj,jj)) = zero
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
        if AbsVal(qraux(jj)) /= zero then
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

end QuadDobl_Complex_QR_Least_Squares;
