with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with QuadDobl_Mathematical_Functions;     use QuadDobl_Mathematical_Functions;

package body Quad_Double_QR_Least_Squares is

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

  function dsign ( a,b : quad_double ) return quad_double is

  -- DESCRIPTION : returns the absolute value of a, set to the same sign as b.

  begin
    if b >= 0.0
     then return abs(a); 
     else return -abs(a);
    end if;
  end dsign;

  procedure dswap ( a : in out Quad_Double_Matrices.Matrix;
                    k1,k2 : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps the columns k1 and k2 in the matrix a.

    tmp : quad_double;

  begin
    for i in a'range(1) loop
      tmp := a(i,k1); a(i,k1) := a(i,k2); a(i,k2) := tmp;
    end loop;
  end dswap;

  procedure dcopy ( n,start : in integer32;
                    x : in Quad_Double_Vectors.Vector;
                    y : out Quad_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Copies the n entries of x to y, beginning at start.

  begin
    for i in start..start+n-1 loop
      y(i) := x(i);
    end loop;
  end dcopy;
 
  function dnrm2 ( a : Quad_Double_Matrices.Matrix;
                   row,col : integer32 ) return Quad_Double is

  -- DESCRIPTION :
  --   Computes the 2-norm of the vector in the column col of the matrix,
  --   starting at the given row.

    sum : quad_double := create(0.0);

  begin
    for i in row..a'last(1) loop
      sum := sum + a(i,col)*a(i,col);
    end loop;
    return SQRT(sum);
  end dnrm2;

  function ddot ( a : Quad_Double_Matrices.Matrix;
                  row,k1,k2 : integer32 ) return Quad_Double is

  -- DESCRIPTION :
  --   Returns the inner product of the vectors in the columns k1 and k2,
  --   starting at the given row.

    res : quad_double := create(0.0);

  begin
    for i in row..a'last(1) loop
      res := res + a(i,k1)*a(i,k2);
    end loop;
    return res;
  end ddot;

  function ddot ( row : integer32;
                  x : Quad_Double_Matrices.Matrix;
                  y : Quad_Double_Vectors.Vector )
                return Quad_Double is

  -- DESCRIPTION :
  --   Dot product of two vectors : x(row..x'last(1),row)*y(row..y'last).

    res : quad_double := create(0.0);

  begin
    for i in row..y'last loop
      res := res + x(i,row)*y(i);
    end loop;
    return res;
  end ddot;

  procedure daxpy ( a : in out Quad_Double_Matrices.Matrix; 
                    f : in quad_double; row,k1,k2 : in integer32 ) is

  -- DESCRIPTION :
  --   The column k2 is added with f times the column k1, starting at row.

  begin
    for i in row..a'last(1) loop
      a(i,k2) := a(i,k2) + f*a(i,k1);
    end loop;
  end daxpy;

  procedure daxpy ( n,row,col : in integer32; f : in Quad_Double;
                    x : in Quad_Double_Matrices.Matrix;
                    y : in out Quad_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   y(i) := y(i) + f*x(i,col) for i in row..row+n-1.

  begin
    for i in row..row+n-1 loop
      y(i) := y(i) + f*x(i,col);
    end loop;
  end daxpy;

  procedure dscal ( a : in out Quad_Double_Matrices.Matrix;
                    f : in quad_double; row,col : in integer32 ) is

  -- DESCRIPTION :
  --   Multiplies the column col of the matrix with f, starting at row.

  begin
    for i in row..a'last(1) loop
      a(i,col) := f*a(i,col);
    end loop;
  end dscal;

  procedure QRD ( x : in out Quad_Double_Matrices.Matrix;
                  qraux : in out Quad_Double_Vectors.Vector;
                  jpvt : in out Standard_Integer_Vectors.Vector;
                  piv : in boolean ) is

    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);
    n : constant integer32 := x'length(1);  -- number of rows
    p : constant integer32 := x'length(2);  -- number of columns
    work : Quad_Double_Vectors.Vector(x'range(2));
    jj,jp,lp1,lup,maxj,pl,pu : integer32;
    maxnrm,tt,nrmxl,t : Quad_Double;
    negj,swapj : boolean;

  begin
    pl := 1; pu := 0;
    if piv then
      for j in x'range(2) loop         -- rearrange columns according to jpvt
        swapj := (jpvt(j) > 0);
        negj := (jpvt(j) < 0);
        jpvt(j) := j;
        if negj
         then jpvt(j) := -j;
        end if;
        if (swapj and then (j /= pl)) then
          dswap(x,pl,j);
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
            dswap(x,pu,jj);
            jp := jpvt(pu);
            jpvt(pu) := jpvt(jj);
            jpvt(jj) := jp;
          end if;
          pu := pu - 1;
        end if;
      end loop; 
    end if;
    for j in pl..pu loop                  -- compute norms of the free columns
      qraux(j) := dnrm2(x,1,j);
      work(j) := qraux(j);
    end loop;
    lup := min0(n,p);                -- perform the householder reduction of x
    for l in 1..lup loop
      if (l >= pl and l < pu) then
        maxnrm := zero;                  -- locate column with largest norm and
        maxj := l;                           -- bring it in the pivot position
        for j in l..pu loop
          if qraux(j) > maxnrm
           then maxnrm := qraux(j); maxj := j;
          end if;
        end loop;
        if maxj /= l then
          dswap(x,l,maxj);
          qraux(maxj) := qraux(l);
          work(maxj) := work(l);
          jp := jpvt(maxj);
          jpvt(maxj) := jpvt(l);
          jpvt(l) := jp;
        end if;
      end if;
      qraux(l) := zero;
      if l /= n then
        nrmxl := dnrm2(x,l,l);   -- householder transformation for column l
        if nrmxl /= zero then
          if x(l,l) /= zero
           then nrmxl := dsign(nrmxl,x(l,l));
          end if;
          dscal(x,one/nrmxl,l,l);
          x(l,l) := one + x(l,l);
          lp1 := l + 1;   --  apply the transformation to the remaining
          for j in lp1..p loop          --  columns, updating the norms
            t := -ddot(x,l,l,j)/x(l,l);
            daxpy(x,t,l,l,j);
            if (j >= pl) and (j <= pu) and (qraux(j) /= zero) then
              tt := one - sqr(abs(x(l,j))/qraux(j));
              tt := dmax1(tt,zero);
              t := tt;
              tt := one + 0.05*tt*sqr(qraux(j)/work(j));
              if tt /= one
               then qraux(j) := qraux(j)*SQRT(t);
               else qraux(j) := dnrm2(x,l+1,j); work(j) := qraux(j);
              end if;
            end if;
          end loop;
          qraux(l) := x(l,l);                -- save the transformation
          x(l,l) := -nrmxl;
        end if;
      end if;
    end loop;
  end QRD;

  procedure Permute_Columns ( x : in out Quad_Double_Matrices.Matrix;
                              jpvt : in Standard_Integer_Vectors.Vector ) is

    res : Quad_Double_Matrices.Matrix(x'range(1),x'range(2));

  begin
    for k in jpvt'range loop
      for i in res'range(1) loop
        res(i,k) := x(i,jpvt(k));
      end loop;
    end loop;
    x := res;
  end Permute_Columns;

  procedure Permute ( x : in out Quad_Double_Vectors.Vector;
                      jpvt : in Standard_Integer_Vectors.Vector ) is

    res : Quad_Double_Vectors.Vector(x'range);

  begin
    for k in jpvt'range loop
      res(k) := x(jpvt(k));
    end loop;
    x := res;
  end Permute;

  procedure Basis ( qr : in out Quad_Double_Matrices.Matrix;
                    x : in Quad_Double_Matrices.Matrix ) is

    sum : quad_double;
    wrk : Quad_Double_Vectors.Vector(qr'range(1));

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

  procedure QRLS ( x : in out Quad_Double_Matrices.Matrix;
                   ldx,n,k : in integer32;
                   qraux,y : in Quad_Double_Vectors.Vector;
                   qy,qty,b,rsd,xb : out Quad_Double_Vectors.Vector;
                   job : in integer32; info : out integer32 ) is

    zero : constant quad_double := create(0.0);
    cb,cqy,cqty,cr,cxb : boolean;
    jj,ju,kp1 : integer32;
    t,temp : quad_double;

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
      if x(1,1) = zero
       then info := 1;
       else b(1) := y(1)/x(1,1);
      end if;
      end if;
      if cr then rsd(1) := zero; end if;
      return;
    end if;
    if cqy                                     -- set up to compute qy or qty
     then dcopy(n,y'first,y,qy);
    end if;
    if cqty
     then dcopy(n,y'first,y,qty);
    end if;
    if cqy then                                                 -- compute qy
      for j in 1..ju loop
        jj := ju - j + 1;
        if qraux(jj) /= zero then
          temp := x(jj,jj);
          x(jj,jj) := qraux(jj);
          t := -ddot(jj,x,qy)/x(jj,jj);
          daxpy(n-jj+1,jj,jj,t,x,qy);
          x(jj,jj) := temp;
        end if;
      end loop;
    end if;
    if cqty then                                       -- compute trans(q)*y
      for j in 1..ju loop
        if qraux(j) /= zero then
          temp := x(j,j);
          x(j,j) := qraux(j);
          t := -ddot(j,x,qty)/x(j,j);
          daxpy(n-j+1,j,j,t,x,qty);
          x(j,j) := temp;
        end if;
      end loop;
    end if;
    if cb                                   -- set up to compute b,rsd, or xb
     then dcopy(k,qty'first,qty,b);
    end if;
    kp1 := k + 1;
    if cxb then dcopy(k,qty'first,qty,xb); end if;
    if (cr and (k < n))
     then dcopy(n-k,kp1,qty,rsd);
    end if;
    if (cxb and (kp1 <= n)) then
      for i in kp1..n loop
        xb(i) := zero;
      end loop;
    end if;
    if cr then
      for i in 1..k loop
        rsd(i) := zero;
      end loop;
    end if;
    if cb then                                                   -- compute b
      for j in 1..k loop
        jj := k - j + 1;
        if x(jj,jj) = zero 
         then info := jj; exit;
        end if;
        b(jj) := b(jj)/x(jj,jj);
        if jj /= 1
         then t := -b(jj); daxpy(jj-1,1,jj,t,x,b);
        end if;
      end loop;
    end if;
    if cr or cxb then                       -- compute rsd or xb as requested
      for j in 1..ju loop
        jj := ju - j + 1;
        if qraux(jj) /= zero then
          temp := x(jj,jj);
          x(jj,jj) := qraux(jj);
          if cr then
            t := -ddot(jj,x,rsd)/x(jj,jj);
            daxpy(n-jj+1,jj,jj,t,x,rsd);
          end if;
          if cxb then
            t := -ddot(jj,x,xb)/x(jj,jj);
            daxpy(n-jj+1,jj,jj,t,x,xb);
          end if;
          x(jj,jj) := temp;
        end if;
      end loop;
    end if;
  end QRLS;

end Quad_Double_QR_Least_Squares;
