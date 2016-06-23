with Multprec_Floating_Numbers;           use Multprec_Floating_Numbers;
with Multprec_Mathematical_Functions;     use Multprec_Mathematical_Functions;

package body Multprec_Floating_QR_Least_Squares is

-- AUXILIARIES :

  function min0 ( a,b : integer32 ) return integer32 is

  -- DESCRIPTION : returns the minimum of a and b.

  begin
    if a <= b
     then return a;
     else return b;
    end if;
  end min0;

  function dmax1 ( a,b : Floating_Number ) return Floating_Number is

  -- DESCRIPTION : returns the maximum of a and b.

  begin
    if a > b
     then return a;
     else return b;
    end if;
  end dmax1;

  function dsign ( a,b : Floating_Number ) return Floating_Number is

  -- DESCRIPTION : returns the absolute value of a, set to the same sign as b.

  begin
    if b < 0.0
     then return -AbsVal(a);
     else return AbsVal(a); 
    end if;
  end dsign;

  procedure dswap ( a : in out Multprec_Floating_Matrices.Matrix;
                    k1,k2 : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps the columns k1 and k2 in the matrix a.

    tmp : Floating_Number;

  begin
    for i in a'range(1) loop
      Copy(a(i,k1),tmp);      -- tmp := a(i,k1);
      Copy(a(i,k2),a(i,k1));  -- a(i,k1) := a(i,k2); 
      Copy(tmp,a(i,k2));      -- a(i,k2) := tmp;
    end loop;
    Clear(tmp);
  end dswap;

  procedure dcopy ( n,start : in integer32;
                    x : in Multprec_Floating_Vectors.Vector;
                    y : out Multprec_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Copies the n entries of x to y, beginning at start.

  begin
    for i in start..start+n-1 loop
      Copy(x(i),y(i));  -- y(i) := x(i);
    end loop;
  end dcopy;
 
  function dnrm2 ( a : Multprec_Floating_Matrices.Matrix;
                   row,col : integer32 ) return Floating_Number is

  -- DESCRIPTION :
  --   Computes the 2-norm of the vector in the column col of the matrix,
  --   starting at the given row.

    sum : Floating_Number := create(0.0);
    acc : Floating_Number;

  begin
    for i in row..a'last(1) loop
     -- sum := sum + a(i,col)*a(i,col);
      acc :=  a(i,col)*a(i,col);
      Add(sum,acc);
      Clear(acc);
    end loop;
    return SQRT(sum);
  end dnrm2;

  function ddot ( a : Multprec_Floating_Matrices.Matrix;
                  row,k1,k2 : integer32 ) return Floating_Number is

  -- DESCRIPTION :
  --   Returns the inner product of the vectors in the columns k1 and k2,
  --   starting at the given row.

    res : Floating_Number := create(0.0);
    acc : Floating_Number;

  begin
    for i in row..a'last(1) loop
     -- res := res + a(i,k1)*a(i,k2);
      acc := a(i,k1)*a(i,k2);
      Add(res,acc);
      Clear(acc);
    end loop;
    return res;
  end ddot;

  function ddot ( row : integer32;
                  x : Multprec_Floating_Matrices.Matrix;
                  y : Multprec_Floating_Vectors.Vector )
                return Floating_Number is

  -- DESCRIPTION :
  --   Dot product of two vectors : x(row..x'last(1),row)*y(row..y'last).

    res : Floating_Number := create(0.0);
    acc : Floating_Number;

  begin
    for i in row..y'last loop
     -- res := res + x(i,row)*y(i);
      acc := x(i,row)*y(i);
      Add(res,acc);
      Clear(acc);
    end loop;
    return res;
  end ddot;

  procedure daxpy ( a : in out Multprec_Floating_Matrices.Matrix; 
                    f : in Floating_Number; row,k1,k2 : in integer32 ) is

  -- DESCRIPTION :
  --   The column k2 is added with f times the column k1, starting at row.

    acc : Floating_Number;

  begin
    for i in row..a'last(1) loop
     -- a(i,k2) := a(i,k2) + f*a(i,k1);
      acc := f*a(i,k1);
      Add(a(i,k2),acc);
      Clear(acc);
    end loop;
  end daxpy;

  procedure daxpy ( n,row,col : in integer32; f : in Floating_Number;
                    x : in Multprec_Floating_Matrices.Matrix;
                    y : in out Multprec_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   y(i) := y(i) + f*x(i,col) for i in row..row+n-1.

    acc : Floating_Number;

  begin
    for i in row..row+n-1 loop
     -- y(i) := y(i) + f*x(i,col);
      acc := f*x(i,col);
      Add(y(i),acc);
      Clear(acc);
    end loop;
  end daxpy;

  procedure dscal ( a : in out Multprec_Floating_Matrices.Matrix;
                    f : in Floating_Number; row,col : in integer32 ) is

  -- DESCRIPTION :
  --   Multiplies the column col of the matrix with f, starting at row.

  begin
    for i in row..a'last(1) loop
      Mul(a(i,col),f);  -- a(i,col) := f*a(i,col);
    end loop;
  end dscal;

  procedure QRD ( x : in out Multprec_Floating_Matrices.Matrix;
                  qraux : in out Multprec_Floating_Vectors.Vector;
                  jpvt : in out Standard_Integer_Vectors.Vector;
                  piv : in boolean ) is

    zero : Floating_Number := create(0.0);
    one : Floating_Number := create(1.0);
    n : constant integer32 := x'length(1);  -- number of rows
    p : constant integer32 := x'length(2);  -- number of columns
    work : Multprec_Floating_Vectors.Vector(x'range(2));
    jj,jp,lp1,lup,maxj,pl,pu : integer32;
    maxnrm,tt,nrmxl,t,acc1,acc2 : Floating_Number;
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
    for j in pl..pu loop               -- compute norms of the free columns
      qraux(j) := dnrm2(x,1,j);
      Copy(qraux(j),work(j));
    end loop;
    lup := min0(n,p);             -- perform the householder reduction of x
    for l in 1..lup loop
      if (l >= pl and l < pu) then
        Copy(zero,maxnrm);           -- locate column with largest norm and
        maxj := l;                        -- bring it in the pivot position
        for j in l..pu loop
          if qraux(j) > maxnrm
           then Copy(qraux(j),maxnrm); maxj := j;
          end if;
        end loop;
        if maxj /= l then
          dswap(x,l,maxj);
          Copy(qraux(l),qraux(maxj));
          Copy(work(l),work(maxj));
          jp := jpvt(maxj);
          jpvt(maxj) := jpvt(l);
          jpvt(l) := jp;
        end if;
      end if;
      Copy(zero,qraux(l));
      if l /= n then
        nrmxl := dnrm2(x,l,l);   -- householder transformation for column l
        if not Equal(nrmxl,0.0) then
          if not Equal(x(l,l),0.0) then
            acc1 := dsign(nrmxl,x(l,l));
            Copy(acc1,nrmxl); Clear(acc1);
          end if;
          acc1 := one/nrmxl;
          dscal(x,acc1,l,l); Clear(acc1);
          Add(x(l,l),1.0);
          lp1 := l + 1;   --  apply the transformation to the remaining
          for j in lp1..p loop          --  columns, updating the norms
            acc1 := ddot(x,l,l,j);
            Div(acc1,x(l,l));
            t := -acc1; Clear(acc1);
            daxpy(x,t,l,l,j);
            if (j >= pl) and (j <= pu) and (qraux(j) /= zero) then
             -- tt := one - (AbsVal(x(l,j))/qraux(j))**2;
              acc1 := AbsVal(x(l,j));
              Div(acc1,qraux(j));
              acc2 := acc1*acc1; Clear(acc1);
              tt := one - acc2; Clear(acc2);
              acc1 := dmax1(tt,zero);
              Copy(acc1,tt); Clear(acc1);
              Copy(tt,t);
             -- tt := one + 0.05*tt*(qraux(j)/work(j))**2;
              acc1 := qraux(j)/work(j);
              acc2 := acc1*acc1; Clear(acc1);
              Mul(acc2,tt);
              Mul(acc2,0.05);
              Clear(tt);
              tt := one + acc2; Clear(acc2);
              if not Equal(tt,1.0) then
                acc1 := SQRT(t);
                Mul(qraux(j),acc1); Clear(acc1);
              else
                Copy(dnrm2(x,l+1,j),qraux(j));
                Copy(qraux(j),work(j));
              end if;
            end if;
            Clear(t); Clear(tt);
          end loop;
          Copy(x(l,l),qraux(l));             -- save the transformation
          Copy(nrmxl,x(l,l));
          Min(x(l,l));
        end if;
        Clear(nrmxl);
      end if;
    end loop;
    Clear(zero); Clear(one);
  end QRD;

  procedure Permute_Columns ( x : in out Multprec_Floating_Matrices.Matrix;
                              jpvt : in Standard_Integer_Vectors.Vector ) is

    res : Multprec_Floating_Matrices.Matrix(x'range(1),x'range(2));

  begin
    for k in jpvt'range loop
      for i in res'range(1) loop
        res(i,k) := x(i,jpvt(k));
      end loop;
    end loop;
    x := res;
  end Permute_Columns;

  procedure Permute ( x : in out Multprec_Floating_Vectors.Vector;
                      jpvt : in Standard_Integer_Vectors.Vector ) is

    res : Multprec_Floating_Vectors.Vector(x'range);

  begin
    for k in jpvt'range loop
      res(k) := x(jpvt(k));
    end loop;
    x := res;
  end Permute;

  procedure Basis ( qr : in out Multprec_Floating_Matrices.Matrix;
                    x : in Multprec_Floating_Matrices.Matrix ) is

    sum,acc : Floating_Number;
    wrk : Multprec_Floating_Vectors.Vector(qr'range(1));

  begin
    for j in x'range(2) loop               -- compute jth column of q
      for i in qr'range(1) loop
        Copy(x(i,j),sum);
        for k in qr'first(2)..(j-1) loop
         -- sum := sum - qr(i,k)*qr(k,j);
          acc := qr(i,k)*qr(k,j);
          Sub(sum,acc);
          Clear(acc);
        end loop;
        wrk(i) := sum/qr(j,j);
        Clear(sum);
      end loop;
      for i in qr'range(1) loop
        Copy(wrk(i),qr(i,j));
        Clear(wrk(i));
      end loop;
    end loop;
  end Basis;

  procedure QRLS ( x : in out Multprec_Floating_Matrices.Matrix;
                   ldx,n,k : in integer32;
                   qraux,y : in Multprec_Floating_Vectors.Vector;
                   qy,qty,b,rsd,xb : out Multprec_Floating_Vectors.Vector;
                   job : in integer32; info : out integer32 ) is

    zero : Floating_Number := create(0.0);
    cb,cqy,cqty,cr,cxb : boolean;
    jj,ju,kp1 : integer32;
    t,temp,acc : Floating_Number;

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
      if Equal(x(1,1),0.0)
       then info := 1;
       else Clear(b(1)); b(1) := y(1)/x(1,1);
      end if;
      end if;
      if cr then Clear(rsd(1)); rsd(1) := Create(0.0); end if;
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
          Copy(x(jj,jj),temp);
          Copy(qraux(jj),x(jj,jj));
          acc := ddot(jj,x,qy);
          t := acc/x(jj,jj); Clear(acc);
          Min(t); 
          daxpy(n-jj+1,jj,jj,t,x,qy);
          Copy(temp,x(jj,jj)); Clear(temp);
        end if;
      end loop;
    end if;
    if cqty then                                       -- compute trans(q)*y
      for j in 1..ju loop
        if qraux(j) /= zero then
          Copy(x(j,j),temp);
          Copy(qraux(j),x(j,j));
          acc := ddot(j,x,qty);
          t := acc/x(j,j); Clear(acc);
          Min(t);
          daxpy(n-j+1,j,j,t,x,qty);
          Copy(temp,x(j,j)); Clear(temp);
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
        Copy(zero,xb(i));
      end loop;
    end if;
    if cr then
      for i in 1..k loop
        Copy(zero,rsd(i));
      end loop;
    end if;
    if cb then                                                   -- compute b
      for j in 1..k loop
        jj := k - j + 1;
        if x(jj,jj) = zero 
         then info := jj; exit;
        end if;
        Div(b(jj),x(jj,jj));
        if jj /= 1
         then Clear(t); t := -b(jj); daxpy(jj-1,1,jj,t,x,b);
        end if;
      end loop;
    end if;
    if cr or cxb then                       -- compute rsd or xb as requested
      for j in 1..ju loop
        jj := ju - j + 1;
        if not Equal(qraux(jj),0.0) then
          Copy(x(jj,jj),temp);
          Copy(qraux(jj),x(jj,jj));
          if cr then
            acc := ddot(jj,x,rsd);
            t := acc/x(jj,jj); Clear(acc);
            Min(t);
            daxpy(n-jj+1,jj,jj,t,x,rsd);
          end if;
          if cxb then
            acc := ddot(jj,x,xb);
            t := acc/x(jj,jj); Clear(acc);
            Min(t);
            daxpy(n-jj+1,jj,jj,t,x,xb);
          end if;
          Copy(temp,x(jj,jj)); Clear(temp);
        end if;
      end loop;
    end if;
  end QRLS;

end Multprec_Floating_QR_Least_Squares;
