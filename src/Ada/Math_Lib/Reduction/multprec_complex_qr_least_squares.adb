with Multprec_Floating_Numbers;           use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;            use Multprec_Complex_Numbers;
with Multprec_Mathematical_Functions;     use Multprec_Mathematical_Functions;

package body Multprec_Complex_QR_Least_Squares is

-- AUXILIARIES :

  function min0 ( a,b : integer32 ) return integer32 is

  -- DESCRIPTION : returns the minimum of a and b.

  begin
    if a <= b
     then return a;
     else return b;
    end if;
  end min0;

  function cdabs ( a : Complex_Number ) return Floating_Number is

  -- DESCRIPTION :
  --   Computes the modulus of the complex number as
  --     SQRT(REAL_PART(a)**2 + IMAG_PART(a)**2);

    res,rep,re2,ima,im2 : Floating_Number;

  begin
    rep := REAL_PART(a);
    re2 := rep*rep;
    Clear(rep);
    ima := IMAG_PART(a);
    im2 := ima*ima;
    Clear(ima);
    Add(re2,im2);
    Clear(im2);
    res := SQRT(re2);
    Clear(re2);
    return res;
  end cdabs;

  function csign ( a,b : Complex_Number ) return Complex_Number is

  -- DESCRIPTION : translated from
  --       csign(zdum1,zdum2) = cdabs(zdum1)*(zdum2/cdabs(zdum2)) 

    res : Complex_Number;
    nrma,nrmb : Floating_Number;

  begin
    nrma := cdabs(a);
    nrmb := cdabs(b);
    Div(nrma,nrmb);
    res := Create(nrma);
    Mul(res,b);
    Clear(nrma); Clear(nrmb);
    return res;
  end csign;

  procedure zswap ( a : in out Multprec_Complex_Matrices.Matrix;
                    k1,k2 : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps the columns k1 and k2 in the matrix a.

    tmp : Complex_Number;

  begin
    for i in a'range(1) loop
      tmp := a(i,k1);
      a(i,k1) := a(i,k2);
      a(i,k2) := tmp;
    end loop;
  end zswap;

  procedure zcopy ( n,start : in integer32;
                    x : in Multprec_Complex_Vectors.Vector;
                    y : out Multprec_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Copies the n entries of x to y, beginning at start.

  begin
    for i in start..start+n-1 loop
      Copy(x(i),y(i));
    end loop;
  end zcopy;
 
  function znrm2 ( a : Multprec_Complex_Matrices.Matrix;
                   row,col : integer32 ) return Floating_Number is

  -- DESCRIPTION :
  --   Computes the 2-norm of the vector in the column col of the matrix,
  --   starting at the given row.

    res : Floating_Number;
    sum,acc : Complex_Number;
    rep_sum : Floating_Number;

  begin
    sum := Create(integer(0));
    for i in row..a'last(1) loop
      acc := Conjugate(a(i,col));
      Mul(acc,a(i,col));
      Add(sum,acc);
      Clear(acc);
    end loop;
    rep_sum := REAL_PART(sum);
    res := SQRT(rep_sum);
    Clear(rep_sum);
    Clear(sum);
    return res;
  end znrm2;

  function zdot ( a : Multprec_Complex_Matrices.Matrix;
                  row,k1,k2 : integer32 ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the inner product of the vectors in the columns k1 and k2,
  --   starting at the given row.

    res : Complex_Number := Create(integer(0));
    acc : Complex_Number;

  begin
    for i in row..a'last(1) loop
      acc := Conjugate(a(i,k1));
      Mul(acc,a(i,k2));
      Add(res,acc);
      Clear(acc);
    end loop;
    return res;
  end zdot;

  function zdotc ( row : integer32;
                   x : Multprec_Complex_Matrices.Matrix;
                   y : Multprec_Complex_Vectors.Vector )
                 return Complex_Number is

  -- DESCRIPTION :
  --   Dot product of two vectors : x(row..x'last(1),row)*y(row..y'last).

    res : Complex_Number := Create(integer(0));
    acc : Complex_Number;

  begin
    for i in row..y'last loop
      acc := Conjugate(x(i,row));
      Mul(acc,y(i));
      Add(res,acc);
      Clear(acc);
    end loop;
    return res;
  end zdotc;

  procedure zaxpy ( a : in out Multprec_Complex_Matrices.Matrix; 
                    f : in Complex_Number; row,k1,k2 : in integer32 ) is

  -- DESCRIPTION :
  --   The column k2 is added with f times the column k1, starting at row.

    acc : Complex_Number;

  begin
    for i in row..a'last(1) loop
      acc := f*a(i,k1);
      Add(a(i,k2),acc);
      Clear(acc);
    end loop;
  end zaxpy;

  procedure zscal ( a : in out Multprec_Complex_Matrices.Matrix;
                    f : in Complex_Number; row,col : in integer32 ) is

  -- DESCRIPTION :
  --   Multiplies the column col of the matrix with f, starting at row.

  begin
    for i in row..a'last(1) loop
      Mul(a(i,col),f);
    end loop;
  end zscal;

  procedure zaxpy ( n,row,col : in integer32; f : in Complex_Number;
                    x : in Multprec_Complex_Matrices.Matrix;
                    y : in out Multprec_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   y(i) := y(i) + f*x(i,col) for i in row..row+n-1.

    acc : Complex_Number;

  begin
    for i in row..row+n-1 loop
      acc := f*x(i,col);
      Add(y(i),acc);
      Clear(acc);
    end loop;
  end zaxpy;

  procedure QRD ( x : in out Multprec_Complex_Matrices.Matrix;
                  qraux : in out Multprec_Complex_Vectors.Vector;
                  jpvt : in out Standard_Integer_Vectors.Vector;
                  piv : in boolean ) is

    n : constant integer32 := x'length(1);  -- number of rows
    p : constant integer32 := x'length(2);  -- number of columns
    work : Multprec_Complex_Vectors.Vector(x'range(2));
    jj,jp,lp1,lup,maxj,pl,pu : integer32;
    maxnrm,tt,fltacc1,fltacc2,fltacc3,fltacc4 : Floating_Number;
    nrmxl,t,one,cmpacc : Complex_Number;
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
    for j in pl..pu loop                   -- compute norms of the free columns
      fltacc1 := znrm2(x,1,j);
      qraux(j) := Create(fltacc1);
      Clear(fltacc1);
      Copy(qraux(j),work(j));
    end loop;
    lup := min0(n,p);                 -- perform the householder reduction of x
    for l in 1..lup loop
      if (l >= pl and l < pu) then
        maxnrm := Create(0.0);       -- locate column with largest norm and
        maxj := l;                        -- bring it in the pivot position
        for j in l..pu loop
          fltacc1 := REAL_PART(qraux(j));
          if fltacc1 > maxnrm
           then Copy(fltacc1,maxnrm); maxj := j;
          end if;
          Clear(fltacc1);
        end loop;
        if maxj /= l then
          zswap(x,l,maxj);
          Copy(qraux(l),qraux(maxj));
          Copy(work(l),work(maxj));
          jp := jpvt(maxj);
          jpvt(maxj) := jpvt(l);
          jpvt(l) := jp;
        end if;
        Clear(maxnrm);
      end if;
      Clear(qraux(l));
      qraux(l) := Create(integer(0));
      if l /= n then
        fltacc1 := znrm2(x,l,l);
        nrmxl := Create(fltacc1);               -- householder transformation
        if not Equal(fltacc1,0.0) then                        -- for column l
          fltacc2 := cdabs(x(l,l));
          if not Equal(fltacc2,0.0) then
            cmpacc := csign(nrmxl,x(l,l));
            Copy(cmpacc,nrmxl); Clear(cmpacc);
          end if;
          Clear(fltacc2);
          one := Create(integer(1));
          cmpacc := one/nrmxl;
          zscal(x,cmpacc,l,l);
          Clear(cmpacc);
          Add(x(l,l),one);
          Clear(one);
          lp1 := l + 1;        --  apply the transformation to the remaining
          for j in lp1..p loop               --  columns, updating the norms
            cmpacc := zdot(x,l,l,j);
            Min(cmpacc);
            t := cmpacc/x(l,l);
            Clear(cmpacc);
            zaxpy(x,t,l,l,j);
            Clear(t);
            fltacc2 := AbsVal(qraux(j));
            if (j >= pl) and (j <= pu) and (not Equal(fltacc2,0.0)) then
              fltacc3 := cdabs(x(l,j));
              fltacc4 := REAL_PART(qraux(j));
              Div(fltacc3,fltacc4);
              Clear(fltacc4);
              fltacc4 := fltacc3*fltacc3;
              tt := 1.0 - fltacc4;
              Clear(fltacc3); Clear(fltacc4);
              if tt < 0.0
               then Clear(tt); tt := Create(0.0);
              end if;
              t := Create(tt);
              fltacc3 := REAL_PART(qraux(j));
              fltacc4 := REAL_PART(work(j));
              Div(fltacc3,fltacc4);
              Clear(fltacc4);
              fltacc4 := fltacc3*fltacc3;
              Mul(tt,fltacc4);
              Clear(fltacc3); Clear(fltacc4);
              Mul(tt,0.05);
              Add(tt,1.0);
              if not Equal(tt,1.0) then
                fltacc3 := REAL_PART(t);
                fltacc4 := SQRT(fltacc3);
                cmpacc := Create(fltacc4);
                Mul(qraux(j),cmpacc);
                Clear(fltacc3); Clear(fltacc4);
                Clear(cmpacc);
              else 
                fltacc3 := znrm2(x,l+1,j);
                qraux(j) := Create(fltacc3);
                Clear(fltacc3);
                Copy(qraux(j),work(j));
              end if;
              Clear(t); Clear(tt);
            end if;
            Clear(fltacc2);
          end loop;
          Copy(x(l,l),qraux(l));                 -- save the transformation
          Copy(nrmxl,x(l,l));
          Min(x(l,l));
        end if;
        Clear(nrmxl);
        Clear(fltacc1);
      end if;
    end loop;
    Multprec_Complex_Vectors.Clear(work);
  end QRD;

  procedure Basis ( qr : in out Multprec_Complex_Matrices.Matrix;
                    x : in Multprec_Complex_Matrices.Matrix ) is

    sum,acc : Complex_Number;
    wrk : Multprec_Complex_Vectors.Vector(qr'range(1));

  begin
    for j in x'range(2) loop               -- compute jth column of q
      for i in qr'range(1) loop
        Copy(x(i,j),sum);
        for k in qr'first(2)..(j-1) loop
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

  procedure QRLS ( x : in out Multprec_Complex_Matrices.Matrix;
                   ldx,n,k : in integer32;
                   qraux,y : in Multprec_Complex_Vectors.Vector;
                   qy,qty,b,rsd,xb : out Multprec_Complex_Vectors.Vector;
                   job : in integer32; info : out integer32 ) is

    cb,cqy,cqty,cr,cxb : boolean;
    jj,ju,kp1 : integer32;
    fltacc : Floating_Number;
    t,temp,cmpacc : Complex_Number;

  begin
    info := 0;                                               -- set info flag
    cqy := (job/10000 /= 0);                     -- determine what to compute
    cqty := (job mod 10000 /= 0);
    cb := ((job mod 1000)/100 /= 0);
    cr := ((job mod 100)/10 /= 0);
    cxb := ((job mod 10) /= 0);
    ju := min0(k,n-1);
    if ju = 0 then                               -- special action when n = 1
      if cqy then Copy(y(1),qy(1)); end if;
      if cqty then Copy(y(1),qty(1)); end if;
      if cxb then Copy(y(1),xb(1)); end if;
      if cb then
        fltacc := AbsVal(x(1,1));
        if Equal(fltacc,0.0)
         then info := 1;
         else Clear(b(1)); b(1) := y(1)/x(1,1);
        end if;
        Clear(fltacc);
      end if;
      if cr
       then Clear(rsd(1)); rsd(1) := Create(integer(0));
      end if;
      return;
    end if;
    if cqy                                     -- set up to compute qy or qty
     then zcopy(n,y'first,y,qy);
    end if;
    if cqty
     then zcopy(n,y'first,y,qty);
    end if;
    if cqy then                                                 -- compute qy
      for j in 1..ju loop
        jj := ju - j + 1;
        fltacc := AbsVal(qraux(jj));
        if not Equal(fltacc,0.0) then
          Copy(x(jj,jj),temp);
          Copy(qraux(jj),x(jj,jj));
          cmpacc := zdotc(jj,x,qy);
          Min(cmpacc);
          t := cmpacc/x(jj,jj);
          Clear(cmpacc);
          zaxpy(n-jj+1,jj,jj,t,x,qy);
          Clear(t);
          Copy(temp,x(jj,jj));
          Clear(temp);
        end if;
        Clear(fltacc);
      end loop;
    end if;
    if cqty then                                        -- compute trans(q)*y
      for j in 1..ju loop
        fltacc := AbsVal(qraux(j));
        if not Equal(fltacc,0.0) then
          Copy(x(j,j),temp);
          Copy(qraux(j),x(j,j));
          cmpacc := zdotc(j,x,qty);
          Min(cmpacc);
          t := cmpacc/x(j,j);
          Clear(cmpacc);
          zaxpy(n-j+1,j,j,t,x,qty);
          Clear(t);
          Copy(temp,x(j,j));
          Clear(temp);
        end if;
        Clear(fltacc);
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
        Clear(xb(i));
        xb(i) := Create(integer(0));
      end loop;
    end if;
    if cr then
      for i in 1..k loop
        Clear(rsd(i));
        rsd(i) := Create(integer(0));
      end loop;
    end if;
    if cb then                                                   -- compute b
      for j in 1..k loop
        jj := k - j + 1;
        fltacc := AbsVal(x(jj,jj));
        if Equal(fltacc,0.0)
         then info := jj;
        end if;
        Clear(fltacc);
        exit when (info /= 0);
        Div(b(jj),x(jj,jj));
        if jj /= 1 then
          t := -b(jj);
          zaxpy(jj-1,1,jj,t,x,b);
          Clear(t);
        end if;
      end loop;
    end if;
    if cr or cxb then                        -- compute rsd or xb as requested
      for j in 1..ju loop
        jj := ju - j + 1;
        fltacc := AbsVal(qraux(jj));
        if not Equal(fltacc,0.0) then
          Copy(x(jj,jj),temp);
          Copy(qraux(jj),x(jj,jj));
          if cr then
            cmpacc := zdotc(jj,x,rsd);
            Min(cmpacc);
            t := cmpacc/x(jj,jj);
            Clear(cmpacc);
            zaxpy(n-jj+1,jj,jj,t,x,rsd);
            Clear(t);
          end if;
          if cxb then
            cmpacc := zdotc(jj,x,xb);
            Min(cmpacc);
            t := cmpacc/x(jj,jj);
            Clear(cmpacc);
            zaxpy(n-jj+1,jj,jj,t,x,xb);
            Clear(t);
          end if;
          Copy(temp,x(jj,jj));
          Clear(temp);
        end if;
        Clear(fltacc);
      end loop;
    end if;
  end QRLS;

end Multprec_Complex_QR_Least_Squares;
