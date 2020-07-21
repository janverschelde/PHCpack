--with text_io,integer_io;                use text_io,integer_io;
--with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
--with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
--with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
--with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;

with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Complex_BLAS_Helpers;     use Standard_Complex_BLAS_Helpers;

package body Standard_Complex_Singular_Values is

  function Min0 ( a,b : integer32 ) return integer32 is

  -- DESCRIPTION : returns the minimum of a and b.
  --   Note that this should correspond to the fortran min0 function.

  begin
    if a < b
     then return a;
     else return b;
    end if;
  end Min0;

  procedure SVD ( x : in out Matrix; n,p : in integer32;
                  s,e : out Vector; u : out Matrix; v : out Matrix; 
                  job : in integer32; info : out integer32 ) is

    work : Vector(1..n);

  begin
    SVD(x,n,p,s,e,u,v,job,info,work);
  end SVD;

  procedure SVD ( x : in out Matrix; n,p : in integer32;
                  s,e : out Vector; u : out Matrix; v : out Matrix; 
                  job : in integer32; info : out integer32;
                  work : in out Vector ) is

    iter,jobu,kase,kk,ll,lm1,lp1,ls,lu,m,maxit : integer32;
    mm,mm1,mp1,nct,nctp1,ncu,nrt,nrtp1 : integer32;
    t,r : Complex_Number;
    b,c,cs,el,emm1,f,g,scale,shift,sl,sm,sn : double_float;
    smm1,t1,test,ztest : double_float;
    wantu,wantv : boolean;

    cntkase2 : integer32 := 0; -- patch to fix with NaN inputs...
 
  begin
   -- put("entering SVD ...");
   -- put_line("the input matrix x : "); put(x);
    maxit := 30;     -- maximal number of iterations
    wantu := false;  -- determine what is to be computed
    wantv := false;
    jobu := (job mod 100)/10;
    ncu := n;
    if jobu > 1 then ncu := min0(n,p); end if;
    if jobu /= 0 then wantu := true; end if;
    if job mod 10 /= 0 then wantv := true; end if;
   -- reduce x to bidiagonal form, storing the diagonal elements in s
   -- and the super diagonal elements in e
    info := 0;
    nct := min0(n-1,p);
    nrt := max0(0,min0(p-2,n));
    lu := max0(nct,nrt);
    if lu >= 1 then
      for l in 1..lu loop
        lp1 := l+1;
        if l <= nct then
          -- compute the transformation for the l-th column
          -- and place the l-th diagonal in s(l)
          s(l) := Create(dznrm2(n-l+1,x,l,l,1),0.0);
          -- put("s("); put(l,1); put(") = "); put(s(l)); new_line;
          if cabs1(s(l)) /= 0.0 then
            if cdabs(x(l,l)) /= 0.0
             then s(l) := csign(s(l),x(l,l));
            end if;
            -- put("s("); put(l,1); put(") = "); put(s(l)); new_line;
            -- put("inverse : "); put(Create(1.0)/s(l)); new_line;
            zscal(n-l+1,Create(1.0)/s(l),x,l,l,1);
            x(l,l) := Create(1.0) + x(l,l);
            -- put("x("); put(l,1); put(","); put(l,1);
            -- put(") = "); put(x(l,l)); new_line;
          end if;
          s(l) := -s(l);
        end if;
        -- put_line("The matrix x : "); put(x);
        if p >= lp1 then
          for j in lp1..p loop
            if l <= nct then
              if (cabs1(s(l)) /= 0.0)
               then -- apply the transformation
                    t := -zdotc(n-l+1,x,l,l,1,x,l,j,1)/x(l,l);
                    zaxpy(n-l+1,t,x,l,l,1,x,l,j,1);
              end if;
            end if;
            -- place the l-th row of x into e for the subsequent
            -- calculation of the row transformation
            e(j) := Conjugate(x(l,j));
          end loop;
        end if;
        -- put_line("The matrix x : "); put(x);
        if wantu and l <= nct then
          -- place the transformation in u for subsequent
          -- back multiplication
          for i in l..n loop
            u(i,l) := x(i,l);
          end loop;
        end if;
        if l <= nrt then
           -- compute the l-th row transformation
           -- and place the l-th super diagonal in e(l)
          e(l) := Create(dznrm2(p-l,e,lp1,1),0.0);
          if cabs1(e(l)) /= 0.0 then
            if cdabs(e(lp1)) /= 0.0
             then e(l) := csign(e(l),e(lp1));
            end if;
            zscal(p-l,Create(1.0)/e(l),e,lp1,1);
            e(lp1) := Create(1.0) + e(lp1);
          end if;
          e(l) := -Conjugate(e(l));
          if lp1 <= n and cabs1(e(l)) /= 0.0 then
             -- apply the transformation
            for i in lp1..n loop
              work(i) := Create(0.0);
            end loop;
            for j in lp1..p loop
              zaxpy(n-l,e(j),x,lp1,j,1,work,lp1,1);
            end loop;
            for j in lp1..p loop
              zaxpy(n-l,Conjugate(-e(j)/e(lp1)),work,lp1,1,x,lp1,j,1);
            end loop;
          end if;
          if wantv then
            -- place the transformation in v
            -- for subsequent back multiplication
            for i in lp1..p loop
              v(i,l) := e(i);
            end loop;
          end if;
        end if;
      end loop;
    end if;
   -- set up the final bidiagonal matrix of order m
    m := min0(p,n+1);
    nctp1 := nct+1;
    nrtp1 := nrt+1;
    if nct < p then s(nctp1) := x(nctp1,nctp1); end if;
    if n < m then s(m) := Create(0.0); end if;
    if nrtp1 < m then e(nrtp1) := x(nrtp1,m); end if;
    e(m) := Create(0.0);
   -- if required, generate u
    if wantu then
      if ncu >= nctp1 then
        for j in nctp1..ncu loop
          for i in 1..n loop
            u(i,j) := Create(0.0);
          end loop;
          u(j,j) := Create(1.0);
        end loop;
      end if;
      if nct >= 1 then
        for l in 1..nct loop
          ll := nct - l + 1;
          if cabs1(s(ll)) = 0.0 then
            for i in 1..n loop
              u(i,ll) := Create(0.0);
            end loop;
            u(ll,ll) := Create(1.0);
          else
            lp1 := ll + 1;
            if ncu >= lp1 then
              for j in lp1..ncu loop
                t := -zdotc(n-ll+1,u,ll,ll,1,u,ll,j,1)/u(ll,ll);
                zaxpy(n-ll+1,t,u,ll,ll,1,u,ll,j,1);
              end loop;
            end if;
            zscal(n-ll+1,Create(-1.0),u,ll,ll,1);
            u(ll,ll) := Create(1.0) + u(ll,ll);
            lm1 := ll - 1;
            if lm1 >= 1 then
              for i in 1..lm1 loop
                u(i,ll) := Create(0.0);
              end loop;
            end if;
          end if;
        end loop;
      end if;
    end if;
   -- put_line("The matrix u : "); put(u);
   -- if required, generate v
    if wantv then
      for l in 1..p loop
        ll := p - l + 1;
        lp1 := ll + 1;
        if ll <= nrt then
          if cabs1(e(ll)) /= 0.0 then
            for j in lp1..p loop
              t := -zdotc(p-ll,v,lp1,ll,1,v,lp1,j,1)/v(lp1,ll);
              zaxpy(p-ll,t,v,lp1,ll,1,v,lp1,j,1);
            end loop;
          end if;
        end if;
        for i in 1..p loop
          v(i,ll) := Create(0.0);
        end loop;
        v(ll,ll) := Create(1.0);
      end loop;
    end if;
   -- put_line("The matrix v = "); put(v);
   -- transform s and e so that they are double precision
    for i in 1..m loop
      if cabs1(s(i)) /= 0.0 then
        t := Create(cdabs(s(i)),0.0);
        r := s(i)/t;
        s(i) := t;
        if i < m then e(i) := e(i)/r; end if;
        if wantu then zscal(n,r,u,1,i,1); end if;
      end if;
      exit when (i = m);
      if cabs1(e(i)) /= 0.0 then
        t := Create(cdabs(e(i)),0.0);
        r := t/e(i);
        e(i) := t;
        s(i+1) := s(i+1)*r;
        if wantv then zscal(p,r,v,1,i+1,1); end if;
      end if;
    end loop;
   -- main iteration loop for the singular values
    mm := m;
    iter := 0;
    loop
     -- put("inside SVD loop with m = "); put(m,1);
     -- put(" and iter = "); put(iter,1); new_line;
      exit when m = 0; -- quit if all the singular values have been found
      if iter > maxit  -- too many iterations have been performed
       then info := m; -- set the flag
            exit;      -- return
      end if;
     -- This section of the program inspects for negligible elements in 
     -- the s and e arrays.  On completion the variables kase and l are
     -- set as follows:
     --   kase = 1     if s(m) and e(l-1) are negligible and l.lt.m
     --   kase = 2     if s(l) is negligible and l.lt.m
     --   kase = 3     if e(l-1) is negligible, l.lt.m, and
     --                  s(l), ..., s(m) are not negligible (qr step).
     --   kase = 4     if e(m-1) is negligible (convergence).
      for l in 1..m loop
        ll := m - l;
        exit when (ll = 0);
        test := cdabs(s(ll)) + cdabs(s(ll+1));
        ztest := test + cdabs(e(ll));
        if ztest = test
         then e(ll) := Create(0.0); exit;
        end if;
      end loop;
      if ll = m - 1 then
        kase := 4;
      else
        lp1 := ll + 1;
        mp1 := m + 1;
        for lls in lp1..mp1 loop
          ls := m - lls + lp1;
          exit when (ls = ll);
          test := 0.0;
          if ls /= n then test := test + cdabs(e(ls)); end if;
          if ls /= ll+1 then test := test + cdabs(e(ls-1)); end if;
          ztest := test + cdabs(s(ls));
          if ztest = test 
           then s(ls) := Create(0.0); exit;
          end if;
        end loop;
        if ls = ll then
          kase := 3;
        elsif ls = m then
          kase := 1;
        else
          kase := 2;
          ll := ls;
        end if;
      end if;
      ll := ll + 1;
     -- perform the task indicated by kase
     -- put(" kase = "); put(kase,1); 
      case kase is
        when 1 => -- deflate negligible s(m)
          mm1 := m-1;
          f := REAL_PART(e(m-1));
          e(m-1) := Create(0.0);
          for k in ll..mm1 loop
            kk := mm1 - k + ll;
            t1 := REAL_PART(s(kk));
            cs := 0.0; sn := 0.0;
            drotg(t1,f,cs,sn);
            s(kk) := Create(t1);
            if kk /= ll
             then f := -sn*REAL_PART(e(kk-1));
                  e(kk-1) := cs*e(kk-1);
            end if;
            if wantv
             then zdrot(p,v,1,kk,1,v,1,m,1,cs,sn);
            end if;
          end loop;
        when 2 => -- split at negligible s(ll)
          f := REAL_PART(e(ll-1));
          e(ll-1) := Create(0.0);
          for k in ll..m loop
            t1 := REAL_PART(s(k));
            drotg(t1,f,cs,sn);
            s(k) := Create(t1);
            f := -sn*REAL_PART(e(k));
            e(k) := cs*e(k);
            if wantu
             then zdrot(n,u,1,k,1,u,1,ll-1,1,cs,sn);
            end if;
          end loop;
          cntkase2 := cntkase2 + 1;  -- patch for fixing NaN inputs ...
          if cntkase2 > 100 then return; end if;
        when 3 => -- perform one qr step
                  -- 1) calculate the shift
           scale := dmax1(cdabs(s(m)),cdabs(s(m-1)),cdabs(e(m-1)),
                          cdabs(s(ll)),cdabs(e(ll)));
           sm := REAL_PART(s(m))/scale;
           smm1 := REAL_PART(s(m-1))/scale;
           emm1 := REAL_PART(e(m-1))/scale;
           sl := REAL_PART(s(ll))/scale;
           el := REAL_PART(e(ll))/scale;
           b := ((smm1 + sm)*(smm1 - sm) + emm1**2)/2.0;
           c := (sm*emm1)**2;
           shift := 0.0;
           if b = 0.0 or c = 0.0 then
             shift := SQRT(b**2+c);
             if b < 0.0 then shift := -shift; end if;
             shift := c/(b + shift);
           end if;
           f := (sl + sm)*(sl - sm) + shift;
           g := sl*el;
          -- 2) chase zeros
           mm1 := m - 1;
           for k in ll..mm1 loop
             drotg(f,g,cs,sn);
             if k /= ll then e(k-1) := Create(f); end if;
             f := cs*REAL_PART(s(k)) + sn*REAL_PART(e(k));
             e(k) := cs*e(k) - sn*s(k);
             g := sn*REAL_PART(s(k+1));
             s(k+1) := cs*s(k+1);
             if wantv
              then zdrot(p,v,1,k,1,v,1,k+1,1,cs,sn);
             end if;
             drotg(f,g,cs,sn);
             s(k) := Create(f);
             f := cs*REAL_PART(e(k)) + sn*REAL_PART(s(k+1));
             s(k+1) := -sn*e(k) + cs*s(k+1);
             g := sn*REAL_PART(e(k+1));
             e(k+1) := cs*e(k+1);
             if wantu and k < n
              then zdrot(n,u,1,k,1,u,1,k+1,1,cs,sn);
             end if;
           end loop;
           e(m-1) := Create(f);
           iter := iter + 1;
        when 4 => -- convergence
                  -- 1) make the singular value positive
          if REAL_PART(s(ll)) < 0.0 then
            s(ll) := -s(ll);
            if wantv 
             then zscal(p,Create(-1.0),v,1,ll,1);
            end if;
          end if;
         -- 2) order the singular values
          while ll /= mm loop
            exit when REAL_PART(s(ll)) >= REAL_PART(s(ll+1));
            t := s(ll);
            s(ll) := s(ll+1);
            s(ll+1) := t;
            if wantv and ll < p
             then zswap(p,v,1,ll,1,v,1,ll+1,1);
            end if;
            if wantu and ll < n
             then zswap(n,u,1,ll,1,u,1,ll+1,1);
            end if;
            ll := ll+ 1;
          end loop;
          iter := 0;
          m := m-1;
        when others => null;
      end case;
    end loop;
   -- put_line("... leaving SVD");
--  exception
--    when others 
--      => put_line("exception caught by SVD");
--         raise;
  end SVD;

  function Rank ( s : Vector ) return integer32 is
  begin
    for i in s'range loop 
      if AbsVal(s(i)) + 1.0 = 1.0
       then return i-1;
      end if;
    end loop;
    return s'length;
  end Rank;

  function Rank ( s : Vector; tol : double_float ) return integer32 is
  begin
    for i in s'range loop 
      if AbsVal(s(i)) < tol
       then return i-1;
      end if;
    end loop;
    return s'length;
  end Rank;

  function Inverse_Condition_Number
             ( s : Vector ) return double_float is

    smax : constant double_float := AbsVal(s(s'first));
    smin,val : double_float;

  begin
    if smax + 1.0 = 1.0 then
      return 0.0;
    else
      smin := smax;
      for i in s'first+1..s'last loop
        val := AbsVal(s(i));
        exit when (val + 1.0 = 1.0);
        smin := val;
      end loop;
      return smin/smax;
    end if;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( s : Vector; tol : double_float ) return double_float is

    smax : constant double_float := AbsVal(s(s'first));
    smin,val : double_float;

  begin
    if smax < tol then
      return 0.0;
    else
      smin := smax;
      for i in s'first+1..s'last loop
        val := AbsVal(s(i));
        exit when (val < tol);
        smin := val;
      end loop;
      return smin/smax;
    end if;
  end Inverse_Condition_Number;

  function Conjugate_Transpose ( z : Matrix ) return Matrix is

    res : Matrix(z'range(2),z'range(1));

  begin
    for i in z'range(1) loop
      for j in z'range(2) loop
        res(j,i) := Conjugate(z(i,j));
      end loop;
    end loop;
    return res;
  end Conjugate_Transpose;

  function Inverse ( u,v : Matrix; s : Vector ) return Matrix is

    ut : constant Matrix(u'range(2),u'range(1)) := Conjugate_Transpose(u);
    su : Matrix(u'range(2),u'range(1));

  begin
    for i in s'range loop
      exit when (AbsVal(s(i)) + 1.0 = 1.0);
      for j in ut'range(2) loop
        su(i,j) := ut(i,j)/s(i);
      end loop;
    end loop;
    return v*su;
  end Inverse;

  function Inverse ( u,v : Matrix; s : Vector; tol : double_float )
                   return Matrix is

    ut : constant Matrix(u'range(2),u'range(1)) := Conjugate_Transpose(u);
    su : Matrix(u'range(2),u'range(1));

  begin
    for i in s'range loop
      exit when (AbsVal(s(i)) < tol);
      for j in ut'range(2) loop
        su(i,j) := ut(i,j)/s(i);
      end loop;
    end loop;
    return v*su;
  end Inverse;

  function Solve ( u,v : Matrix; s,b : Vector ) return Vector is

    ut : constant Matrix(u'range(2),u'range(1)) := Conjugate_Transpose(u);
    utb : constant Vector(u'range(2)) := ut*b;
    sub : Vector(v'range(1)) := (v'range(1) => Create(0.0));

  begin
    for i in s'range loop
      exit when (AbsVal(s(i)) + 1.0 = 1.0);
      exit when ((i > sub'last) or (i > utb'last));
      sub(i) := utb(i)/s(i);
    end loop;
    return v*sub;
  end Solve;

  procedure Solve ( ut,v : in Matrix; s,b : in Vector;
                    utb,sub : in out Vector; sol : out Vector ) is
  begin
    utb := ut*b;
    sub := (v'range(1) => Create(0.0));
    for i in s'range loop
      exit when (AbsVal(s(i)) + 1.0 = 1.0);
      exit when ((i > sub'last) or (i > utb'last));
      sub(i) := utb(i)/s(i);
    end loop;
    sol := v*sub;
  end Solve;

  function Solve ( u,v : Matrix; s,b : Vector; tol : double_float )
                 return Vector is

    ut : constant Matrix(u'range(2),u'range(1)) := Conjugate_Transpose(u);
    utb : constant Vector(u'range(2)) := ut*b;
    sub : Vector(v'range(1)) := (v'range(1) => Create(0.0));

  begin
    for i in s'range loop
      exit when (AbsVal(s(i)) < tol);
      sub(i) := utb(i)/s(i);
    end loop;
    return v*sub;
  end Solve;

end Standard_Complex_Singular_Values;
