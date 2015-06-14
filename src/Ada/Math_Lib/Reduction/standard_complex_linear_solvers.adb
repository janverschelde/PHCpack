--with text_io; use text_io;
--with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;

with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package body Standard_Complex_Linear_Solvers is

-- lufac_count : natural32 := 0;

-- AUXLILIARIES :

  function cabs ( c : Complex_Number ) return double_float is
  begin
    return (ABS(REAL_PART(c)) + ABS(IMAG_PART(c)));
  end cabs;

  function dconjg ( x : Complex_Number ) return Complex_Number is
  begin
    return Create(REAL_PART(x),-IMAG_PART(x));
  end dconjg;

  function csign ( x,y : Complex_Number ) return Complex_Number is
  begin
    return (Create(cabs(x)) * y / Create(cabs(y)));
  end csign;

  procedure Swap ( x,y : in out Complex_Number ) is

  -- DESCRIPTION :
  --   Swaps x with y.

    z : constant Complex_Number := y;

  begin
    y := x;
    x := z;
  end Swap;

-- TARGET ROUTINES :

  procedure Scale ( a : in out Standard_Complex_Matrices.Matrix;
                    b : in out Standard_Complex_Vectors.Vector ) is

    fac : Complex_Number;

    function Maximum ( a : in Standard_Complex_Matrices.Matrix;
                       i : in integer32 ) return Complex_Number is

      res : integer32 := a'first(2);
      max : double_float := cabs(a(i,res));
      tmp : double_float;

    begin
      for j in a'first(2)+1..a'last(2) loop
        tmp := cabs(a(i,j));
        if tmp > max
         then max := tmp; res := j;
        end if;
      end loop;
      return a(i,res);
    end Maximum;

    procedure Divide ( a : in out Standard_Complex_Matrices.Matrix;
                       b : in out Standard_Complex_Vectors.Vector;
                       i : in integer32; fac : in Complex_Number ) is
    begin
      for j in a'range(2) loop
        a(i,j) := a(i,j)/fac;
      end loop;
      b(i) := b(i)/fac;
    end Divide;
  
  begin
    for i in a'range(1) loop
      fac := Maximum(a,i);
      Divide(a,b,i,fac);
    end loop;
  end Scale;

  function Norm1 ( a : Standard_Complex_Matrices.Matrix )
                 return double_float is

    res : double_float := 0.0;
    sum : double_float;

  begin
    for j in a'range(2) loop
      sum := 0.0;
      for i in a'range(1) loop
        sum := sum + cabs(a(i,j));
      end loop;
      if sum > res
       then res := sum;
      end if; 
    end loop;
    return res;
  end Norm1;

  function Norm1 ( a : Standard_Complex_VecVecs.VecVec )
                 return double_float is

    res : double_float := 0.0;
    sum : double_float;
    aj : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for j in a'range loop
      sum := 0.0;
      aj := a(j);
      for i in aj'range loop
        sum := sum + cabs(aj(i));
      end loop;
      if sum > res
       then res := sum;
      end if; 
    end loop;
    return res;
  end Norm1;

-- TARGET ROUTINES :

  procedure lufac ( a : in out Standard_Complex_Matrices.Matrix;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    info : out integer32 ) is

    kp1,L,nm1 : integer32;
    smax,ikabs : double_float;
    temp : Complex_Number;

  begin
    info := 0;
    nm1 := n - 1;
    if nm1 >= 1 then
      for k in 1..nm1 loop
        kp1 := k + 1;                              -- find the pivot index L
        L := k; smax := cabs(a(k,k));
        for i in kp1..n loop
          ikabs := cabs(a(i,k));
          if ikabs > smax then
            L := i;
            smax := ikabs;
          end if;
        end loop;
        ipvt(k) := L;
        if smax = 0.0 then          -- this column is already triangularized
          info := k;
        else
          if L /= k then                         -- interchange if necessary
            temp := a(L,k);
            a(L,k) := a(k,k);
            a(k,k) := temp;
          end if;
          temp := -Create(1.0)/a(k,k);                -- compute multipliers
          for i in kp1..n loop
            a(i,k) := temp*a(i,k);
          end loop;
          for j in kp1..n loop       -- row elimination with column indexing
            temp := a(L,j);
            if L /= k then
              a(L,j) := a(k,j);
              a(k,j) := temp;
            end if;
            for i in kp1..n loop
              a(i,j) := a(i,j) + temp*a(i,k);
            end loop;
          end loop;
        end if;
      end loop;
    end if;
    ipvt(n) := n;
    if AbsVal(a(n,n)) = 0.0
     then info := n;
    end if;
   -- lufac_count := lufac_count + 1;
   -- put("lufac count : "); put(lufac_count,1); new_line;
  end lufac;

  procedure lufac ( a : in out Standard_Complex_VecVecs.VecVec;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    info : out integer32 ) is

    kp1,ell,nm1 : integer32;
    smax,amax : double_float;
    acc : Complex_Number;
    ak,aj : Standard_Complex_Vectors.Link_to_Vector;
    zero : constant Complex_Number := Create(0.0);

  begin
    info := 0;
    nm1 := integer32(n - 1);
    if nm1 >= 1 then
      for k in 1..nm1 loop
        ak := a(k);
        kp1 := k + 1;
        ell := k;                               -- find the pivot index ell
        smax := cabs(ak(k));
        for i in kp1..n loop
          amax := cabs(ak(i));
          if amax > smax 
           then ell := i; smax := amax;
          end if;
        end loop;
        ipvt(k) := ell;
        if smax = 0.0 then
          info := k;               -- this column is already triangulated
        else
          if ell /= k then                     -- interchange if necessary
            acc := ak(ell);
            ak(ell) := ak(k);
            ak(k) := acc;
          end if;                                  -- compute multipliers
          acc := -Create(1.0)/ak(k); 
          for i in kp1..n loop
            ak(i) := ak(i)*acc;
          end loop;
          for j in kp1..n loop                         -- row elimination
            aj := a(j);
            acc := aj(ell);
            if ell /= k then
              aj(ell) := aj(k);
              aj(k) := acc;
            end if;
            for i in kp1..n loop
              aj(i) := aj(i) + acc*ak(i);
            end loop;
          end loop;
        end if;
      end loop;
    end if;
    ipvt(n) := n;
    if AbsVal(a(n)(n)) = 0.0
     then info := n;
    end if;
  end lufac;

  procedure estco ( a : in Standard_Complex_Matrices.Matrix;
                    n : in integer32;
                    ipvt : in Standard_Integer_Vectors.Vector;
                    anorm : in double_float; rcond : out double_float ) is

    z : Standard_Complex_Vectors.Vector(1..n);
    kb,kp1,l : integer32;
    s,sm,sum,ynorm : double_float;
    ek,t,wk,wkm : Complex_Number;

    use Standard_Complex_Vectors;

  begin
    ek := Create(1.0);                              -- solve ctrans(u)*w = e
    for j in 1..n loop
      z(j) := Create(0.0);
    end loop;
    for k in 1..n loop
      if cabs(z(k)) /= 0.0
       then ek := csign(ek,-z(k));
      end if;
      if cabs(ek-z(k)) > cabs(a(k,k)) then
        s := cabs(a(k,k))/cabs(ek-z(k));
        z := Create(s) * z;
        ek := Create(s) * ek;
      end if;
      wk := ek - z(k);
      wkm := -ek - z(k);
      s := cabs(wk);
      sm := cabs(wkm);
      if cabs(a(k,k)) = 0.0  then
        wk := Create(1.0);
        wkm := Create(1.0);
      else
        wk := wk / dconjg(a(k,k));
        wkm := wkm / dconjg(a(k,k));
      end if;
      kp1 := k + 1;
      if kp1 <= n then
        for j in kp1..n loop
          sm := sm + cabs(z(j)+wkm*dconjg(a(k,j)));
          z(j) := z(j) + wk*dconjg(a(k,j));
          s := s + cabs(z(j));
        end loop;
        if s < sm then
          t := wkm - wk;
          wk := wkm;
          for j in kp1..n loop
            z(j) := z(j) + t*dconjg(a(k,j));
          end loop;
        end if;
      end if;
      z(k) := wk;
    end loop;
    sum := 0.0;
    for i in 1..n loop
      sum := sum + cabs(z(i));
    end loop;
    s := 1.0 / sum;
    z := Create(s) * z;
    for k in 1..n loop                           -- solve ctrans(l)*y = w
      kb := n+1-k;
      if kb < n then
        t := Create(0.0);
        for i in (kb+1)..n loop
          t := t + dconjg(a(i,kb))*z(i);
        end loop;
        z(kb) := z(kb) + t;
      end if;
      if cabs(z(kb)) > 1.0 then
        s := 1.0 / cabs(z(kb));
        z := Create(s) * z;
      end if;
      l := ipvt(kb);
      t := z(l);
      z(l) := z(kb);
      z(kb)  := t;
    end loop;
    sum := 0.0;
    for i in 1..n loop
      sum := sum + cabs(z(i));
    end loop;
    s := 1.0 / sum;
    z := Create(s) * z;
    ynorm := 1.0;
    for k in 1..n loop                                    -- solve l*v = y
      l := ipvt(k);
      t := z(l);
      z(l) := z(k);
      z(k) := t;
      if k < n then
        for i in (k+1)..n loop
          z(i) := z(i) + t * a(i,k);
        end loop;
      end if;
      if cabs(z(k)) > 1.0 then
        s := 1.0 / cabs(z(k));
        z := Create(s) * z;
        ynorm := s * ynorm;
      end if;
    end loop;
    sum := 0.0;
    for i in 1..n loop
      sum := sum + cabs(z(i));
    end loop;
    s := 1.0 / sum;
    z := Create(s) * z;
    ynorm := s * ynorm;
    for k in 1..n loop                                    -- solve u*z = v
      kb := n+1-k;
      if cabs(z(kb)) > cabs(a(kb,kb)) then
        s := cabs(a(kb,kb)) / cabs(z(kb));
        z := Create(s) * z;
        ynorm := s * ynorm;
      end if;
      if cabs(a(kb,kb)) = 0.0
       then z(kb) := Create(1.0);
       else z(kb) := z(kb) / a(kb,kb);
      end if;
      t := -z(kb);
      for i in 1..(kb-1) loop
        z(i) := z(i) + t * a(i,kb);
      end loop;
    end loop;
    sum := 0.0;                                       -- make znorm = 1.0
    for i in 1..n loop
      sum := sum + cabs(z(i));
    end loop;
    s := 1.0 / sum;
   -- z := Create(s) * z; -- deemed useless by GNAT GPL 2009 compiler
    ynorm := s * ynorm;
    if anorm = 0.0
     then rcond := 0.0;
     else rcond := ynorm/anorm;
    end if;
  end estco;

  procedure estco ( a : in Standard_Complex_VecVecs.VecVec;
                    n : in integer32;
                    ipvt : in Standard_Integer_Vectors.Vector;
                    anorm : in double_float; rcond : out double_float ) is

    z : Standard_Complex_Vectors.Vector(1..n);
    ak,aj : Standard_Complex_Vectors.Link_to_Vector;
    kb,kp1,l : integer32;
    s,sm,sum,ynorm : double_float;
    ek,t,wk,wkm : Complex_Number;

    use Standard_Complex_Vectors;

  begin
    ek := Create(1.0);                              -- solve ctrans(u)*w = e
    for j in 1..n loop
      z(j) := Create(0.0);
    end loop;
    for k in 1..n loop
      ak := a(k);
      if cabs(z(k)) /= 0.0
       then ek := csign(ek,-z(k));
      end if;
      if cabs(ek-z(k)) > cabs(ak(k)) then
        s := cabs(ak(k))/cabs(ek-z(k));
        z := Create(s) * z;
        ek := Create(s) * ek;
      end if;
      wk := ek - z(k);
      wkm := -ek - z(k);
      s := cabs(wk);
      sm := cabs(wkm);
      if cabs(ak(k)) = 0.0  then
        wk := Create(1.0);
        wkm := Create(1.0);
      else
        wk := wk / dconjg(ak(k));
        wkm := wkm / dconjg(ak(k));
      end if;
      kp1 := k + 1;
      if kp1 <= n then
        for j in kp1..n loop
          aj := a(j);
          sm := sm + cabs(z(j)+wkm*dconjg(aj(k)));
          z(j) := z(j) + wk*dconjg(aj(k));
          s := s + cabs(z(j));
        end loop;
        if s < sm then
          t := wkm - wk;
          wk := wkm;
          for j in kp1..n loop
            aj := a(j);
            z(j) := z(j) + t*dconjg(aj(k));
          end loop;
        end if;
      end if;
      z(k) := wk;
    end loop;
    sum := 0.0;
    for i in 1..n loop
      sum := sum + cabs(z(i));
    end loop;
    s := 1.0 / sum;
    z := Create(s) * z;
    for k in 1..n loop                           -- solve ctrans(l)*y = w
      kb := n+1-k;
      if kb < n then
        ak := a(kb);
        t := Create(0.0);
        for i in (kb+1)..n loop
          t := t + dconjg(ak(i))*z(i);
        end loop;
        z(kb) := z(kb) + t;
      end if;
      if cabs(z(kb)) > 1.0 then
        s := 1.0 / cabs(z(kb));
        z := Create(s) * z;
      end if;
      l := ipvt(kb);
      t := z(l);
      z(l) := z(kb);
      z(kb)  := t;
    end loop;
    sum := 0.0;
    for i in 1..n loop
      sum := sum + cabs(z(i));
    end loop;
    s := 1.0 / sum;
    z := Create(s) * z;
    ynorm := 1.0;
    for k in 1..n loop                                    -- solve l*v = y
      l := ipvt(k);
      t := z(l);
      z(l) := z(k);
      z(k) := t;
      if k < n then
        ak := a(k);
        for i in (k+1)..n loop
          z(i) := z(i) + t * ak(i);
        end loop;
      end if;
      if cabs(z(k)) > 1.0 then
        s := 1.0 / cabs(z(k));
        z := Create(s) * z;
        ynorm := s * ynorm;
      end if;
    end loop;
    sum := 0.0;
    for i in 1..n loop
      sum := sum + cabs(z(i));
    end loop;
    s := 1.0 / sum;
    z := Create(s) * z;
    ynorm := s * ynorm;
    for k in 1..n loop                                    -- solve u*z = v
      kb := n+1-k;
      ak := a(kb);
      if cabs(z(kb)) > cabs(ak(kb)) then
        s := cabs(ak(kb)) / cabs(z(kb));
        z := Create(s) * z;
        ynorm := s * ynorm;
      end if;
      if cabs(ak(kb)) = 0.0
       then z(kb) := Create(1.0);
       else z(kb) := z(kb) / ak(kb);
      end if;
      t := -z(kb);
      for i in 1..(kb-1) loop
        z(i) := z(i) + t * ak(i);
      end loop;
    end loop;
    sum := 0.0;                                       -- make znorm = 1.0
    for i in 1..n loop
      sum := sum + cabs(z(i));
    end loop;
    s := 1.0 / sum;
   -- z := Create(s) * z; -- deemed useless by GNAT GPL 2009 compiler
    ynorm := s * ynorm;
    if anorm = 0.0
     then rcond := 0.0;
     else rcond := ynorm/anorm;
    end if;
  end estco;

  procedure lufco ( a : in out Standard_Complex_Matrices.Matrix;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    rcond : out double_float ) is

    anorm : constant double_float := Norm1(a);
    info : integer32;

  begin
    lufac(a,n,ipvt,info);
    if info = 0
     then estco(a,n,ipvt,anorm,rcond);
     else rcond := 0.0;
    end if;
  end lufco;

  procedure lufco ( a : in out Standard_Complex_VecVecs.VecVec;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    rcond : out double_float ) is

    anorm : constant double_float := Norm1(a);
    info : integer32;

  begin
    lufac(a,n,ipvt,info);
    if info = 0
     then estco(a,n,ipvt,anorm,rcond);
     else rcond := 0.0;
    end if;
  end lufco;

  procedure lusolve ( a : in Standard_Complex_Matrices.Matrix;
                      n : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      b : in out Standard_Complex_Vectors.Vector ) is

    ell,nm1,kb : integer32;
    temp : Complex_Number;
 
  begin
    nm1 := n-1;
    if nm1 >= 1 then                                       -- solve l*y = b
      for k in 1..nm1 loop
        ell := ipvt(k);
        temp := b(ell);
        if ell /= k then
          b(ell) := b(k);
          b(k) := temp;
        end if;
        for i in (k+1)..n loop
          b(i) := b(i) + temp*a(i,k);
        end loop;
      end loop;
    end if;
    for k in 1..n loop                                     -- solve u*x = y
      kb := n+1-k;
      b(kb) := b(kb) / a(kb,kb);
      temp := -b(kb);
      for j in 1..(kb-1) loop
        b(j) := b(j) + temp*a(j,kb);
      end loop;
    end loop;
  end lusolve;

  procedure lusolve ( a : in Standard_Complex_VecVecs.VecVec;
                      n : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      b : in out Standard_Complex_Vectors.Vector ) is

    ell,nm1,kb : integer32;
    temp : Complex_Number;
    ak : Standard_Complex_Vectors.Link_to_Vector;
 
  begin
    nm1 := n-1;
    if nm1 >= 1 then                                       -- solve l*y = b
      for k in 1..nm1 loop
        ell := ipvt(k);
        temp := b(ell);
        if ell /= k then
          b(ell) := b(k);
          b(k) := temp;
        end if;
        ak := a(k);
        for i in (k+1)..n loop
          b(i) := b(i) + temp*ak(i);
        end loop;
      end loop;
    end if;
    for k in 1..n loop                                     -- solve u*x = y
      kb := n+1-k;
      ak := a(kb);
      b(kb) := b(kb)/ak(kb);
      temp := -b(kb);
      for j in 1..(kb-1) loop
        b(j) := b(j) + temp*ak(j);
      end loop;
    end loop;
  end lusolve;

  procedure Triangulate ( a : in out Standard_Complex_Matrices.Matrix;
                          tol : in double_float;
                          n,m : in integer32 ) is

    max,cbs : double_float;
    temp : Complex_Number;
    pivot,k,kcolumn : integer32;

  begin
    k := 1;
    kcolumn := 1;
    while (k <= n) and (kcolumn <= m) loop
      max := 0.0;                                             -- find pivot
      pivot := 0;
      for l in k..n loop
        cbs := cabs(a(l,kcolumn));
        if (cbs > tol) and then (cbs > max) then
          max := cbs;
          pivot := l;
        end if;
      end loop;
      if pivot = 0 then
        kcolumn := kcolumn + 1;
      else
        if pivot /= k then                      -- interchange if necessary
          for i in 1..m loop
            temp := a(pivot,i);
            a(pivot,i) := a(k,i);
            a(k,i) := temp;
          end loop;
        end if;
        for j in (kcolumn+1)..m loop                       -- triangulate a
          a(k,j) := a(k,j) / a(k,kcolumn);
        end loop;
        a(k,kcolumn) := Create(1.0);
        for i in (k+1)..n loop
          for j in (kcolumn+1)..m loop
            a(i,j) := a(i,j) - a(i,kcolumn) * a(k,j);
          end loop;
          a(i,kcolumn) := Create(0.0);
        end loop;
        k := k + 1;
        kcolumn := kcolumn + 1;
      end if;
    end loop;
  end Triangulate;

  procedure Diagonalize ( a : in out Standard_Complex_Matrices.Matrix;
                          n,m : in integer32 ) is

    max : double_float;
    temp : Complex_Number;
    pivot,k,kcolumn : integer32;

  begin
    k := 1;
    kcolumn := 1;
    while (k <= n) and (kcolumn <= m) loop
      max := 0.0;                                             -- find pivot
      for l in k..n loop
        if cabs(a(l,kcolumn)) > max then
          max := cabs(a(l,kcolumn));
          pivot := l;
        end if;
      end loop;
      if max = 0.0 then
        kcolumn := kcolumn + 1;
      else
        if pivot /= k then                      -- interchange if necessary
          for i in 1..m loop
            temp := a(pivot,i);
            a(pivot,i) := a(k,i);
            a(k,i) := temp;
          end loop;
        end if;
        for j in (kcolumn+1)..m loop                       -- diagonalize a
          a(k,j) := a(k,j) / a(k,kcolumn);
        end loop;
        a(k,kcolumn) := Create(1.0);
        for i in 1..(k-1) loop
          for j in (kcolumn+1)..m loop
            a(i,j) := a(i,j) - a(i,kcolumn) * a(k,j);
          end loop;
        end loop;
        for i in (k+1)..n loop
          for j in (kcolumn+1)..m loop
            a(i,j) := a(i,j) - a(i,kcolumn) * a(k,j);
          end loop;
        end loop;
        for j in 1..(k-1) loop
          a(j,kcolumn) := Create(0.0);
        end loop;
        for j in (k+1)..n loop
          a(j,kcolumn) := Create(0.0);
        end loop;
        k := k + 1;
        kcolumn := kcolumn + 1;
      end if;
    end loop;
  end Diagonalize;

-- TO TEST THE LU FACTORIZATION :

  function Permutation_Matrix
              ( ipvt : Standard_Integer_Vectors.Vector )
              return Standard_Natural_Matrices.Matrix is

    res : Standard_Natural_Matrices.Matrix(ipvt'range,ipvt'range);
    tmp : natural32;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if i = j
         then res(i,j) := 1;
         else res(i,j) := 0;
        end if;
      end loop;
    end loop;
    for i in ipvt'range loop
      if ipvt(i) /= i then
        for j in res'range(2) loop
          tmp := res(ipvt(i),j);
          res(ipvt(i),j) := res(i,j);
          res(i,j) := tmp;
        end loop;
      end if;
    end loop;
    return res;
  end Permutation_Matrix;

  function Permute ( P : Standard_Natural_Matrices.Matrix;
                     A : Standard_Complex_Matrices.Matrix )
                   return Standard_Complex_Matrices.Matrix is

    fP : Standard_Complex_Matrices.Matrix(P'range(1),P'range(2));

    use Standard_Complex_Matrices;

  begin
    for i in P'range(1) loop
      for j in P'range(2) loop
        fP(i,j) := Create(P(i,j));
      end loop;
    end loop;
    return fP*A;
  end Permute;

  function Lower_Diagonal ( A : Standard_Complex_Matrices.Matrix )
                          return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if i < j then
          res(i,j) := Create(integer(0));
        elsif i = j then
          res(i,j) := Create(integer(1));
        else
          res(i,j) := -A(i,j);
        end if;
      end loop;
    end loop;
    return res;
  end Lower_Diagonal;

  function Upper_Diagonal ( A : Standard_Complex_Matrices.Matrix )
                          return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if i <= j 
         then res(i,j) := A(i,j);
         else res(i,j) := Create(integer(0));
        end if;
      end loop;
    end loop;
    return res;
  end Upper_Diagonal;

  procedure Permute_Lower
              ( L : in out Standard_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector ) is

    tmp : Complex_Number;

  begin
    for i in ipvt'range loop
      if ipvt(i) /= i then
        for j in 1..(i-1) loop
          tmp := L(i,j);
          L(i,j) := L(ipvt(i),j);
          L(ipvt(i),j) := tmp;
        end loop;
      end if;
    end loop;
  end Permute_Lower;

--begin
--  lufac_count := 0;
end Standard_Complex_Linear_Solvers;
