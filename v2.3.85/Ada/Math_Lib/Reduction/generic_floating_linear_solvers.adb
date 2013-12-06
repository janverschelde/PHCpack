with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package body Generic_Floating_Linear_Solvers is

  use Field,Vectors;

-- AUXILIARIES :

  function csign ( x,y : number ) return number is

  -- DESCRIPTION :
  --   return AbsVal(x) * y / AbsVal(y);

    res : number;

  begin
    if y > zero
     then res := AbsVal(x);
     else res := AbsVal(x); Min(res);
    end if;
    return res;
  end csign;

  function Inverse_Abs_Sum ( z : Vectors.Vector ) return number is

  -- DESCRIPTION :
  --   Returns the reciprocal of the sum of the absolute values in z.

    res,sum,acc : number;

  begin
    sum := Create(0);
    for i in z'range loop
      acc := AbsVal(z(i));
      Add(sum,acc);
      Clear(acc);
    end loop;
    res := one/sum;
    Clear(sum);
    return res;
  end Inverse_Abs_Sum;

  procedure Swap ( a,b : in out number ) is

  -- DESCRIPTION :
  --   Swaps the two numbers by switching pointers.

    temp : number;

  begin
    temp := a;
    a := b;
    b := temp;
  end Swap;

-- STATIC TRIANGULATORS :

  procedure lufac ( a : in out matrix; n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    info : out integer32 ) is

    kp1,l,nm1 : integer32;
    smax,acc : number;

  begin
    info := 0;
    nm1 := integer32(n - 1);
    if nm1 >= 1 then
      for k in 1..nm1 loop
        kp1 := k + 1;
        l := k;                                  -- find the pivot index l
        smax := AbsVal(a(k,k));
        for i in kp1..n loop
          acc := AbsVal(a(i,k));
          if acc > smax 
           then l := i; Copy(acc,smax);
          end if;
          Clear(acc);
        end loop;
        ipvt(k) := l;
        if Equal(smax,zero) then
          info := k;               -- this column is already triangulated
        else
          if l /= k                           -- interchange if necessary
           then Swap(a(l,k),a(k,k));
          end if;                                  -- compute multipliers
          acc := one/a(k,k); 
          Min(acc);
          for i in kp1..n loop
            Mul(a(i,k),acc);
          end loop;
          Clear(acc);
          for j in kp1..n loop                         -- row elimination
            if l /= k
             then Swap(a(l,j),a(k,j));
            end if;
            for i in kp1..n loop
              acc := a(k,j)*a(i,k);
              Add(a(i,j),acc);
              Clear(acc);
            end loop;
          end loop;
        end if;
        Clear(smax);
      end loop;
    end if;
    ipvt(n) := n;
    if Equal(a(n,n),zero)
     then info := n;
    end if;
  end lufac;

  procedure lufco ( a : in out Matrix; n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    rcond : out number ) is

  -- NOTE :
  --   rcond = 1/(norm(a)*(estimate of norm(inverse(a))))
  --   estimate = norm(z)/norm(y) where a*z = y and ctrans(a)*y = e.
  --   ctrans(a) is the conjugate transpose of a.
  --   The components of e are chosen to cause maximum local
  --   growth in teh elements of w where ctrans(u)*w = e.
  --   The vectors are frequently rescaled to avoid overflow.

    z : Vectors.Vector(1..n);
    info,kb,kp1,l : integer32;
    s,sm,sum,anorm,ynorm,ek,t,wk,wkm,acc,absacc1,absacc2 : number;
    ipvtt : Standard_Integer_Vectors.Vector(1..n);

  begin
    anorm := Create(0);
    for j in 1..n loop                                -- compute 1-norm of a
      sum := Create(0);
      for i in 1..n loop
        acc := AbsVal(a(i,j));
        Add(sum,acc);
        Clear(acc);
      end loop;
      if sum > anorm
       then Copy(sum,anorm);
      end if; 
      Clear(sum);
    end loop;
    lufac(a,n,ipvtt,info);                                         -- factor
    for i in 1..n loop
      ipvt(i) := ipvtt(i);
    end loop;
    ek := Create(1);                                 -- solve ctrans(u)*w = e
    for j in 1..n loop
      z(j) := Create(0);
    end loop;
    for k in 1..n loop
      acc := AbsVal(z(k));
      if not Equal(acc,zero) then
        Clear(acc);
        acc := -z(k);
        absacc1 := csign(ek,acc);
        Copy(absacc1,ek);
        Clear(absacc1);
      end if;
      Clear(acc);
      acc := ek-z(k);
      absacc1 := AbsVal(acc);
      absacc2 := AbsVal(a(k,k));
      if absacc1 > absacc2 then
        s := absacc2/absacc1;
        Mul(z,s);
        Mul(ek,s);
        Clear(s);
      end if;
      Clear(absacc1); Clear(absacc2);
      Clear(acc);
      wk := ek - z(k);
      wkm := -ek;
      Sub(wkm,z(k));
      s := AbsVal(wk);
      sm := AbsVal(wkm);
      acc := AbsVal(a(k,k));
      if Equal(acc,zero) then
        Copy(one,wk);
        Copy(one,wkm);
      else
        Div(wk,a(k,k)); --wk := wk / a(k,k);
        Div(wkm,a(k,k)); --wkm := wkm / a(k,k);
      end if;
      Clear(acc);
      kp1 := k + 1;
      if kp1 <= n then
        for j in kp1..n loop
          acc := wkm*a(k,j);
          Add(acc,z(j));
          absacc1 := AbsVal(acc);
          Add(sm,absacc1);
          Clear(acc); Clear(absacc1);
          acc := wk*a(k,j);
          Add(z(j),acc);
          Clear(acc);
          absacc1 := AbsVal(z(j));
          Add(s,absacc1);
          Clear(absacc1);
        end loop;
        if s < sm then
          t := wkm - wk;
          Copy(wkm,wk);
          for j in kp1..n loop
            acc := t*a(k,j);
            Add(z(j),acc);
            Clear(acc);
          end loop;
          Clear(t);
        end if;
      end if;
      Copy(wk,z(k)); Clear(wk); Clear(wkm);
      Clear(s); Clear(sm);
    end loop;
    Clear(ek);
    s := Inverse_Abs_Sum(z);
    Mul(z,s); Clear(s);
    for k in 1..n loop                              -- solve ctrans(l)*y = w
      kb := n+1-k;
      if kb < n then
        Copy(zero,t);
        for i in (kb+1)..n loop
          acc := a(i,kb)*z(i);
          Add(t,acc);
          Clear(acc);
        end loop;
        Add(z(kb),t); Clear(t);
      end if;
      acc := AbsVal(z(kb));
      if acc > one 
       then s := one / acc; Mul(z,s); Clear(s);
      end if;
      Clear(acc);
      l := ipvtt(kb);
      if l /= kb then
        Copy(z(l),t);
        Copy(z(kb),z(l));
        Copy(t,z(kb)); Clear(t);
      end if;
    end loop;
    s := Inverse_Abs_Sum(z);
    Mul(z,s); Clear(s);
    ynorm := Create(1);
    for k in 1..n loop                             -- solve l*v = y
      l := ipvtt(k);
      if l /= k 
       then Copy(z(l),t); Copy(z(k),z(l)); Copy(t,z(k));
       else Copy(z(l),t);
      end if;
      if k < n then
        for i in (k+1)..n loop
          acc := t*a(i,k);
          Add(z(i),acc);
          Clear(acc);
        end loop;
      end if;
      Clear(t);
      acc := AbsVal(z(k));
      if acc > one
       then s := one / acc; Mul(z,s); Mul(ynorm,s); Clear(s);
      end if;
      Clear(acc);
    end loop;
    s := Inverse_Abs_Sum(z);
    Mul(z,s);
    Mul(ynorm,s); Clear(s);
    for k in 1..n loop                            -- solve u*z = v
      kb := n+1-k;
      absacc1 := AbsVal(z(kb));
      absacc2 := AbsVal(a(kb,kb));
      if absacc1 > absacc2 then
        s := absacc2 / absacc1;
        Mul(z,s);
        Mul(ynorm,s); Clear(s);
      end if;
      Clear(absacc1); 
      if Equal(absacc2,zero)
       then Copy(one,z(kb));
       else Div(z(kb),a(kb,kb));
      end if;
      Clear(absacc2);
      t := -z(kb);
      for i in 1..(kb-1) loop
        acc := t*a(i,kb);
        Add(z(i),acc);
        Clear(acc);
      end loop;
      Clear(t);
    end loop;
    s := Inverse_Abs_Sum(z);                             -- make znorm = 1.0
    Mul(z,s);
    Mul(ynorm,s); Clear(s);
    if Equal(anorm,zero)
     then rcond := Create(0);
     else rcond := ynorm/anorm;
    end if;
    Clear(z); Clear(ynorm); Clear(anorm);
  end lufco;

  procedure lusolve ( a : in matrix; n : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      b : in out Vectors.vector ) is

    l,nm1,kb : integer32;
    temp,acc : number;
 
  begin
    temp := Create(0);
    nm1 := n-1;
    if nm1 >= 1 then                                        -- solve l*y = b
      for k in 1..nm1 loop
        l := ipvt(k);
        Copy(b(l),temp);
        if l /= k 
         then Copy(b(k),b(l)); Copy(temp,b(k));
        end if;
        for i in (k+1)..n loop
          acc := temp*a(i,k);
          Add(b(i),acc);
          Clear(acc);
        end loop;
        Clear(temp);
      end loop;
    end if;
    for k in 1..n loop                                      -- solve u*x = y
      kb := n+1-k;
      Div(b(kb),a(kb,kb));
      temp := -b(kb);
      for j in 1..(kb-1) loop
        acc := temp*a(j,kb);
        Add(b(j),acc);
        Clear(acc);
      end loop;
      Clear(temp);
    end loop;
  end lusolve;

  procedure Triangulate ( a : in out Matrix; n,m : in integer32 ) is

    pivot,k,kcolumn : integer32;
    max,temp,acc : number;

  begin
    k := 1;
    kcolumn := 1;
    while (k <= n) and (kcolumn <= m) loop
      max := zero;                                            -- find pivot
      for l in k..n loop
        if AbsVal(a(l,kcolumn)) > max
         then max := AbsVal(a(l,kcolumn)); pivot := l;
        end if;
      end loop;
      if Equal(max,zero) then
        kcolumn := kcolumn + 1;
      else
        if pivot /= k then                      -- interchange if necessary
          for i in 1..m loop
            temp := a(pivot,i);
            a(pivot,i) := a(k,i);
            a(k,i) := temp;
          end loop;
        end if;
        for j in (kcolumn+1)..m loop                      -- triangulate a
          Div(a(k,j),a(k,kcolumn));
        end loop;
        Copy(one,a(k,kcolumn));
        for i in (k+1)..n loop
          for j in (kcolumn+1)..m loop
            acc := a(i,kcolumn)*a(k,j);
            Sub(a(i,j),acc);
            Clear(acc);
          end loop;
        end loop;
        for j in (k+1)..n loop
        Copy(zero,a(j,kcolumn));
        end loop;
        k := k + 1;
        kcolumn := kcolumn + 1;
      end if;
    end loop;
  end Triangulate;

  procedure Diagonalize ( a : in out Matrix; n,m : in integer32 ) is

    max,temp,acc : number;
    pivot,k,kcolumn : integer32;

  begin
    k := 1;
    kcolumn := 1;
    while (k <= n) and (kcolumn <= m) loop
      max := zero;                                              -- find pivot
      for l in k..n loop
        if AbsVal(a(l,kcolumn)) > max
         then max := AbsVal(a(l,kcolumn)); pivot := l;
        end if;
      end loop;
      if Equal(max,zero) then
        kcolumn := kcolumn + 1;
      else
        if pivot /= k then                        -- interchange if necessary
          for i in 1..m loop
            temp := a(pivot,i);
            a(pivot,i) := a(k,i);
            a(k,i) := temp;
          end loop;
        end if;
        for j in (kcolumn+1)..m loop                        -- diagonalize a
          Div(a(k,j),a(k,kcolumn));
        end loop;
        Copy(one,a(k,kcolumn));
        for i in 1..(k-1) loop
          for j in (kcolumn+1)..m loop
            acc := a(i,kcolumn)*a(k,j);
            Sub(a(i,j),acc);
            Clear(acc);
          end loop;
        end loop;
        for i in (k+1)..n loop
          for j in (kcolumn+1)..m loop
            acc := a(i,kcolumn)*a(k,j);
            Sub(a(i,j),acc);
            Clear(acc);
          end loop;
        end loop;
        for j in 1..(k-1) loop
          Copy(zero,a(j,kcolumn));
        end loop;
        for j in (k+1)..n loop
          Copy(zero,a(j,kcolumn));
        end loop;
        k := k + 1;
        kcolumn := kcolumn + 1;
      end if;
    end loop;
  end Diagonalize;

-- DYNAMIC TRIANGULATORS :

  procedure Upper_Triangulate
               ( row : in integer32; mat : in out Matrix; tol : in number;
                 ipvt : in out Standard_Integer_Vectors.Vector;
                 pivot : out integer32 ) is

   factor,tmp,max,acc : number;
   piv : integer32 := 0;
   tpi : integer32 := 0;

  begin
    for j in mat'first(1)..(row-1) loop
      if AbsVal(mat(row,j)) > tol then           -- make mat(row,j) zero
        factor := mat(row,j)/mat(j,j);
        for k in j..mat'last(2) loop
          acc := factor*mat(j,k);
          Sub(mat(row,k),acc);
          Clear(acc);
        end loop;
      end if;
    end loop;
    for j in row..ipvt'last loop                -- search pivot
      tmp := AbsVal(mat(row,j));
      if tmp > tol then
        if piv = 0 then
          max := tmp; piv := j;
        elsif tmp > max then
          max := tmp; piv := j;
        end if;
      end if;
    end loop;
    pivot := piv;
    if piv /= 0 then                            -- zero row
      if piv /= row  then                       -- interchange columns
        for k in mat'range(1) loop
          tmp := mat(k,row); mat(k,row) := mat(k,piv);
          mat(k,piv) := tmp;
        end loop;
        tpi := ipvt(row);
        ipvt(row) := ipvt(piv);
        ipvt(piv) := tpi;
      end if;
    end if;
  end Upper_Triangulate;

  procedure Upper_Triangulate
               ( roweli : in integer32; elim : in Matrix; tol : in number;
                 rowmat : in integer32; mat : in out Matrix ) is

    factor,acc : number;

  begin
    if AbsVal(mat(rowmat,roweli)) > tol then
      factor := mat(rowmat,roweli)/elim(roweli,roweli);
      for i in roweli..mat'last(2) loop
        acc := factor*elim(roweli,i);
        Sub(mat(rowmat,i),acc);
        Clear(acc);
      end loop;
    end if;
  end Upper_Triangulate;

  procedure Upper_Triangulate
               ( roweli : in integer32; elim : in Matrix; tol : in number;
                 firstrow,lastrow : in integer32; mat : in out Matrix ) is
  begin
    for i in firstrow..lastrow loop
      Upper_Triangulate(roweli,elim,tol,i,mat);
    end loop;
  end Upper_Triangulate;

  procedure Switch ( ipvt : in Standard_Integer_Vectors.Vector;
                     row : in integer32; mat : in out Matrix ) is

    tmp : Vectors.Vector(mat'range(2));

  begin
    for k in tmp'range loop
      tmp(k) := mat(row,k);
    end loop;
    for k in ipvt'range loop
      mat(row,k) := tmp(ipvt(k));
    end loop;
    for k in ipvt'last+1..mat'last(2) loop
      mat(row,k) := tmp(k);
    end loop;
  end Switch;

  procedure Switch ( k,pivot,first,last : in integer32; 
                     mat : in out Matrix ) is

    tmp : number;

  begin
    if k /= pivot then
      for i in first..last loop
        tmp := mat(i,k);
        mat(i,k) := mat(i,pivot);
        mat(i,pivot) := tmp;
      end loop;
    end if;
  end Switch;

  function Solve ( mat : Matrix; tol : number;
                   ipvt : Standard_Integer_Vectors.Vector )
                 return Vectors.Vector is

    res,x : Vectors.Vector(mat'range(2)) := (mat'range(2) => zero);
    index : integer32;
    acc : number;

  begin
    for i in mat'range(1) loop
      index := i;
      exit when i > mat'last(2);
      exit when AbsVal(mat(i,i)) < tol;
    end loop;
    if (AbsVal(mat(index,index)) > tol) and then (index < mat'last(2))
     then index := index + 1;
    end if;
    Copy(one,x(index));
    for i in reverse mat'first(1)..(index-1) loop
      x(i) := -mat(i,index);
      for j in i+1..index-1 loop
        acc := mat(i,j)*x(j);
        Sub(x(i),acc);
        Clear(acc);
      end loop;
      Div(x(i),mat(i,i));
    end loop;
    for k in ipvt'range loop
      res(ipvt(k)) := x(k);
    end loop;
    for k in ipvt'last+1..res'last loop
      res(k) := x(k);
    end loop;
    return res;
  end Solve;

  function Solve ( n,col : integer32; mat : Matrix )
                 return Vectors.Vector is

    res : Vectors.Vector(1..(n+1));
    acc : number;

  begin
    res(n+1) := Create(1);
    for i in reverse 1..n loop
      res(i) := -mat(i,col);
      for j in i+1..n loop
        acc := mat(i,j)*res(j);
        Sub(res(i),acc);
        Clear(acc);
      end loop;
      Div(res(i),mat(i,i));
    end loop;
    return res;
  end Solve;

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
                     A : Matrix ) return Matrix is

    fP : Matrix(P'range(1),P'range(2));

  begin
    for i in P'range(1) loop
      for j in P'range(2) loop
        fP(i,j) := Create(integer(P(i,j)));
      end loop;
    end loop;
    return fP*A;
  end Permute;

  function Lower_Diagonal ( A : Matrix ) return Matrix is

    res : Matrix(A'range(1),A'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if i < j then
          res(i,j) := Create(0);
        elsif i = j then
          res(i,j) := Create(1);
        else
          res(i,j) := -A(i,j);
        end if;
      end loop;
    end loop;
    return res;
  end Lower_Diagonal;

  function Upper_Diagonal ( A : Matrix ) return Matrix is

    res : Matrix(A'range(1),A'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if i <= j 
         then res(i,j) := A(i,j);
         else res(i,j) := Create(0);
        end if;
      end loop;
    end loop;
    return res;
  end Upper_Diagonal;

  procedure Permute_Lower
              ( L : in out Matrix;
                ipvt : in Standard_Integer_Vectors.Vector ) is

    tmp : number;

  begin
    for i in ipvt'range loop
      if ipvt(i) /= i then
        for j in 1..(i-1) loop
          Copy(L(i,j),tmp);
          Copy(L(ipvt(i),j),L(i,j));
          Copy(tmp,L(ipvt(i),j));
          Clear(tmp);
        end loop;
      end if;
    end loop;
  end Permute_Lower;

end Generic_Floating_Linear_Solvers;
