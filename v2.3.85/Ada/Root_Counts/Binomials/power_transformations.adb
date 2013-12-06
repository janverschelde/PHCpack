with Standard_Common_Divisors;          use Standard_Common_Divisors;

package body Power_Transformations is

  function Identity_Matrix ( n : natural ) return Matrix is

    res : Matrix(1..n,1..n);

  begin
    for i in 1..n loop
      for j in 1..n loop
        if i = j
         then res(i,j) := 1;
         else res(i,j) := 0;
        end if;
      end loop;
    end loop;
    return res;
  end Identity_Matrix;

  function Pivot ( v : Vector ) return integer is
  begin
    for i in v'range loop
      if v(i) /= 0
       then return i;
      end if;
    end loop;
    return v'first-1; 
  end Pivot;

  function Pivot ( v : Vector; i : integer ) return integer is
  begin
    for k in i..v'last loop
      if v(k) /= 0
       then return k;
      end if;
    end loop;
    return v'last + 1;
  end Pivot;

  function Rotate ( v : Vector; i : integer ) return Matrix is

    n : constant natural := v'last;
    t : Matrix(1..n,1..n) := Identity_Matrix(n);
    res : Matrix(1..n,1..n);
    w : Vector(v'range);
    j,a,b,ka,lb,d : integer;

  begin
    if v(i) = 0 then
      return t;
    else
      w := v; res := t;
      loop
        j := 0;
        for k in w'range loop
          if (w(k) /= 0) and (k /= i)
           then j := k;
          end if;
	  exit when j /= 0;
        end loop;
        if j /= 0 then
          a := w(i); b := w(j);
          gcd(a,b,ka,lb,d);
          a := a/d;  b := b/d;
          t(i,i) := ka; t(i,j) := lb;
          t(j,i) := -b; t(j,j) := a;
          w := t*w;
          res := t*res;
          t(i,i) := 1; t(j,j) := 1;
          t(i,j) := 0; t(j,i) := 0;
        end if;
        exit when j = 0;
      end loop;
      return res;
    end if;
  end Rotate;

  function Eliminate ( v : Vector; i : integer ) return Matrix is

    n : constant natural := v'last;
    t : Matrix(1..n,1..n) := Identity_Matrix(n);
    res : Matrix(1..n,1..n) := t;
    w : Vector(v'range) := v;
    j1,j2,a,b,ka,lb,d : integer;
    ind : integer := Pivot(v,i);
    zeroes : boolean;

  begin
    if ind = v'last + 1 then
      return res;
    else
      if v(i) = 0 then
        res(i,i) := 0;   res(i,ind) := 1;
        res(ind,i) := 1; res(ind,ind) := 0;
        w(i) := w(ind); w(ind) := 0;
      end if;
      for k in v'first..i loop
        if w(k) /= 0
         then j1 := k; exit;
        end if;
      end loop;
      if j1 /= i then
        zeroes := false;
        loop
          for k in (j1+1)..i loop
            if w(k) /= 0
             then j2 := k; exit;
            end if;
          end loop;
          a := w(j1); b := w(j2);
          gcd(a,b,ka,lb,d);
          a := a/d;  b := b/d;
          t(j1,j1) := lb; t(j1,j2) := -ka;
          t(j2,j1) := a;  t(j2,j2) := b;
          w(j2) := d;
          res := t*res;
          t(j1,j1) := 1; t(j2,j2) := 1;
          t(j1,j2) := 0; t(j2,j1) := 0;
          if j2 < i
           then j1 := j2;
           else exit;
          end if;
        end loop;
      else
        zeroes := true;
      end if;
      for k in reverse i..v'last loop
        if w(k) /= 0
         then j2 := k; exit;
        end if;
      end loop;
      if j2 /= i then
        loop
          for k in reverse i..(j2-1) loop
            if w(k) /= 0
             then j1 := k; exit;
            end if;
          end loop;
          a := w(j1); b := w(j2);
          gcd(a,b,ka,lb,d);
          a := a/d;  b := b/d;
          t(j1,j1) := a;   t(j1,j2) := b;
          t(j2,j1) := -lb; t(j2,j2) := ka;
          w(j1) := d;
          res := t*res;
          t(j1,j1) := 1; t(j2,j2) := 1;
          t(j1,j2) := 0; t(j2,j1) := 0;
          if j1 > i
	   then j2 := j1;
	   else exit;
          end if;
        end loop;
      elsif zeroes then
        if w(i) < 0 then
          t(i,i) := -1;
          res := t*res;
        end if;
      end if;
      return res;
    end if;
  end Eliminate;

end Power_Transformations;
