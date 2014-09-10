-- with text_io;                            use text_io;
-- with Multprec_Integer_Numbers_io;        use Multprec_Integer_Numbers_io;

with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Common_Divisors;           use Multprec_Common_Divisors;

package body Multprec_Power_Transformations is

  function Identity_Matrix
              ( n : natural32 ) return Multprec_Integer_Matrices.Matrix is

    res : Multprec_Integer_Matrices.Matrix(1..integer32(n),1..integer32(n));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if i = j
         then res(i,j) := Multprec_Integer_Numbers.Create(integer32(1));
         else res(i,j) := Multprec_Integer_Numbers.Create(integer32(0));
        end if;
      end loop;
    end loop;
    return res;
  end Identity_Matrix;

  function Pivot ( v : Multprec_Integer_Vectors.Vector ) return integer32 is
  begin
    for i in v'range loop
      if not Equal(v(i),0)
       then return i;
      end if;
    end loop;
    return v'first-1; 
  end Pivot;

  function Pivot ( v : Multprec_Integer_Vectors.Vector;
                   i : integer32 ) return integer32 is
  begin
    for k in i..v'last loop
      if not Equal(v(k),0)
       then return k;
      end if;
    end loop;
    return v'last + 1;
  end Pivot;

  function Rotate ( v : Multprec_Integer_Vectors.Vector; i : integer32 )
                  return Multprec_Integer_Matrices.Matrix is

    n : constant integer32 := v'last;
    t : Multprec_Integer_Matrices.Matrix(1..n,1..n)
      := Identity_Matrix(natural32(n));
    res : Multprec_Integer_Matrices.Matrix(1..n,1..n);
    w : Multprec_Integer_Vectors.Vector(v'range);
    j : integer32;
    a,b,ka,lb,d : Integer_Number;

  begin
    if Equal(v(i),0) then
      return t;
    else
      Multprec_Integer_Matrices.Copy(t,res);
      Multprec_Integer_Vectors.Copy(v,w);
      loop
        j := 0;
        for k in w'range loop
          if (not Equal(w(k),0)) and (k /= i)
           then j := k;
          end if;
	  exit when j /= 0;
        end loop;
        if j /= 0 then
          Copy(w(i),a); Copy(w(j),b);
          Multprec_Common_Divisors.gcd(a,b,ka,lb,d);
          Div(a,d);  Div(b,d);
          Copy(ka,t(i,i)); Copy(lb,t(i,j));
          Min(b); Copy(b,t(j,i)); Copy(a,t(j,j));
          Clear(a); Clear(b); Clear(ka); Clear(lb); Clear(d);
          Multprec_Integer_Matrices.Mul(t,w);
          Multprec_Integer_Matrices.Mul2(t,res);
          Clear(t(i,i)); Clear(t(j,j));
          Clear(t(i,j)); Clear(t(j,i));
          t(i,i) := Create(integer32(1));
          t(j,j) := Create(integer32(1));
          t(i,j) := Create(integer32(0));
          t(j,i) := Create(integer32(0));
        end if;
        exit when j = 0;
      end loop;
      Multprec_Integer_Vectors.Clear(w);
      Multprec_Integer_Matrices.Clear(t);
      return res;
    end if;
  end Rotate;

  function Eliminate ( v : Multprec_Integer_Vectors.Vector; i : integer32 )
                     return Multprec_Integer_Matrices.Matrix is

    n : constant integer32 := v'last;
    t : Multprec_Integer_Matrices.Matrix(1..n,1..n)
      := Identity_Matrix(natural32(n));
    res : Multprec_Integer_Matrices.Matrix(1..n,1..n);
    w : Multprec_Integer_Vectors.Vector(v'range);
    j1,j2 : integer32;
    a,b,ka,lb,d : Integer_Number;
    ind : constant integer32 := Pivot(v,i);
    zeroes : boolean;

  begin
    if ind < v'first or ind > v'last then
      return t;
    else
      Multprec_Integer_Matrices.Copy(t,res);
      Multprec_Integer_Vectors.Copy(v,w);
      if Equal(v(i),0) then
        Clear(res(i,i));   Clear(res(i,ind));
        Clear(res(ind,i)); Clear(res(ind,ind));
        res(i,i) := Create(integer32(0));
        res(i,ind) := Create(integer32(1));
        res(ind,i) := Create(integer32(1));
        res(ind,ind) := Create(integer32(0));
        Copy(w(ind),w(i)); Clear(w(ind));
        w(ind) := Create(integer32(0));
      end if;
      for k in v'first..i loop
        if not Equal(w(k),0)
         then j1 := k; exit;
        end if;
      end loop;
      if j1 /= i then
        zeroes := false;
        loop
          for k in (j1+1)..i loop
            if not Equal(w(k),0)
             then j2 := k; exit;
            end if;
          end loop;
          Copy(w(j1),a); Copy(w(j2),b);
          Multprec_Common_Divisors.gcd(a,b,ka,lb,d);
         -- put("gcd("); put(a,1); put(","); put(b,1); put(") = "); put(d,1); 
         -- put("  k = "); put(ka,1); put("  l = "); put(lb,1); new_line;
          Div(a,d); Div(b,d);
          Copy(lb,t(j1,j1)); Min(ka); Copy(ka,t(j1,j2));
          Copy(a,t(j2,j1)); Copy(b,t(j2,j2));
          Copy(d,w(j2));
          Multprec_Integer_Matrices.Mul2(t,res);
          Clear(a); Clear(b); Clear(ka); Clear(lb); Clear(d);
          Clear(t(j1,j1)); Clear(t(j2,j2));
          Clear(t(j1,j2)); Clear(t(j2,j1));
          t(j1,j1) := Create(integer32(1)); t(j2,j2) := Create(integer32(1));
          t(j1,j2) := Create(integer32(0)); t(j2,j1) := Create(integer32(0));
          if j2 < i
           then j1 := j2;
           else exit;
          end if;
        end loop;
      else
        zeroes := true;
      end if;
      for k in reverse i..v'last loop
        if not Equal(w(k),0)
         then j2 := k; exit;
        end if;
      end loop;
      if j2 /= i then
        loop
          for k in reverse i..(j2-1) loop
            if not Equal(w(k),0)
             then j1 := k; exit;
            end if;
          end loop;
          Copy(w(j1),a); Copy(w(j2),b);
          Multprec_Common_Divisors.gcd(a,b,ka,lb,d);
         -- put("gcd("); put(a,1); put(","); put(b,1); put(") = "); put(d,1); 
         -- put("  k = "); put(ka,1); put("  l = "); put(lb,1); new_line;
          Div(a,d); Div(b,d);
          Copy(a,t(j1,j1)); Copy(b,t(j1,j2));
          Min(lb); Copy(lb,t(j2,j1)); Copy(ka,t(j2,j2));
          Copy(d,w(j1));
          Multprec_Integer_Matrices.Mul2(t,res);
          Clear(a); Clear(b); Clear(ka); Clear(lb); Clear(d);
          Clear(t(j1,j1)); Clear(t(j2,j2));
          Clear(t(j1,j2)); Clear(t(j2,j1));
          t(j1,j1) := Create(integer32(1)); t(j2,j2) := Create(integer32(1));
          t(j1,j2) := Create(integer32(0)); t(j2,j1) := Create(integer32(0));
          if j1 > i
	   then j2 := j1;
	   else exit;
          end if;
        end loop;
      elsif zeroes then
        if not Equal(w(i),0) then
          if w(i) < 0 then
            Clear(t(i,i)); t(i,i) := Create(integer32(1)); Min(t(i,i));
            Multprec_Integer_Matrices.Mul2(t,res);
          end if;
        end if;
      end if;
      Multprec_Integer_Vectors.Clear(w);
      Multprec_Integer_Matrices.Clear(t);
      return res;
    end if;
  end Eliminate;

end Multprec_Power_Transformations;
