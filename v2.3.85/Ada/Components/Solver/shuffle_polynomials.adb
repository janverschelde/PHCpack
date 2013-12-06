with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Permute_Operations;                 use Permute_Operations;

package body Shuffle_Polynomials is

  function Degrees ( f : Poly_Sys ) return Permutation is

  -- DESCRIPTION :
  --   Returns the degree vector of the equations in f.

    res : Permutation(f'range);

  begin
    for i in res'range loop
      res(i) := Degree(f(i));
    end loop;
    return res;
  end Degrees;

  function Identity ( n : integer32 ) return Permutation is

  -- DESCRIPTION :
  --   Returns the identity permutation for n elements.

    res : Permutation(1..n);

  begin
    for i in 1..n loop
      res(i) := i;
    end loop;
    return res;
  end Identity;

  function Increasing_Sort ( v : Permutation ) return Permutation is

  -- DESCRIPTION :
  --   Returns the permutation to sort the elements of v
  --   in increasing order.

    n : constant integer32 := v'last;
    res : Permutation(v'range) := Identity(n);
    wrk : Permutation(v'range) := v;
    tmp : integer32;

  begin
    for i in 1..n-1 loop
      for j in i+1..n loop
        if wrk(j) < wrk(i) then
          tmp := wrk(i); wrk(i) := wrk(j); wrk(j) := tmp;
          tmp := res(i); res(i) := res(j); res(j) := tmp;
        end if;
      end loop;
    end loop;
    return res;
  end Increasing_Sort;

  function Decreasing_Sort ( v : Permutation ) return Permutation is

  -- DESCRIPTION :
  --   Returns the permutation to sort the elements of v
  --   in decreasing order.

    n : constant integer32 := v'last;
    res : Permutation(v'range) := Identity(n);
    wrk : Permutation(v'range) := v;
    tmp : integer32;

  begin
    for i in 1..n-1 loop
      for j in i+1..n loop
        if wrk(j) > wrk(i) then
          tmp := wrk(i); wrk(i) := wrk(j); wrk(j) := tmp;
          tmp := res(i); res(i) := res(j); res(j) := tmp;
        end if;
      end loop;
    end loop;
    return res;
  end Decreasing_Sort;

-- TARGET FUNCTIONS :

  function Increasing_Degree_Sort ( f : Poly_Sys ) return Poly_Sys is

    d : Permutation(f'range) := Degrees(f);
    p : Permutation(f'range) := Increasing_Sort(d);
  
  begin
    return p*f;
  end Increasing_Degree_Sort;

  function Decreasing_Degree_Sort ( f : Poly_Sys ) return Poly_Sys is

    d : Permutation(f'range) := Degrees(f);
    p : Permutation(f'range) := Decreasing_Sort(d);

  begin
    return p*f;
  end Decreasing_Degree_Sort;

  function Permute ( f : Poly_Sys; p : Permutation ) return Poly_Sys is
  begin
    return p*f;
  end Permute;

  function Random_Permutation ( n : integer32 ) return Permutation is

    res : Permutation(1..n);
    pool : Permutation(1..n) := Identity(n);
    r : integer32;

  begin
    for k in reverse 2..n loop
      r := Random(1,k);
      res(k) := pool(r);
      pool(r) := pool(k);
    end loop;
    res(1) := pool(1);
    return res;
  end Random_Permutation;

  function Random_Permute ( f : Poly_Sys ) return Poly_Sys is

    p : Permutation(f'range) := Random_Permutation(f'last);

  begin
    return p*f;
  end Random_Permute;

end Shuffle_Polynomials;
