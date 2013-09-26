with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package body Complex_Osculating_Planes is

  function Standard_Basis
              ( n,d : natural32; s : Complex_Number ) return Matrix is

    res : Matrix(1..integer32(n),1..integer32(d));
    acc : Complex_Number;
    fac,lim : integer32;

  begin
    for i in 1..integer32(d) loop                -- set the zeros and ones
      res(i,i) := Create(1.0);
      for j in (i+1)..integer32(d) loop
        res(i,j) := Create(0.0);
      end loop;
    end loop;
    for j in 1..integer32(d) loop        -- set the powers of the s-values
      acc := s;
      for i in (j+1)..integer32(n) loop
        res(i,j) := acc;
        acc := acc*s;
      end loop;
    end loop;
    for i in 3..integer32(n) loop    -- compute the factors from derivation
      fac := 1;
      if i-1 > integer32(d)
       then lim := integer32(d);
       else lim := i-1;
      end if;
      for j in 2..lim loop
        fac := fac*(i+1-j);
        res(i,j) := Create(fac)*res(i,j);
      end loop;
      if i <= integer32(d)
       then res(i,i) := Create(fac);
      end if;
    end loop;
    for j in 3..integer32(d) loop                     -- scale the columns
      for i in (j+1)..integer32(n) loop
        res(i,j) := res(i,j)/res(j,j);
      end loop;
      res(j,j) := Create(1.0);
    end loop;
    return res;
  end Standard_Basis;

 -- function Chebychev_Basis
 --               ( n,d : natural; s : Complex_Number ) return Matrix is
--
--    res : Matrix(1..n,1..d);
--
 -- begin
  --  return res;
  --end Chebychev_Basis;

  --function Orthogonal_Basis
  --            ( n,d : natural; s : Complex_Number ) return Matrix is
--
 --   res : Matrix(1..n,1..d);
--
 -- begin
  --  return res;
  --end Orthogonal_Basis;

end Complex_Osculating_Planes;
