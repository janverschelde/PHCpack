with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package body Cyclic_Laurent_System is

  function Monomial_Support
             ( i,j,n : natural32 ) return Standard_Integer_Vectors.Vector is

    dim : constant integer32 := integer32(n)-1;
    res : Standard_Integer_Vectors.Vector(1..dim) := (1..dim => 0);
    idx : integer32;

  begin
    if j = 1 then        -- first monomial is x(i)/x(0) = x(i)
      res(integer32(i)) := 1;
    else
      res(integer32(j)-1) := -1;    -- divide by x(j-1)
      idx := integer32(i + j) - 1;
      if idx >= integer32(n)
       then idx := idx - integer32(n);
      end if;
      if idx > 0
       then res(idx) := 1;
      end if;
    end if;
    return res;
  end Monomial_Support;

  function Cyclic_Monomial ( i,j,n : natural32 ) return Term is

    res : Term;

  begin
    res.cf := Create(1.0);
    res.dg := new Standard_Integer_Vectors.Vector'(Monomial_Support(i,j,n));
    return res;
  end Cyclic_Monomial;

  function Polynomial_Support ( i,n : natural32 ) return List is

    res,last : List;

  begin
    for j in 1..n loop
      Append(res,last,Monomial_Support(i,j,n));
    end loop;
    return res;
  end Polynomial_Support;

  function Cyclic_Polynomial ( i,n : natural32 ) return Poly is

    res : Poly := Null_Poly;
    trm : Term;

  begin
    for j in 1..n loop
      trm := Cyclic_Monomial(i,j,n);
      Add(res,trm);
      Clear(trm);
    end loop;
    return res;
  end Cyclic_Polynomial;

  function System_Support ( n : natural32 ) return Array_of_Lists is

    dim : constant integer32 := integer32(n)-1;
    res : Array_of_Lists(1..dim);

  begin
    for i in 1..n-1 loop
      res(integer32(i)) := Polynomial_Support(i,n);
    end loop;
    return res;
  end System_Support;

  function Cyclic_System ( n : natural32 ) return Laur_Sys is

    dim : constant integer32 := integer32(n)-1;
    res : Laur_Sys(1..dim);

  begin
    for i in 1..n-1 loop
      res(integer32(i)) := Cyclic_Polynomial(i,n);
    end loop;
    return res;
  end Cyclic_System;

end Cyclic_Laurent_System;
