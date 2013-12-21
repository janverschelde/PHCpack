with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural64_Vectors;

-- with text_io,integer_io;  use text_io,integer_io;

package body Monomial_Hashing is

  procedure Enumerate_Monomials ( d,n : in natural32 ) is

    accu : Standard_Natural_Vectors.Vector(1..integer32(n))
         := (1..integer32(n) => 0);
    sum : natural32 := 0;

    procedure Enumerate ( i,dg : natural32 ) is

    -- DESCRIPTION :
    --   Enumerates the distribution of degree dg over entries i to n.

    begin
      if sum = d then
        Monomial(accu);
      elsif i = n then
        accu(integer32(n)) := dg;
        Monomial(accu);
        accu(integer32(n)) := 0;
      else
        for j in reverse 0..dg loop
          accu(integer32(i)) := j;
          sum := sum + j;
          Enumerate(i+1,d-sum);
          sum := sum - j;
          accu(integer32(i)) := 0;
        end loop;
      end if;
    end Enumerate;

  begin
    Enumerate(1,d);
  end Enumerate_Monomials;

  procedure Enumerate_Leaves_of_Monomial_Tree ( d,n : in natural32 ) is

    accu : Standard_Natural_Vectors.Vector(1..integer32(n))
         := (1..integer32(n) => 0);

    procedure Enumerate ( level : in natural32 ) is

    -- DESCRIPTION :
    --   Enumerates all monomials at the current level <= d.

    begin
      for i in 1..integer32(n) loop
        accu(i) := accu(i) + 1;
        if level = d
         then Monomial(accu);
         else Enumerate(level+1);
        end if;
        accu(i) := accu(i) - 1;
      end loop;
    end Enumerate;

  begin
    if d = 0
     then Monomial(accu);
     else Enumerate(1);
    end if;
  end Enumerate_Leaves_of_Monomial_Tree;

  function Monomial_Count ( d,n : natural32 ) return natural32 is

    res : natural32 := 0;

    procedure Count ( m : Standard_Natural_Vectors.Vector ) is
    begin
      res := res + 1;
    end Count;
    procedure Count_Monomials is new Enumerate_Monomials(Count);

  begin
    Count_Monomials(d,n);
    return res;
  end Monomial_Count;

  function Monomial_Code
             ( d : natural32; m : Standard_Natural_Vectors.Vector )
             return natural64 is

    res : natural64 := 0;

  begin
    for i in m'range loop
      res := res*natural64(d) + natural64(m(i));
    end loop;
    return res;
--  exception
--    when others => put_line("exception raised in Monomial code");
--      put("  m = "); 
--      for i in m'range loop
--        put(" "); put(m(i),1);
--      end loop;
--      new_line;
--      put("  d = "); put(d,1); new_line;
--     -- put("  res = "); put(res,1); new_line;
--      raise;
  end Monomial_Code;

  function Monomial_Keys
             ( k,n : natural32 ) return Standard_Natural64_VecVecs.VecVec is

    res : Standard_Natural64_VecVecs.VecVec(1..integer32(k));
    ind,code_ind : integer32 := 0;

    procedure Label ( m : in Standard_Natural_Vectors.Vector ) is
    begin
      code_ind := code_ind + 1;
      res(ind)(code_ind) := Monomial_Code(k+1,m);
    end Label;
    procedure Enumerate is new Enumerate_Monomials(Label);

  begin
    for i in 1..integer32(k) loop
      res(i) := new Standard_Natural64_Vectors.Vector
                      (1..integer32(Monomial_Count(natural32(i),n)));
      code_ind := 0; ind := i;
      Enumerate(natural32(i),n);
    end loop;
    return res;
  end Monomial_Keys;

  function Search ( monkeys : Standard_Natural64_VecVecs.VecVec;
                    c : natural64; s : natural32 ) return natural32 is

    keys : Standard_Natural64_Vectors.Link_to_Vector;

  begin
    if s <= natural32(monkeys'last) then
      keys := monkeys(integer32(s));
      for i in keys'range loop   -- linear search is not most efficient!
        if keys(i) = c
         then return natural32(i);
        end if;
      end loop;
    end if;
    return 0;
  end Search;

  function Search ( monkeys : Standard_Natural64_VecVecs.VecVec;
                    m : Standard_Natural_Vectors.Vector; s : natural32 )
                  return natural32 is

    k : constant natural32 := natural32(monkeys'last);
    c : constant natural64 := Monomial_Code(k+1,m);

  begin
    return Search(monkeys,c,s);
  end Search;

  function Search ( monkeys : Standard_Natural64_VecVecs.VecVec;
                    m : Standard_Natural_Vectors.Vector ) return natural32 is

    s : constant natural32 := Standard_Natural_Vectors.Sum(m);

  begin
    return Search(monkeys,m,s);
  end Search;

end Monomial_Hashing;
