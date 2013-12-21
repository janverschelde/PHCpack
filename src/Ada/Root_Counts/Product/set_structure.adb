with unchecked_deallocation;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Generate_Unions;

package body Set_Structure is

-- DATASTRUCTURES :

  type set is array ( natural32 range <> ) of boolean;
  type link_to_set is access set;
  procedure free is new unchecked_deallocation(set,link_to_set);

  type set_equations is array ( natural32 range <> ) of link_to_set;
  type link_to_set_equations is access set_equations;
  procedure free is 
    new unchecked_deallocation(set_equations,link_to_set_equations);

  type set_system is array ( natural32 range <> ) of link_to_set_equations;
  type link_to_set_system is access set_system;
  procedure free is
    new unchecked_deallocation(set_system,link_to_set_system);

-- INTERNAL DATA :

  n : natural32 := 0;  -- the number of unknowns and equations

  ls : link_to_set_system := null;

-- CONSTRUCTORS :

  procedure Init ( ns : in Standard_Natural_Vectors.Vector ) is
  begin
    n := natural32(ns'length);
    ls := new set_system(1..n);
    for i in ls'range loop
      ls(i) := new set_equations(1..ns(integer32(i)));
      for j in ls(i)'range loop
	ls(i).all(j) := new set'(1..n => false);
      end loop;
    end loop;
  end Init;

  procedure Add ( i,j,k : in natural32 ) is

    s : set renames ls(i).all(j).all;

  begin
    s(k) := true;
  end Add;

  procedure Remove ( i,j,k : in natural32 ) is

    s : set renames ls(i).all(j).all;

  begin
    s(k) := false;
  end Remove;

-- SELECTORS :

  function Empty return boolean is
  begin
    return (ls = null);
  end Empty;

  function Dimension return natural32 is
  begin
    return n;
  end Dimension;

  function Number_of_Sets ( i : natural32 ) return natural32 is
  begin
    return ls(i)'last;
  end Number_of_Sets;

  function Is_In ( i,j,k : natural32 ) return boolean is 

    s : set renames ls(i).all(j).all;

  begin
    return s(k);
  end Is_In;

-- COMPUTING THE UPPER BOUND :

  function Extent_Of ( s : in set ) return natural32 is

   -- DESCRIPTION : return the number of elements in s

    sum : natural32 := 0;

  begin
    for i in s'range loop
      if s(i)
       then sum := sum + 1;
      end if;
    end loop;
    return sum;
  end Extent_Of;

  procedure Union ( s : in set; u : in out set ) is

   -- DESCRIPTION : u = u U s

  begin
    for i in s'range loop
      if s(i)
       then u(i) := true;
      end if;
    end loop;
  end Union;

  function acceptable ( lset_eq : link_to_set_equations;
                        k,n : natural32; lset : link_to_set )
                      return boolean is

    type arr is array ( natural32 range <> ) of boolean;
    accep : boolean := true;

    procedure check ( a : in arr; continue : out boolean ) is

      u : set(lset'range);

    begin
      u := lset.all;
      for i in a'range loop
        if a(i)
         then Union(lset_eq(i).all,u);
        end if;
      end loop;
      accep := ( Extent_Of(u) >= k+1 );
      continue := accep;
    end check;

    procedure gen is new Generate_Unions(arr,check);

  begin
    gen(k,1,n);  -- generates all possible unions of k sets
		 -- out of the range 1..n
    return accep;
  end acceptable;

  function acceptable ( lset_eq : link_to_set_equations; 
		        n : natural32; lset : link_to_set ) return boolean is

   -- DESCRIPTION :
   --   if acceptable(lset_eq,n)
   --    then verify if acceptable(lset_eq + lset,n+1)

  begin
    for k in 1..n loop
      if not acceptable(lset_eq,k,n,lset)
       then return false;
      end if;
    end loop;
    return true;
  end acceptable;

  procedure Compute ( i,n,sum : in natural32; res : in out natural32;
                      lset_eq : in out link_to_set_equations ) is
  begin
    if i > n then
      res := res + sum;
    else -- Pick out a set and check if it is allowed :
      for j in ls(i)'range loop
        if acceptable(lset_eq,i-1,ls(i).all(j)) then
          lset_eq(i) := ls(i).all(j);
          Compute(i+1,n,sum,res,lset_eq);
        end if;
      end loop;
    end if;
  end Compute;

  function B return natural32 is

    res : natural32 := 0;
    lset_eq : link_to_set_equations := new set_equations(1..n);

  begin
    for i in lset_eq'range loop
      lset_eq(i) := new set'(1..n => false);
    end loop;
    Compute(1,n,1,res,lset_eq);
    return res;
  end B;

  procedure Compute ( i,n,sum : in natural32; res : in out natural32;
                      lset_eq : in out link_to_set_equations;
		      pos : in out Standard_Integer_Vectors.Vector;
		      first,last : in out List ) is
  begin
    if i > n then
      res := res + sum;
      Append(first,last,pos);
    else -- Pick out a set and check if it is allowed :
      for j in ls(i)'range loop
        pos(integer32(i)) := integer32(j);
        if acceptable(lset_eq,i-1,ls(i).all(j)) then
          lset_eq(i) := ls(i).all(j);
          Compute(i+1,n,sum,res,lset_eq,pos,first,last);
        end if;
      end loop;
    end if;
  end Compute;

  procedure B ( bn : out natural32; L : in out List ) is

    res : natural32 := 0;
    lset_eq : link_to_set_equations := new set_equations(1..n);
    pos : Standard_Integer_Vectors.Vector(1..integer32(n))
        := (1..integer32(n) => 1);
    last : List;

  begin
    for i in lset_eq'range loop
      lset_eq(i) := new set'(1..n => false);
    end loop;
    Compute(1,n,1,res,lset_eq,pos,L,last);
    bn := res;
  end B;

-- DESTRUCTOR :

  procedure Clear is
  begin
    if ls /= null then
      for i in ls'range loop
        for j in ls(i)'range loop
          free(ls(i).all(j));
        end loop;
        free(ls(i));
      end loop;
      free(ls);
    end if;
    n := 0; ls := null;
  end Clear;

end Set_Structure;
