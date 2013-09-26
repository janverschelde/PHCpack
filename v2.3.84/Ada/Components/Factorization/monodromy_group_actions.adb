with unchecked_deallocation;
with Standard_Natural_VecVecs;          use Standard_Natural_VecVecs;

package body Monodromy_Group_Actions is

-- DATA STRUCTURES TO MANAGE THE DECOMPOSITION :

  type Boolean_Array is array ( integer32 range <> ) of boolean;

  type Irreducible_Components_Rep ( n : integer32 ) is record
    sets : Standard_Natural_VecVecs.VecVec(1..n);
    actv : Boolean_Array(1..n);
  end record;   

-- AUXILIARIES TO THE CONSTRUCTOR Act :

  procedure Sort ( s : in out Vector ) is

  -- DESCRIPTION :
  --   Sorts the elements in s in increasing order.

    n : constant integer32 := s'last;
    min : natural32;
    ind : integer32;

  begin
    for i in 1..n loop
      exit when (s(i) = 0);
      min := s(i); ind := i;
      for j in i+1..n loop
        exit when (s(j) = 0);
        if s(j) < min
         then ind := j; min := s(j);
        end if;
      end loop;
      if ind /= i 
       then s(ind) := s(i); s(i) := min;
      end if;
    end loop;
  end Sort;

  function Active ( ic : Irreducible_Components; k : integer32 )
                  return integer32 is

  -- DESCRIPTION :
  --   Returns the entry of the leading set where k belongs to.

  begin
    if ic.actv(k) then
      return k;
    else
      for i in 1..ic.n loop
        if ic.actv(i) then
          if Is_In(ic,i,natural32(k))
           then return i;
          end if;
        end if;
      end loop;
    end if;
    return 0;
  end Active; 

-- CONSTRUCTORS :

  function Create ( n : integer32 ) return Irreducible_Components is

    res : Irreducible_Components;
    res_rep : Irreducible_Components_Rep(n);

  begin
    for i in 1..n loop
      res_rep.sets(i) := new Standard_Natural_Vectors.Vector'(1..n => 0);
      res_rep.sets(i)(1) := natural32(i);
      res_rep.actv(i) := true;
    end loop;           
    res := new Irreducible_Components_Rep'(res_rep);
    return res;
  end Create;

  procedure Add ( ic : in out Irreducible_Components;
                  i : in integer32; j : in natural32 ) is

    found : boolean := false;

  begin
    for k in 1..ic.n loop
      if ic.sets(i)(k) = j then
        found := true;
      elsif ic.sets(i)(k) = 0 then
        ic.sets(i)(k) := j;
        Sort(ic.sets(i).all);
        found := true;
      end if;
      exit when found;
    end loop;
  end Add;

  procedure Merge ( ic : in out Irreducible_Components;
                    i,j : in integer32 ) is

    actvi,actvj : integer32;

  begin
    if ic.actv(i)
     then actvi := i;
     else actvi := Active(ic,i);
    end if;
    if ic.actv(j)
     then actvj := j;
     else actvj := Active(ic,j);
    end if;
    for k in 1..integer32(Cardinality(ic,actvj)) loop
      Add(ic,actvi,ic.sets(actvj)(k));
    end loop;
    ic.actv(actvj) := false;
  end Merge;

  procedure Act ( ic : in out Irreducible_Components; map : in Vector ) is
  begin
    for i in map'range loop
      if ic.actv(i) then
        for j in map'range loop
          if j = i then
            if map(i) /= natural32(i) then
              if not Is_In(ic,i,map(i))
               then Merge(ic,i,integer32(map(i)));
              end if;
            end if;
          elsif ic.actv(j) then
            if map(j) = natural32(i)
             then Merge(ic,i,j);
            end if;
          end if;
        end loop;
      end if;
    end loop;
  end Act;

-- SELECTORS :

  function Empty ( ic : Irreducible_Components;
                   i : integer32 ) return boolean is
  begin
    if ic = null then
      return true;
    elsif i <= ic.n then
      return not ic.actv(i);
    else
      return true;
    end if;
  end Empty;

  function Cardinality ( ic : Irreducible_Components; i : integer32 )
                       return natural32 is

    res : natural32 := 0;

  begin
    if not Empty(ic,i) then
      for j in ic.sets(i)'range loop
        exit when (ic.sets(i)(j) = 0);
        res := res + 1;
      end loop;
    end if;
    return res;
  end Cardinality;

  function Component ( ic : Irreducible_Components; i : integer32 )
                     return Vector is

    crd : constant natural32 := Cardinality(ic,i);
    res : Vector(1..integer32(crd));

  begin
    if crd > 0
     then res := ic.sets(i)(res'range);
    end if;
    return res;
  end Component;

  function Is_In ( ic : Irreducible_Components;
                   i : integer32; j : natural32 ) return boolean is
  begin
    if not Empty(ic,i) then
      for k in ic.sets(i)'range loop
        exit when (ic.sets(i)(k) = 0);
        if ic.sets(i)(k) = j
         then return true;
        end if;
      end loop;
    end if;
    return false;
  end Is_In;

  function Empty ( ic : Irreducible_Components ) return boolean is
  begin
    return (ic = null);
  end Empty;

  function Cardinality ( ic : Irreducible_Components ) return natural32 is

    res : natural32 := 0;

  begin
    if not Empty(ic) then
      for i in 1..ic.n loop
        if ic.actv(i)
         then res := res + 1;
        end if;
      end loop;
    end if;
    return res;
  end Cardinality;

  function Sum_of_Degrees ( ic : Irreducible_Components ) return integer32 is
  begin
    if Empty(ic)
     then return 0;
     else return ic.n;
    end if;
  end Sum_of_Degrees;

  function Degrees ( ic : Irreducible_Components ) return Vector is

    res : Vector(1..ic.n);
    cnt : integer32 := 0;

  begin
    if not Empty(ic) then
      for i in 1..ic.n loop
        if ic.actv(i) then
          cnt := cnt + 1;
          res(cnt) := Cardinality(ic,i);
        end if;
      end loop;
      Sort(res(1..cnt));
    end if;
    return res(1..cnt);
  end Degrees;

  function Nonempty_Sets ( ic : Irreducible_Components ) return Vector is

    res : Vector(1..ic.n);
    cnt : integer32 := 0;

  begin
    if not Empty(ic) then
      for i in 1..ic.n loop
        if ic.actv(i) then
          cnt := cnt + 1;
          res(cnt) := natural32(i);
        end if;
      end loop;
    end if;
    return res(1..cnt);
  end Nonempty_Sets;

  function Representatives ( ic : Irreducible_Components; k : natural32 )
                           return Vector is

    res : Vector(1..ic.n);
    cnt : integer32 := 0;
    ind : natural32;

  begin
    if not Empty(ic) then
      for i in 1..ic.n loop
        if ic.actv(i) then
          cnt := cnt + 1;
          ind := k mod Cardinality(ic,i);
          if ind = 0
           then res(cnt) := ic.sets(i)(1);
           else res(cnt) := ic.sets(i)(integer32(ind));
          end if;
        end if;
      end loop;
    end if;
    Sort(res(1..cnt));
    return res(1..cnt);
  end Representatives;

  function Representatives ( ic : Irreducible_Components; k,j : integer32 )
                           return Vector is

    res : Vector(1..ic.n);
    cnt : integer32 := 0;
    start,ind,crd : integer32;

  begin
    if not Empty(ic) then
      for i in 1..ic.n loop
        if ic.actv(i) then
          cnt := cnt + 1;
          crd := integer32(Cardinality(ic,i));
          start := k mod crd;
          if start = 0
           then start := 1;
          end if;
          res(cnt) := ic.sets(i)(start);
          if j > 0 then
            ind := start + j;
            while ind <= crd loop
              cnt := cnt + 1;
              res(cnt) := ic.sets(i)(ind);
              ind := ind + j;
            end loop;
            ind := ind mod crd;
            if ind = 0
             then ind := 1;
            end if;
            while ind < start loop
              cnt := cnt + 1;
              res(cnt) := ic.sets(i)(ind);
              ind := ind + j;
            end loop;
          end if;
        end if;
      end loop;
    end if;
    Sort(res(1..cnt));
    return res(1..cnt);
  end Representatives;  

-- DESTRUCTOR :

  procedure Clear ( ic : in out Irreducible_Components ) is

    procedure free is new unchecked_deallocation
      (Irreducible_Components_Rep,Irreducible_Components);

  begin
    if ic /= null then
      for i in ic.sets'range loop
        Clear(ic.sets(i));
      end loop;
      free(ic);
    end if;
  end Clear;

end Monodromy_Group_Actions;
