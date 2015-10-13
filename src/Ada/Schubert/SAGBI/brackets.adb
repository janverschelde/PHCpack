with unchecked_deallocation;

package body Brackets is

-- AUXILIARY OPERATION :

  procedure Swap ( v : in out Standard_Natural_Vectors.Vector;
                   i,j : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps the i-th and j-th entry in the vector v.

    tmp : constant natural32 := v(i);

  begin
    v(i) := v(j); v(j) := tmp;
  end Swap;

-- CONSTRUCTORS :

  procedure Create ( v : in Standard_Natural_Vectors.Vector;
                     b : out Bracket; sign : out integer32 ) is

    sig : integer32 := +1;
    min : natural32;
    ind : integer32;
    bb : Bracket(v'range) := Bracket(v);
  
  begin
    for i in bb'first..bb'last-1 loop
      min := bb(i);
      ind := i;
      for j in i+1..bb'last loop
        if bb(j) < min then
          ind := j;
          min := bb(j);
        end if;
      end loop;
      if ind /= i then
        Swap(Standard_Natural_Vectors.Vector(bb),i,ind);
        sig := -sig;
      end if;
    end loop;
    b := bb;
    sign := sig;
  end Create;

  procedure Create ( v : in Standard_Natural_Vectors.Vector;
                     perm : out Standard_Natural_Vectors.Vector;
                     b : out Bracket; sign : out integer32 ) is

    sig : integer32 := +1;
    min : natural32;
    ind : integer32;
    bb : Bracket(v'range) := Bracket(v);
    pp : Standard_Natural_Vectors.Vector(v'range);
  
  begin
    for i in pp'range loop
      pp(i) := natural32(i);
    end loop;
    for i in bb'first..bb'last-1 loop
      min := bb(i);
      ind := i;
      for j in i+1..bb'last loop
        if bb(j) < min then
          ind := j;
          min := bb(j);
        end if;
      end loop;
      if ind /= i then
        Swap(Standard_Natural_Vectors.Vector(bb),i,ind);
        Swap(pp,i,ind);
        sig := -sig;
      end if;
    end loop;
    perm := pp;
    b := bb;
    sign := sig;
  end Create;

  function Modulo ( b : Bracket; n : natural32 ) return Bracket is

    res : Bracket(b'range);
    modvec : Standard_Natural_Vectors.Vector(b'range);
    sig : integer32;

  begin
    for i in b'range loop
      modvec(i) := b(i) mod n;
      if modvec(i) = 0
       then modvec(i) := n;
      end if;
    end loop;
    Create(modvec,res,sig);
    return res;
  end Modulo;

  procedure Modulo ( b : in Bracket; n : in natural32;
                     perm : out Standard_Natural_Vectors.Vector;
                     mb : out Bracket ) is

    res : Bracket(b'range);
    modvec : Standard_Natural_Vectors.Vector(b'range);
    sig : integer32;

  begin
    for i in b'range loop
      modvec(i) := b(i) mod n;
      if modvec(i) = 0
       then modvec(i) := n;
      end if;
    end loop;
    Create(modvec,perm,res,sig);
    mb := res;
  end Modulo;

-- SELECTORS :

  function Is_Zero ( b : Bracket ) return boolean is

  begin
    for i in b'first..b'last-1 loop
      if b(i) = b(i+1)
       then return true;
      end if;
    end loop;
    return false;
  end Is_Zero;

  function Is_Equal ( b1,b2 : Bracket ) return boolean is

    use Standard_Natural_Vectors;

  begin
    if b1'length /= b2'length
     then return false;
     else return Equal(Vector(b1),Vector(b2));
    end if;
  end Is_Equal;

  function "<" ( b1,b2 : Bracket ) return boolean is
  begin
    for i in b1'range loop
      if b1(i) < b2(i) then
        return true;
      elsif b1(i) > b2(i) then
        return false;
      end if;
    end loop;
    return false;
  end "<";

  function ">" ( b1,b2 : Bracket ) return boolean is
  begin
    for i in b1'range loop
      if b1(i) > b2(i) then
        return true;
      elsif b1(i) < b2(i) then
        return false;
      end if;
    end loop;
    return false;
  end ">";

  function Is_Standard ( b1,b2 : Bracket ) return natural32 is
  begin
    for i in b1'range loop
      if b1(i) > b2(i)
       then return natural32(i);
      end if;
    end loop;
    return 0;
  end Is_Standard;

  function "<=" ( alpha,beta : Bracket ) return boolean is
  begin
    for i in alpha'range loop
      if alpha(i) > beta(i)
       then return false;
      end if;
    end loop;
    return true;
  end "<=";

-- DESTRUCTORS :

  procedure Clear ( b : in out Link_to_Bracket ) is

    procedure free is new unchecked_deallocation(Bracket,Link_to_Bracket);
 
  begin
    free(b);
  end Clear;

  procedure Clear ( b : in out Array_of_Brackets ) is
  begin
    for i in b'range loop
      Clear(b(i));
    end loop;
  end Clear;

end Brackets;
