package body Lexicographical_Supports is

  function LexLess ( a,b : Standard_Integer_Vectors.Vector ) return boolean is
  begin
    for i in a'range loop
      if a(i) < b(i) then
        return true;
      elsif a(i) > b(i) then
        return false;
      end if;
    end loop;
    return false;
  end LexLess;

  function DegLess ( a,b : Standard_Integer_Vectors.Vector ) return boolean is

    deg_a : constant integer32 := Standard_Integer_Vectors.Sum(a);
    deg_b : constant integer32 := Standard_Integer_Vectors.Sum(b);

  begin
    if deg_a < deg_b then
      return true;
    elsif deg_a > deg_b then
      return false;
    else
      return LexLess(a,b);
    end if;
  end DegLess;

  procedure Swap ( a,b : in out Standard_Integer_Vectors.Link_to_Vector ) is

    tmp : integer32;

  begin
    for i in a'range loop
      tmp := a(i);
      a(i) := b(i);
      b(i) := tmp;
    end loop;
  end Swap;

  function Sort ( s : List ) return List is

    res,first,second : List;
    A,B : Standard_Integer_Vectors.Link_to_Vector;

  begin
    Copy(s,res);
    if Length_Of(res) >= 2 then
      first := res;
      while not Is_Null(first) loop
        A := Head_Of(first);
        second := Tail_Of(first);
        while not Is_Null(second) loop
          B := Head_Of(second);
         -- if LexLess(B.all,A.all)
          if DegLess(B.all,A.all)
           then Swap(A,B);
          end if;
          second := Tail_Of(second);
        end loop;
        first := Tail_Of(first);
      end loop;
    end if;
    return res;
  end Sort;

  function Index ( s : Standard_Integer_VecVecs.VecVec;
                   v : Standard_Integer_Vectors.Vector ) return integer32 is

    sfirst : integer32 := s'first;
    slast : integer32 := s'last;
    middle : integer32;

  begin
    while sfirst <= slast loop
      if sfirst = slast then
        if Standard_Integer_Vectors.Equal(s(sfirst).all,v)
         then return sfirst;
         else return 0;
        end if;
      end if;
      middle := (sfirst+slast)/2;
      if Standard_Integer_Vectors.Equal(s(middle).all,v) then
        return middle;
     -- elsif LexLess(s(middle).all,v) then
      elsif DegLess(s(middle).all,v) then
        sfirst := middle + 1;
      else
        slast := middle - 1;
      end if;
    end loop;
    return 0;
  end Index;

-- FACTORED REPRESENTATION :

  function First_Positive ( s : List ) return natural32 is

    res : natural32 := 0;
    tmp : List := s;
    ptr : Standard_Integer_Vectors.Link_to_Vector;
    neg,pos : boolean;

  begin
    while not Is_Null(tmp) loop
      ptr := Head_Of(tmp);
      res := res + 1;
      neg := false; pos := false;
      for i in ptr'range loop
        if ptr(i) < 0 then
          neg := true;
        elsif ptr(i) > 0 then
          pos := true;
        end if;
        exit when neg;
      end loop;
      if pos and not neg
       then return res;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return 0;
  end First_Positive;

  procedure Factor ( v : in out Standard_Integer_VecVecs.VecVec;
                     org : in out Standard_Integer_VecVecs.VecVec;
                     start,k : in integer32 ) is

    w : constant Standard_Integer_Vectors.Link_to_Vector := v(k);
    p : Standard_Integer_Vectors.Link_to_Vector;
    divides : boolean;
    d : Standard_Integer_Vectors.Vector(w'range);

  begin
    for i in reverse start..(k-1) loop -- reverse: find largest factor
      p := org(i); d(0) := i;
      divides := true;
      for j in 1..p'last loop
        d(j) := w(j) - p(j);
        if d(j) < 0
         then divides := false; exit;
        end if;
      end loop;
      if divides
       then w.all := d; exit;
      end if;
    end loop;
  end Factor;

  function Factor ( s : List ) return Standard_Integer_VecVecs.VecVec is

    res : Standard_Integer_VecVecs.VecVec(1..integer32(Length_Of(s)));
    org : Standard_Integer_VecVecs.VecVec(res'range);
    tmp : List := s;
    ind : integer32 := 0;
    ptr : Standard_Integer_Vectors.Link_to_Vector;
    start : constant integer32 := integer32(First_Positive(s));

  begin
    while not Is_Null(tmp) loop
      ptr := Head_Of(tmp);
      ind := ind + 1;
      declare
        v : Standard_Integer_Vectors.Vector(0..ptr'last);
      begin
        v(0) := 0;
        v(ptr'range) := ptr.all;
        org(ind) := ptr;
        res(ind) := new Standard_Integer_Vectors.Vector'(v);
        Factor(res,org,start,ind);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Factor;

-- COMPRESSION :

  function Compress ( v : Standard_Integer_Vectors.Vector )
                    return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(v'first..v'last*2);
    ind : integer32 := 1;

  begin
    if v'first = 0 then
      res(res'first) := v(v'first);   -- copy factorization
    end if;
    for i in 1..v'last loop
      if v(i) > 0 then                -- nonzero exponent
        res(ind) := i;
        res(ind+1) := v(i);
        ind := ind + 2;
      end if;
    end loop;
    return res(res'first..ind-1);
  end Compress;

  function Compress ( v : Standard_Integer_VecVecs.VecVec )
                    return Standard_Integer_VecVecs.VecVec is

    res : Standard_Integer_VecVecs.VecVec(v'range);

  begin
    for i in v'range loop
      declare
        p : Standard_Integer_Vectors.Link_to_Vector := v(i);
        c : constant Standard_Integer_Vectors.Vector := Compress(p.all);
      begin
        res(i) := new Standard_Integer_Vectors.Vector'(c);
      end;
    end loop;
    return res;
  end Compress;

-- PRE-FACTORIZATION :

  function Decrement ( v : Standard_Integer_Vectors.Vector ) 
                     return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      if v(i) <= 1 
       then res(i) := 0;
       else res(i) := v(i) - 1;
      end if;
    end loop;
    return res;
  end Decrement;

  function Nodes ( s : Standard_Integer_VecVecs.VecVec ) return List is

    res,res_last : List;
    ptr : Standard_Integer_Vectors.Link_to_Vector;

  begin
    for i in s'range loop
      ptr := s(i);
      declare
        d : constant Standard_Integer_Vectors.Vector := Decrement(ptr.all);
      begin
        if Standard_Integer_Vectors.Sum(d) > 1 then
          if Index(s,d) = 0 then
            if not Is_In(res,d)
             then Append(res,res_last,d);
            end if;
          end if;
        end if;
      end;
    end loop;
    return res;
  end Nodes;

end Lexicographical_Supports;
