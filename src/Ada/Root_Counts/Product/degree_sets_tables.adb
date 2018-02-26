with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Boolean_Matrices;
with Set_Structure;
with Pivot_Selection;

package body Degree_Sets_Tables is

-- AUXILIAIRIES :

  function Number_of_Sets return natural32 is

    res : natural32 := 0;

  begin
    for i in 1..Set_Structure.Dimension loop
      res := res + Set_Structure.Number_of_Sets(i);
    end loop;
    return res;
  end Number_of_Sets;

  function Create_Set ( n,i,j : natural32 ) return Set is

  -- DESCRIPTION :
  --   Returns the jth set for the ith equation in the set structure.

    res : Set := Create(n);

  begin
    for k in 1..n loop
      if Set_Structure.Is_In(i,j,k)
       then Add(res,k);
      end if;
    end loop;
    return res;
  end Create_Set;

  function Is_In ( ase : Array_of_Sets; s : Set ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the given set s occurs in the array of sets.

  begin
    for i in ase'range loop
      if Is_Equal(ase(i),s)
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Different_Sets return Array_of_Sets is

  -- DESCRIPTION :
  --   Returns the array of different sets of the set structure.

    n : constant natural32 := Set_Structure.Dimension;
    nbs : constant natural32 := Number_Of_Sets;
    res : Array_of_Sets(1..integer32(nbs));
    cnt : integer32 := 0;

  begin
    for i in 1..n loop
      for j in 1..Set_Structure.Number_of_Sets(i) loop
        declare
          s : Set := Create_Set(n,i,j);
        begin
          if not Is_In(res(1..cnt),s)
           then cnt := cnt + 1;
                res(cnt) := s;
           else Clear(s);
          end if;
        end;
      end loop;
    end loop;
    return res(1..cnt);
  end Different_Sets;

  function Index ( ase : Array_of_Sets; s : Set ) return integer32 is

  -- DESCRIPTION :
  --   Returns the index of the given set in the array of sets.
  --   If the set does not occur in ase, then ase'last+1 will be returned.

  begin
    for i in ase'range loop
      if Is_Equal(ase(i),s)
       then return i;
      end if;
    end loop;
    return ase'last+1;
  end Index;

-- CONSTRUCTOR :

  function Create return Degree_Sets_Table is

    n : constant natural32 := Set_Structure.Dimension;
    ase : constant Array_of_Sets := Different_Sets;
    res : Degree_Sets_Table(integer32(n),ase'length);

  begin
    res.s := ase;
    for i in res.a'range(1) loop
      for j in res.a'range(2) loop
        res.a(i,j) := 0;
      end loop;
    end loop;
    for i in 1..n loop
      for j in 1..Set_Structure.Number_of_Sets(i) loop
        declare
          s : Set := Create_Set(n,i,j);
          k : constant integer32 := Index(res.s,s);
        begin
          res.a(integer32(i),k) := res.a(integer32(i),k) + 1;
          Clear(s);
        end;
      end loop;
    end loop;
    return res;
  end Create;

-- PERMANENT COMPUTATIONS :

  function Union_Acceptable ( s : Array_of_Sets ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the union of all sets in s contains at least
  --   as many elements as the length of s, returns false otherwise.

    res : boolean;
    uni : Set := Create(s(s'first));

  begin
    for i in s'first+1..s'last loop
      Union(uni,s(i));
    end loop;
    res := (Extent_Of(uni) >= s'length);
    Clear(uni);
    return res;
  end Union_Acceptable;

  function Partial_Acceptable ( s : Array_of_Sets; k : natural32 )
                              return boolean is

  -- DESCRIPTION :
  --   Checks whether any union of k sets out of s(s'first)..s(s'last-1),
  --   together with s(s'last) forms an acceptable tuple.

    res : boolean := true;
    accu : Set := Create(s(s'last));

    function Partial_Acceptable ( s : Array_of_Sets; k,l,start : natural32;
                                  uni : Set ) return boolean is

    -- DESCRIPTION : recursive enumeration of all candidates.

    -- ON ENTRY :
    --   l         the number of sets still to choose;
    --   start     choose out of s(start..s'last-1);
    --   uni       partial union.

      res : boolean := true;

    begin
      if l = 0 then
        res := (Extent_Of(uni) >= k+1);
      else
        for ll in integer32(start)..(s'last-integer32(l)) loop
          declare
            newuni : Set := Create(uni);
          begin
            Union(newuni,s(ll));
            res := Partial_Acceptable(s,k,l-1,natural32(ll)+1,newuni);
            exit when not res;
            Clear(newuni);
          end;
        end loop;
      end if;
    --  if not res
    --   then put("Not acceptable with "); put(uni); put(" for k = ");
    --        put(k,1); new_line;
    --  end if;
      return res;
    end Partial_Acceptable;

  begin
    res := Partial_Acceptable(s,k,k,natural32(s'first),accu);
    Clear(accu);
    return res;
  end Partial_Acceptable;

  function Acceptable ( s : Array_of_Sets ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the array of sets is an acceptable tuple.
  --   The first s'last-1 sets form already an acceptable tuple and
  --   are ordered according to the cardinality of their union with
  --   the last set, from low to high.

    extlast : constant natural32 := Extent_Of(s(s'last));
    res : boolean;

  begin
   -- put_line("The array of sets "); put(s); new_line;
    if not Union_Acceptable(s) then
      res := false;
    else
      res := true;
      for k in extlast..natural32(s'last)-2 loop
        res := Partial_Acceptable(s,k);
        exit when not res;
      end loop;
    end if;
   -- if res
   --  then put_line("is an acceptable tuple.");
   --  else put_line("is not an acceptable tuple.");
   -- end if;
    return res;
  end Acceptable;

  function Acceptable ( s : Array_of_Sets; v : Vector; i : integer32 )
                      return boolean is
  -- DESCRIPTION :
  --   Returns true if the choice of sets { s(v(j)) }, j=1,2,..,i, is 
  --   acceptable.  The first i-1 sets form already an acceptable tuple.

  begin 
    if (i = v'first) or (integer32(Extent_Of(s(integer32(v(i))))) = i) then
      return true;
    else
      declare
        sv,osv : Array_of_Sets(1..i);
        min,extset : natural32;
        minind : integer32;
        u : Set;
      begin
        for j in 1..i loop                   -- create tuple of sets
          sv(j) := s(integer32(v(j)));
        end loop;
        for j in 1..(i-1) loop               -- order tuple of sets
          u := Union(sv(j),sv(i));
          min := Extent_Of(u); Clear(u);
          minind := j;
          for k in j+1..(i-1) loop
            u := Union(sv(k),sv(i));
            extset := Extent_Of(u); Clear(u);
            if extset < min
             then min := extset; minind := k;
            end if;
          end loop;
          osv(j) := sv(minind);
          if j /= minind
           then sv(minind) := sv(j);
          end if;
        end loop;
        osv(i) := sv(i);
        return Acceptable(osv);
      end;
    end if;
  end Acceptable;

  function Permanent ( a : matrix; s : Array_of_Sets; v : Vector;
                       i,n : integer32 ) return integer32 is

  -- ALGORITHM : Row expansion without memory.

  begin
    if i = n+1 then
      return 1;
    else
      declare
        res : integer32 := 0;
        vv : Vector(v'range) := v;
      begin
        for j in a'range(2) loop
          if a(i,j) /= 0 then
            vv(i) := natural32(j);
            if Acceptable(s,vv,i)
             then res := res + a(i,j)*Permanent(a,s,vv,i+1,n);
            end if;
          end if;
        end loop;
        return res;
      end;
    end if;
  end Permanent;

  function Permanent ( dst : Degree_Sets_Table ) return integer32 is

    v : constant Vector(1..dst.n) := (1..dst.n => 0);

  begin
    return Permanent(dst.a,dst.s,v,1,dst.n);
  end Permanent;

  function Adjacency_Matrix
             ( s : Array_of_Sets; v : Vector; i,n : integer32 )
             return Boolean_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the adjacency matrix for the selection of sets
  --   defined by s and the vector v, with i rows and n columns.

    res : Boolean_Matrices.Matrix(1..i,1..n);
    vrs : Set;

  begin
    for k in 1..i loop
      vrs := s(integer32(v(k)));
      for j in 1..n loop
        res(i,j) := Is_In(vrs,natural32(j));
      end loop;
    end loop;
    return res;
  end Adjacency_Matrix;
 
  function Consistent
             ( s : Array_of_Sets; v : Vector; i,n : integer32 )
             return boolean is

  -- DESCRIPTION :
  --   Returns true if the selection of the sets will lead to a linear
  --   system that is consistent, i.e.: that has a solution.

    mat : constant Boolean_Matrices.Matrix(1..i,1..n)
        := Adjacency_Matrix(s,v,i,n);
    prm : Standard_Natural_Vectors.Vector(1..i);
    prm_size : integer32;
    match : Standard_Integer_Vectors.Vector(1..i);
    match_size : natural32;

  begin
    Pivot_Selection.Greedy_Search(mat,prm,prm_size);
    if prm_size = i then
      return true;
    else
      Pivot_Selection.maxBPM(mat,match,match_size);
      return (integer32(match_size) = i);
    end if;
  end Consistent;

  function Matching_Permanent
             ( a : matrix; s : Array_of_Sets; v : Vector;
               i,n : integer32 ) return integer32 is

  -- ALGORITHM : Row expansion without memory.

  begin
    if i = n+1 then
      return 1;
    else
      declare
        res : integer32 := 0;
        vv : Vector(v'range) := v;
      begin
        for j in a'range(2) loop
          if a(i,j) /= 0 then
            vv(i) := natural32(j);
            if Consistent(s,vv,i,n)
             then res := res + a(i,j)*Matching_Permanent(a,s,vv,i+1,n);
            end if;
          end if;
        end loop;
        return res;
      end;
    end if;
  end Matching_Permanent;

  function Matching_Permanent
             ( dst : Degree_Sets_Table ) return integer32 is

    v : constant Vector(1..dst.n) := (1..dst.n => 0);

  begin
    return Matching_Permanent(dst.a,dst.s,v,1,dst.n);
  end Matching_Permanent;

-- DESTRUCTOR :

  procedure Clear ( ase : in out Array_of_Sets ) is
  begin
    for i in ase'range loop
      Clear(ase(i));
    end loop;
  end Clear;

  procedure Clear ( dst : in out Degree_Sets_Table ) is
  begin
    Clear(dst.s);
  end Clear;

end Degree_Sets_Tables;
