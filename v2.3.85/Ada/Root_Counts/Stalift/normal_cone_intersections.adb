with Inner_Normal_Cones;                 use Inner_Normal_Cones;

package body Normal_Cone_Intersections is

-- AUXILIARY :

  function Get ( L : List; i : integer32 ) return Vector is

  -- DESCRIPTION :
  --   Returns the ith point vector in the list.

    tmp : List := L;
    cnt : integer32 := 0;
    nll : constant Vector(0..0) := (0..0 => 0);

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      if cnt = i
       then return Head_Of(tmp).all;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return nll;
  end Get;

-- CONSTRUCTORS :

  function Number_of_Cones ( L : Array_of_Lists; i : integer32 )
                           return integer32 is

    res : integer32 := 0;

  begin
    for j in L'range loop
      if j /= i
       then res := res + integer32(Length_Of(L(j)));
      end if;
    end loop;
    return res;
  end Number_of_Cones;

  function Lengths ( L : Array_of_Lists; i : integer32 ) return Vector is

    res : Vector(L'range);

  begin
    res(res'first) := 1;
    for j in l'first..(i-1) loop
      res(j+1) := res(j) + integer32(Length_Of(l(j)));
    end loop;
    for j in (i+1)..l'last loop
      res(j) := res(j-1) + integer32(Length_Of(l(j)));
    end loop;
    return res;
  end Lengths;

  function Create ( L : Array_of_Lists; g : List; i : integer32 ) 
                  return Intersection_Matrix is

    n : constant integer32 := L'length - 1;
    m : constant integer32 := integer32(Length_Of(g));
    ll : constant Vector := Lengths(L,i);
    nc : constant integer32 := ll(ll'last)-1;
    res : Intersection_Matrix(n,m,nc);

  begin
    res.sv := ll(ll'first..ll'last-1);
    for j in l'range loop
      if j /= i then
        declare
          ind : integer32 := j;
          tmpl : List := L(j);
          cntl : integer32 := 0;
        begin
          if ind > i
           then ind := ind - 1;
          end if;
          while not Is_Null(tmpl) loop
            declare
              cone : constant Matrix
                   := Inner_Normal_Cone(l(j),Head_Of(tmpl).all);
              tmpg : List := g;
              cntg,sum : integer32 := 0;
            begin
             -- put_line("The inequalities of the normal cone : "); put(cone);
              while not Is_Null(tmpg) loop
                cntg := cntg + 1;
               -- put(" "); put(Head_Of(tmpg).all); 
                if Satisfies(cone,Head_Of(tmpg).all)
                 then res.im(cntg,res.sv(ind)+cntl) := 1; sum := sum + 1;
                     -- put_line(" satisfies.");
                 else res.im(cntg,res.sv(ind)+cntl) := 0;
                     -- put_line(" does not satisfy.");
                end if;
                tmpg := Tail_Of(tmpg);
              end loop;
              res.im(0,res.sv(ind)+cntl) := sum;
            end;
            tmpl := Tail_Of(tmpl);
            cntl := cntl + 1;
          end loop;
        end;
      end if;
    end loop;
    return res;
  end Create;

-- ELEMENTARY SELECTORS :

  function Is_In ( ima : Intersection_Matrix; i,j,k : integer32 )
                 return boolean is
  begin
    if ima.im(i,ima.sv(j)+k-1) = 1
     then return true;
     else return false;
    end if;
  end Is_In;

  function Maximal_Column ( ima : Intersection_Matrix ) return integer32 is

    res : integer32 := ima.im'first(2);
    max : integer32 := ima.im(0,ima.im'first(2));

  begin
    for j in ima.im'first(2)+1..ima.im'last(2) loop
      if ima.im(0,j) > max
       then max := ima.im(0,j); res := j;
      end if;
    end loop;
    return res;
  end Maximal_Column;

  function Component ( ima : Intersection_Matrix; column : integer32 )
                     return integer32 is
  begin
    for i in ima.sv'range loop
      if ima.sv(i) > column
       then return i-1;
      end if;
    end loop;
    return ima.sv'last;
  end Component;

  function Length ( ima : Intersection_Matrix; i : integer32 )
                  return integer32 is
  begin
    if i < ima.sv'last
     then return (ima.sv(i+1) - ima.sv(i));
     else return (ima.im'last(2) - ima.sv(i) + 1);
    end if;
  end Length;

  function Row_Sum ( ima : Intersection_Matrix; i,j : integer32 )
                   return integer32 is

    res : integer32 := 0;
    lst : integer32;

  begin
    if j < ima.sv'last 
     then lst := ima.sv(j+1)-1;
     else lst := ima.im'last(2);
    end if;
    for k in ima.sv(j)..lst loop
      res := res + ima.im(i,k);
    end loop;
    return res;
  end Row_Sum;

-- ENUMERATING COMPLEMENTARY COLUMNS :

  procedure Complementary_Columns ( ima : in Intersection_Matrix ) is

    acc : Standard_Integer_Vectors.Vector(ima.sv'range) := (ima.sv'range => 0);
     -- acc(j) = 0 if no cone from jth component has been chosen yet,
     --        = k if kth cone from jth component is selected.

    continue : boolean := true;

    function Is_In ( acc : in Vector; i : in integer32 ) return boolean is

    -- DESCRIPTION :
    --   Returns true if the ith generator already belongs to one of the
    --   chosen cones in acc.

    begin
      for j in acc'range loop               -- enumerate the components
        if acc(j) /= 0 then                 -- acc(j) = cone selected
          if Is_In(ima,i,j,acc(j))          -- is in selected cone ?
           then return true;
          end if;
        end if;
      end loop;
      return false;
    end Is_In;

    procedure Select_Columns ( i : in integer32 ) is

    -- DESCRIPTION :
    --   Selects all columns such that the ith generator belongs to the
    --   collection of columns.

      lst : integer32;

    begin
      if i > ima.im'last(1) then
        Process(acc,continue);
      else
        if Is_In(acc,i) then
          Select_Columns(i+1);
        else
          for j in acc'range loop
            if acc(j) = 0 then
              if j < ima.sv'last
               then lst := ima.sv(j+1)-1;
               else lst := ima.im'last(2);
              end if;
              for k in ima.sv(j)..lst loop
                if ima.im(i,k) = 1 then
                  acc(j) := k-ima.sv(j)+1;
                  Select_Columns(i+1);
                  acc(j) := 0;
                end if;
                exit when not continue;
              end loop;
            end if;
            exit when not continue;
          end loop;
        end if;
      end if;
    end Select_Columns;

  begin
    Select_Columns(1);
  end Complementary_Columns;

  function Partition ( ima : Intersection_Matrix; cols : Vector; g : List )
                     return Array_of_Lists is

    res,res_last : Array_of_Lists(cols'range);
    tmp : List := g;

    procedure Search_and_Update ( v : in Vector; i : in integer32 ) is

    -- DESCRIPTION :
    --   Given the ith generator v from the list g, this procedures searches
    --   for the first cone that contains it and updates the partition.

      found : boolean := false;

    begin
      for j in cols'range loop
        if (cols(j) /= 0) and then Is_In(ima,i,j,cols(j)) then
          found := true;
          Append(res(j),res_last(j),v);
        end if;
        exit when found;
      end loop;
    end Search_and_Update;

  begin
    for i in ima.im'first(1)+1..ima.im'last(1) loop
      Search_and_Update(Head_Of(tmp).all,i);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Partition;

  function Partition_in_Union
             ( partg,points : Array_of_Lists; i : integer32;
               cols : Vector ) return boolean is

  -- ALGORITHM : lexicographic enumeration of all couples of lists in partg,
  --             with each time a check whether it belongs to the union of
  -- the normal cones as given by the set of complementary columns.

    function Index ( j : integer32 ) return integer32 is
    begin
      if j < i
       then return j;
       else return j+1;
      end if;
    end Index;

  begin
    for k1 in partg'range loop
      if not Is_Null(partg(k1)) then
        for k2 in (k1+1)..partg'last loop
          if not Is_Null(partg(k2)) then
            declare
              ind1 : constant integer32 := Index(k1);
              ind2 : constant integer32 := Index(k2);
              x1 : constant Vector := Get(points(ind1),cols(k1));
              x2 : constant Vector := Get(points(ind2),cols(k2));
              ic1 : constant Matrix := Inner_Normal_Cone(points(ind1),x1);
              ic2 : constant Matrix := Inner_Normal_Cone(points(ind2),x2);
            begin
              if not In_Union(partg(k1),partg(k2),ic1,ic2)
               then return false;
              end if;
            end;
          end if;
        end loop;
      end if;
    end loop;
    return true;
  end Partition_in_Union;

  function Contained_in_Union
             ( l : Array_of_Lists; i : integer32; g : List;
               ima : Intersection_Matrix; cols : Vector ) return boolean is

    p : Array_of_Lists(cols'range) := Partition(ima,cols,g);
    res : constant boolean := Partition_in_Union(p,l,i,cols);

  begin
    Deep_Clear(p);
    return res;
  end Contained_in_Union;

-- FINAL TARGET ROUTINE :

  function Contained_in_Union
             ( l : Array_of_Lists; i : integer32; g : List;
               ima : Intersection_Matrix ) return boolean is

    res : boolean := false;

    procedure Examin_Selection ( cols : in Vector; continue : out boolean ) is
    begin
      res := Contained_in_Union(l,i,g,ima,cols);
      continue := not res;
    end Examin_Selection;
    procedure Enumerate_Complementary_Columns is
      new Complementary_Columns(Examin_Selection);

  begin
    Enumerate_Complementary_Columns(ima);
    return res;
  end Contained_in_Union;

end Normal_Cone_Intersections;
