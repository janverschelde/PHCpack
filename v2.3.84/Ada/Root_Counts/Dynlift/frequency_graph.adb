with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package body Frequency_Graph is

-- CREATORS :

  function Occurrences ( i : integer32; l : List ) return integer32 is

    res : integer32 := 0;
    tmp : List := l;
    pt : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if pt(i) /= 0
       then res := res + 1;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Occurrences;

  function Graph ( n : integer32; supports : Array_of_Lists ) return Matrix is

    res : Matrix(1..n,supports'range);

  begin
    for i in 1..n loop
      for j in supports'range loop
        res(i,j) := Occurrences(i,supports(j));
      end loop;
    end loop;
    return res;
  end Graph;

-- MODIFIER :

  procedure Ignore ( m : in out Matrix; point : in Vector ) is
  begin
    for i in point'range loop
      if point(i) /= 0 then
        for j in m'range(2) loop
          m(i,j) := 1;
        end loop;
      end if;
    end loop;
  end Ignore;

-- SELECTORS :

  function Occurrence ( i : integer32; m : Matrix ) return integer32 is

    res : integer32 := 0;

  begin
    for j in m'range(2) loop
      if m(i,j) /= 0
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Occurrence;

  function Occurrence
             ( i : integer32; m : Matrix; col : integer32; perm : Vector )
             return integer32 is

    res : integer32 := 0;

  begin
    for j in col+1..m'last(2) loop
      if m(i,perm(j)) /= 0
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Occurrence;

  function Occurrence ( v : Vector; m : Matrix ) return integer32 is

    min : integer32 := 1000000;
    occ : integer32;

  begin
    for i in v'range loop
      if v(i) /= 0 then
        occ := Occurrence(i,m);
        if occ < min
         then min := occ;
        end if;
      end if;
    end loop;
    return min;
  end Occurrence;

  function Occurrence
             ( v : Vector; m : Matrix; col : integer32; perm : Vector )
             return integer32 is

    min : integer32 := 1000000;
    occ : integer32;

  begin
    for i in v'range loop
      if v(i) /= 0
       then occ := Occurrence(i,m,col,perm);
            if occ < min
             then min := occ;
            end if;
      end if;
    end loop;
    return min;
  end Occurrence;

  function Lowest_Occurrence
             ( Vec : VecVec; start : integer32; m : Matrix )
             return integer32 is

    res : integer32 := start;
    min : integer32 := Occurrence(vec(start).all,m);
    occ : integer32;

  begin
    for i in start+1..vec'last loop
      occ := Occurrence(vec(i).all,m);
      if occ < min
       then min := occ; res := i;
      end if;
    end loop;
    return res;
  end Lowest_Occurrence;

  function Lowest_Occurrence
               ( vec : VecVec; start : integer32; m : Matrix;
                 col : integer32; perm : Vector ) return integer32 is

    res : integer32 := start;
    min : integer32 := Occurrence(vec(start).all,m,col,perm);
    occ : integer32;

  begin
    for i in start+1..vec'last loop
      occ := Occurrence(vec(i).all,m,col,perm);
      if occ < min
       then min := occ; res := i;
      end if;
    end loop;
    return res;
  end Lowest_Occurrence;

-- CONSTRUCTORS :

  function Sort ( L : List; m : Matrix ) return List is

    res : List;
    vec : VecVec(1..integer32(Length_Of(L))) := Deep_Create(L);
    tmp : Link_to_Vector;
    low : integer32;

  begin
    if Length_Of(L) <= 1 then
      Copy(L,res);
    else
      for i in vec'first..vec'last-1 loop
        low := Lowest_Occurrence(vec,i,m);
        if low /= i then
          tmp := vec(i);
          vec(i) := vec(low);
          vec(low) := tmp;
        end if;
      end loop;
      res := Deep_Create(vec);
      Clear(vec);
    end if;
    return res;
  end Sort;

  procedure Sort ( L : in out List; m : in Matrix ) is

    res : List := Sort(L,m);

  begin
    Copy(res,L); Deep_Clear(res);
  end Sort;

  function Sort ( L : List; m : Matrix; col : integer32; perm : Vector )
                return List is

    res : List;
    vec : VecVec(1..integer32(Length_Of(L))) := Deep_Create(l);
    tmp : Link_to_Vector;
    low : integer32;

  begin
    if Length_Of(L) <= 1 then
      Copy(l,res);
    else
      for i in vec'first..vec'last-1 loop
        low := Lowest_Occurrence(vec,i,m,col,perm);
        if low /= i then
          tmp := vec(i);
          vec(i) := vec(low);
          vec(low) := tmp;
        end if;
      end loop;
      res := Deep_Create(vec);
      Clear(vec);
    end if;
    return res;
  end Sort;

  procedure Sort ( L : in out List; m : in Matrix;
                   col : in integer32; perm : in Vector ) is

    res : List := Sort(L,m,col,perm);

  begin
    Copy(res,L); Deep_Clear(res);
  end Sort;

end Frequency_Graph;
