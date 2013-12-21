with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Monomial_Maps_io;         use Standard_Monomial_Maps_io;

package body Standard_Monomial_Map_Filters is

  function Is_Pure_Dimensional
              ( maps : Monomial_Map_List ) return boolean is

    dim : integer32;
    tmp : Monomial_Map_List;

  begin
    if not Is_Null(maps) then
      dim := Head_Of(maps).d;
      tmp := Tail_Of(maps);
      while not Is_Null(tmp) loop
        if Head_Of(tmp).d /= dim 
         then return false;
         else tmp := Tail_Of(tmp);
        end if;
     end loop;
    end if;
    return true;
  end Is_Pure_Dimensional;

  function Pure_Dimensional_Maps
              ( maps : Monomial_Map_List; dim : natural32 )
              return Monomial_Map_List is

    res,res_last : Monomial_Map_List;
    tmp : Monomial_Map_List := maps;
    link_to_map : Link_to_Monomial_Map;

  begin
    while not Is_Null(tmp) loop
      link_to_map := Head_Of(tmp);
      if link_to_map.d = integer32(dim)
       then Append(res,res_last,link_to_map.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Pure_Dimensional_Maps;

  function Pure_Dimensional_Maps
              ( maps : Monomial_Map_List )
              return Array_of_Monomial_Map_Lists is

   res : Array_of_Monomial_Map_Lists(0..integer32(Top_Dimension(maps)));

  begin
    for i in res'range loop
      res(i) := Pure_Dimensional_Maps(maps,natural32(i));
    end loop;
    return res;
  end Pure_Dimensional_Maps;

  function Remove_Duplicates
               ( maps : Monomial_Map_List ) return Monomial_Map_List is

    res,res_last : Monomial_Map_List;
    tmp : Monomial_Map_List := maps;
    link_to_map : Link_to_Monomial_Map;

  begin
    while not Is_Null(tmp) loop
      link_to_map := Head_Of(tmp);
      if not Is_In(res,link_to_map.all)
       then Append(res,res_last,link_to_map.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Duplicates;

  procedure Filter_Duplicates ( maps : in out Monomial_Map_List ) is

    f : constant Monomial_Map_List := Remove_Duplicates(maps);

  begin
    Clear(maps);
    maps := f;
  end Filter_Duplicates;

  procedure Filter_Duplicates ( maps : in out Array_of_Monomial_Map_Lists ) is
  begin
    for i in maps'range loop
      Filter_Duplicates(maps(i));
    end loop;
  end Filter_Duplicates;

  function Is_Free ( map : Monomial_Map ) return boolean is

    cnt_free : integer32 := 0;
    v : Standard_Integer_Vectors.Link_to_Vector;

  begin
    for i in map.c'range loop
      if not Is_Zero(map.c(i)) then
        if not Is_One(map.c(i)) then
          return false;
        else
          cnt_free := cnt_free + 1;
          v := map.v(i);
          for j in v'range loop
            if j = cnt_free then
              if v(j) /= 1
               then return false; 
              end if;
            else
              if v(j) /= 0
               then return false;
              end if;
            end if;
          end loop;
        end if;
      end if;
    end loop;
    return true;
  end Is_Free;

  function Has_Zeroes ( map : Monomial_Map ) return boolean is
  begin
    for i in map.c'range loop
      if Is_Zero(map.c(i))
       then return true;
      end if;
    end loop;
    return false;
  end Has_Zeroes;

  function Is_Free_Submap ( m1,m2 : Monomial_Map ) return boolean is
  begin
    if m1.d > m2.d then
      return false;
    else
      for i in m1.c'range loop
        if not Is_Zero(m1.c(i)) then
          if Is_Zero(m2.c(i))
           then return false;
          end if;
        end if;
      end loop;
    end if;
    return true;
  end Is_Free_Submap;

  function Remove_Free_Submaps
             ( m1 : Monomial_Map_List; m2 : Monomial_Map )
             return Monomial_Map_List is

    res,res_last : Monomial_Map_List;
    tmp : Monomial_Map_List := m1;
    map : Link_to_Monomial_Map;

  begin
    while not Is_Null(tmp) loop
      map := Head_Of(tmp);
      if not Is_Free(map.all) then
        Append(res,res_last,map.all);
      elsif not Is_Free_Submap(map.all,m2) then
        Append(res,res_last,map.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Free_Submaps;

  function Remove_Free_Submaps
             ( m1,m2 : Monomial_Map_List ) return Monomial_Map_List is

    res,res_last : Monomial_Map_List;
    tmp : Monomial_Map_List := m2;
    map : Link_to_Monomial_Map;

  begin
    Concatenate(m1,res,res_last);
    while not Is_Null(tmp) loop
      map := Head_Of(tmp);
      if Is_Free(map.all) then
        declare
          res_new : constant Monomial_Map_List
                  := Remove_Free_Submaps(res,map.all);
        begin
          Clear(res);
          res := res_new;
        end;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Free_Submaps;

  procedure Filter_Free_Submaps
              ( m1 : in out Monomial_Map_List;
                m2 : in Monomial_Map_List ) is

    m1filtered : constant Monomial_Map_List := Remove_Free_Submaps(m1,m2);
 
  begin
    Clear(m1);
    m1 := m1filtered;
  end Filter_Free_Submaps;

  procedure Filter_Free_Submaps
              ( maps : in out Array_of_Monomial_Map_Lists ) is
  begin
    for i in reverse maps'range loop
      if not Is_Null(maps(i)) then
        for j in maps'first..(i-1) loop
          if not Is_Null(maps(j)) 
           then Filter_Free_Submaps(maps(j),maps(i));
          end if;
        end loop;
      end if;
    end loop;
  end Filter_Free_Submaps;

  function Is_Zero_Submap ( m1,m2 : Monomial_Map ) return boolean is
  begin
    for i in m2.c'range loop
      if Is_Zero(m2.c(i)) then
        if not Is_Zero(m1.c(i))
         then return false;
        end if;
      end if;
    end loop;
    return true;
  end Is_Zero_Submap;

  function Filter ( p : Poly; map : Monomial_Map ) return Poly is

    res : Poly := Null_Poly;
    tol : constant double_float := 1.0E-8;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is

      is_zero : boolean := false;

    begin
      for i in t.dg'range loop
        if t.dg(i) > 0 then
          if AbsVal(map.c(i)) < tol
           then is_zero := true;
          end if;
        end if;
        exit when is_zero;
      end loop;
      if not is_zero
       then Add(res,t);
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Filter;

  function Filter ( p : Laur_Sys; map : Monomial_Map )
                  return Link_to_Laur_Sys is

    res : Link_to_Laur_Sys;
    sys : Laur_Sys(p'range);
    cnt : integer32 := sys'first-1;

  begin
    for i in p'range loop
      declare
        f : constant Poly := Filter(p(i),map);
      begin
        if f /= Null_Poly then
          cnt := cnt + 1;
          sys(cnt) := f;
        end if;
      end;
    end loop;
    res := new Laur_Sys'(sys(sys'first..cnt));
    return res;
  end Filter;

  function Is_Generated_by_Monomials
              ( p : Poly; map : Monomial_Map ) return boolean is

    res : boolean := true;

    procedure Visit ( t : in Term; continue : out boolean ) is

    -- DESCRIPTION :
    --   For t to be generated by the free variables in the map,
    --   there must be at least one index i for which t.dg(i) > 0
    --   corresponding to a zero monomial in the map.

      found : boolean := false;

    begin
      for i in t.dg'range loop
        if t.dg(i) > 0
         then found := Is_Zero(map.c(i));
        end if;
        exit when found;
      end loop;
      if found
       then continue := true;
       else continue := false; res := false;
      end if;
    end Visit;
    procedure Visit_Terms is new Visiting_Iterator(Visit);

  begin
    Visit_Terms(p);
    return res;
  end Is_Generated_by_Monomials;

  function Is_Generated_by_Monomials
              ( p : Laur_Sys; m1,m2 : Monomial_Map ) return boolean is

    res : boolean := true;
    fp : Link_to_Laur_Sys := Filter(p,m2);

  begin
    for i in fp'range loop
      if not Is_Generated_by_Monomials(fp(i),m1)
       then res := false;
      end if;
      exit when (not res);
    end loop;
    Clear(fp);
    return res;
  end Is_Generated_by_Monomials;

  function Is_Free_of_Affine_Submap
              ( p : Laur_Sys; m1,m2 : Monomial_Map ) return boolean is
  begin
    if not Is_Zero_Submap(m1,m2)
     then return false;
     else return Is_Generated_by_Monomials(p,m1,m2);
    end if;
  end Is_Free_of_Affine_Submap;

  function Remove_Free_of_Affine_Submaps
              ( p : Laur_Sys; m1 : Monomial_Map_List; m2 : Monomial_Map )
              return Monomial_Map_List is

    res,res_last : Monomial_Map_List;
    tmp : Monomial_Map_List := m1;
    link_to_map : Link_to_Monomial_Map;

  begin
    while not Is_Null(tmp) loop
      link_to_map := Head_Of(tmp);
      if not Is_Free(link_to_map.all) then
        Append(res,res_last,link_to_map.all);
      elsif not Is_Free_of_Affine_Submap(p,link_to_map.all,m2) then
        Append(res,res_last,link_to_map.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Free_of_Affine_Submaps;

  function Remove_Free_of_Affine_Submaps
             ( p : Laur_Sys; m1,m2 : Monomial_Map_List )
             return Monomial_Map_List is

    res,res_last : Monomial_Map_List;
    tmp : Monomial_Map_List := m2;
    map : Link_to_Monomial_Map;

  begin
    Concatenate(m1,res,res_last);
    while not Is_Null(tmp) loop
      map := Head_Of(tmp);
      if not Is_Free(map.all) then
        declare
          res_new : constant Monomial_Map_List
                  := Remove_Free_of_Affine_Submaps(p,res,map.all);
        begin
          Clear(res);
          res := res_new;
        end;
      end if;
      exit when Is_Null(res);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Free_of_Affine_Submaps;

  procedure Filter_Free_of_Affine_Submaps
              ( p : Laur_Sys; m1 : in out Monomial_Map_List;
                m2 : in Monomial_Map_List ) is

    m1filtered : constant Monomial_Map_List
               := Remove_Free_of_Affine_Submaps(p,m1,m2);
 
  begin
    Clear(m1);
    m1 := m1filtered;
  end Filter_Free_of_Affine_Submaps;

  procedure Filter_Free_of_Affine_Submaps
              ( p : Laur_Sys; maps : in out Array_of_Monomial_Map_Lists ) is
  begin
    for i in reverse maps'range loop
      if not Is_Null(maps(i)) then
        for j in maps'first..(i-1) loop
          if not Is_Null(maps(j)) 
           then Filter_Free_of_Affine_Submaps(p,maps(j),maps(i));
          end if;
        end loop;
      end if;
    end loop;
  end Filter_Free_of_Affine_Submaps;

  function Is_Affine_Submap
              ( p : Laur_Sys; m1,m2 : Monomial_Map ) return boolean is
  begin
    if not Is_Zero_Submap(m1,m2) then
      return false;
    else
      declare
        res : boolean := true;
        f1 : Link_to_Laur_Sys := Filter(p,m1);
        f2 : Link_to_Laur_Sys := Filter(p,m2);
        found : boolean;
      begin
        for j in f2'range loop
          found := false;
          for i in f1'range loop
            found := Equal(f1(i),f2(j));
            exit when found;
          end loop;
          if not found then
            if not Is_Generated_by_Monomials(f2(j),m1)
             then res := false;
            end if;
          end if;
          exit when not res;
        end loop;
        Clear(f1); Clear(f2);
        return res;
      end;
    end if;
  end Is_Affine_Submap;

  function Remove_Affine_Submaps
              ( p : Laur_Sys; m1 : Monomial_Map_List; m2 : Monomial_Map )
              return Monomial_Map_List is

    res,res_last : Monomial_Map_List;
    tmp : Monomial_Map_List := m1;
    link_to_map : Link_to_Monomial_Map;

  begin
    while not Is_Null(tmp) loop
      link_to_map := Head_Of(tmp);
      if not Is_Affine_Submap(p,link_to_map.all,m2)
       then Append(res,res_last,link_to_map.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Affine_Submaps;

  function Remove_Affine_Submaps
              ( p : Laur_Sys; m1,m2 : Monomial_Map_List )
              return Monomial_Map_List is

    res,res_last : Monomial_Map_List;
    tmp : Monomial_Map_List := m2;
    map : Link_to_Monomial_Map;

  begin
    Concatenate(m1,res,res_last);
    while not Is_Null(tmp) loop
      map := Head_Of(tmp);
      if not Is_Free(map.all) then
        declare
          res_new : constant Monomial_Map_List
                  := Remove_Affine_Submaps(p,res,map.all);
        begin
          Clear(res);
          res := res_new;
        end;
      end if;
      exit when Is_Null(res);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Affine_Submaps;

  procedure Filter_Affine_Submaps
              ( p : Laur_Sys; m1 : in out Monomial_Map_List;
                m2 : in Monomial_Map_List ) is

    m1filtered : constant Monomial_Map_List
               := Remove_Affine_Submaps(p,m1,m2);
 
  begin
    Clear(m1);
    m1 := m1filtered;
  end Filter_Affine_Submaps;

  procedure Filter_Affine_Submaps
              ( p : Laur_Sys; maps : in out Array_of_Monomial_Map_Lists ) is
  begin
    for i in reverse maps'range loop
      if not Is_Null(maps(i)) then
        for j in maps'first..(i-1) loop
          if not Is_Null(maps(j)) 
           then Filter_Affine_Submaps(p,maps(j),maps(i));
          end if;
        end loop;
      end if;
    end loop;
  end Filter_Affine_Submaps;

-- MAIN DRIVERS :

  procedure Silent_Filter
               ( p : in Laur_Sys; maps : in Monomial_Map_List;
                 c : out Link_to_Array_of_Monomial_Map_Lists ) is

    td : constant natural32 := Top_Dimension(maps);
    components : Array_of_Monomial_Map_Lists(0..integer32(td))
               := Pure_Dimensional_Maps(maps);

  begin
    Filter_Duplicates(components);
    Filter_Free_Submaps(components);
    Filter_Free_of_Affine_Submaps(p,components);
    Filter_Affine_Submaps(p,components);
    c := new Array_of_Monomial_Map_Lists'(components);
  end Silent_Filter;

  procedure Reporting_Filter
              ( file : in file_type;
                p : in Laur_Sys; maps : in Monomial_Map_List;
                c : out Link_to_Array_of_Monomial_Map_Lists ) is

    td : constant natural32 := Top_Dimension(maps);
    components : Array_of_Monomial_Map_Lists(0..integer32(td))
               := Pure_Dimensional_Maps(maps);

  begin
    put(file,"The top dimension : "); put(file,td,1); new_line(file);
    if Is_Pure_Dimensional(maps)
     then put_line(file,"The list of maps is pure dimensional.");
     else put_line(file,"The list of maps is not pure dimensional.");
    end if;
    Write_Lengths(file,components);
    new_line(file);
    put_line(file,"filtering the duplicate maps out ...");
    Filter_Duplicates(components);
    Write_Lengths(file,components);
    put_line(file,"the maps after filtering duplicates : ");
    put(file,components);
    new_line(file);
    put_line(file,"filtering the free submaps out ...");
    Filter_Free_Submaps(components);
    Write_Lengths(file,components);
    put_line(file,"the maps after filtering free submaps : ");
    put(file,components);
    new_line(file);
    put_line(file,"filtering free maps as submaps of affine maps ...");
    Filter_Free_of_Affine_Submaps(p,components);
    Write_Lengths(file,components);
    put_line(file,"the maps after filtering free submaps : ");
    put(file,components);
    new_line(file);
    put_line(file,"filtering maps as affine submaps of affine maps ...");
    Filter_Affine_Submaps(p,components);
    Write_Lengths(file,components);
    put_line(file,"the maps after filtering affine submaps : ");
    put(file,components);
    c := new Array_of_Monomial_Map_Lists'(components);
  end Reporting_Filter;

end Standard_Monomial_Map_Filters;
