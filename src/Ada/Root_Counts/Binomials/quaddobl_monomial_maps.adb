with unchecked_deallocation;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Binomial_Varieties;

package body QuadDobl_Monomial_Maps is

-- CREATORS :

  function Create ( n,d : integer32;
                    c : QuadDobl_Complex_Vectors.Vector;
                    v : Standard_Integer_VecVecs.VecVec )
                  return Monomial_Map is

    res : Monomial_Map(n);

  begin
    res.d := d;
    res.c := c;
    for i in res.v'range loop
      res.v(i) := new Standard_Integer_Vectors.Vector'(v(i).all);
    end loop;
    return res;
  end Create;

  function Create ( maps : Monomial_Map_Array ) return Monomial_Map_List is

    res,res_last : Monomial_Map_List;

  begin
    for i in maps'range loop
      Append(res,res_last,maps(i));
    end loop;
    return res;
  end Create;

  function Create ( maps : Monomial_Map_List ) return Monomial_Map_Array is

    res : Monomial_Map_Array(1..integer32(Length_Of(maps)));
    tmp : Monomial_Map_List := maps;
    link_to_map : Link_to_Monomial_Map;

  begin
    for i in res'range loop
      link_to_map := Head_Of(tmp);
      Copy(link_to_map,res(i));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

-- SELECTORS :

  function Is_Zero ( c : Complex_Number ) return boolean is

    zero : constant quad_double := create(0.0);

  begin
    if not (REAL_PART(c) = zero) then
      return false;
    elsif not (IMAG_PART(c) = zero) then
      return false;
    else
      return true;
    end if;
  end Is_Zero;

  function Is_One ( c : Complex_Number ) return boolean is

    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);

  begin
    if not (REAL_PART(c) = one) then
      return false;
    elsif not (IMAG_PART(c) = zero) then
      return false;
    else
      return true;
    end if;
  end Is_One;

  function Is_Equal ( m1,m2 : Monomial_Map ) return boolean is

    d : quad_double;
    one : constant quad_double := create(1.0);

  begin
    if m1.d /= m2.d then
      return false;
    elsif m1.n /= m2.n then
      return false;
    else
      for i in m1.v'range loop
        for j in m1.v(i)'range loop
          if m1.v(i)(j) /= m2.v(i)(j)
           then return false;
          end if;
        end loop;
        d := (REAL_PART(m1.c(i)) - REAL_PART(m2.c(i)));
        if abs(d) + one /= one
         then return false;
        end if;
        d := (IMAG_PART(m1.c(i)) - IMAG_PART(m2.c(i)));
        if abs(d) + one /= one 
         then return false;
        end if;
      end loop;
      return true;
    end if;
  end Is_Equal;

  function Is_In ( maps : Monomial_Map_List; m : Monomial_Map )
                 return boolean is

    tmp : Monomial_Map_List := maps;
 
  begin
    while not Is_Null(tmp) loop
      if Is_Equal(Head_Of(tmp).all,m)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  function Tropisms
             ( map : Monomial_Map ) return Standard_Integer_VecVecs.VecVec is

    res : Standard_Integer_VecVecs.VecVec(1..map.d);

  begin
    for i in res'range loop
      res(i) := new Standard_Integer_Vectors.Vector(1..map.n);
    end loop;
    for i in map.v'range loop
      for j in map.v(i)'range loop
        res(j)(i) := map.v(i)(j);
      end loop;
    end loop;
    return res;
  end Tropisms;

  function Tropism_Configuration
             ( map : Monomial_Map ) return Standard_Integer_Matrices.Matrix is

    res : Standard_Integer_Matrices.Matrix(1..map.d,1..map.n);
    lv : Standard_Integer_Vectors.Link_to_Vector;

  begin
    for j in map.v'range loop
      lv := map.v(j);
      for i in lv'range loop
        res(i,j) := lv(i);
      end loop;
    end loop;
    return res;
  end Tropism_Configuration;

  function Degree ( map : Monomial_Map ) return natural32 is

    T : Standard_Integer_VecVecs.VecVec(1..map.d) := Tropisms(map);
    A : Standard_Integer_Matrices.Matrix(1..map.n,1..map.d);

  begin
    for j in 1..map.d loop
      for i in 1..map.n loop
        A(i,j) := T(j)(i);
      end loop;
    end loop;
    Standard_Integer_VecVecs.Clear(T);
    return Standard_Binomial_Varieties.Degree(A);
  end Degree;

  function Degrees ( maps : Monomial_Map_List )
                   return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(1..integer32(Length_Of(maps)));
    tmp : Monomial_Map_List := maps;
    link_to_map : Link_to_Monomial_Map;

  begin
    for i in res'range loop
      link_to_map := Head_Of(tmp);
      res(i) := Degree(link_to_map.all);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Degrees;

  function Top_Dimension ( maps : Monomial_Map_Array ) return natural32 is

    res : integer32 := maps(maps'first).d;

  begin
    for i in maps'first+1..maps'last loop
      if maps(i).d > res
       then res := maps(i).d;
      end if;
    end loop;
    return natural32(res);
  end Top_Dimension;

  function Top_Dimension ( maps : Monomial_Map_List ) return natural32 is

    res : integer32 := 0;
    tmp : Monomial_Map_List;

  begin
    if not Is_Null(maps) then
      res := Head_Of(maps).d;
      tmp := Tail_Of(maps);
      while not Is_Null(tmp) loop
        if Head_Of(tmp).d > res
         then res := Head_Of(tmp).d;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    return natural32(res);
  end Top_Dimension;

  function Lengths ( maps : Array_of_Monomial_Map_Lists )
                   return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(maps'range);

  begin
    for i in res'range loop
      res(i) := Length_Of(maps(i));
    end loop;
    return res;
  end Lengths;

-- CONSTRUCTORS :

  procedure Copy ( m1 : in Monomial_Map; m2 : out Monomial_Map ) is
  begin
    m2.d := m1.d;
    m2.c := m1.c;
    for i in m2.v'range loop
      m2.v(i) := new Standard_Integer_Vectors.Vector'(m1.v(i).all);
    end loop;
  end Copy;
  
  procedure Copy ( m1 : in Link_to_Monomial_Map;
                   m2 : out Link_to_Monomial_Map ) is

    map : Monomial_Map(m1.n);

  begin
    Copy(m1.all,map);
    m2 := new Monomial_Map'(map);
  end Copy;

  procedure Append ( first,last : in out Monomial_Map_List;
                     map : in Link_to_Monomial_Map ) is
  begin
    Append(first,last,map.all);
  end Append;

  procedure Append ( first,last : in out Monomial_Map_List;
                     map : in Monomial_Map ) is

    copy_map : Monomial_Map(map.n);
    link_to_copy_map : Link_to_Monomial_Map;

  begin
    Copy(map,copy_map);
    link_to_copy_map := new Monomial_Map'(copy_map);
    if Is_Null(first) then
      Construct(link_to_copy_map,first);
      last := first;
    else
      declare
        tmp : Monomial_Map_List;
      begin
        Construct(link_to_copy_map,tmp);
        Swap_Tail(last,tmp);
        last := Tail_Of(last);
      end;
    end if;
  end Append;

  procedure Concatenate ( from : in Monomial_Map_List;
                          first,last : in out Monomial_Map_List ) is

    tmp : Monomial_Map_List := from;
    link_to_map : Link_to_Monomial_Map;

  begin
    while not Is_Null(tmp) loop
      link_to_map := Head_Of(tmp);
      Append(first,last,link_to_map.all);
      tmp := Tail_Of(tmp);
    end loop;
  end Concatenate;

  procedure Concatenate ( from : in Monomial_Map_Array;
                          first,last : in out Monomial_Map_List ) is
  begin
    for i in from'range loop
      Append(first,last,from(i).all);
    end loop;
  end Concatenate;

-- DESTRUCTORS :

  procedure Clear ( lmm : in out Link_to_Monomial_Map ) is

    procedure free is
      new unchecked_deallocation(Monomial_Map,Link_to_Monomial_Map);

  begin
    for i in lmm.v'range loop
      Standard_Integer_Vectors.Clear(lmm.v(i));
    end loop;
    free(lmm);
  end Clear;

  procedure Clear ( amm : in out Monomial_Map_Array ) is
  begin
    for i in amm'range loop
      Clear(amm(i));
    end loop;
  end Clear;

  procedure Clear ( amm : in out Link_to_Monomial_Map_Array ) is

    procedure free is
      new unchecked_deallocation(Monomial_Map_Array,
                                 Link_to_Monomial_Map_Array);

  begin
    for i in amm'range loop
      Clear(amm(i));
    end loop;
    free(amm);
  end Clear;

  procedure Clear ( lm : in out Monomial_Map_List ) is

    tmp : Monomial_Map_List := lm;
    link_to_map : Link_to_Monomial_Map;

  begin
    while not Is_Null(tmp) loop
      link_to_map := Head_Of(tmp);
      Clear(link_to_map);
      tmp := Tail_Of(tmp);
    end loop;
    List_of_Monomial_Maps.Clear(List_of_Monomial_Maps.List(lm));
  end Clear;

  procedure Clear ( maps : in out Array_of_Monomial_Map_Lists ) is
  begin
    for i in maps'range loop
      Clear(maps(i));
    end loop;
  end Clear;

  procedure Clear ( maps : in out Link_to_Array_of_Monomial_Map_Lists ) is

    procedure free is 
      new unchecked_deallocation(Array_of_Monomial_Map_Lists,
                                 Link_to_Array_of_Monomial_Map_Lists);

  begin
    if maps /= null
     then Clear(maps.all);
    end if;
    free(maps);
  end Clear;

end QuadDobl_Monomial_Maps;
