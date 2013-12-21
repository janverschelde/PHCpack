with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package body Monomial_Maps_Container is

  the_maps : Link_to_Array_of_Monomial_Map_Lists;

  procedure Initialize ( sols : in Array_of_Monomial_Map_Lists ) is
  begin
    Clear(the_maps);
    the_maps := new Array_of_Monomial_Map_Lists'(sols);
  end Initialize;

  function Retrieve return Link_to_Array_of_Monomial_Map_Lists is
  begin
    return the_maps;
  end Retrieve;

  function Top_Dimension return integer32 is
  begin
    if the_maps = null
     then return -1;
     else return the_maps'last;
    end if;
  end Top_Dimension;

  function Number_of_Maps ( d : integer32 ) return integer32 is
  begin
    if the_maps = null then
      return -1;
    elsif d > the_maps'last then
      return 0;
    elsif d < the_maps'first then
      return 0;
    else
      return integer32(Length_Of(the_maps(d)));
    end if;
  end Number_of_Maps;

  function Retrieve_Map ( d,k : integer32 ) return Link_to_Monomial_Map is
  begin
    if the_maps = null then
      return null;
    elsif d > the_maps'last then
      return null;
    elsif d < the_maps'first then
      return null;
    else
      declare
        maps : Monomial_Map_List := the_maps(d);
      begin
        for i in 1..(k-1) loop
          exit when Is_Null(maps);
          maps := Tail_Of(maps);
        end loop;
        if Is_Null(maps)
         then return null;
         else return Head_Of(maps);
        end if;
      end;
    end if;
  end Retrieve_Map;

  function Degree ( d,k : integer32 ) return integer32 is

    map : constant Link_to_Monomial_Map := Retrieve_Map(d,k);

  begin
    if map = null
     then return -1;
     else return integer32(Degree(map.all));
    end if;
  end Degree;

  function Coefficients
             ( d,k : integer32 ) return Standard_Complex_Vectors.Vector is

    map : constant Link_to_Monomial_Map := Retrieve_Map(d,k);
    nullvec : Standard_Complex_Vectors.Vector(0..0);

  begin
    if map = null then
      nullvec(0) := Create(-1.0);
      return nullvec;
    else
      return map.c;
    end if;
  end Coefficients;

  function Exponents
             ( d,k : integer32 ) return Standard_Integer_VecVecs.VecVec is

    map : constant Link_to_Monomial_Map := Retrieve_Map(d,k);
    nullvec : Standard_Integer_VecVecs.VecVec(0..0);

  begin
    if map = null then
      nullvec(0) := null;
      return nullvec;
    else
      return map.v;
    end if;
  end Exponents;

  procedure Coefficients_and_Exponents
              ( d,k : in integer32;
                cff : out Standard_Complex_Vectors.Vector;
                exp : out Standard_Integer_VecVecs.VecVec ) is

    map : constant Link_to_Monomial_Map := Retrieve_Map(d,k);

  begin
    if map /= null then
      cff := map.c;
      exp := map.v;
    end if;
  end Coefficients_and_Exponents;

  procedure Clear is
  begin
    Clear(the_maps);
  end Clear;

begin
  the_maps := null;
end Monomial_Maps_Container;
