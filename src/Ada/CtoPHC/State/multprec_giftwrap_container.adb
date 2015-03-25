with String_Splitters;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Multprec_Integer_Vectors;
with Multprec_Lattice_Polytopes;

package body Multprec_Giftwrap_Container is

-- INTERNAL DATA :

  f3d : Multprec_Lattice_3d_Facets.Facet_3d_List;
  f4d : Multprec_Lattice_4d_Facets.Facet_4d_List;
  pts3 : Multprec_Integer_Matrices.Link_to_Matrix;
  pts4 : Multprec_Integer_Matrices.Link_to_Matrix;
  supp : String_Splitters.Link_to_String := null;

-- CONSTRUCTOR :

  procedure Create ( A : in Matrix ) is
  begin
    if A'last(1) = 3 then
      f3d := Multprec_Lattice_3d_Facets.Convex_Hull_3D(A);
      pts3 := new Multprec_Integer_Matrices.Matrix'(A);
    end if;
    if A'last(1) = 4 then
      declare
        r : natural32;
        v : Multprec_Integer_Vectors.Vector(A'range(1));
      begin
        Multprec_Lattice_Polytopes.Initial_Facet_Normal(A,r,v);
        f4d := Multprec_Lattice_4d_Facets.Convex_Hull_4D(A,v);
      end;
      pts4 := new Multprec_Integer_Matrices.Matrix'(A);
    end if;
  end Create;

  procedure Store_String ( s : in string ) is
  begin
    supp := new string'(s);
  end Store_String;

-- SELECTORS :

  function Point_Configuration_in_3d return Link_to_Matrix is
  begin
    return pts3;
  end Point_Configuration_in_3d;

  function Point_Configuration_in_4d return Link_to_Matrix is
  begin
    return pts4;
  end Point_Configuration_in_4d;

  function Number_of_3d_Facets return natural32 is
  begin
    return Multprec_Lattice_3d_Facets.Length_Of(f3d);
  end Number_of_3d_Facets;

  function Number_of_4d_Facets return natural32 is
  begin
    return Multprec_Lattice_4d_Facets.Length_Of(f4d);
  end Number_of_4d_Facets;

  function Facet_3d_Data
             ( k : natural32 )
             return Multprec_Lattice_3d_Facets.Link_to_3d_Facet is

    use Multprec_Lattice_3d_Facets;

    res : Link_to_3d_Facet;
    tmp : Facet_3d_List := f3d;
    lft : Link_to_3d_Facet;

  begin
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      if lft.label = integer32(k) 
       then res := lft; exit;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return res;
  end Facet_3d_Data;

  function Facet_4d_Data
             ( k : natural32 )
             return Multprec_Lattice_4d_Facets.Link_to_4d_Facet is

    use Multprec_Lattice_4d_Facets;

    res : Link_to_4d_Facet;
    tmp : Facet_4d_List := f4d;
    lft : Link_to_4d_Facet;

  begin
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      if lft.label = integer32(k)
       then res := lft; exit;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return res;
  end Facet_4d_Data;

  function Retrieve_String return string is

    use String_Splitters;

  begin
    if supp = null
     then return "";
     else return supp.all;
    end if;
  end Retrieve_String;

-- DESTRUCTOR : 

  procedure Clear_3d is
  begin
    Multprec_Lattice_3d_Facets.Clear(f3d);
    pts3 := null; -- memory for A managed outside container
  end Clear_3d;

  procedure Clear_4d is
  begin
    Multprec_Lattice_4d_Facets.Clear(f4d);
    pts4 := null; -- memory for A managed outside container
  end Clear_4d;

  procedure Clear_String is
  begin
    String_Splitters.Clear(supp);
    supp := null;
  end Clear_String;

  procedure Clear is
  begin
    Clear_3d;
    Clear_4d;
    Clear_String;
  end Clear;

end Multprec_Giftwrap_Container;
