with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Multprec_Integer_Vectors;
with Multprec_Lattice_Polytopes;

package body Multprec_Giftwrap_Container is

-- INTERNAL DATA :

  f3d : Multprec_Lattice_3d_Facets.Facet_3d_List;
  f4d : Multprec_Lattice_4d_Facets.Facet_4d_List;

-- CONSTRUCTOR :

  procedure Create ( A : in Matrix ) is
  begin
    if A'last(1) = 3 
     then f3d := Multprec_Lattice_3d_Facets.Convex_Hull_3D(A);
    end if;
    if A'last(1) = 4 then
      declare
        r : natural32;
        v : Multprec_Integer_Vectors.Vector(A'range(1));
      begin
        Multprec_Lattice_Polytopes.Initial_Facet_Normal(A,r,v);
        f4d := Multprec_Lattice_4d_Facets.Convex_Hull_4D(A,v);
      end;
    end if;
  end Create;

-- SELECTORS :

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

-- DESTRUCTOR : 

  procedure Clear_3d is
  begin
    Multprec_Lattice_3d_Facets.Clear(f3d);
  end Clear_3d;

  procedure Clear_4d is
  begin
    Multprec_Lattice_4d_Facets.Clear(f4d);
  end Clear_4d;

  procedure Clear is
  begin
    Clear_3d;
    Clear_4d;
  end Clear;

end Multprec_Giftwrap_Container;
