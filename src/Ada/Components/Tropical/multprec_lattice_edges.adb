with unchecked_deallocation;
with Standard_Lattice_Supports;

package body Multprec_Lattice_Edges is

-- CONSTRUCTORS :

  function Edges_of_3D_Hull ( A : Matrix ) return Edge_List is

    f : constant Facet_3d_List := Convex_Hull_3D(A);
    m : constant integer32 := A'last(2);

  begin
    return Edges_of_3D_Hull(m,f);
  end Edges_of_3D_Hull;

  function Edges_of_3D_Hull
              ( m : integer32; f : Facet_3d_List ) return Edge_List is

    res,res_last : Edge_List;
    cnt : integer32 := 0;        -- defines the labels of the edges

    procedure Append_New_Edge
                  ( p,q : in integer32; pf,qf : in Link_to_3d_Facet ) is

    -- DESCRIPTION :
    --   Appends the edge defined by vertices p and q
    --   in the intersection of the facets pf and qf to the list res.

      e : Edge;
      lk : Link_to_Edge;

    begin
      e.label := cnt;
      e.a := p; e.b := q;
      e.f := pf; e.g := qf;
      lk := new Edge'(e);
      Append(res,res_last,lk);
      cnt := cnt + 1;
    end Append_New_Edge;

    procedure Edge_Crawl ( g : in Facet_3d_List; v : in integer32;
                           b : in out Boolean_Matrix ) is

    -- DESCRIPTION :
    --   The adjusted Multprec_Lattice_Facets.Crawl appends edges
    --   to the resulting list res.

      tmp : Facet_3d_List := g;
      lft,z : Link_to_3d_Facet;
      ind,w : integer32;

    begin
      while not Is_Null(tmp) loop
        lft := Head_Of(tmp);
        ind := Standard_Lattice_Supports.Member(lft.points,v);
        if ind >= lft.points'first then   -- v is incident to lft
          if ind > lft.points'first then  -- w will be previous neighbor
            w := lft.points(ind-1);
            z := lft.neighbors(ind-1);
          else
            w := lft.points(lft.points'last);
            z := lft.neighbors(lft.points'last);
          end if;
          if not b(v,w) then
            Append_New_Edge(v,w,lft,z); 
            b(v,w) := true; b(w,v) := true; Edge_Crawl(g,w,b);
          end if;
          z := lft.neighbors(ind);
          if ind < lft.points'last       -- w will be next neighbor
           then w := lft.points(ind+1);
           else w := lft.points(lft.points'first);
          end if;
          if not b(v,w) then
            Append_New_Edge(v,w,lft,z);
            b(v,w) := true; b(w,v) := true; Edge_Crawl(g,w,b);
          end if;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end Edge_Crawl;

    procedure Main is

      b : Boolean_Matrix(1..m,1..m);
      v : constant integer32 := First_Incident_Vertex(f);

    begin
      for i in b'range(1) loop
        for j in b'range(2) loop
          b(i,j) := false;
        end loop;
      end loop;
      Edge_Crawl(f,v,b);
    end Main;

  begin
    Main;
    return res;
  end Edges_of_3D_Hull;

-- DESTRUCTORS :

  procedure Clear ( e : in out Link_to_Edge ) is

    procedure free is new unchecked_deallocation(Edge,Link_to_Edge);

  begin
    if e /= null
     then free(e);
    end if;
  end Clear;

  procedure Clear ( e : in out Edge_List ) is

    tmp : Edge_List := e;

  begin
    while not Is_Null(tmp) loop
      declare
        lk : Link_to_Edge := Head_Of(tmp);
      begin
        Clear(lk);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_Edges.Clear(Lists_of_Edges.List(e));
  end Clear;

end Multprec_Lattice_Edges;
