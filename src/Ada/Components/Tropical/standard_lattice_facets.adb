with text_io;                          use text_io;
--with Standard_Integer64_Vectors_io; use Standard_Integer64_Vectors_io;
--with Standard_Integer64_Matrices_io; use Standard_Integer64_Matrices_io;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with unchecked_deallocation;
with Standard_Lattice_Supports;
with Standard_Lattice_Polytopes;
with Standard_Power_Transformations;
with Standard_Lattice_4d_Facets;

package body Standard_Lattice_Facets is

-- CONSTRUCTORS :

  function Create_4d_Facet
             ( f : Standard_Lattice_4d_Facets.Facet_in_4d ) return Facet is

  -- DESCRIPTION :
  --   Returns the equivalent facet representation of the 4d given facet.

  -- REQUIRED : f.d = 4.

    res : Facet(f.d,f.n,f.m);

  begin
    res.label := f.label;
    res.normal := f.normal;
    res.points := f.points;
    res.ridges := f.ridges;
    for i in res.neighbors'range loop
      declare
        use Standard_Lattice_4d_Facets;
        nf : constant Link_to_4d_Facet := f.neighbors(i);
      begin
        if nf = null
         then res.neighbors(i) := null;
         else res.neighbors(i) := new Facet'(Create_4d_Facet(nf.all));
        end if;
      end;
    end loop;
    return res;
  end Create_4d_Facet;

  procedure Convex_Hull_of_Ridge
              ( A : in Matrix; v : in Standard_Integer64_Vectors.Vector;
                s : in Standard_Integer_Vectors.Vector;
                p : out integer32; M : out Matrix; F : out Facet_List ) is

    B : Matrix(A'range(1),s'range);
    MB : Matrix(A'range(1),s'range);
    B2 : Matrix(A'first(1)..A'last(1)-1,s'range);

  begin
    B := Standard_Lattice_Supports.Support_Points(A,s);
    p := Standard_Power_Transformations.Pivot(v);
   -- put_line("calling eliminate ...");
    M := Standard_Power_Transformations.Eliminate(v,p);
   -- put_line("performing matrix multiplication ...");
    MB := M*B;
   -- put_line("dropping coordinate ...");
    B2 := Standard_Lattice_3d_Facets.Drop(MB,p);
    F := Convex_Hull(B2);
 -- exception
 --   when others => put("exception in Convex_Hull_of_Ridge"); 
 --     put("  v = "); put(v); put(", p = "); put(p,1); new_line;
 --     put_line("B2 = "); put(B2); raise;
  end Convex_Hull_of_Ridge;

  function Filter_non_Vertices
              ( f : Facet; v : Standard_Integer_Vectors.Vector )
              return Facet is

    res : Facet(f.d,v'last,f.m);
    ind : integer32 := res.points'first-1;

  begin
    res.label := f.label;
    res.normal := f.normal;
    for i in f.points'range loop
      if Standard_Lattice_Supports.Member(v,f.points(i)) >= v'first
       then ind := ind + 1; res.points(ind) := f.points(i);
      end if;
    end loop;
    res.faces := f.faces;
    res.neighbors := f.neighbors;
    return res;
  end Filter_non_Vertices;

  function Faces_of_Facet
              ( A : Matrix; v : Standard_Integer64_Vectors.Vector )
              return Facet is

    s : constant Standard_Integer_Vectors.Vector
      := Standard_Lattice_Supports.Support(A,v);
    p : integer32;
    M : Matrix(A'range(1),A'range(1));
    Cf : Facet_List;
    d : constant integer32 := v'last;
    n : constant integer32 := s'last;

  begin
    Convex_Hull_of_Ridge(A,v,s,p,M,Cf);
    declare
      mc : constant integer32 := integer32(Length_Of(Cf));
      res : Facet(d,n,mc);
      ptr : Facet_List;
    begin
      res.label := 0; res.normal := v; res.points := s;
      ptr := Cf;
      for i in res.faces'range loop
        declare
          lft : constant Link_to_Facet := Head_Of(ptr);
          fct : Facet(lft.d+1,lft.n,lft.m);
        begin
          fct.label := lft.label;
          fct.normal := Standard_Lattice_4d_Facets.Insert_Zero(lft.normal,p);
          Standard_Lattice_4d_Facets.Multiply_with_Transpose(fct.normal,M);
          fct.points := Standard_Lattice_4d_Facets.Relabel(lft.points,s);
          fct.neighbors := lft.neighbors;
          res.faces(i) := new Facet'(fct);
        end;
        ptr := Tail_Of(ptr);
      end loop;
      res.neighbors := (1..mc => null);
      Clear(Cf);
      declare
        vtp : constant Standard_Integer_Vectors.Vector
            := Vertices(A'last(2),res.faces);
      begin
        if vtp'last /= res.points'last
         then return Filter_non_Vertices(res,vtp);
         else return res;
        end if;
      end;
    end;
  end Faces_of_Facet;

  function Ridges_of_Facet
              ( A : Matrix; v : Standard_Integer64_Vectors.Vector )
              return Facet is

    d : constant integer32 := A'last(1);

  begin
   -- put("in Ridges_of_Facet with d = "); put(d,1); new_line;
    if d = 4 then
      declare
        use Standard_Lattice_4d_Facets;
        f : Link_to_4d_Facet;
      begin
        f := new Facet_in_4d'(Standard_Lattice_4d_Facets.Ridges_of_Facet(A,v));
        return Create_4d_Facet(f.all);
     -- exception
     --   when others => put_line("exception in Ridges_of_Facet with A = ");
     --     put(A); raise;
      end;
    else
      return Faces_of_Facet(A,v);
    end if;
  end Ridges_of_Facet;

  procedure Connect_along_Ridge
               ( f : in Link_to_Facet; r : in integer32 ) is

  -- DESCRIPTION :
  --   Connects f.neighbor(r) to f along the corresponding ridge.

  -- REQUIRED : f.d = 4.

    rdg : constant Standard_Lattice_3d_Facets.Link_to_3d_Facet := f.ridges(r);
    g : constant Link_to_Facet := f.neighbors(r);

  begin
    for i in g.ridges'range loop
      declare
        lft : constant Standard_Lattice_3d_Facets.Link_to_3d_Facet
            := g.ridges(i);
      begin
        if Standard_Lattice_4d_Facets.Is_Same(lft.points,rdg.points)
         then g.neighbors(i) := f; exit;
        end if;
      end;
    end loop;
  end Connect_along_Ridge;

  procedure Connect_along_Face
               ( f : in Link_to_Facet; r : in integer32 ) is

  -- DESCRIPTION :
  --   Connects f.neighbor(r) to f along the corresponding face.

  -- REQUIRED : f.d > 4.

    rfg : constant Link_to_Facet := f.faces(r);
    g : constant Link_to_Facet := f.neighbors(r);

  begin
    for i in g.faces'range loop
      declare
        lft : constant Link_to_Facet := g.faces(i);
      begin
        if Standard_Lattice_4d_Facets.Is_Same(lft.points,rfg.points)
         then g.neighbors(i) := f; exit;
        end if;
      end;
    end loop;
  end Connect_along_Face;

  procedure Connect_along_Ridge
              ( f,g : in Link_to_Facet; connected : out boolean ) is

  -- DESCRIPTION :
  --   Connects f and g along a common ridge.

  -- REQUIRED : f.d = 4

  begin
    connected := false;
    for r in f.neighbors'range loop
      if f.neighbors(r) = null then
        for i in g.neighbors'range loop
          if g.neighbors(i) = null then
            if Standard_Lattice_4d_Facets.Is_Same
                 (f.ridges(r).points,g.ridges(i).points) then
              f.neighbors(r) := g;
              g.neighbors(i) := f;
              connected := true;
            end if;
          end if;
          exit when connected;
        end loop;
      end if;
      exit when connected;
    end loop;
  end Connect_along_Ridge;

  procedure Connect_along_Face
              ( f,g : in Link_to_Facet; connected : out boolean ) is

  -- DESCRIPTION :
  --   Connects f and g along a common face.

  -- REQUIRED : f.d > 4

  begin
    connected := false;
    for r in f.neighbors'range loop
      if f.neighbors(r) = null then
        for i in g.neighbors'range loop
          if g.neighbors(i) = null then
            if Standard_Lattice_4d_Facets.Is_Same
                 (f.faces(r).points,g.faces(i).points) then
              f.neighbors(r) := g;
              g.neighbors(i) := f;
              connected := true;
            end if;
          end if;
          exit when connected;
        end loop;
      end if;
      exit when connected;
    end loop;
  end Connect_along_Face;

  procedure Neighbors ( A : in Matrix; f : in Link_to_Facet;
                        idcnt : in out integer32 ) is

    k : integer32;
    v : Standard_Integer64_Vectors.Vector(A'range(1));
    connected : boolean;

    use Standard_Lattice_4d_Facets;
    
  begin
    if f.d = 4 then
      for r in f.ridges'range loop
        if f.neighbors(r) = null then
          k := Extreme_Angle(A,f.ridges(r).points,f.normal,f.ridges(r).normal);
          v := Inner_Normal(A,f.ridges(r).points,k);
          f.neighbors(r) := new Facet'(Ridges_of_Facet(A,v));
          idcnt := idcnt + 1;
          f.neighbors(r).label := idcnt;
          Connect_along_Ridge(f,r);
          for i in f.neighbors'range loop
            if i /= r and f.neighbors(i) /= null then
              Connect_along_Ridge(f.neighbors(r),f.neighbors(i),connected);
            end if;
          end loop;
        end if;
      end loop;
    else
      for r in f.faces'range loop
        if f.neighbors(r) = null then
          k := Extreme_Angle(A,f.faces(r).points,f.normal,f.faces(r).normal);
          v := Inner_Normal(A,f.faces(r).points,k);
          f.neighbors(r) := new Facet'(Ridges_of_Facet(A,v));
          idcnt := idcnt + 1;
          f.neighbors(r).label := idcnt;
          Connect_along_Face(f,r);
          for i in f.neighbors'range loop
            if i /= r and f.neighbors(i) /= null then
              Connect_along_Face(f.neighbors(r),f.neighbors(i),connected);
            end if;
          end loop;
        end if;
      end loop;
    end if;
  end Neighbors;

  procedure Connect ( f : in Facet_List; lft : in Link_to_Facet ) is

  -- DESCRIPTION :
  --   Connects the facet lft to the facets in the list f.

    tmp : Facet_List := f;
    ptr : Link_to_Facet;
    connected : boolean;

  begin
    if lft.d = 4 then
      while not Is_Null(tmp) loop
        ptr := Head_Of(tmp);
        if not Is_Connected(ptr)
         then Connect_along_Ridge(ptr,lft,connected);
        end if;
        exit when Is_Connected(lft); -- we stop if lft is connected
        tmp := Tail_Of(tmp);
      end loop;
    else
      while not Is_Null(tmp) loop
        ptr := Head_Of(tmp);
        if not Is_Connected(ptr)
         then Connect_along_Face(ptr,lft,connected);
        end if;
        exit when Is_Connected(lft); -- we stop if lft is connected
        tmp := Tail_Of(tmp);
      end loop;
    end if;
  end Connect;

  function Convex_Hull
              ( A : Matrix; v : Standard_Integer64_Vectors.Vector )
              return Facet_List is

    res : Facet_List;
    first : constant Link_to_Facet := new Facet'(Ridges_of_Facet(A,v));
    lft : Link_to_Facet;
    cnt,inc,idcnt : integer32 := 0;

  begin
    Construct(first,res); cnt := 1;
    loop
      lft := Pop(res);
      exit when (lft = null);
      Neighbors(A,lft,idcnt);
      inc := 0;
      for i in lft.neighbors'range loop
        if lft.neighbors(i).label >= cnt then
          if not Is_Connected(lft.neighbors(i))
           then Connect(res,lft.neighbors(i));
          end if;
          Construct(lft.neighbors(i),res);
          inc := inc + 1;
        end if;
      end loop;
      cnt := cnt + inc;
    end loop;
    return res;
  end Convex_Hull;

  function Convex_Hull ( A : Matrix ) return Facet_List is

    res : Facet_List;
    r : natural32;
    v : Standard_Integer64_Vectors.Vector(A'range(1));

  begin
    Standard_Lattice_Polytopes.Initial_Facet_Normal(A,r,v);
    res := Convex_Hull(A,v);
    return res;
 -- exception
 --   when others => put("exception raised in Convex_Hull ");
 --      put("  v = "); put(v); new_line; raise;
  end Convex_Hull;

-- SELECTORS :

  function Is_Connected ( f : Link_to_Facet ) return boolean is
  begin
    for i in f.neighbors'range loop
      if f.neighbors(i) = null
       then return false;
      end if;
    end loop;
    return true;
  end Is_Connected;

  function Pop ( f : Facet_List ) return Link_to_Facet is

    tmp : Facet_List := f;
    lft : Link_to_Facet;

  begin
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      if not Is_Connected(lft)
       then return lft;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return null;
  end Pop;

  function Vertices ( n : integer32; f : Facet_List )
                    return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..n) := (1..n => 0);
    ind : integer32 := 0;
    ptr : Facet_List := f;
    lft : Link_to_Facet;
    found : boolean;

  begin
    while not Is_Null(ptr) loop
      lft := Head_Of(ptr);
      for i in lft.points'range loop
        found := false;
        for j in 1..ind loop
          if res(j) = lft.points(i)
           then found := true;
          end if;
          exit when found;
        end loop;
        if not found then
          ind := ind + 1;
          res(ind) := lft.points(i);
        end if;
      end loop;
      ptr := Tail_Of(ptr);
    end loop;
    return res(1..ind);
  end Vertices;

  function Vertices ( n : integer32; f : Array_of_Facets )
                    return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..n) := (1..n => 0);
    ind : integer32 := 0;
    lft : Link_to_Facet;
    found : boolean;

  begin
    for i in f'range loop
      lft := f(i);
      for i in lft.points'range loop
        found := false;
        for j in 1..ind loop
          if res(j) = lft.points(i)
           then found := true;
          end if;
          exit when found;
        end loop;
        if not found then
          ind := ind + 1;
          res(ind) := lft.points(i);
        end if;
      end loop;
    end loop;
    return res(1..ind);
  end Vertices;

  procedure Concatenate 
              ( first,last : in out Lists_of_Integer_Vectors.List;
                e : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Appends all elements in e to the list first, without duplicates.

    p : Lists_of_Integer_Vectors.List := e;
    x : Standard_Integer_Vectors.Link_to_Vector;

  begin
    while not Lists_of_Integer_Vectors.Is_Null(p) loop
      x := Lists_of_Integer_Vectors.Head_Of(p);
      if not Lists_of_Integer_Vectors.Is_In(first,x.all)
       then Lists_of_Integer_Vectors.Append(first,last,x.all);
      end if;
      p := Lists_of_Integer_Vectors.Tail_Of(p);
    end loop;
  end Concatenate;

  function Edges ( f : Facet_List )
                  return Lists_of_Integer_Vectors.List is

    res,res_last,e : Lists_of_Integer_Vectors.List;
    ptr : Facet_List := f;
    lft : Link_to_Facet;

  begin
    while not Is_Null(ptr) loop
      lft := Head_Of(ptr);
      put("In Edges, d = "); put(lft.d,1); new_line;
      if lft.d = 4
       then e := Standard_Lattice_3D_Facets.Edges(lft.ridges);
       else e := Edges(lft.faces);
      end if;
      Concatenate(res,res_last,e);
      Lists_of_Integer_Vectors.Clear(e);
      ptr := Tail_Of(ptr);
    end loop;
    return res;
  end Edges;

  function Edges ( f : Array_of_Facets )
                 return Lists_of_Integer_Vectors.List is

    res,res_last,e : Lists_of_Integer_Vectors.List;
    lft : Link_to_Facet;

  begin
    for i in f'range loop
      put("in Edges with i = "); put(i,1); new_line;
      if f(i) /= null then
        lft := f(i);
        put("face "); put(i,1); put(" has dimension "); put(lft.d,1); new_line;
        if lft.d = 4
         then e := Standard_Lattice_3D_Facets.Edges(lft.ridges);
         else e := Edges(lft.faces);
        end if;
        Concatenate(res,res_last,e);
        Lists_of_Integer_Vectors.Clear(e);
      else
        put("face "); put(i,1); put_line(" is null");
      end if;
    end loop;
    return res;
  end Edges;

  function Ridges ( f : Facet_List )
                  return Lists_of_Integer_Vectors.List is

    res,res_last : Lists_of_Integer_Vectors.List;
    ptr : Facet_List := f;
    lft : Link_to_Facet;

  begin
    while not Is_Null(ptr) loop
      lft := Head_Of(ptr);
      if lft.d > 4 then
        for j in lft.faces'range loop
          declare
            r : constant Standard_Integer_Vectors.Vector
              := lft.faces(j).points;
          begin      
            if not Standard_Lattice_4D_Facets.Belongs_To(res,r)
             then Lists_of_Integer_Vectors.Append(res,res_last,r);
            end if;
          end;
        end loop;
      else
        for j in lft.ridges'range loop
          declare
            r : constant Standard_Integer_Vectors.Vector
              := lft.ridges(j).points;
          begin      
            if not Standard_Lattice_4D_Facets.Belongs_To(res,r)
             then Lists_of_Integer_Vectors.Append(res,res_last,r);
            end if;
          end;
        end loop;
      end if;
      ptr := Tail_Of(ptr);
    end loop;
    return res;
  end Ridges;

  function Ridges ( f : Array_of_Facets )
                  return Lists_of_Integer_Vectors.List is

    res,res_last : Lists_of_Integer_Vectors.List;
    lft : Link_to_Facet;

  begin
    for i in f'range loop
      lft := f(i);
      if lft.d > 4 then
        for j in lft.faces'range loop
          declare
            r : constant Standard_Integer_Vectors.Vector
              := lft.faces(j).points;
          begin      
            if not Standard_Lattice_4D_Facets.Belongs_To(res,r)
             then Lists_of_Integer_Vectors.Append(res,res_last,r);
            end if;
          end;
        end loop;
      else
        for j in lft.ridges'range loop
          declare
            r : constant Standard_Integer_Vectors.Vector
              := lft.ridges(j).points;
          begin      
            if not Standard_Lattice_4D_Facets.Belongs_To(res,r)
             then Lists_of_Integer_Vectors.Append(res,res_last,r);
            end if;
          end;
        end loop;
      end if;
    end loop;
    return res;
  end Ridges;

  function Faces ( f : Facet_List; d : integer32 )
                 return Lists_of_Integer_Vectors.List is

    res,res_last : Lists_of_Integer_Vectors.List;
    ptr : Facet_List := f;
    lft : Link_to_Facet;

  begin
    while not Is_Null(ptr) loop
      lft := Head_Of(ptr);
     -- if lft.normal'last = d+1 then
      if lft.d = d+1 then
        if not Standard_Lattice_4D_Facets.Belongs_To(res,lft.points)
         then Lists_of_Integer_Vectors.Append(res,res_last,lft.points);
        end if;
      elsif lft.d > 4 then
        declare
          rdf : constant Lists_of_Integer_Vectors.List := Faces(lft.faces,d);
        begin
          Lists_of_Integer_Vectors.Concat(res,res_last,rdf);
        end;
      end if;
      ptr := Tail_Of(ptr);
    end loop;
    return res;
  end Faces;

  function Faces ( f : Array_of_Facets; d : integer32 )
                 return Lists_of_Integer_Vectors.List is

    res,res_last : Lists_of_Integer_Vectors.List;
    lft : Link_to_Facet;

  begin
    for i in f'range loop
      lft := f(i);
     -- put("the normal : "); put(lft.normal); new_line;
     -- if lft.normal'last = d+1 then
      if lft.d = d+1 then
        if not Standard_Lattice_4D_Facets.Belongs_To(res,lft.points)
         then Lists_of_Integer_Vectors.Append(res,res_last,lft.points);
        end if;
      elsif lft.d > 4 then
       -- put("calling for lft.d = "); put(lft.d,1); new_line;
        declare
          rdf : constant Lists_of_Integer_Vectors.List := Faces(lft.faces,d);
        begin
          Lists_of_Integer_Vectors.Concat(res,res_last,rdf);
        end;
      end if;
    end loop;
    return res;
  end Faces;

  function fvector ( n,d : integer32; f : Facet_List )
                   return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(0..d-1) := (0..d-1 => 0);

  begin
    res(d-1) := Length_Of(f);
    declare
      v : constant Standard_Integer_Vectors.Vector := Vertices(n,f);
      r : Lists_of_Integer_Vectors.List := Ridges(f);
      e : Lists_of_Integer_Vectors.List := Edges(f);
    begin
      res(0) := natural32(v'last);
      res(1) := Lists_of_Integer_Vectors.Length_Of(e);
      Lists_of_Integer_Vectors.Clear(e);
      res(d-2) := Lists_of_Integer_Vectors.Length_Of(r);
      Lists_of_Integer_Vectors.Clear(r);
    end;
   -- for i in reverse d-2..d-2 loop
   --   declare
   --     lfs : Lists_of_Integer_Vectors.List := Faces(f,i);
   --   begin
   --     res(i) := Lists_of_Integer_Vectors.Length_Of(lfs);
   --     Lists_of_Integer_Vectors.Clear(lfs);
   --   end;
   -- end loop;
    return res;
  end fvector;

-- DESTRUCTORS :

  procedure Clear ( f : in out Facet ) is
  begin
    if f.d = 4 then
      for i in f.ridges'range loop
        Standard_Lattice_3d_Facets.Clear(f.ridges(i));
      end loop;
    else
      for i in f.faces'range loop
        Clear(f.faces(i));
      end loop;
    end if;
  end Clear;

  procedure Clear ( f : in out Link_to_Facet ) is

    procedure free is 
      new unchecked_deallocation(Facet,Link_to_Facet);

  begin
    if f /= null
     then free(f);
    end if;
  end Clear;

  procedure Clear ( f : in out Array_of_Facets ) is
  begin
    for i in f'range loop
      Clear(f(i));
    end loop;
  end Clear;

end Standard_Lattice_Facets;
