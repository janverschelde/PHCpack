--with text_io; use text_io;
--with Standard_Integer_Numbers_io; use Standard_Integer_Numbers_io;
--with Standard_Integer_Vectors_io; use Standard_Integer_Vectors_io;
--with Standard_Integer64_Vectors_io; use Standard_Integer64_Vectors_io;
--with Standard_Integer64_Matrices_io; use Standard_Integer64_Matrices_io;
--with Standard_Lattice_4D_Facets_io; use Standard_Lattice_4D_Facets_io;

with unchecked_deallocation;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Lattice_Supports;
with Standard_Integer_Orthogonals;
with Standard_Power_Transformations;

package body Standard_Lattice_4d_Facets is

-- CONSTRUCTORS :

  function Insert_Zero
              ( v : Standard_Integer64_Vectors.Vector; p : integer32 )
              return Standard_Integer64_Vectors.Vector is

    res : Standard_Integer64_Vectors.Vector(v'first..v'last+1);

  begin
    for i in v'first..(p-1) loop
      res(i) := v(i);
    end loop;
    res(p) := 0;
    for i in p..v'last loop
      res(i+1) := v(i);
    end loop;
    return res;
  end Insert_Zero;

  procedure Multiply_with_Transpose
              ( v : in out Standard_Integer64_Vectors.Vector;
                M : in Matrix ) is

    w : Standard_Integer64_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      w(k) := 0;
      for i in M'range(1) loop
        w(k) := w(k) + M(i,k)*v(i);
      end loop;
    end loop;
    v := w;
  end Multiply_with_Transpose;

  function Relabel ( p,s : Standard_Integer_Vectors.Vector )
                   return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(p'range);

  begin
    for i in p'range loop
      res(i) := s(p(i));
    end loop;
    return res;
  end Relabel;

  procedure Convex_Hull_of_Ridge
              ( A : in Matrix; v : in Standard_Integer64_Vectors.Vector;
                s : in Standard_Integer_Vectors.Vector;
                p : out integer32; M : out Matrix;
                F : out Standard_Lattice_3d_Facets.Facet_3d_List ) is

    B : Matrix(A'range(1),s'range);
    MB : Matrix(A'range(1),s'range);
    B2 : Matrix(1..3,s'range);

  begin
   -- put_line("In Convex_Hull_of_Ridge ...");
   -- put_line("The matrix A :"); put(A);
   -- put("with v = " ); put(v); new_line;
    B := Standard_Lattice_Supports.Support_Points(A,s);
   -- put_line(" ridge is spanned by the points :"); put(B);
    p := Standard_Power_Transformations.Pivot(v);
   -- put("the pivot p = "); put(p,1); new_line;
    M := Standard_Power_Transformations.Eliminate(v,p);
   -- put("The matrix M : "); put(M);
    MB := M*B;
   -- put("The matrix MB : "); put(MB);
    B2 := Standard_Lattice_3d_Facets.Drop(MB,p);
   -- put_line("B2 after dropping one coordinate :"); put(B2);
    F := Standard_Lattice_3d_Facets.Convex_Hull_3D(B2);
 -- exception
 --   when others => put_line("exception in Convex_Hull_of_Ridge B2 =");
 --     put(B2);  
 --     put_line("A = "); put(A);
 --     put("v = "); put(v); new_line;
 --     raise;
  end Convex_Hull_of_Ridge;

  function Filter_non_Vertices
              ( f : Facet_in_4d; v : Standard_Integer_Vectors.Vector )
              return Facet_in_4d is

    res : Facet_in_4d(f.d,v'last,f.m);
    ind : integer32 := res.points'first-1;

  begin
    res.label := f.label;
    res.normal := f.normal;
    for i in f.points'range loop
      if Standard_Lattice_Supports.Member(v,f.points(i)) >= v'first
       then ind := ind + 1; res.points(ind) := f.points(i);
      end if;
    end loop;
    res.ridges := f.ridges;
    res.neighbors := f.neighbors;
   -- put("res.points after adjustments : "); put(res.points); new_line;
    return res;
  end Filter_non_Vertices;

  function Ridges_of_Facet
              ( A : Matrix; v : Standard_Integer64_Vectors.Vector )
              return Facet_in_4d is

    s : constant Standard_Integer_Vectors.Vector
      := Standard_Lattice_Supports.Support(A,v);
    p : integer32;
    M : Matrix(A'range(1),A'range(1));
    Cf : Standard_Lattice_3d_Facets.Facet_3d_List;
    d : constant integer32 := v'last;
    n : constant integer32 := s'last;

  begin
    Convex_Hull_of_Ridge(A,v,s,p,M,Cf);
    declare
      mc : constant integer32
         := integer32(Standard_Lattice_3d_Facets.Length_Of(Cf));
      res : Facet_in_4d(d,n,mc);
      ptr : Standard_Lattice_3d_Facets.Facet_3d_List;
    begin
     -- put("in Ridges_of_Facet with mc = "); put(mc,1); new_line;
      res.label := 0; res.normal := v; res.points := s;
      ptr := Cf;
      for i in res.ridges'range loop
        declare
          lft : constant Standard_Lattice_3d_Facets.Link_to_3d_Facet
              := Standard_Lattice_3d_Facets.Head_Of(ptr);
          fct : Standard_Lattice_3d_Facets.Facet_in_3d(lft.d+1,lft.n);
        begin
          fct.label := lft.label;
         -- put("normal to 3d facet : "); put(lft.normal);
         -- put(" inserting 0 at p = "); put(p,1); 
          fct.normal := Insert_Zero(lft.normal,p);
         -- put(" : "); put(fct.normal); new_line;
          Multiply_with_Transpose(fct.normal,M);
         -- put("after multiplying with M^T : "); put(fct.normal); new_line;
          fct.points := Relabel(lft.points,s);
          fct.neighbors := lft.neighbors;
          res.ridges(i) := new Standard_Lattice_3d_Facets.Facet_in_3d'(fct);
        end;
        ptr := Standard_Lattice_3d_Facets.Tail_Of(ptr);
      end loop;
     -- put_line("outside the loop in Ridges_of_Facet");
      res.neighbors := (1..mc => null);
     -- put_line("clearing the facets list Cf");
      Standard_Lattice_3d_Facets.Clear(Cf);
     -- put_line("The new facet : ");
     -- Write_4D_Facet(A,res);
     -- put_line("checking the vertices");
      declare
        vtp : constant Standard_Integer_Vectors.Vector
            := Standard_Lattice_3d_Facets.Vertices(A'last(2),res.ridges);
      begin
       -- put("vertex points : "); put(vtp);
       -- put("  points in f : "); put(res.points); new_line;
        if vtp'last /= res.points'last
         then return Filter_non_Vertices(res,vtp);
         else return res;
        end if;
      end;
    end;
 -- exception
 --   when others => put_line("problem in Ridges_of_Facet"); 
 --                  put("v = "); put(v); raise;
  end Ridges_of_Facet;

  function Extreme_Angle
              ( A : Matrix; f : Standard_Integer_Vectors.Vector;
                v,w : Standard_Integer64_Vectors.Vector ) return integer32 is

    res : integer32;
    ind : integer32 := A'first(2) - 1;
    x1,y1,x2,y2 : integer64;
    d : Standard_Integer64_Vectors.Vector(A'range(1));
 
  begin
   -- put("inside Extreme_Angle ...");
   -- put(" f = "); put(f); put(" v = "); put(v);
   -- put(" w = "); put(w); new_line;
    y1 := 0;
    for k in A'range(2) loop
      if Standard_Lattice_Supports.Member(f,k) < f'first then
       -- put("calling Shift with k = "); put(k,1); new_line;
        d := Standard_Lattice_3d_Facets.Shift(A,f(f'first),k);
        y1 := Standard_Lattice_Supports.Inner_Product(d,v);
        x1 := Standard_Lattice_Supports.Inner_Product(d,w);
        if y1 < 0 then y1 := -y1; end if;
        ind := k; res := k;
      end if;
      exit when ((y1 /= 0) and (ind >= A'first(2)));
    end loop;
   -- put("  x1 = "); put(x1,1); put("  y1 = "); put(y1,1);
   -- put("  ind = "); put(ind,1); new_line;
    for k in (ind+1)..A'last(2) loop
      if Standard_Lattice_Supports.Member(f,k) < f'first then
       -- put("calling Shift with k = "); put(k,1); new_line;
        d := Standard_Lattice_3d_Facets.Shift(A,f(f'first),k);
        y2 := Standard_Lattice_Supports.Inner_Product(d,v);
        x2 := Standard_Lattice_Supports.Inner_Product(d,w);
        --put("  x2 = "); put(x2,1); put("  y2 = "); put(y2,1); new_line;
       -- also added flipping sign of x2 ???
        if y2 < 0 then y2 := -y2; end if;
       -- code copied from Standard_Lattice_3d_Facets.Extreme_Angle
        if y2 = 0 then
          null; -- put_line("  y2 = 0, should not happend ?");
        else
          if x2 = 0 then
            if x1 > 0
             then x1 := x2; y1 := y2; res := k; --put("  res = "); put(res,1);
            end if;
            --new_line;
          elsif x2 < 0 and x1 >= 0 then
            x1 := x2; y1 := y2; res := k; 
            --put(" res = "); put(res,1); new_line;
          elsif x1 > 0 and x2 > 0 then
            if x1*y2 > x2*y1
             then x1 := x2; y1 := y2; res := k; --put("  res = "); put(res,1);
            end if;
            --new_line;
         -- else case x1 < 0 and x2 > 0: no new extreme
          elsif x1 < 0 and x2 < 0 then
            if x1*y2 > x2*y1
             then x1 := x2; y1 := y2; res := k; --put("  res = "); put(res,1);
            end if;
            --new_line;
         -- else
         --   put_line(" no case for new extreme");
          end if;
        end if;
      end if;
    end loop;
   -- put("Extreme_Angle returns res = "); put(res,1); new_line;
    return res;
 -- exception
 --   when others => put_line("exception happens in Extreme_Angle"); raise;
  end Extreme_Angle;

  function Inner_Normal
              ( A : Matrix; p : Standard_Integer_Vectors.Vector;
                k : integer32 ) return Standard_Integer64_Vectors.Vector is

    res : Standard_Integer64_Vectors.Vector(A'range(1));
    B : Matrix(A'range(1),p'range);
    s,t : integer64;

  begin
   -- put_line("inside Inner_Normal :");
    for j in p'range loop
      for i in A'range(1) loop
        B(i,j) := A(i,p(j)) - A(i,k);
      end loop;
    end loop;
   -- put_line("calling Complement to B = "); put(B);
    res := Standard_Integer_Orthogonals.Complement(B);
   -- put("returning res = "); put(res); new_line;
    s := Standard_Lattice_Supports.Inner_Product(A,k,res);
   -- put("computed normal "); put(res);
   -- put(" with IP s = "); put(s,1);
    t := s;
    for i in A'range(2) loop
      if i /= k and Standard_Lattice_Supports.Member(p,i) < p'first then
        t := Standard_Lattice_Supports.Inner_Product(A,i,res);
      end if;
      exit when (t /= s);
    end loop;
   -- put(" t = "); put(t,1); new_line;
    if s > t
     then Standard_Integer64_Vectors.Min(res);
    end if;
    return res;
  end Inner_Normal;

  function Is_Same ( a,b : Standard_Integer_Vectors.Vector ) return boolean is
  begin
    if a'last /= b'last then
      return false;
    else
      for i in a'range loop
        if Standard_Lattice_Supports.Member(b,a(i)) < b'first
         then return false;
        end if;
      end loop;
      return true;
    end if;
  end Is_Same;

  procedure Connect_along_Ridge
               ( f : in Link_to_4d_Facet; r : in integer32 ) is

  -- DESCRIPTION :
  --   Connects f.neighbor(r) to f along the corresponding ridge.

    rdg : constant Standard_Lattice_3d_Facets.Link_to_3d_Facet := f.ridges(r);
    g : constant Link_to_4d_Facet := f.neighbors(r);

  begin
   -- put("connecting f.neighbor("); put(r,1);
   -- put_line(") to f along ridge");
    for i in g.ridges'range loop
      declare
        lft : constant Standard_Lattice_3d_Facets.Link_to_3d_Facet
            := g.ridges(i);
      begin
        if Is_Same(lft.points,rdg.points)
         then g.neighbors(i) := f; exit;
        end if;
      end;
    end loop;
   -- put_line("-> leaving connect_along_ridge");
  end Connect_along_Ridge;

  procedure Connect_along_Ridge
              ( f,g : in Link_to_4d_Facet; connected : out boolean ) is

  -- DESCRIPTION :
  --   Connects f and g along a common ridge.

  begin
    connected := false;
    for r in f.neighbors'range loop
      if f.neighbors(r) = null then
        for i in g.neighbors'range loop
          if g.neighbors(i) = null then
            if Is_Same(f.ridges(r).points,g.ridges(i).points) then
              f.neighbors(r) := g;
              g.neighbors(i) := f;
             -- put("connecting ridge "); put(r,1);
             -- put(" to ridge "); put(i,1); new_line;
              connected := true;
            end if;
          end if;
          exit when connected;
        end loop;
      end if;
      exit when connected;
    end loop;
   -- put_line("-> leaving second connect_along_ridge");
  end Connect_along_Ridge;

  procedure Neighbors ( A : in Matrix; f : in Link_to_4d_Facet;
                        idcnt : in out integer32 ) is

    k : integer32;
    v : Standard_Integer64_Vectors.Vector(A'range(1));
    connected : boolean;

  begin
    for r in f.ridges'range loop
      if f.neighbors(r) = null then
        k := Extreme_Angle(A,f.ridges(r).points,f.normal,f.ridges(r).normal);
       -- k := Extreme_Angle(A,f.points,f.normal,f.ridges(r).normal);
       -- put("extreme angle at column "); put(k,1); new_line;
        v := Inner_Normal(A,f.ridges(r).points,k);
       -- declare -- sanity check
       --   s : constant Standard_Integer_Vectors.Vector
       --     := Standard_Lattice_Supports.Support(A,v);
       -- begin
       --   put("the inner normal v : "); put(v); put(" supports ");
       --   put(s); new_line;
       --   put(" IP :");
       --   put(Standard_Lattice_Supports.Inner_Products(A,v));
       --   if s'length < A'last(1)
       --    then put_line(" BUG !!!");
       --    else put_line(" ok");
       --   end if;
       -- end;
       -- put("r = "); put(r,1); new_line;
       -- put_line("Calling Ridges_of_Facet ...");
        f.neighbors(r) := new Facet_in_4d'(Ridges_of_Facet(A,v));
       -- put_line(" done with Ridges_of_Facet");
        idcnt := idcnt + 1;
        f.neighbors(r).label := idcnt;
        Connect_along_Ridge(f,r);
       -- for i in f.neighbors'first..(r-1) loop
       --   if f.neighbors(i) /= null then
        for i in f.neighbors'range loop
          if i /= r and f.neighbors(i) /= null then
           -- put("making connections for new neighbor "); put(i,1);
           -- put_line(" ...");
            Connect_along_Ridge(f.neighbors(r),f.neighbors(i),connected);
          end if;
        end loop;
      end if;
    end loop;
   -- put_line("leaving Neighbors...");
 -- exception
 --   when others => put_line("exception happens in Neighbors"); raise;
  end Neighbors;

  procedure Connect ( f : in Facet_4d_List; lft : in Link_to_4d_Facet ) is

  -- DESCRIPTION :
  --   Connects the facet lft to the facets in the list f.

    tmp : Facet_4d_List := f;
    ptr : Link_to_4d_Facet;
    connected : boolean;

  begin
    while not Is_Null(tmp) loop
      ptr := Head_Of(tmp);
      if not Is_Connected(ptr)
       then Connect_along_Ridge(ptr,lft,connected);
      end if;
      exit when Is_Connected(lft); -- we stop if lft is connected
      tmp := Tail_Of(tmp);
    end loop;
  end Connect;

  function Convex_Hull_4D
              ( A : Matrix; v : Standard_Integer64_Vectors.Vector )
              return Facet_4d_List is

    res : Facet_4d_List;
    first : constant Link_to_4d_Facet
          := new Facet_in_4d'(Ridges_of_Facet(A,v));
    lft : Link_to_4d_Facet;
    cnt,inc,idcnt : integer32 := 0;

  begin
    Construct(first,res); cnt := 1;
    loop
      lft := Pop(res);
      exit when (lft = null);
     -- put("Popped facet with label : "); put(lft.label,1); new_line;
     -- put("constructing the neigbhors of facet cnt = "); put(cnt,1); new_line;
     -- put("-> before neighbors, idcnt = "); put(idcnt,1); new_line;
      Neighbors(A,lft,idcnt);
     -- put("-> after neighbors, idcnt = "); put(idcnt,1); new_line;
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
     -- put("  inc = "); put(inc,1);
     -- put("  cnt = "); put(cnt,1); new_line;
    end loop;
    return res;
  end Convex_Hull_4D;

-- SELECTORS :

  function Is_Connected ( f : Link_to_4d_Facet ) return boolean is
  begin
    for i in f.neighbors'range loop
      if f.neighbors(i) = null
       then return false;
      end if;
    end loop;
    return true;
  end Is_Connected;

  function Pop ( f : Facet_4d_List ) return Link_to_4d_Facet is

    tmp : Facet_4d_List := f;
    lft : Link_to_4d_Facet;

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

-- for plain check of Euler characteristic

  function Vertices ( n : integer32; f : Facet_4d_List )
                    return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..n) := (1..n => 0);
    ind : integer32 := 0;
    tmp : Facet_4d_List := f;
    lft : Link_to_4d_Facet;
    found : boolean;

  begin
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
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
      tmp := Tail_Of(tmp);
    end loop;
    return res(1..ind);
  end Vertices;

  function Belongs_To ( R : Lists_of_Integer_Vectors.List;
                        p : Standard_Integer_Vectors.Vector ) return boolean is

    t : Lists_of_Integer_Vectors.List := R;
    v : Standard_Integer_Vectors.Link_to_Vector;

  begin
    while not Lists_of_Integer_Vectors.Is_Null(t) loop
      v := Lists_of_Integer_Vectors.Head_Of(t);
      if Is_Same(v.all,p)
       then return true;
       else t := Lists_of_Integer_Vectors.Tail_Of(t);
      end if;
    end loop;
    return false;
  end Belongs_To;

  function Edges ( f : Facet_4d_List ) return Lists_of_Integer_Vectors.List is

    res,res_last : Lists_of_Integer_Vectors.List;
    ptr : Facet_4d_List := f;
    lft : Link_to_4d_Facet;

  begin
    while not Is_Null(ptr) loop
      lft := Head_Of(ptr);
      for i in lft.ridges'range loop
        declare
          r : constant Standard_Integer_Vectors.Vector := lft.ridges(i).points;
          e : Standard_Integer_Vectors.Vector(1..2);
        begin
          for i in r'range loop
            e(1) := r(i);
            if i < r'last
             then e(2) := r(i+1);
             else e(2) := r(r'first);
            end if;
            if not Belongs_To(res,e)
             then Lists_of_Integer_Vectors.Append(res,res_last,e);
            end if;
          end loop;
        end;
      end loop;
      ptr := Tail_Of(ptr);
    end loop;
    return res;
  end Edges;

  function Ridges ( f : Facet_4d_List ) return Lists_of_Integer_Vectors.List is

    res,res_last : Lists_of_Integer_Vectors.List;
    ptr : Facet_4d_List := f;
    lft : Link_to_4d_Facet;

  begin
    while not Is_Null(ptr) loop
      lft := Head_Of(ptr);
      for i in lft.ridges'range loop
        declare
          r : constant Standard_Integer_Vectors.Vector := lft.ridges(i).points;
        begin
          if not Belongs_To(res,r)
           then Lists_of_Integer_Vectors.Append(res,res_last,r);
          end if;
        end;
      end loop;
      ptr := Tail_Of(ptr);
    end loop;
    return res;
  end Ridges;

-- DESTRUCTORS :

  procedure Clear ( f : in out Facet_in_4d ) is
  begin
    for i in f.ridges'range loop
      Standard_Lattice_3d_Facets.Clear(f.ridges(i));
    end loop;
  end Clear;

  procedure Clear ( f : in out Link_to_4d_Facet ) is

    procedure free is 
      new unchecked_deallocation(Facet_in_4d,Link_to_4d_Facet);

  begin
    if f /= null
     then free(f);
    end if;
  end Clear;

  procedure Clear ( f : in out Array_of_4d_Facets ) is
  begin
    for i in f'range loop
      Clear(f(i));
    end loop;
  end Clear;

end Standard_Lattice_4d_Facets;
