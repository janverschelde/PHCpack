-- with text_io; use text_io;
-- with Standard_Integer_Numbers_io; use Standard_Integer_Numbers_io;
-- with Multprec_Integer_Matrices_io; use Multprec_Integer_Matrices_io;
-- with Multprec_Integer_Numbers_io; use Multprec_Integer_Numbers_io;
-- with Standard_Integer_Vectors_io; use Standard_Integer_Vectors_io;
-- with Multprec_Integer_Vectors_io; use Multprec_Integer_Vectors_io;

with unchecked_deallocation;
with Standard_Lattice_Supports;
with Multprec_Lattice_Supports;
with Multprec_Integer_Orthogonals;
with Multprec_Power_Transformations;
with Multprec_Lattice_Polytopes;

package body Multprec_Lattice_4d_Facets is

-- CONSTRUCTORS :

  function Insert_Zero
              ( v : Multprec_Integer_Vectors.Vector; p : integer32 )
              return Multprec_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Inserts a zero at position p of the inner normal v.

    res : Multprec_Integer_Vectors.Vector(v'first..v'last+1);

  begin
    for i in v'first..(p-1) loop
      Copy(v(i),res(i));
    end loop;
    res(p) := Multprec_Integer_Numbers.Create(integer32(0));
    for i in p..v'last loop
      Copy(v(i),res(i+1));
    end loop;
    return res;
  end Insert_Zero;

  procedure Multiply_with_Transpose
              ( v : in out Multprec_Integer_Vectors.Vector;
                M : in Matrix ) is

  -- DESCRIPTION :
  --   Multiplies v with the transpose of M.

    w : Multprec_Integer_Vectors.Vector(v'range);
    acc : Integer_Number;

  begin
    for k in v'range loop
      w(k) := Multprec_Integer_Numbers.Create(integer32(0));
      for i in M'range(1) loop
        if not Equal(v(i),0) then
          if not Equal(M(i,k),0) then
            acc := M(i,k)*v(i);
            Add(w(k),acc);
            Clear(acc);
          end if;
        end if;
      end loop;
    end loop;
   -- Multprec_Integer_Vectors.Clear(v); v := w;
    Multprec_Integer_Vectors.Copy(w,v);
    Multprec_Integer_Vectors.Clear(w);
  end Multiply_with_Transpose;

  function Relabel ( p,s : Standard_Integer_Vectors.Vector )
                   return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the relabeling of p (points supported on a ridge),
  --   using the indices of the supported points on the facet.

    res : Standard_Integer_Vectors.Vector(p'range);

  begin
    for i in p'range loop
      res(i) := s(p(i));
    end loop;
    return res;
  end Relabel;

  procedure Convex_Hull_of_Ridge
              ( A : in Matrix; v : in Multprec_Integer_Vectors.Vector;
                s : in Standard_Integer_Vectors.Vector;
                p : out integer32; M : out Matrix;
                F : out Multprec_Lattice_3d_Facets.Facet_3d_List ) is

    B : Matrix(A'range(1),s'range);
    MB : Matrix(A'range(1),s'range);
    B2 : Matrix(1..3,s'range);
    w : Multprec_Integer_Vectors.Vector(v'range);

  begin
   -- put_line("In Convex_Hull_of_Ridge ...");
   -- put_line("The matrix A :"); put(A);
   -- put("with v = " ); put(v); new_line;
    B := Multprec_Lattice_Supports.Support_Points(A,s);
   -- put_line("ridge is spanned by the poins : "); put(B);
    Multprec_Lattice_Polytopes.Normalize(B);
   -- put_line("The matrix B after normalization : "); put(B);
   -- put_line("The vector v before normalize : "); put(v);
    Multprec_Integer_Vectors.Copy(v,w);
    Multprec_Lattice_Polytopes.Normalize(w);
   -- put_line("The vector v after normalize : "); put(w);
    p := Multprec_Power_Transformations.Pivot(w);
   -- put("the pivot p = "); put(p,1); new_line;
    M := Multprec_Power_Transformations.Eliminate(w,p);
   -- put_line("The matrix M : "); put(M);
    MB := M*B;
   -- put_line("The matrix MB : "); put(MB);
    B2 := Multprec_Lattice_3d_Facets.Drop(MB,p);
   -- put_line("The matrix B2 before normalization : "); put(B2);
    Multprec_Lattice_Polytopes.Normalize(B2);
   -- put_line("B2 after dropping one coordinate : "); put(B2);
    F := Multprec_Lattice_3d_Facets.Convex_Hull_3D(B2);
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
              ( A : Matrix; v : Multprec_Integer_Vectors.Vector )
              return Facet_in_4d is

    s : constant Standard_Integer_Vectors.Vector
      := Multprec_Lattice_Supports.Support(A,v);
    p : integer32;
    M : Matrix(A'range(1),A'range(1));
    Cf : Multprec_Lattice_3d_Facets.Facet_3d_List;
    d : constant integer32 := v'last;
    n : constant integer32 := s'last;

  begin
    Convex_Hull_of_Ridge(A,v,s,p,M,Cf);
    declare
      mc : constant integer32
         := integer32(Multprec_Lattice_3d_Facets.Length_Of(Cf));
      res : Facet_in_4d(d,n,mc);
      ptr : Multprec_Lattice_3d_Facets.Facet_3d_List;
    begin
     -- put("in Ridges_of_Facet with mc = "); put(mc,1); new_line;
      res.label := 0; res.normal := v; res.points := s;
      ptr := Cf;
      for i in res.ridges'range loop
        declare
          lft : constant Multprec_Lattice_3d_Facets.Link_to_3d_Facet
              := Multprec_Lattice_3d_Facets.Head_Of(ptr);
          fct : Multprec_Lattice_3d_Facets.Facet_in_3d(lft.d+1,lft.n);
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
          res.ridges(i) := new Multprec_Lattice_3d_Facets.Facet_in_3d'(fct);
        end;
        ptr := Multprec_Lattice_3d_Facets.Tail_Of(ptr);
      end loop;
      res.neighbors := (1..mc => null);
      Multprec_Lattice_3d_Facets.Clear(Cf);
      declare
        vtp : constant Standard_Integer_Vectors.Vector
            := Multprec_Lattice_3d_Facets.Vertices(A'last(2),res.ridges);
      begin
       -- put("vertex points : "); put(vtp);
       -- put("  points in f : "); put(res.points); new_line;
        if vtp'last /= res.points'last
         then return Filter_non_Vertices(res,vtp);
         else return res;
        end if;
      end;
    end;
  end Ridges_of_Facet;

  function Extreme_Angle
              ( A : Matrix; f : Standard_Integer_Vectors.Vector;
                v,w : Multprec_Integer_Vectors.Vector ) return integer32 is

    res : integer32;
    ind : integer32 := A'first(2) - 1;
    x1,y1,x2,y2,p,q : Integer_Number;
    d : Multprec_Integer_Vectors.Vector(A'range(1));
 
  begin
   -- put("inside Extreme_Angle ...");
   -- put(" f = "); put(f); put(" v = "); put(v);
   -- put(" w = "); put(w); new_line;
    y1 := Create(integer32(0));
    for k in A'range(2) loop
      if Standard_Lattice_Supports.Member(f,k) < f'first then
        d := Multprec_Lattice_3d_Facets.Shift(A,f(f'first),k);
        y1 := Multprec_Lattice_Supports.Inner_Product(d,v);
        x1 := Multprec_Lattice_Supports.Inner_Product(d,w);
        if y1 < 0 then Min(y1); end if;
        ind := k; res := k;
        Multprec_Integer_Vectors.Clear(d);
      end if;
      exit when (not Equal(y1,0) and (ind >= A'first(2)));
    end loop;
   -- put("  x1 = "); put(x1); put("  y1 = "); put(y1);
   -- put("  ind = "); put(ind,1); new_line;
    for k in (ind+1)..A'last(2) loop
      if Standard_Lattice_Supports.Member(f,k) < f'first then
        d := Multprec_Lattice_3d_Facets.Shift(A,f(f'first),k);
        y2 := Multprec_Lattice_Supports.Inner_Product(d,v);
        x2 := Multprec_Lattice_Supports.Inner_Product(d,w);
       -- put("  x2 = "); put(x2); put("  y2 = "); put(y2); new_line;
        if y2 < 0 then Min(y2); end if;
       -- code copied from Standard_Lattice_3d_Facets.Extreme_Angle
        if not Equal(y2,0) then
          if Equal(x2,0) then
            if x1 > 0
             then Copy(x2,x1); Copy(y2,y1); res := k;
            end if;
          elsif x2 < 0 and (x1 > 0 or Equal(x1,0)) then
            Copy(x2,x1); Copy(y2,y1); res := k;
          elsif x1 > 0 and x2 > 0 then
            p := x1*y2; q := x2*y1;
            if p > q
             then Copy(x2,x1); Copy(y2,y1); res := k;
            end if;
            Clear(p); Clear(q);
         -- else case x1 < 0 and x2 > 0: no new extreme
          elsif x1 < 0 and x2 < 0 then
            p := x1*y2; q := x2*y1;
            if p > q
             then Copy(x2,x1); Copy(y2,y1); res := k;
            end if;
          end if;
        end if;
        Clear(x2); Clear(y2);
        Multprec_Integer_Vectors.Clear(d);
      end if;
    end loop;
    Clear(x1); Clear(y1);
   -- put("Extreme_Angle returns res = "); put(res,1); new_line;
    return res;
  end Extreme_Angle;

  function Inner_Normal
              ( A : Matrix; p : Standard_Integer_Vectors.Vector;
                k : integer32 ) return Multprec_Integer_Vectors.Vector is

    res : Multprec_Integer_Vectors.Vector(A'range(1));
    B : Matrix(A'range(1),p'range);
    s,t : Integer_Number;

  begin
    for j in p'range loop
      for i in A'range(1) loop
        B(i,j) := A(i,p(j)) - A(i,k);
      end loop;
    end loop;
    res := Multprec_Integer_Orthogonals.Complement(B);
    Clear(B);
    s := Multprec_Lattice_Supports.Inner_Product(A,k,res);
    Copy(s,t);
    for i in A'range(2) loop
      if i /= k and Standard_Lattice_Supports.Member(p,i) < p'first then
        Clear(t);
        t := Multprec_Lattice_Supports.Inner_Product(A,i,res);
      end if;
      exit when not Multprec_Integer_Numbers.Equal(t,s);
    end loop;
    if s > t
     then Multprec_Integer_Vectors.Min(res);
    end if;
    Clear(s); Clear(t);
    return res;
  end Inner_Normal;

  function Is_Same ( a,b : Standard_Integer_Vectors.Vector ) return boolean is

  -- DESCRIPTION :
  --   The index sets a and b are the same if they have the same length
  --   and every element of a occurs as well in b.

  begin
    if a'length /= b'length then
      return false;
    else
      for i in a'range loop
        if Standard_Lattice_Supports.Member(b,a(i)) < b'first
         then return false;
        end if;
      end loop;
    end if;
    return true;
  end Is_Same;

  procedure Connect_along_Ridge
               ( f : in Link_to_4d_Facet; r : in integer32 ) is

  -- DESCRIPTION :
  --   Connects f.neighbor(r) to f along the corresponding ridge.

    rdg : constant Multprec_Lattice_3d_Facets.Link_to_3d_Facet := f.ridges(r);
    g : constant Link_to_4d_Facet := f.neighbors(r);

  begin
    for i in g.ridges'range loop
      declare
        lft : constant Multprec_Lattice_3d_Facets.Link_to_3d_Facet
            := g.ridges(i);
      begin
        if Is_Same(lft.points,rdg.points)
         then g.neighbors(i) := f; exit;
        end if;
      end;
    end loop;
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
  end Connect_along_Ridge;

  procedure Neighbors ( A : in Matrix; f : in Link_to_4d_Facet;
                        idcnt : in out integer32 ) is

    k : integer32;
    v : Multprec_Integer_Vectors.Vector(A'range(1));
    connected : boolean;

  begin
    for r in f.ridges'range loop
      if f.neighbors(r) = null then
        k := Extreme_Angle(A,f.ridges(r).points,f.normal,f.ridges(r).normal);
       -- put("extreme angle at column "); put(k,1); new_line;
        v := Inner_Normal(A,f.ridges(r).points,k);
       -- put("the inner normal v : "); put(v); put(" supports ");
       -- put(Multprec_Lattice_Supports.Support(A,v)); new_line;
       -- put(" IP :");
       -- put(Multprec_Lattice_Supports.Inner_Products(A,v)); new_line;
        f.neighbors(r) := new Facet_in_4d'(Ridges_of_Facet(A,v));
        idcnt := idcnt + 1;
        f.neighbors(r).label := idcnt;
        Connect_along_Ridge(f,r);
        for i in f.neighbors'range loop
          if i /= r and f.neighbors(i) /= null then
           -- put("making connections for new neighbor "); put(i,1);
           -- put_line(" ...");
            Connect_along_Ridge(f.neighbors(r),f.neighbors(i),connected);
          end if;
        end loop;
       -- Multprec_Integer_Vectors.Clear(v);
       -- the link to v is assigned in Ridges_of_Facet(A,v) !!!!
      end if;
    end loop;
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
      exit when Is_Connected(lft);
      tmp := Tail_Of(tmp);
    end loop;
  end Connect;

  function Convex_Hull_4D
              ( A : Matrix; v : Multprec_Integer_Vectors.Vector )
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
     -- put("constructing the neigbhors of facet cnt = ");
     -- put(cnt,1); new_line;
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
  end Convex_Hull_4D;

-- SELECTORS :

  function Support_Value_of_Facet
             ( A : Matrix; f : Facet_in_4d ) return Integer_Number is

    use Multprec_Lattice_Supports;

    ind : constant integer32 := f.points(f.points'first);
    res : Integer_Number := Inner_Product(A,ind,f.normal);

  begin
    return res;
  end Support_Value_of_Facet;

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

  -- DESCRIPTION :
  --   Returns true if p is the same as one of the vectors in R.

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
      Multprec_Lattice_3d_Facets.Clear(f.ridges(i));
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

end Multprec_Lattice_4d_Facets;
