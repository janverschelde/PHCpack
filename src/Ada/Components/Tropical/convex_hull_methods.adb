with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer64_Vectors_io;      use Standard_Integer64_Vectors_io;
with Standard_Integer64_Matrices_io;     use Standard_Integer64_Matrices_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Standard_Lattice_Supports;
with Standard_Lattice_Polygons;
with Standard_Lattice_3d_Facets_io;
with Standard_Lattice_4d_Facets_io;      use Standard_Lattice_4d_Facets_io;
with Standard_Lattice_Polytopes;
with Standard_Lattice_Facets;
with Standard_Lattice_Facets_io;         use Standard_Lattice_Facets_io;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Numbers_io;        use Multprec_Integer_Numbers_io;
with Multprec_Integer_Vectors_io;        use Multprec_Integer_Vectors_io;
with Multprec_Integer_Matrices_io;       use Multprec_Integer_Matrices_io;
with Multprec_Lattice_Polygons;
with Multprec_Lattice_Supports;
with Multprec_Lattice_3d_Facets_io;
with Multprec_Lattice_4d_Facets_io;      use Multprec_Lattice_4d_Facets_io;
with Multprec_Lattice_Polytopes;

package body Convex_Hull_Methods is

  function Is_In ( A : Standard_Integer64_Matrices.Matrix;
                   k : integer32 ) return boolean is

    found : boolean := false;

  begin
    for j in A'first(2)..(k-1) loop
      found := true;
      for i in A'range(1) loop
        if A(i,j) /= A(i,k)
         then found := false; exit;
        end if;
      end loop;
      if found
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Is_In ( A : Multprec_Integer_Matrices.Matrix;
                   k : integer32 ) return boolean is

    found : boolean := false;

  begin
    for j in A'first(2)..(k-1) loop
      found := true;
      for i in A'range(1) loop
        if not Equal(A(i,j),A(i,k))
         then found := false; exit;
        end if;
      end loop;
      if found
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Filter_Duplicates
             ( A : Standard_Integer64_Matrices.Matrix )
             return Standard_Integer64_Matrices.Matrix is

    B : Standard_Integer64_Matrices.Matrix(A'range(1),A'range(2));
    ind : integer32 := A'first(2);

  begin
    for i in A'range(1) loop
      B(i,ind) := A(i,ind);
    end loop;
    for j in A'first(2)+1..A'last(2) loop
      if not Is_In(A,j) then
        ind := ind + 1;
        for i in A'range(1) loop
          B(i,ind) := A(i,j);
        end loop;
      end if;
    end loop;
    declare
      res : Standard_Integer64_Matrices.Matrix(A'range(1),B'first(2)..ind);
    begin
      for i in res'range(1) loop
        for j in res'range(2) loop
          res(i,j) := B(i,j);
        end loop;
      end loop;
      return res;
    end;
  end Filter_Duplicates;

  function Filter_Duplicates
             ( A : Multprec_Integer_Matrices.Matrix )
             return Multprec_Integer_Matrices.Matrix is

    B : Multprec_Integer_Matrices.Matrix(A'range(1),A'range(2));
    ind : integer32 := A'first(2);

  begin
    for i in A'range(1) loop
      B(i,ind) := A(i,ind);
    end loop;
    for j in A'first(2)+1..A'last(2) loop
      if not Is_In(A,j) then
        ind := ind + 1;
        for i in A'range(1) loop
          B(i,ind) := A(i,j);
        end loop;
      end if;
    end loop;
    declare
      res : Multprec_Integer_Matrices.Matrix(A'range(1),B'first..ind);
    begin
      for i in res'range(1) loop
        for j in res'range(2) loop
          res(i,j) := B(i,j);
        end loop;
      end loop;
      return res;
    end;
  end Filter_Duplicates;

  function User_Input ( n,m : integer32 )
                      return Standard_Integer64_Matrices.Matrix is

    res : Standard_Integer64_Matrices.Matrix(1..n,1..m);

  begin
    put("Give a "); put(n,1); put("-by-"); put(m,1);
    put_line(" matrix :"); get(res);
    return res;
  end User_Input;

  function User_Input ( n,m : integer32 )
                      return Multprec_Integer_Matrices.Matrix is

    res : Multprec_Integer_Matrices.Matrix(1..n,1..m);

  begin
    put("Give a "); put(n,1); put("-by-"); put(m,1);
    put_line(" matrix :"); get(res);
    return res;
  end User_Input;

  function Random_Data ( n,m : integer32 )
                       return Standard_Integer64_Matrices.Matrix is

    res : Standard_Integer64_Matrices.Matrix(1..n,1..m);
    lower,upper : integer64 := 0;

  begin
    put("Give a lower bound for the coordinates : "); get(lower);
    put("Give an upper bound for the coordinates : "); get(upper);
    res := Random_Matrix(natural32(n),natural32(m),lower,upper);
    return res;
  end Random_Data;

  function Convert ( A : Standard_Integer64_Matrices.Matrix )
                   return Multprec_Integer_Matrices.Matrix is

    res : Multprec_Integer_Matrices.Matrix(A'range(1),A'range(2));
    int : integer32;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        int := integer32(A(i,j));
        res(i,j) := Multprec_Integer_Numbers.Create(int);
      end loop;
    end loop;
    return res;
  end Convert;

  function Random_Data ( n,m : integer32 )
                       return Multprec_Integer_Matrices.Matrix is

    res : Multprec_Integer_Matrices.Matrix(1..n,1..m);
    A : Standard_Integer64_Matrices.Matrix(1..n,1..m);
   -- size : natural := 0;

  begin
   -- put("Give the size of the random numbers : "); get(size);
   -- res := Random_Matrix(n,m,size);
    A := Random_Data(n,m);
    res := Convert(A);
    return res;
  end Random_Data;

  function Cyclic_Polytope
              ( d : integer32; t : Standard_Integer64_Vectors.Vector )
              return Standard_Integer64_Matrices.Matrix is

    res : Standard_Integer64_Matrices.Matrix(1..d,t'range);

  begin
    for j in t'range loop
      res(1,j) := t(j);
      for i in 2..d loop
        res(i,j) := res(i-1,j)*t(j);
      end loop;
    end loop;
    return res;
  end Cyclic_Polytope;

  function Cyclic_Polytope
              ( d : integer32; t : Multprec_Integer_Vectors.Vector )
              return Multprec_Integer_Matrices.Matrix is

    res : Multprec_Integer_Matrices.Matrix(1..d,t'range);

  begin
    for j in t'range loop
      Copy(t(j),res(1,j));
      for i in 2..d loop
        res(i,j) := res(i-1,j)*t(j);
      end loop;
    end loop;
    return res;
  end Cyclic_Polytope;

  function Cyclic_Polytope 
              ( d,n : integer32 )
              return Standard_Integer64_Matrices.Matrix is

    res : Standard_Integer64_Matrices.Matrix(1..d,1..n);
    t : Standard_Integer64_Vectors.Vector(1..n);
    ans : character;
    lower,upper : integer64 := 0;

  begin
    put("-> generate random points at the moment curve ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give a lower bound for the coordinates : "); get(lower);
      put("Give an upper bound for the coordinates : "); get(upper);
      t := Random_Vector(1,n,lower,upper);
    else
      put("Give "); put(n,1); put(" integers : "); get(t);
    end if;
    res := Cyclic_Polytope(d,t);
    return res;
  end Cyclic_Polytope;

  function Cyclic_Polytope 
              ( d,n : integer32 )
              return Multprec_Integer_Matrices.Matrix is

    res : Multprec_Integer_Matrices.Matrix(1..d,1..n);
    t : Multprec_Integer_Vectors.Vector(1..n);
    rt : Standard_Integer64_Vectors.Vector(1..n);
    ans : character;
    lower,upper : integer64 := 0;

  begin
    put("-> generate random points at the moment curve ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give a lower bound for the coordinates : "); get(lower);
      put("Give an upper bound for the coordinates : "); get(upper);
      rt := Random_Vector(1,n,lower,upper);
      for i in rt'range loop
        t(i) := Multprec_Integer_Numbers.Create(integer(rt(i)));
      end loop;
    else
      put("Give "); put(n,1); put(" integers : "); get(t);
    end if;
    res := Cyclic_Polytope(d,t);
    return res;
  end Cyclic_Polytope;

  procedure Standard_Planar_Hull
              ( A : in Standard_Integer64_Matrices.Matrix ) is

    use Standard_Integer64_Matrices,Standard_Lattice_Polygons;
    B : Standard_Integer64_Matrices.Matrix(A'range(1),A'range(2)) := A;

  begin
    put_line("The support set : "); put(A);
    Lexicographic_Decreasing_Sort(B);
    put_line("After lexicographic sort "); put(B);
    declare
      V : constant Matrix := Convex_Hull_2D(B);
      N : constant Matrix := Inner_Normals(V);
      bug : boolean;
    begin
      put_line("The vertices : "); put(V);
      put_line("The inner normals : "); put(N);
      Check(B,V,N,true,bug);
      if bug
       then put_line("Reported a bug...");
       else put_line("Passed rank checks.");
      end if;
    end;
  end Standard_Planar_Hull;

  procedure Multprec_Planar_Hull
              ( A : in Multprec_Integer_Matrices.Matrix ) is

    use Multprec_Integer_Matrices,Multprec_Lattice_Polygons;
    B : Multprec_Integer_Matrices.Matrix(A'range(1),A'range(2));

  begin
    Copy(A,B);
    put_line("The support set : "); put(A);
    Lexicographic_Decreasing_Sort(B);
    put_line("After lexicographic sort : "); put(B);
    declare
      V : constant Matrix := Convex_Hull_2D(B);
      N : constant Matrix := Inner_Normals(V);
      bug : boolean;
    begin
      put_line("The vertices : "); put(V);
      put_line("The inner normals : "); put(N);
      Check(B,V,N,true,bug);
      if bug
       then put_line("Reported a bug...");
       else put_line("Passed rank checks.");
      end if;
    end;
  end Multprec_Planar_Hull;

  procedure Standard_Show_Initial_Facet
               ( A : in Standard_Integer64_Matrices.Matrix ) is

    use Standard_Lattice_3d_Facets,Standard_Lattice_3d_Facets_io;

    k,e,t,f : integer32;
    s,h : integer64;
    w,u,g : Standard_Integer64_Vectors.Vector(A'range(1));

  begin
    k := Lowest(A);
    put("The lexicographically lowest column : "); put(k,1); new_line;
    e := Initial_Edge(A,k);
    put("The initial edge is spanned by"); Write_Coordinates(A,k);
    put(" and"); Write_Coordinates(A,e); put_line(".");
    w := Edge_Normal(A,k,e);
    put("-> inner normal to the initial edge :"); put(w); 
    s := Standard_Lattice_Supports.Minimum(A,w);
    put(" minimum value is "); put(s,1); new_line;
    put("with inner products :");
    put(Standard_Lattice_Supports.Inner_Products(A,w)); put_line(".");
    t := Third_Point(A,k,e,s,w);
    if t = 0
     then put_line("There is no facet with inner normal to initial edge.");
     else put("Third linearly independent point :");
          Write_Coordinates(A,t); new_line;
    end if;
    u := Normal(A,k,e,w);
    put("The normal to edge and its edge normal :"); put(u); new_line;
    put("-> its product with edge normal : ");
    put(u(1)*w(1) + u(2)*w(2) + u(3)*w(3),1);
    put("  values at edge points : ");  
    put(Standard_Lattice_Supports.Inner_Product(A,k,u),1); put(" ");
    put(Standard_Lattice_Supports.Inner_Product(A,e,u),1); put_line(".");
    f := Largest_Angle(A,k,e,w,u);
    put("Index of point with largest angle : "); put(f,1); new_line;
    put("The initial facet is spanned by"); Write_Coordinates(A,k);
    put(","); Write_Coordinates(A,e); put(", and"); Write_Coordinates(A,f);
    new_line;
    g := Inner_Normal(A,k,e,f);
    put("The facet inner normal :"); put(g); 
    h := Standard_Lattice_Supports.Minimum(A,g);
    put(" with minimum value : "); put(h,1); new_line;
    put("  with inner products :");
    put(Standard_Lattice_Supports.Inner_Products(A,g)); new_line;
    put("  and support :");
    put(Standard_Lattice_Supports.Support(A,g)); put_line(".");
  end Standard_Show_Initial_Facet; 

  procedure Multprec_Show_Initial_Facet
               ( A : in Multprec_Integer_Matrices.Matrix ) is

    use Multprec_Lattice_3d_Facets,Multprec_Lattice_3d_Facets_io;

    k,e,t,f : integer32;
    s,h : Integer_Number;
    w,u,g : Multprec_Integer_Vectors.Vector(A'range(1));

  begin
    k := Lowest(A);
    put("The lexicographically lowest column : "); put(k,1); new_line;
    e := Initial_Edge(A,k);
    put("The initial edge is spanned by"); Write_Coordinates(A,k);
    put(" and"); Write_Coordinates(A,e); put_line(".");
    w := Edge_Normal(A,k,e);
    put("-> inner normal to the initial edge :"); put(w); 
    s := Multprec_Lattice_Supports.Minimum(A,w);
    put(" minimum value is "); put(s,1); new_line;
    put("with inner products :");
    put(Multprec_Lattice_Supports.Inner_Products(A,w)); put_line(".");
    t := Third_Point(A,k,e,s,w);
    if t = 0
     then put_line("There is no facet with inner normal to initial edge.");
     else put("Third linearly independent point :");
          Write_Coordinates(A,t); new_line;
    end if;
    u := Normal(A,k,e,w);
    put("The normal to edge and its edge normal :"); put(u); new_line;
    put("-> its product with edge normal : ");
    put(u(1)*w(1) + u(2)*w(2) + u(3)*w(3),1);
    put("  values at edge points : ");  
    put(Multprec_Lattice_Supports.Inner_Product(A,k,u)); put(" ");
    put(Multprec_Lattice_Supports.Inner_Product(A,e,u)); put_line(".");
    f := Largest_Angle(A,k,e,w,u);
    put("Index of point with largest angle : "); put(f,1); new_line;
    put("The initial facet is spanned by"); Write_Coordinates(A,k);
    put(","); Write_Coordinates(A,e); put(", and"); Write_Coordinates(A,f);
    new_line;
    g := Inner_Normal(A,k,e,f);
    put("The facet inner normal :"); put(g); 
    h := Multprec_Lattice_Supports.Minimum(A,g);
    put(" with minimum value : "); put(h,1); new_line;
    put("  with inner products :");
    put(Multprec_Lattice_Supports.Inner_Products(A,g)); new_line;
    put("  and support :");
    put(Multprec_Lattice_Supports.Support(A,g)); put_line(".");
  end Multprec_Show_Initial_Facet; 

  procedure Standard_Show_Neighbors
               ( A : in Standard_Integer64_Matrices.Matrix;
                 f : in Standard_Lattice_3d_Facets.Link_to_3d_Facet ) is

    use Standard_Lattice_3d_Facets,Standard_Lattice_3d_Facets_io;

    u,v,w : Standard_Integer64_Vectors.Vector(A'range(1));
    i,j,k : integer32;
    lbl : integer32 := f.label;

  begin
    put("Indices of points spanning facet : "); put(f.points); new_line;
    for p in f.points'range loop
      put_line("*********************************************************");
      i := f.points(p);
      if p < f.points'last
       then j := f.points(p+1);
       else j := f.points(f.points'first);
      end if;
      u := Shift(A,i,j);
      w := Normal(u,f.normal);
      Standard_Lattice_Supports.Inner(A,i,j,f.points,w);
      put("Normal to edge and facet : "); put(w); new_line;
      put("  product with edge : "); 
      put(Standard_Lattice_Supports.Inner_product(w,u),1);
      put("  product with facet : "); 
      put(Standard_Lattice_Supports.Inner_product(w,f.normal),1);
      new_line;
      put("  IP :"); put(Standard_Lattice_Supports.Inner_Products(A,w));
      put("  support : "); put(Standard_Lattice_Supports.Support(A,w));
      new_line;
      k := Extreme_Angle(A,i,j,f.points,w,f.normal); 
      put("index of extreme angle : "); put(k,1); new_line;
      v := Inner_Normal(A,i,j,k);
      f.neighbors(p) := new Facet_in_3d'(Edges_of_Facet(A,v));
      lbl := lbl + 1;
      f.neighbors(p).label := lbl;
      Connect(f,p,i,j);
      Write_Facet(A,f.neighbors(p).all);
    end loop;
    put_line("*********************************************************");
  end Standard_Show_Neighbors;

  procedure Multprec_Show_Neighbors
               ( A : in Multprec_Integer_Matrices.Matrix;
                 f : in Multprec_Lattice_3d_Facets.Link_to_3d_Facet ) is

    use Multprec_Lattice_3d_Facets,Multprec_Lattice_3d_Facets_io;

    u,v,w : Multprec_Integer_Vectors.Vector(A'range(1));
    i,j,k : integer32;
    lbl : integer32 := f.label;

  begin
    put("Indices of points spanning facet : "); put(f.points); new_line;
    for p in f.points'range loop
      put_line("*********************************************************");
      i := f.points(p);
      if p < f.points'last
       then j := f.points(p+1);
       else j := f.points(f.points'first);
      end if;
      u := Shift(A,i,j);
      w := Normal(u,f.normal);
      Multprec_Lattice_Supports.Inner(A,i,j,f.points,w);
      put("Normal to edge and facet : "); put(w); new_line;
      put("  product with edge : "); 
      put(Multprec_Lattice_Supports.Inner_product(w,u),1);
      put("  product with facet : "); 
      put(Multprec_Lattice_Supports.Inner_product(w,f.normal),1);
      new_line;
      put("  IP :"); put(Multprec_Lattice_Supports.Inner_Products(A,w));
      put("  support : "); put(Multprec_Lattice_Supports.Support(A,w));
      new_line;
      k := Extreme_Angle(A,i,j,f.points,w,f.normal); 
      put("index of extreme angle : "); put(k,1); new_line;
      v := Inner_Normal(A,i,j,k);
      f.neighbors(p) := new Facet_in_3d'(Edges_of_Facet(A,v));
      lbl := lbl + 1;
      f.neighbors(p).label := lbl;
      Connect(f,p,i,j);
      Write_Facet(A,f.neighbors(p).all);
    end loop;
    put_line("*********************************************************");
  end Multprec_Show_Neighbors;

  procedure Standard_Start_3D_Giftwrapping
              ( A : in Standard_Integer64_Matrices.Matrix ) is

    use Standard_Lattice_3d_Facets,Standard_Lattice_3d_Facets_io;

    ans : character;
    f : constant Link_to_3d_Facet := new Facet_in_3d'(Initial_Facet(A));
    lbl : integer32 := 0;

  begin
    put("Extra output for initial facet construction ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Standard_Show_Initial_Facet(A);
    end if;
    put_line("The initial facet is spanned by ");
    put(Standard_Lattice_Supports.Support_Points(A,f.points));
    put_line("Computing neighboring facets ...");
    put("Do you want extra output ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Standard_Show_Neighbors(A,f);
     else Neighbors(A,f,lbl);
    end if;
    Write_Facet(A,f.all);
    for i in f.neighbors'range loop
      Write_Facet(A,f.neighbors(i).all);
    end loop;
  end Standard_Start_3D_Giftwrapping;

  procedure Multprec_Start_3D_Giftwrapping
              ( A : in Multprec_Integer_Matrices.Matrix ) is

    use Multprec_Lattice_3d_Facets,Multprec_Lattice_3d_Facets_io;

    ans : character;
    f : constant Link_to_3d_Facet := new Facet_in_3d'(Initial_Facet(A));
    lbl : integer32 := 0;

  begin
    put("Extra output for initial facet construction ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Multprec_Show_Initial_Facet(A);
    end if;
    put_line("The initial facet is spanned by ");
    put(Multprec_Lattice_Supports.Support_Points(A,f.points));
    put_line("Computing neighboring facets ...");
    put("Do you want extra output ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Multprec_Show_Neighbors(A,f);
     else Neighbors(A,f,lbl);
    end if;
    Write_Facet(A,f.all);
    for i in f.neighbors'range loop
      Write_Facet(A,f.neighbors(i).all);
    end loop;
  end Multprec_Start_3D_Giftwrapping;

  procedure Walk_Vertices
               ( m : in integer32;
                 f : in Standard_Lattice_3d_Facets.Facet_3d_List ) is

    b : Standard_Lattice_3d_Facets.Boolean_Array(1..m) := (1..m => false);
    v : constant integer32
      := Standard_Lattice_3d_Facets.First_Incident_Vertex(f);
    cnt : integer32 := 1;

    procedure Write ( w : in integer32; continue : out boolean ) is
    begin
      cnt := cnt + 1;
      put(" "); put(w,1);
      continue := true;
    end Write;
    procedure Do_Walk is new Standard_Lattice_3d_Facets.Walk(Write);

  begin
    put("Vertices : "); put(v,1);
    Do_Walk(f,v,b);
    put("  count : "); put(cnt,1); new_line;
  end Walk_Vertices;

  procedure Walk_Edges
               ( m : in integer32;
                 f : in Standard_Lattice_3d_Facets.Facet_3d_List ) is

    b : Standard_Lattice_3d_Facets.Boolean_Matrix(1..m,1..m);
    v : constant integer32
      := Standard_Lattice_3d_Facets.First_Incident_Vertex(f);
    cnt : integer32 := 0;

    procedure Write ( p,q : in integer32; continue : out boolean ) is
    begin
      cnt := cnt + 1;
      if cnt mod 8 = 0
       then new_line; put("  ");
      end if;
      put(" ("); put(p,1); put(","); put(q,1); put(")");
      continue := true;
    end Write;
    procedure Do_Walk is new Standard_Lattice_3d_Facets.Crawl(Write);

  begin
    for i in b'range(1) loop
      for j in b'range(2) loop
        b(i,j) := false;
      end loop;
    end loop;
    put("Edges : ");
    Do_Walk(f,v,b);
    put("  count : "); put(cnt,1); new_line;
  end Walk_Edges;

  procedure Check_Edge
              ( e : in Standard_Lattice_Edges.Link_to_Edge;
                bug : out boolean ) is
  begin
    bug := false;
    put("edge "); put(e.label,1); put(" is intersected by");
    put(e.f.normal); put(" and"); put(e.g.normal); new_line;
    put("  vertex "); put(e.a,1);
    if Standard_Lattice_Supports.Member(e.f.points,e.a) > 0 and
       Standard_Lattice_Supports.Member(e.g.points,e.a) > 0
     then put(" belongs");
     else put(" does NOT belong, BUG!"); bug := true;
    end if;
    put("  vertex "); put(e.b,1);
    if Standard_Lattice_Supports.Member(e.f.points,e.b) > 0 and
       Standard_Lattice_Supports.Member(e.g.points,e.b) > 0
     then put_line(" belongs");
     else put_line(" does NOT belong, BUG!"); bug := true;
    end if;
    if bug then
      put(e.f.normal); put(" supports :"); put(e.f.points); new_line;
      put(e.g.normal); put(" supports :"); put(e.g.points); new_line;
    end if;
  end Check_Edge;

  procedure Check_Edge
              ( e : in Multprec_Lattice_Edges.Link_to_Edge;
                bug : out boolean ) is
  begin
    bug := false;
    put("edge "); put(e.label,1); put(" is intersected by");
    put(e.f.normal); put(" and"); put(e.g.normal); new_line;
    put("  vertex "); put(e.a,1);
    if Standard_Lattice_Supports.Member(e.f.points,e.a) > 0 and
       Standard_Lattice_Supports.Member(e.g.points,e.a) > 0
     then put(" belongs");
     else put(" does NOT belong, BUG!"); bug := true;
    end if;
    put("  vertex "); put(e.b,1);
    if Standard_Lattice_Supports.Member(e.f.points,e.b) > 0 and
       Standard_Lattice_Supports.Member(e.g.points,e.b) > 0
     then put_line(" belongs");
     else put_line(" does NOT belong, BUG!"); bug := true;
    end if;
    if bug then
      put(e.f.normal); put(" supports :"); put(e.f.points); new_line;
      put(e.g.normal); put(" supports :"); put(e.g.points); new_line;
    end if;
  end Check_Edge;

  procedure Check_Edges
              ( A : in Standard_Integer64_Matrices.Matrix;
                f : in Standard_Lattice_3d_Facets.Facet_3d_List ) is

    use Standard_Lattice_Edges;

    e : constant Edge_List := Edges_of_3D_Hull(A'last(2),f);
    tmp : Edge_List := e;
    lnk : Link_to_Edge;
    bug : boolean := false;
   
  begin
    put("Checking "); put(Length_Of(e),1); put_line(" edges ...");
    while not Is_Null(tmp) loop
      lnk := Head_Of(tmp);
      Check_Edge(lnk,bug);
      exit when bug;
      tmp := Tail_Of(tmp);
    end loop;
    if bug
     then put_line("BUG detected in edge connectivity!");
     else put_line("No bugs found in edge connectivity.");
    end if;
  end Check_Edges;

  procedure Check_Edges
              ( A : in Multprec_Integer_Matrices.Matrix;
                f : in Multprec_Lattice_3d_Facets.Facet_3d_List ) is

    use Multprec_Lattice_Edges;

    e : constant Edge_List := Edges_of_3D_Hull(A'last(2),f);
    tmp : Edge_List := e;
    lnk : Link_to_Edge;
    bug : boolean := false;
   
  begin
    put("Checking "); put(Length_Of(e),1); put_line(" edges ...");
    while not Is_Null(tmp) loop
      lnk := Head_Of(tmp);
      Check_Edge(lnk,bug);
      exit when bug;
      tmp := Tail_Of(tmp);
    end loop;
    if bug
     then put_line("BUG detected in edge connectivity!");
     else put_line("No bugs found in edge connectivity.");
    end if;
  end Check_Edges;

  procedure Standard_3D_Giftwrapping
                ( A : in Standard_Integer64_Matrices.Matrix ) is

    use Standard_Lattice_3d_Facets,Standard_Lattice_3d_Facets_io;

    f : Facet_3d_List;
    e : Lists_of_Integer_Vectors.List;
    nf,ne,nv : integer32;

  begin
    put_line("The support : "); put(A);
    f := Convex_Hull_3D(A); nf := integer32(Length_Of(f));
    put("Found "); put(nf,1); put_line(" facets :");
    Write_Facets(A,f);
    e := Edges(f); ne := integer32(Lists_of_Integer_Vectors.Length_Of(e));
    put("Found "); put(ne,1); put_line(" edges :"); put(e);
    declare
      v : constant Standard_Integer_Vectors.Vector := Vertices(A'last(2),f);
    begin
      nv := v'last;
      put("Found "); put(v'last,1); put(" vertices : "); put(v);
      new_line;
    end;
    put("#facets : "); put(nf,1);
    put("  #edges : "); put(ne/2,1);
    put("  #vertices : "); put(nv,1); new_line;
    put("Euler characteristic : "); put(nf-ne/2+nv,1); new_line;
    Walk_Vertices(A'last(2),f);
    Walk_Edges(A'last(2),f);
    Check_Edges(A,f);
  end Standard_3D_Giftwrapping;

  procedure Multprec_3D_Giftwrapping
                ( A : in Multprec_Integer_Matrices.Matrix ) is

    use Multprec_Lattice_3d_Facets,Multprec_Lattice_3d_Facets_io;

    f : Facet_3d_List;
    nf : natural32;

  begin
    put_line("The support : "); put(A);
    f := Convex_Hull_3D(A); nf := Length_Of(f);
    put("Found "); put(nf,1); put_line(" facets :");
    Write_Facets(A,f);
    Check_Euler_Characteristic(A'last(2),f);
    Check_Edges(A,f);
  end Multprec_3D_Giftwrapping;

  procedure Standard_3D_Hull
              ( A : in Standard_Integer64_Matrices.Matrix ) is

    ans : character;

  begin
    put_line("The support set : "); put(A);
    put("Show only start of giftwrapping ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Standard_Start_3D_Giftwrapping(A);
     else Standard_3D_Giftwrapping(A);
    end if;
  end Standard_3D_Hull;

  procedure Multprec_3D_Hull
              ( A : in Multprec_Integer_Matrices.Matrix ) is

    ans : character;

  begin
    put_line("The support set : "); put(A);
    put("Show only start of giftwrapping ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Multprec_Start_3D_Giftwrapping(A);
     else Multprec_3D_Giftwrapping(A);
    end if;
  end Multprec_3D_Hull;

  procedure Standard_Start_4D_hull
              ( A : in Standard_Integer64_Matrices.Matrix;
                v : in Standard_Integer64_Vectors.Vector ) is

    use Standard_Lattice_4D_Facets;
    init,g : Link_to_4d_Facet;
    lbl : integer32 := 0;

  begin
    new_line;
    put_line("starting giftwrapping in 4-space ...");
    init := new Facet_in_4d'(Ridges_of_Facet(A,v));
    Write_4D_Facet(A,init);
    put_line("computing the neighbors ...");
    Neighbors(A,init,lbl);
    put_line("-> done computing the neighbors");
    for i in init.neighbors'range loop
      g := init.neighbors(i);
      put("neighbor "); put(i,1); put_line(" :");
      Write_4D_Facet(A,g);
    end loop;
    put_line("Checking on the initial facet : ");
    Write_4D_Facet(A,init);
  end Standard_Start_4D_Hull;

  procedure Multprec_Start_4D_hull
              ( A : in Multprec_Integer_Matrices.Matrix;
                v : in Multprec_Integer_Vectors.Vector ) is

    use Multprec_Lattice_4D_Facets;
    init,g : Link_to_4d_Facet;
    lbl : integer32 := 0;

  begin
    new_line;
    put_line("starting giftwrapping in 4-space ...");
    init := new Facet_in_4d'(Ridges_of_Facet(A,v));
    Write_4D_Facet(A,init);
    put_line("computing the neighbors ...");
    Neighbors(A,init,lbl);
    for i in init.neighbors'range loop
      g := init.neighbors(i);
      put("neighbor "); put(i,1); put_line(" :");
      Write_4D_Facet(A,g);
    end loop;
  end Multprec_Start_4D_Hull;

  procedure Standard_4D_Euler_Characteristic
              ( n : in integer32;
                f : in Standard_Lattice_4d_Facets.Facet_4d_List ) is

    use Standard_Lattice_4d_Facets;

    nf : constant integer32 := integer32(Length_Of(f));
    v : constant Standard_Integer_Vectors.Vector := Vertices(n,f);
    nv : constant integer32 := v'last;
    r : Lists_of_Integer_Vectors.List := Ridges(f);
    nr : constant integer32
       := integer32(Lists_of_Integer_Vectors.Length_Of(r));
    e : Lists_of_Integer_Vectors.List := Edges(f);
    ne : constant integer32
       := integer32(Lists_of_Integer_Vectors.Length_Of(e));
    ch : constant integer32 := nf - nr + ne - nv;

  begin
   -- put_line("Points spanning the edges :"); put(e);
   -- put_line("Points spanning the ridges :"); put(r);
    put("Computed "); put(nf,1); put(" facets");
    put(" - "); put(nr,1); put(" ridges");
    put(" + "); put(ne,1); put(" edges");
    put(" - "); put(nv,1); put(" vertices :"); put(v); new_line;
    put(" = "); put(ch,1); put_line(" Euler characteristic");
    Lists_of_Integer_Vectors.Clear(r);
    Lists_of_Integer_Vectors.Clear(e);
  end Standard_4D_Euler_Characteristic;

  procedure Multprec_4D_Euler_Characteristic
              ( n : in integer32;
                f : in Multprec_Lattice_4d_Facets.Facet_4d_List ) is

    use Multprec_Lattice_4d_Facets;

    nf : constant integer32 := integer32(Length_Of(f));
    v : constant Standard_Integer_Vectors.Vector := Vertices(n,f);
    nv : constant integer32 := v'last;
    r : Lists_of_Integer_Vectors.List := Ridges(f);
    nr : constant integer32
       := integer32(Lists_of_Integer_Vectors.Length_Of(r));
    e : Lists_of_Integer_Vectors.List := Edges(f);
    ne : constant integer32
       := integer32(Lists_of_Integer_Vectors.Length_Of(e));
    ch : constant integer32 := nf - nr + ne - nv;

  begin
   -- put_line("Points spanning the edges :"); put(e);
   -- put_line("Points spanning the ridges :"); put(r);
    put("Computed "); put(nf,1); put(" facets");
    put(" - "); put(nr,1); put(" ridges");
    put(" + "); put(ne,1); put(" edges");
    put(" - "); put(nv,1); put(" vertices :"); put(v); new_line;
    put(" = "); put(ch,1); put_line(" Euler characteristic");
    Lists_of_Integer_Vectors.Clear(r);
    Lists_of_Integer_Vectors.Clear(e);
  end Multprec_4D_Euler_Characteristic;

  procedure Standard_4D_Hull
              ( A : in Standard_Integer64_Matrices.Matrix ) is

    use Standard_Lattice_4D_Facets,Standard_Lattice_Polytopes;

    r : natural32;
    v : Standard_Integer64_Vectors.Vector(A'range(1));
    ans : character;
    ch : Facet_4d_List;

  begin
    new_line;
    put_line("The support set :"); put(A);
    put("Do only start of giftwrapping ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Initial_Facet_Normal(A,r,v);
      if integer32(r) < A'last(1) then
        put("The point configuration is of rank "); put(r,1);
        put(" < "); put(A'last(1),1); new_line;
      else
        put_line("The point configuration has full rank.");
        Standard_Start_4D_Hull(A,v);
      end if;
    else
      Initial_Facet_Normal(A,r,v);
      ch := Convex_Hull_4D(A,v);
      Write_4D_Facets(A,ch);
      Standard_4D_Euler_Characteristic(A'last(2),ch);
    end if;
  end Standard_4D_Hull;

  procedure Multprec_4D_Hull
              ( A : in Multprec_Integer_Matrices.Matrix ) is

    use Multprec_Lattice_4d_Facets,Multprec_Lattice_Polytopes;

    r : natural32;
    v : Multprec_Integer_Vectors.Vector(A'range(1));
    ans : character;
    ch : Facet_4D_List;

  begin
    new_line;
    put_line("The support set :"); put(A);
    put("Do only start of giftwrapping ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Initial_Facet_Normal(A,r,v);
      if integer32(r) < A'last(1) then
        put("The point configuration is of rank "); put(r,1);
        put(" < "); put(A'last(1),1); new_line;
      else
        put_line("The point configuration has full rank.");
        Multprec_Start_4D_Hull(A,v);
      end if;
    else
      Initial_Facet_Normal(A,r,v);
      ch := Convex_Hull_4D(A,v);
      Write_4D_Facets(A,ch);
      Multprec_4D_Euler_Characteristic(A'last(2),ch);
    end if;
  end Multprec_4D_Hull;

  procedure Standard_Start_hull
              ( A : in Standard_Integer64_Matrices.Matrix;
                v : in Standard_Integer64_Vectors.Vector ) is

    use Standard_Lattice_Facets;
    d : constant integer32 := v'last;
    init,g : Link_to_Facet;
    lbl : integer32 := 0;

  begin
    new_line;
    put("starting giftwrapping in "); put(d,1); put_line("-space ...");
    init := new Facet'(Ridges_of_Facet(A,v));
    Write_Facet(A,init);
    put_line("computing the neighbors ...");
    Neighbors(A,init,lbl);
    put_line("-> done computing the neighbors");
    for i in init.neighbors'range loop
      g := init.neighbors(i);
      put("neighbor "); put(i,1); put_line(" :");
      Write_Facet(A,g);
    end loop;
  end Standard_Start_Hull;

  procedure Standard_General_Hull
              ( A : in Standard_Integer64_Matrices.Matrix ) is

    use Standard_Lattice_Facets,Standard_Lattice_Polytopes;

    r : natural32;
    v : Standard_Integer64_Vectors.Vector(A'range(1));
    ans : character;
    ch : Facet_List;

  begin
    new_line;
    put_line("The support set :"); put(A);
    put("Do only start of giftwrapping ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Initial_Facet_Normal(A,r,v);
      put("the initial facet normal : "); put(v); new_line;
      if integer32(r) < A'last(1) then
        put("The point configuration is of rank "); put(r,1);
        put(" < "); put(A'last(1),1); new_line;
      else
        put_line("The point configuration has full rank.");
        Standard_Start_Hull(A,v);
      end if;
    else
      Initial_Facet_Normal(A,r,v);
      ch := Convex_Hull(A,v);
      Write_Facets(A,ch);
      declare
        d : constant integer32 := A'last(1);
        n : constant integer32 := A'last(2);
        fv : constant Standard_Natural_Vectors.Vector := fvector(n,d,ch);
      begin
        put("The f-vector : "); put(fv); new_line;
      end;
    end if;
  end Standard_General_Hull;

  procedure Multprec_General_Hull
              ( A : in Multprec_Integer_Matrices.Matrix ) is

    use Multprec_Lattice_Polytopes;

    r : natural32;
    v : Multprec_Integer_Vectors.Vector(A'range(1));
   -- ans : character;

  begin
    new_line;
    put_line("The support set :"); put(A);
   -- put("Do only start of giftwrapping ? (y/n) ");
   -- Ask_Yes_or_No(ans);
   -- if ans = 'y' then
      Initial_Facet_Normal(A,r,v);
      if integer32(r) < A'last(1) then
        put("The point configuration is of rank "); put(r,1);
        put(" < "); put(A'last(1),1); new_line;
      else
        put_line("The point configuration has full rank.");
      end if;
   -- end if;
  end Multprec_General_Hull;

end Convex_Hull_Methods;
