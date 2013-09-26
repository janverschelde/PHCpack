with text_io;                          use text_io;
with Communications_with_User;         use Communications_with_User;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;      use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;     use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;      use Standard_Natural_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;     use Standard_Floating_Vectors_io;
with Standard_Floating_VecVecs;        use Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;     use Standard_Floating_VecVecs_io;
with Standard_Random_Vectors;          use Standard_Random_Vectors;
with Standard_Random_VecVecs;          use Standard_Random_VecVecs;

procedure ts_giftwrap is

-- DESCRIPTION :
--   Development of giftwrapping to enumerate all vertices and edges
--   of a polytope.

  procedure Normal_in_Cone
              ( u,v : in Standard_Floating_Vectors.Vector;
                a,b : out double_float;
                w : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Returns a vector in the cone spanned by u and v
  --   that is perpendicular to v, i.e.: w = a*u + b*v.
  --   For the standard inner product denoted by < . , . >,
  --   solves the system
  --        < w , v > = a < u , v > + b < v, v > = 0
  --                  = a           + b          = 1

  -- REQUIRED : u is not parallel to v.

    use Standard_Floating_Vectors;
    vv : constant double_float := v*v;
    uv : constant double_float := u*v;
    d : constant double_float := uv - vv;

  begin
   -- put("d = "); put(d); new_line;
    if 1.0 + d = 1.0 then
      a := 0.0; b := 0.0;
      w := (v'range => 0.0);
    else
      a := - vv/d;
      b := uv/d;
     -- put("a = "); put(a); new_line;
     -- put("b = "); put(b); new_line;
      w := a*u + b*v;
     -- put_line("The normal in the cone :"); put_line(w);
    end if;
  end Normal_in_Cone;

  function Inner_Products
              ( p : Standard_Floating_VecVecs.VecVec;
                h : Standard_Floating_Vectors.Vector )
              return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range equal to p'range that contains for
  --   each point in p the inner product with the vector h.

    res : Standard_Floating_Vectors.Vector(p'range);
    use Standard_Floating_Vectors;

  begin
    for i in p'range loop
      res(i) := p(i).all*h;
    end loop;
    return res;
  end Inner_Products;

  procedure Sort ( p : in out Standard_Floating_VecVecs.VecVec;
                   v : in out Standard_Floating_Vectors.Vector;
                   L : in out Standard_Natural_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Sorts the points in p in increasing order along the values in v,
  --   with corresponding labels in L.

    min : double_float;
    ind : integer32;
    lbl : natural32;
    tmp : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in v'first..v'last-1 loop
      ind := i; min := v(ind);
      for j in i+1..v'last loop
        if v(j) < min
         then ind := j; min := v(ind);
        end if;
      end loop;
      if ind /= i then
        v(ind) := v(i); v(i) := min;
        lbl := L(ind); L(ind) := L(i); L(i) := lbl;
        tmp := p(i); p(i) := p(ind); p(ind) := tmp;
      end if;
    end loop;
  end Sort;

  function Sign_Edge ( ipw : Standard_Floating_Vectors.Vector;
                       i,j : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns +1 if ipw(i) = ipw(j) is minimal,
  --   returns -1 if ipw(i) = ipw(j) is maximal,
  --   returns 0 otherwise.
  --   Note: the -1 means that the normal needs flipping.

  -- REQUIRED : ipw(i) = ipw(j)

    ismin,ismax : boolean := false;

  begin
    for k in ipw'range loop
      if k /= i and k /= j then
        if ipw(k) < ipw(i) then
          ismax := true;
        elsif ipw(k) > ipw(i) then
          ismin := true;
        end if;
      end if;
      if ismin and ismax
       then return 0;
      end if;
    end loop;
    if ismin and not ismax then
      return +1;
    elsif ismax and not ismin then
      return -1;
    else
      return 0;
    end if;
  end Sign_Edge;

  procedure Walk ( n,m,k : in integer32;
                   p : in Standard_Floating_VecVecs.VecVec;
                   h : in out Standard_Floating_VecVecs.VecVec;
                   v : in Standard_Floating_Vectors.Vector;
                   L : in Standard_Natural_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Walks from the current vertex at p(k) to the sink at p'last.

    e,w : Standard_Floating_Vectors.Vector(1..n);
    pw : Standard_Floating_Vectors.Vector(1..m);
    a,b,ew : double_float;
    s : integer32;
    use Standard_Floating_Vectors;

  begin
    for i in k+1..m loop
      e := p(i).all - p(k).all;
     -- put_line("The edge vector e : "); put_line(e);
      Normal_in_Cone(h(k).all,e,a,b,w);
      ew := e*w;
      put("normal condition : "); put(ew,3); new_line;
      pw := Inner_Products(p,w);
     -- put_line("The inner products : "); put_line(pw);
      put(L(k),1); put(" and "); put(L(i),1);
      s := Sign_Edge(pw,k,i);
      if s = 0 then
        put_line(" do not define an edge");
      elsif s = +1 then
        put_line(" define an edge");
      else
        put_line(" define an edge, flip its normal");
      end if;
    end loop;
  end Walk;

  procedure Rank_Points
              ( n,m : in integer32; 
                p : in out Standard_Floating_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Generates a random vector normal to a hyperplane,
  --   use to rank all the points in p.  
  --   Sorting the inner products for this normal with the points in p
  --   picks up the source and sink vertex.

    h : constant Standard_Floating_Vectors.Vector(1..n) := Random_Vector(1,n);
    v : Standard_Floating_Vectors.Vector(p'range) := Inner_Products(p,h);
    L : Standard_Natural_Vectors.Vector(p'range);
    w : Standard_Floating_VecVecs.VecVec(p'range);

  begin
    for i in p'range loop
      L(i) := natural32(i);
    end loop;
    put_line("The inner products : "); put_line(v);
    Sort(p,v,L);
    put_line("The sorted inner products :"); put_line(v);
    put("The labels : "); put(L); new_line;
    w(1) := new Standard_Floating_Vectors.Vector'(h);
    Walk(n,m,1,p,w,v,L);
  end Rank_Points;

  function Ask_Coordinates
             ( n,m : integer32 ) return Standard_Floating_VecVecs.VecVec is

    v : Standard_Floating_VecVecs.VecVec(1..m);
    point : Standard_Floating_Vectors.Vector(1..n) := (1..n => 0.0);

  begin
    put("Reading "); put(m,1); put_line(" points ...");
    for i in 1..m loop
      put("give "); put(n,1);
      put(" coordinates for point "); put(i,1); put_line(" :");
      for j in 1..n loop
        get(point(j)); 
      end loop;
      v(i) := new Standard_Floating_Vectors.Vector'(point);
    end loop;
    return v;
  end Ask_Coordinates;

  procedure Main is

    n,m : integer32 := 0;
    pts : Link_to_VecVec;
    ans : character;

  begin
    new_line;
    put_line("Gift wrapping to enumerate vertices and edges");
    new_line;
    put("Give the dimension : "); get(n);
    put("Give the number of points : "); get(m);
    put("Random points ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then pts := new VecVec'(Random_VecVec(natural32(n),natural32(m)));
     else pts := new VecVec'(Ask_Coordinates(n,m));
    end if;
    put("The "); put(m,1); put(" input ");
    put(n,1); put_line("-vectors :"); put_line(pts.all);
    Rank_Points(n,m,pts.all);
  end Main;

begin
  Main;
end ts_giftwrap;
