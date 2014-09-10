with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;

-- with Multprec_Integer_Numbers_io; use Multprec_Integer_Numbers_io;
-- with Multprec_Integer_Vectors_io; use Multprec_Integer_Vectors_io;
-- with Multprec_Integer_Matrices_io; use Multprec_Integer_Matrices_io;

with unchecked_deallocation;
with Multprec_Common_Divisors;          use Multprec_Common_Divisors;
with Multprec_Power_Transformations;
with Multprec_Lattice_Polygons;
with Standard_Lattice_Supports;
with Multprec_Lattice_Supports;
with Multprec_Lattice_Polytopes;

package body Multprec_Lattice_3d_Facets is

-- INITIAL FACET :

  function Lower ( A : Matrix; i,j : integer32 ) return boolean is
  begin
    for k in A'range(1) loop
      if A(k,i) < A(k,j) then
        return true;
      elsif A(k,i) > A(k,j) then
        return false;
      end if;
    end loop;
    return (A(A'last(1),i) < A(A'last(1),j));
  end Lower;

  function Lowest ( A : Matrix ) return integer32 is

    res : integer32 := A'first(2);

  begin
    for k in A'first(2)+1..A'last(2) loop
      if Lower(A,k,res)
       then res := k;
      end if;
    end loop;
    return res;
  end Lowest;

  function Second_Lowest ( A : Matrix; k : integer32 ) return integer32 is

    res : integer32;
    min : Integer_Number;

  begin
   -- put("In Second_Lowest with k = "); put(k,1); 
    if k /= A'first(2) 
     then res := A'first(2);
     else res := A'first(2) + 1;
    end if;
    Copy(A(A'first(1),res),min);
    for j in A'range(2) loop
      if j /= k then
        if A(A'first(1),j) < min then
          for i in A'range(1) loop
            if not Equal(A(i,j),A(i,k))
             then res := j; Copy(A(A'first(1),j),min); exit;
            end if;
          end loop;
        end if;
      end if;
    end loop;
    Clear(min);
   -- put(" returns "); put(res,1); new_line;
    return res;
  end Second_Lowest;

  function Largest_Angle ( A : Matrix; k : integer32 ) return integer32 is

    res : integer32;
    v1,v2,w1,w2,p,q : Integer_Number;

  begin
   -- put("inside Largest_Angle for k = "); put(k,1);
    if k = A'first(2)
     then res := k+1;
     else res := A'first(2);
    end if;
    v1 := A(A'first(1),res) - A(A'first(1),k);
    v2 := A(A'first(1)+1,res) - A(A'first(1)+1,k);
    for j in A'range(2) loop
      if j /= k then
        w1 := A(A'first(1),j) - A(A'first(1),k);
        w2 := A(A'first(1)+1,j) - A(A'first(1)+1,k);
        p := v1*w2; q := w1*v2;
        if p > q
         then res := j; Copy(w1,v1); Copy(w2,v2);
        end if;
        Clear(p); Clear(q);
      end if;
      Clear(w1); Clear(w2);
    end loop;
    Clear(v1); Clear(v2);
   -- put(" returns "); put(res,1); new_line;
    return res;
  end Largest_Angle;

  function Initial_Edge ( A : Matrix; k : integer32 ) return integer32 is

    res : integer32 := Second_Lowest(A,k);

  begin
   -- put("A("); put(A'first(1),1); put(","); put(res,1); put(") = ");
   -- put(A(A'first(1),res)); put(" ? = ");
   -- put("A("); put(A'first(1),1); put(","); put(k,1); put(") = ");
   -- put(A(A'first(1),k)); new_line;
    if not Equal(A(A'first(1),res),A(A'first(1),k))
     then res := Largest_Angle(A,k);
    end if;
    return res;
  end Initial_Edge;

  function Edge_Normal ( A : Matrix; i,j : integer32 )
                       return Multprec_Integer_Vectors.Vector is

    res : Multprec_Integer_Vectors.Vector(A'range(1));
    d : Integer_Number;

  begin
    for i in res'range loop
      res(i) := Multprec_Integer_Numbers.Create(integer32(0));
    end loop;
    if Equal(A(A'first(1),i),A(A'first(1),j)) then
      res(res'first) := Multprec_Integer_Numbers.Create(integer32(1));
    elsif Equal(A(A'first(1)+1,i),A(A'first(1)+1,j)) then
      res(res'first+1) := Multprec_Integer_Numbers.Create(integer32(1));
    else
      res(res'first+1) := A(A'first(1),i) - A(A'first(1),j);
      res(res'first) := A(A'first(1)+1,j) - A(A'first(1)+1,i);
      if res(res'first) < 0 then
        for i in res'range loop
          if not Equal(res(i),0)
           then Multprec_Integer_Numbers.Min(res(i));
          end if;
        end loop;
      end if;
      d := gcd(res(res'first),res(res'first+1));
      if not Equal(d,1) then
        Div(res(res'first),d);
        Div(res(res'first+1),d);
      end if;
    end if;
    Multprec_Lattice_Supports.Inner(A,i,j,res);
    return res;
  end Edge_Normal;

  function Independent ( A : Matrix; i,j,k : integer32 ) return boolean is

    res : boolean;
    v,w : Multprec_Integer_Vectors.Vector(A'range(1));
    d,e : Integer_Number;

  begin
   -- put("inside independent with i = "); put(i,1);
   -- put(", j = "); put(j,1); put(", and k = "); put(k,1); 
    for p in A'range(1) loop
      v(p) := A(p,j) - A(p,i);
      w(p) := A(p,k) - A(p,i);
    end loop;
    d := v(v'first)*w(w'first+1);
    e := v(v'first+1)*w(w'first);
    Sub(d,e); Clear(e);
    if not Equal(d,0) then
      res := true;
    else
      Clear(d);
      d := v(v'first)*w(w'first+2);
      e := v(v'first+2)*w(w'first);
      Sub(d,e); Clear(e);
      if not Equal(d,0) then
        res := true;
      else
        Clear(d);
        d := v(v'first+1)*w(w'first+2);
        e := v(v'first+2)*w(w'first+1);
        Sub(d,e); Clear(e);
        if not Equal(d,0)
         then res := true;
         else res := false;
        end if;
      end if;
    end if;
    Multprec_Integer_Vectors.Clear(v);
    Multprec_Integer_Vectors.Clear(w);
   -- if res
   --  then put_line(" returning true");
   --  else put_line(" returning false");
   -- end if;
    return res;
  end Independent;

  function Third_Point
              ( A : Matrix; i,j : integer32; m : Integer_Number;
                v : Multprec_Integer_Vectors.Vector ) return integer32 is

    ip : Integer_Number;

  begin
    for k in A'range(2) loop
      if k /= i and k /= j then
        ip := Multprec_Lattice_Supports.Inner_Product(A,k,v);
        if Equal(ip,m) then
          if Independent(A,i,j,k)
           then Clear(ip); return k;
          end if;
        end if;
        Clear(ip);
      end if;
    end loop;
    return 0;
  end Third_Point;

  procedure Normalize ( v : in out Multprec_Integer_Vectors.Vector ) is

    h : Integer_Number := gcd(v(v'first+1),v(v'first+2));
    h2 : Integer_Number;

  begin
    if not Equal(h,1)
     then h2 := gcd(v(v'first),h); Copy(h2,h); Clear(h2);
    end if;
    if not Equal(h,1) and not Equal(h,0) then
      for k in v'range loop
        Div(v(k),h);
      end loop;
    end if;
    Clear(h);
  end Normalize;

  function Normal ( A : Matrix; i,j : integer32;
                    v : Multprec_Integer_Vectors.Vector )
                  return Multprec_Integer_Vectors.Vector is

    res,d : Multprec_Integer_Vectors.Vector(A'range(1));
    h,h2 : Integer_Number;

  begin
   -- put_line("The matrix A : "); put(A);
   -- put("in normal with i = "); put(i,1);
   -- put(" j = "); put(j,1); put(" and v = "); put(v); new_line;
    for k in d'range loop
      res(k) := Multprec_Integer_Numbers.Create(integer32(0));
      d(k) := A(k,j) - A(k,i);
    end loop;
    Multprec_Lattice_Polytopes.Normalize(d);
   -- put("the d = "); put(d); new_line;
    if Equal(d(d'first+2),0) then
      res(res'first+2) := Multprec_Integer_Numbers.Create(integer32(1)); 
    elsif Equal(v(v'first),0) then
      Copy(d(d'first+2),res(res'first));
      res(res'first+1) := Multprec_Integer_Numbers.Create(integer32(0));
      Copy(d(d'first),res(res'first+2)); Min(res(res'first+2));
      h := gcd(res(res'first),res(res'first+2));
      if not Equal(h,1) then
        Div(res(res'first),h);
        Div(res(res'first+2),h);
      end if;
      Clear(h);
    elsif Equal(v(v'first+1),0) then
      res(res'first) := Multprec_Integer_Numbers.Create(integer32(0));
      Copy(d(d'first+2),res(res'first+1));
      Copy(d(d'first+1),res(res'first+2)); Min(res(res'first+2));
      h := gcd(res(res'first+1),res(res'first+2));
      if not Equal(h,1) then
        Div(res(res'first+1),h);
        Div(res(res'first+2),h);
      end if;
      Clear(h);
    else
      Copy(v(v'first+1),res(res'first)); Min(res(res'first));
      Copy(v(v'first),res(res'first+1));
      h := res(res'first)*d(d'first);
      h2 := res(res'first+1)*d(d'first+1);
      Add(h,h2); Clear(h2);
      Copy(h,res(res'first+2)); Min(res(res'first+2));
      Mul(res(res'first),d(d'first+2));
      Mul(res(res'first+1),d(d'first+2));
      Normalize(res);
      Clear(h);
    end if;
    Multprec_Integer_Vectors.Clear(d);
    if res(res'first) < 0
     then Multprec_Integer_Vectors.Min(res);
    end if;
   -- put("Normal returns res = "); put(res); new_line;
    return res;
  end Normal;

  function Shift ( A : Matrix; i,j : integer32 )
                 return Multprec_Integer_Vectors.Vector is

    res : Multprec_Integer_Vectors.Vector(A'range(1));

  begin
   -- put_line("The matrix before the shift : "); put(A);
    for i in A'range(1) loop
      for j in A'range(2) loop
        if Empty(A(i,j)) then
          put("Element ("); put(i,1); put(","); put(j,1);
          put_line(") is empty!");
        elsif Negative(A(i,j)) then
          if Equal(A(i,j),0) then
            put("Element ("); put(i,1); put(","); put(j,1);
            put_line(") equals -0!");
          end if;
        end if;
      end loop;
    end loop;
    for k in A'range(1) loop
      res(k) := A(k,j) - A(k,i);
      declare
        unsresk : constant Natural_Number := Unsigned(res(k));
      begin
        if Equal(unsresk,0) and then Negative(res(k))
         then Min(res(k));
        end if;
      end;
    end loop;
    Multprec_Lattice_Polytopes.Normalize(res);
   -- put_line("The vector after the shift : "); put(res); new_line;
    return res;
  end Shift;

  function Largest_Angle
              ( A : Matrix; i,j : integer32;
                v,w : Multprec_Integer_Vectors.Vector ) return integer32 is

    res : integer32;
    ind : integer32 := A'first(2)-1;
    x1,x2,y1,y2,p,q : Integer_Number;
    d : Multprec_Integer_Vectors.Vector(A'range(1));

  begin
   -- new_line;
   -- put("inside Largest_Angle with i = "); put(i,1);
   -- put(", j = "); put(j,1); put(", v = "); put(v);
   -- put(" and w = "); put(w);
    y1 := Create(integer32(0));
    for k in A'range(2) loop
      if k /= i and k /= j then
        d := Shift(A,i,k);
        Clear(x1); x1 := Multprec_Lattice_Supports.Inner_Product(d,v);
        Clear(y1); y1 := Multprec_Lattice_Supports.Inner_Product(d,w);
        if y1 < 0
         then Min(y1);
        end if;
        ind := k; res := k;
        Multprec_Integer_Vectors.Clear(d);
      end if;
      exit when ((not Equal(y1,0)) and (ind >= A'first(2)));
    end loop;
   -- put("ind = "); put(ind,1);
   -- put("  x1 = "); put(x1); put("  y1 = "); put(y1); new_line;
    for k in ind+1..A'last(2) loop
      if k /= i and k /= j then
        d := Shift(A,i,k);
        x2 := Multprec_Lattice_Supports.Inner_Product(d,v);
        y2 := Multprec_Lattice_Supports.Inner_Product(d,w);
        if y2 < 0
         then Min(y2);
        end if;
       -- put("at k = "); put(k,1); 
       -- put("  x2 = "); put(x2); put("  y2 = "); put(y2); new_line;
        p := x1*y2; q := x2*y1;
        if not Equal(x2,0) and x2 > 0 then
          if p > q
           then Copy(x2,x1); Copy(y2,y1); res := k;
          end if;
        elsif Equal(x2,0) then
          if not Equal(p,0) and p > 0
           then Copy(x2,x1); Copy(y2,y1); res := k;
          end if;
        else
          if p < q 
           then Copy(x2,x1); Copy(y2,y1); res := k;
          end if;
        end if;
        Clear(x2); Clear(y2); Clear(p); Clear(q);
        Multprec_Integer_Vectors.Clear(d);
      end if;
    end loop;
    Clear(x1); Clear(y1);
   -- put(" returns "); put(res,1); new_line;
    return res;
  end Largest_Angle;

  function Extreme_Angle
              ( A : Matrix; i,j : integer32;
                f : Standard_Integer_Vectors.Vector;
                v,w : Multprec_Integer_Vectors.Vector ) return integer32 is

    res : integer32;
    ind : integer32 := A'first(2)-1;
    x1,x2,y1,y2,p,q : Integer_Number;
    d : Multprec_Integer_Vectors.Vector(A'range(1));

  begin
    for k in A'range(2) loop
      if Standard_Lattice_Supports.Member(f,k) < f'first then
        if Independent(A,i,j,k) then
          d := Shift(A,i,k);
          x1 := Multprec_Lattice_Supports.Inner_Product(d,v);
          y1 := Multprec_Lattice_Supports.Inner_Product(d,w);
          if y1 < 0 then Min(y1); end if;
          ind := k; res := k;
          Multprec_Integer_Vectors.Clear(d);
        end if;
      end if;
      exit when ((not Equal(y1,0)) and (ind >= A'first(2)));
    end loop;
    for k in ind+1..A'last(2) loop
      if Standard_Lattice_Supports.Member(f,k) < f'first then
        if Independent(A,i,j,k) then
          d := Shift(A,i,k);
          x2 := Multprec_Lattice_Supports.Inner_Product(d,v);
          y2 := Multprec_Lattice_Supports.Inner_Product(d,w);
          if y2 < 0 then Min(y2); end if;
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
            Clear(p); Clear(q);
          end if;
          Clear(x2); Clear(y2);
          Multprec_Integer_Vectors.Clear(d);
        end if;
      end if;
    end loop;
    Clear(x1); Clear(y1);
    return res;
  end Extreme_Angle;

  function Normal ( u,v : Multprec_Integer_Vectors.Vector )
                  return Multprec_Integer_Vectors.Vector is

    res : Multprec_Integer_Vectors.Vector(u'range);
    acc : Integer_Number;

  begin
    res(res'first) := u(u'first+1)*v(v'first+2);
    acc := u(u'first+2)*v(v'first+1);
    Sub(res(res'first),acc); Clear(acc);
    res(res'first+1) := u(u'first+2)*v(v'first);
    acc := u(u'first)*v(v'first+2);
    Sub(res(res'first+1),acc); Clear(acc);
    res(res'first+2) := u(u'first)*v(v'first+1);
    acc := u(u'first+1)*v(v'first);
    Sub(res(res'first+2),acc); Clear(acc);
    for i in res'first+3..res'last loop
      res(i) := Create(integer32(0));
    end loop;
    Normalize(res);
    return res;
  end Normal;

  function Normal ( A : Matrix; i,j,k : integer32 )
                  return Multprec_Integer_Vectors.Vector is

    res : Multprec_Integer_Vectors.Vector(A'range(1));
    u : Multprec_Integer_Vectors.Vector(A'range(1)) := Shift(A,i,j);
    v : Multprec_Integer_Vectors.Vector(A'range(1)) := Shift(A,i,k);

  begin
   -- put("inside Normal, i = "); 
   -- put(i,1); put(" j = "); put(j,1); put(" k = "); put(k,1);
   -- put(" with u = "); put(u); 
   -- put(" and v = "); put(v); new_line;
    res := Normal(u,v);
    Multprec_Integer_Vectors.Clear(u);
    Multprec_Integer_Vectors.Clear(v);
    return res;
  end Normal;

  function Inner_Normal ( A : Matrix; i,j,k : integer32 )
                        return Multprec_Integer_Vectors.Vector is

    res : Multprec_Integer_Vectors.Vector(A'range(1)) := Normal(A,i,j,k);

  begin
    Multprec_Lattice_Supports.Inner(A,i,j,k,res);
    return res;
  end Inner_Normal;

  procedure Initial_Facet_Normal
              ( A : in Matrix; i,j,k : out integer32;
                v : out Multprec_Integer_Vectors.Vector ) is

    u,w : Multprec_Integer_Vectors.Vector(A'range(1));
    s : Integer_Number;

  begin
    i := Lowest(A);
    j := Initial_Edge(A,i);
    u := Edge_Normal(A,i,j);
   -- put("edge normal : "); put(u); new_line;
    s := Multprec_Lattice_Supports.Minimum(A,u);
    k := Third_Point(A,i,j,s,u);
   -- put("  i = "); put(i,1);
   -- put("  j = "); put(j,1);
   -- put("  k = "); put(k,1); new_line;
    if k /= 0 then
      Multprec_Lattice_Supports.Inner(A,i,j,k,u);
      Multprec_Integer_Vectors.Copy(u,v); -- v := u;
    else
     -- put("  u = "); put(u);
      w := Normal(A,i,j,u);
     -- put("  w = "); put(w);
      k := Largest_Angle(A,i,j,u,w);
     -- put("  k = "); put(k,1);
      v := Inner_Normal(A,i,j,k);
     -- put("  v = "); put(v); new_line;
    end if;
    Multprec_Lattice_Polytopes.Normalize(v);
   -- put("Initial_Facet_Normal returns v = "); put(v); new_line;
  exception
    when others => put_line("exception in Initial_Facet_Normal"); raise;
  end Initial_Facet_Normal;

  function Drop ( A : Matrix; p : integer32 ) return Matrix is

    res : Matrix(A'first(1)..A'last(1)-1,A'range(2));

  begin
    for j in A'range(2) loop
      for i in A'range(1) loop
        if i < p then
          Copy(A(i,j),res(i,j));
        elsif i > p then
          Copy(A(i,j),res(i-1,j));
        end if;
      end loop;
    end loop;
    return res;
  end Drop;

  function Match ( A,B : Matrix; p,i,j : integer32 ) return boolean is
  begin
    for k in A'range(1) loop
      if k < p then
        if not Equal(A(k,i),B(k,j))
         then return false;
        end if;
      elsif k > p then
        if not Equal(A(k,i),B(k-1,j))
         then return false;
        end if;
      end if;
    end loop;
    return true;
  end Match;

  function Match_Vertices
              ( A,V : Matrix; p : integer32 )
              return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(V'range(2)) := (V'range(2) => 0);

  begin
    for j in V'range(2) loop
      for k in A'range(2) loop
        if Match(A,V,p,k,j)
         then res(j) := k; exit;
        end if;
      end loop;
    end loop;
    return res;
  end Match_Vertices;

  function Edges_of_Facet
              ( A : Matrix; v : Multprec_Integer_Vectors.Vector )
              return Facet_in_3d is

    s : constant Standard_Integer_Vectors.Vector
      := Multprec_Lattice_Supports.Support(A,v);
    B : Matrix(A'range(1),s'range);
    p : integer32;
    M : Matrix(A'range(1),A'range(1));
    MB : Matrix(A'range(1),B'range(2));
    B2 : Matrix(1..2,B'range(2));
    w : Multprec_Integer_Vectors.Vector(v'range);

  begin
    B := Multprec_Lattice_Supports.Support_Points(A,s);
   -- put_line("the matrix B before normalization : "); put(B);
    Multprec_Lattice_Polytopes.Normalize(B);
   -- put_line("the matrix B after normalization : "); put(B);
    Multprec_Integer_Vectors.Copy(v,w);
   -- put_line("The vector v before normalization : "); put(v);
    Multprec_Lattice_Polytopes.Normalize(w);
   -- put_line("The after v before normalization : "); put(w);
    p := Multprec_Power_Transformations.Pivot(w);
    M := Multprec_Power_Transformations.Eliminate(w,p);
    MB := M*B;
    B2 := Drop(MB,p);
    Multprec_Lattice_Polygons.Lexicographic_Decreasing_Sort(B2);
    declare
      C2 : constant Matrix := Multprec_Lattice_Polygons.Convex_Hull_2D(B2);
      m2 : constant Standard_Integer_Vectors.Vector := Match_Vertices(MB,C2,p);
      pts : constant Matrix := Multprec_Lattice_Supports.Support_Points(B,m2);
      res : Facet_in_3d(3,m2'last);
    begin
      res.label := 0;
      Multprec_Integer_Vectors.Copy(w,res.normal);
      res.points := Multprec_Lattice_Supports.Indices(A,pts);
      res.neighbors := (m2'range => null);
      Multprec_Integer_Vectors.Clear(w);
      return res;
    end;
 -- exception
 --   when others => put_line("exception in Edges_of_Facet"); raise;
  end Edges_of_Facet;

  function Initial_Facet ( A : Matrix ) return Facet_in_3d is

    i,j,k : integer32;
    v : Multprec_Integer_Vectors.Vector(A'range(1));

  begin
   -- put_line("Calling Initial_Facet with A = "); put(A);
    Initial_Facet_Normal(A,i,j,k,v);
    return Edges_of_Facet(A,v);
 -- exception
 --   when others => put_line("exception in Initial_Facet"); raise;
  end Initial_Facet;

-- ROTATE facet normal :

  procedure Invert ( p : in out Standard_Integer_Vectors.Vector;
                     k : in integer32 ) is

    up,down : integer32 := k;
    cp : constant Standard_Integer_Vectors.Vector(p'range) := p;

  begin
    loop
      up := up + 1;
      if up > p'last
       then up := p'first;
      end if;
      exit when (up = k);
      down := down - 1;
      if down < p'first
       then down := p'last;
      end if;
      p(down) := cp(up);
    end loop;  
  end Invert;

  procedure Connect ( f : in Link_to_3d_Facet; p,i,j : in integer32 ) is

    fp : constant Link_to_3d_Facet := f.neighbors(p);

  begin
    for k in fp.points'range loop
      if fp.points(k) = j then
        fp.neighbors(k) := f;
        if k < fp.points'last then
          if fp.points(k+1) /= i
           then Invert(fp.points,k);
          end if;
        else
          if fp.points(fp.points'first) /= i
           then Invert(fp.points,k);
          end if;
        end if;
      end if;
    end loop;
  end Connect;

  procedure Previous_and_Next_Edge
              ( f : in Link_to_3d_Facet;
                i : in integer32; p,q : out integer32 ) is
  begin
    if i = f.points'first then
      p := f.points'last;
      q := i + 1;
    elsif i = f.points'last then
      p := i - 1;
      q := f.points'first;
    else
      p := i - 1;
      q := i + 1;
    end if;
  end Previous_and_Next_Edge;

  procedure Connect ( f,g : in Link_to_3d_Facet; i,j : in integer32;
                      connected : out boolean ) is

    f_prev,f_next,g_prev,g_next : integer32;

  begin
    connected := false;
    Previous_and_Next_Edge(f,i,f_prev,f_next);
    Previous_and_Next_Edge(g,j,g_prev,g_next);
    if f.neighbors(f_prev) = null then
      if g.neighbors(j) = null then
        if f.points(f_prev) = g.points(g_next) then
          f.neighbors(f_prev) := g;
          g.neighbors(j) := f;
          connected := true;
        end if;
      end if;
      if g.neighbors(g_next) = null then
        if f.points(f_prev) = g.points(g_prev) then
          f.neighbors(f_prev) := g;
          g.neighbors(g_next) := f;
          connected := true;
        end if;
      end if;
    end if;
    if f.neighbors(i) = null then
      if g.neighbors(g_prev) = null then
        if f.points(f_next) = g.points(g_prev) then
          f.neighbors(i) := g;
          g.neighbors(g_prev) := f;
          connected := true;
        end if;
      end if;
      if g.neighbors(j) = null then
        if f.points(f_next) = g.points(g_next) then
          f.neighbors(i) := g;
          g.neighbors(j) := f;
          connected := true;
        end if;
      end if;
    end if;
  end Connect;

  procedure Connect ( f,g : in Link_to_3d_Facet ) is

    connected : boolean := false;

  begin
    for i in f.points'range loop
      for j in g.points'range loop
        if f.points(i) = g.points(j)
         then Connect(f,g,i,j,connected);
        end if;
        exit when connected;
      end loop;
      exit when connected;
    end loop;
  end Connect;

  function Is_Connected ( f : Link_to_3d_Facet ) return boolean is
  begin
    for i in f.neighbors'range loop
      if f.neighbors(i) = null
       then return false;
      end if;
    end loop;
    return true;
  end Is_Connected;

  procedure Neighbors ( A : in Matrix; f : in Link_to_3d_Facet;
                        idcnt : in out integer32 ) is

    u,v,w : Multprec_Integer_Vectors.Vector(A'range(1));
    i,j,k : integer32;

  begin
    for p in f.points'range loop
      if f.neighbors(p) = null then
        i := f.points(p);
        if p < f.points'last
         then j := f.points(p+1);
         else j := f.points(f.points'first);
        end if;
        u := Shift(A,i,j);
       -- put("u = "); put(u); new_line;
        w := Normal(u,f.normal);
       -- put("w = "); put(w); new_line;
        Multprec_Lattice_Supports.Inner(A,i,j,f.points,w);
        k := Extreme_Angle(A,i,j,f.points,w,f.normal); 
        v := Inner_Normal(A,i,j,k);
       -- put("Inner_Normal returned v = "); put(v); new_line;
        f.neighbors(p) := new Facet_in_3d'(Edges_of_Facet(A,v));
        idcnt := idcnt + 1;
        f.neighbors(p).label := idcnt;
        Connect(f,p,i,j);
        for q in f.neighbors'range loop
          if q /= p and f.neighbors(q) /= null
           then Connect(f.neighbors(q),f.neighbors(p));
          end if;
        end loop;
        Multprec_Integer_Vectors.Clear(u);
        Multprec_Integer_Vectors.Clear(v);
        Multprec_Integer_Vectors.Clear(w);
      end if;
    end loop;
  end Neighbors;

-- MAIN ALGORITHM :

  function Pop ( f : Facet_3d_List ) return Link_to_3d_Facet is

    tmp : Facet_3d_List := f;
    lft : Link_to_3d_Facet;

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

  procedure Connect ( f : in Facet_3d_List; lft : in Link_to_3d_Facet ) is

    tmp : Facet_3d_List := f;
    ptr : Link_to_3d_Facet;

  begin
    while not Is_Null(tmp) loop
      ptr := Head_Of(tmp); 
      if not Is_Connected(ptr)
       then Connect(ptr,lft);
      end if;
      exit when Is_Connected(lft); -- we stop if lft is connected
      tmp := Tail_Of(tmp);
    end loop;
  end Connect;

  function Convex_Hull_3D ( A : Matrix ) return Facet_3d_List is

    res : Facet_3d_List;
    first : constant Link_to_3d_Facet := new Facet_in_3d'(Initial_Facet(A));
    lft : Link_to_3d_Facet;
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
  end Convex_Hull_3D;

-- SELECTORS :

  function Support_Value_of_Facet
             ( A : Matrix; f : Facet_in_3d ) return Integer_Number is

    use Multprec_Lattice_Supports;

    ind : constant integer32 := f.points(f.points'first);
    res : Integer_Number := Inner_Product(A,ind,f.normal);

  begin
    return res;
  end Support_Value_of_Facet;

  function Facet_Normals
             ( f : Facet_3d_List ) return Lists_of_Integer64_Vectors.List is

    res,res_last : Lists_of_Integer64_Vectors.List;
    tmp : Facet_3d_List := f;
    lft : Link_to_3d_Facet;

  begin
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      -- Lists_of_Integer64_Vectors.Append(res,res_last,lft.normal);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Facet_Normals;

  function Is_Facet_Normal
              ( f : Facet_3d_List; v : Multprec_Integer_Vectors.Vector )
              return boolean is

    tmp : Facet_3d_List := f;
    lft : Link_to_3d_Facet;

  begin
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      if Multprec_Integer_Vectors.Equal(lft.normal,v)
       then return true;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return false;
  end Is_Facet_Normal;

  function Select_Facet_Normals
              ( f : Facet_3d_List; v : Lists_of_Integer64_Vectors.List )
              return Lists_of_Integer64_Vectors.List is

    res,res_last : Lists_of_Integer64_Vectors.List;
    tmp : Lists_of_Integer64_Vectors.List := v;
   -- lv : Standard_Integer64_Vectors.Link_to_Vector;

  begin
    while not Lists_of_Integer64_Vectors.Is_Null(tmp) loop
     -- lv := Lists_of_Integer64_Vectors.Head_Of(tmp);
     -- if Is_Facet_Normal(f,lv.all)
     --  then Lists_of_Integer64_Vectors.Append(res,res_last,lv.all);
     -- end if;
      tmp := Lists_of_Integer64_Vectors.Tail_Of(tmp);
    end loop;
    return res;
  end Select_Facet_Normals;

  function Get_Facet
              ( f : Facet_3d_List; v : Multprec_Integer_Vectors.Vector )
              return Link_to_3d_Facet is    

    tmp : Facet_3d_List := f;
    lft : Link_to_3d_Facet;

  begin
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      if Multprec_Integer_Vectors.Equal(lft.normal,v)
       then return lft;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return null;
  end Get_Facet;

  function Edges ( f : Facet_3d_List ) return Lists_of_Integer_Vectors.List is

    res,res_last : Lists_of_Integer_Vectors.List;
    tmp : Facet_3d_List := f;
    lft : Link_to_3d_Facet;
    edge : Standard_Integer_Vectors.Vector(1..2);

  begin
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      for i in lft.points'range loop
        edge(1) := lft.points(i);
        if i < lft.points'last
         then edge(2) := lft.points(i+1);
         else edge(2) := lft.points(lft.points'first);
        end if;
        Lists_of_Integer_Vectors.Append(res,res_last,edge);
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Edges;

  function Vertices ( n : integer32; f : Facet_3d_List )
                    return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..n) := (1..n => 0);
    ind : integer32 := 0;
    tmp : Facet_3d_List := f;
    lft : Link_to_3d_Facet;
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

  function Vertices ( n : integer32; f : Array_of_3d_Facets )
                    return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..n) := (1..n => 0);
    ind : integer32 := 0;
    lft : Link_to_3d_Facet;
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

  function First_Incident_Vertex ( f : Facet_3d_List ) return integer32 is

    lft : Link_to_3d_Facet;

  begin
    if Is_Null(f) then
      return 0;
    else
      lft := Head_Of(f);
      if lft = null
       then return 0;
       else return lft.points(lft.points'first);
      end if;
    end if;
  end First_Incident_Vertex;

  procedure Check_Euler_Characteristic
               ( m : in integer32; f : in Facet_3d_List ) is

    e : Lists_of_Integer_Vectors.List := Edges(f);
    v : constant Standard_Integer_Vectors.Vector := Vertices(m,f);
    nf : constant natural32 := Length_Of(f);
    ne : constant natural32 := Lists_of_Integer_Vectors.Length_Of(e);
    nv : constant natural32 := natural32(v'last);

  begin
    put("#facets : "); put(nf,1);
    put("  #edges : "); put(ne/2,1);
    put("  #vertices : "); put(nv,1);
    put("  Euler characteristic : "); put(nf-ne/2+nv,1); new_line;
    Lists_of_Integer_Vectors.Clear(e);
  end Check_Euler_Characteristic;

-- ENUMERATORS : walk to enumerate vertices, crawl for edges

  procedure Walk ( f : in Facet_3d_List; v : in integer32;
                   b : in out Boolean_Array ) is

    tmp : Facet_3d_List := f;
    lft : Link_to_3d_Facet;
    ind,w : integer32;
    c : boolean := true;

  begin
    b(v) := true;
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      ind := Standard_Lattice_Supports.Member(lft.points,v);
      if ind >= lft.points'first then  -- v is incident to lft
        if ind > lft.points'first      -- w will be previous neighbor
         then w := lft.points(ind-1);
         else w := lft.points(lft.points'last);
        end if;
        if not b(w)
         then Report(w,c); b(w) := true; if c then Walk(f,w,b); end if;
        end if;
        if ind < lft.points'last       -- w will be next neighbor
         then w := lft.points(ind+1);
         else w := lft.points(lft.points'first);
        end if;
        if not b(w)
         then Report(w,c); b(w) := true; if c then Walk(f,w,b); end if;
        end if;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Walk;

  procedure Crawl ( f : in Facet_3d_List; v : in integer32;
                    b : in out Boolean_Matrix ) is

    tmp : Facet_3d_List := f;
    lft : Link_to_3d_Facet;
    ind,w : integer32;
    c : boolean := true;

  begin
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      ind := Standard_Lattice_Supports.Member(lft.points,v);
      if ind >= lft.points'first then  -- v is incident to lft
        if ind > lft.points'first      -- w will be previous neighbor
         then w := lft.points(ind-1);
         else w := lft.points(lft.points'last);
        end if;
        if not b(v,w) then
          Report(v,w,c); b(v,w) := true; b(w,v) := true;
          if c then Crawl(f,w,b); end if;
        end if;
        if ind < lft.points'last       -- w will be next neighbor
         then w := lft.points(ind+1);
         else w := lft.points(lft.points'first);
        end if;
        if not b(v,w) then
          Report(v,w,c); b(v,w) := true; b(w,v) := true;
          if c then Crawl(f,w,b); end if;
        end if;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Crawl;

-- DESTRUCTORS :

  procedure Clear ( f : in out Link_to_3d_Facet ) is

    procedure free is 
      new unchecked_deallocation(Facet_in_3d,Link_to_3d_Facet);

  begin
    if f /= null then
     -- Clear(f.neighbors); -- else infinite loop !
      free(f);
    end if;
  end Clear;

  procedure Clear ( f : in out Array_of_3d_Facets ) is
  begin
    for i in f'range loop
      Clear(f(i));
    end loop;
  end Clear;

  procedure Clear ( f : in out Facet_3d_List ) is

    tmp : Facet_3d_List := f;

  begin
    while not Is_Null(tmp) loop
      declare
        ft : Link_to_3d_Facet := Head_Of(tmp);
      begin
        Clear(ft);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_3d_Facets.Clear(Lists_of_3d_Facets.List(f));
  end Clear;

end Multprec_Lattice_3d_Facets;
