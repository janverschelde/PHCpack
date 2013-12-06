--with Standard_Integer64_Matrices_io; use Standard_Integer64_Matrices_io;
with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
--with Standard_Integer64_Vectors_io;      use Standard_Integer64_Vectors_io;
with Standard_Integer64_Linear_Solvers;  use Standard_Integer64_Linear_Solvers;
with Standard_Lattice_Supports;
with Standard_Lattice_3d_Facets;
--with Standard_Lattice_3d_Facets_io;      use Standard_Lattice_3d_Facets_io;
with Standard_Integer_Orthogonals;

package body Standard_Lattice_Polytopes is

-- INITIAL FACET :

  function Shift ( A : Matrix; k : integer32 ) return Matrix is

  -- DESCRIPTION :
  --   The matrix on return has the k-th column removed, and the coordinates
  --   all other columns have been subtracted by the k-th column.

    res : Matrix(A'range(1),A'first(2)..A'last(2)-1);

  begin
    for j in A'range(2) loop
      if j < k then
        for i in A'range(1) loop
          res(i,j) := A(i,j) - A(i,k);
        end loop;
      elsif j > k then
        for i in A'range(1) loop
          res(i,j-1) := A(i,j) - A(i,k);
        end loop;
      end if;
    end loop;
    return res;
  end Shift;

  function Rank_of_Upper ( A : Matrix ) return natural32 is

  -- DESCRIPTION :
  --   Returns the rank of an upper triangular matrix.

    res : natural32 := 0;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        if A(i,j) /= 0
         then res := res + 1; exit;
        end if;
      end loop;
    end loop;
    return res;
  end Rank_of_Upper;

  function Rank ( A : Matrix ) return natural32 is

    res : natural32;
    k : constant integer32 := Standard_Lattice_3d_Facets.Lowest(A);
    V : Matrix(A'range(1),A'first(2)..A'last(2)-1) := Shift(A,k);

  begin
   -- put_line("V before Upper_Triangulate : "); put(V); new_line;
    Upper_Triangulate(V);
   -- put_line("V after Upper_Triangulate : "); put(V); new_line;
    res := Rank_of_Upper(V);
    return res;
  end Rank;

  function Normal ( A : Matrix; v : Standard_Integer64_Vectors.Vector )
                  return Standard_Integer64_Vectors.Vector is

    res : Standard_Integer64_Vectors.Vector(v'range);
    s : constant Standard_Integer_Vectors.Vector
      := Standard_Lattice_Supports.Support(A,v);
    B : constant Matrix(A'range(1),s'range)
      := Standard_Lattice_Supports.Support_Points(A,s);
    C : Matrix(A'range(1),B'range(2));

  begin
    for j in C'first(2)..C'last(2)-1 loop
      for i in C'range(1) loop
        C(i,j) := B(i,j+1) - B(i,j);
      end loop;
    end loop;
    for i in C'range(1) loop
      C(i,C'last(2)) := v(i);
    end loop;
    res := Standard_Integer_Orthogonals.Complement(C);
    return res;
  end Normal;

  procedure Inner ( A : in Matrix; k : in integer32;
                    v : in out Standard_Integer64_Vectors.Vector ) is

    m : constant integer64 := Standard_Lattice_Supports.Minimum(A,v);
    p : constant integer64 := Standard_Lattice_Supports.Inner_Product(A,k,v);

  begin
    if p > m
     then Standard_Integer64_Vectors.Min(v);
    end if;
  end Inner;

  function Normal ( A : Matrix; v : Standard_Integer64_Vectors.Vector;
                    k : integer32 ) return Standard_Integer64_Vectors.Vector is

    res : Standard_Integer64_Vectors.Vector(v'range);
    s : constant Standard_Integer_Vectors.Vector
      := Standard_Lattice_Supports.Support(A,v);
    B : constant Matrix(A'range(1),s'range)
      := Standard_Lattice_Supports.Support_Points(A,s);
    C : Matrix(A'range(1),B'range(2));

  begin
    for j in C'first(2)..C'last(2)-1 loop
      for i in C'range(1) loop
        C(i,j) := B(i,j+1) - B(i,B'first(2));
      end loop;
    end loop;
    for i in C'range(1) loop
      C(i,C'last(2)) := A(i,k) - B(i,B'first(2));
    end loop;
    res := Standard_Integer_Orthogonals.Complement(C);
    Inner(A,k,res);
    return res;
  end Normal;

  function Largest_Angle
              ( A : Matrix; f,g : Standard_Integer64_Vectors.Vector )
              return integer32 is

    res,ind : integer32 := 0;
    s : constant Standard_Integer_Vectors.Vector
      := Standard_Lattice_Supports.Support(A,f);
    w : Standard_Integer64_Vectors.Vector(A'range(1));
    x1,y1,x2,y2 : integer64;

  begin
    y1 := 0;
    for k in A'range(2) loop
      if Standard_Lattice_Supports.Member(s,k) < s'first then
        w := Standard_Lattice_3d_Facets.Shift(A,s(s'first),k);
        x1 := Standard_Lattice_Supports.Inner_Product(w,f);
        y1 := Standard_Lattice_Supports.Inner_Product(w,g);
        if y1 < 0 then y1 := -y1; end if;
        ind := k; res := k;
      end if;
      exit when ((y1 /= 0) and (ind > 0));
    end loop;
    put("inside Largest_Angle with ind = "); put(ind,1); new_line;
    for k in ind+1..A'last(2) loop
      if Standard_Lattice_Supports.Member(s,k) < s'first then
        w := Standard_Lattice_3d_Facets.Shift(A,s(s'first),k);
        x2 := Standard_Lattice_Supports.Inner_Product(w,f);
        y2 := Standard_Lattice_Supports.Inner_Product(w,g);
        put("  x1 = "); put(x1,1); put("  y1 = "); put(y1,1);
        put("  x2 = "); put(x2,1); put("  y2 = "); put(y2,1);
        put("  res = "); put(res,1); new_line;
        if y2 /= 0 then
          if y2 < 0 then y2 := -y2; end if;
          if x1*y2 > x2*y1
           then x1 := x2; y1 := y2; res := k;
          end if;
        end if;
      end if;
    end loop;
    return res;
  end Largest_Angle;

  function Is_Zero ( v : Standard_Integer64_Vectors.Vector ) return boolean is
  begin
    for i in v'range loop
      if v(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  function Rank_of_Supported_Points
              ( A : Matrix; v : Standard_Integer64_Vectors.Vector )
              return natural32 is

  -- DESCRIPTON :
  --   Returns the rank of affine plane spanned by those points of A
  --   supported in the direction of v.

    s : constant Standard_Integer_Vectors.Vector
      := Standard_Lattice_Supports.Support(A,v);
    B : constant Matrix(A'range(1),s'range) 
      := Standard_Lattice_Supports.Support_Points(A,s);

  begin
    return Rank(B);
  end Rank_of_Supported_Points;

  procedure Initial_Facet_Normal
              ( A : in Matrix; rnk : out natural32;
                v : out Standard_Integer64_Vectors.Vector ) is

    k : constant integer32 := Standard_Lattice_3d_Facets.Lowest(A);
    e : constant integer32 := Standard_Lattice_3d_Facets.Initial_Edge(A,k);
    w : constant Standard_Integer64_Vectors.Vector(A'range(1))
      := Standard_Lattice_3d_Facets.Edge_Normal(A,k,e);
    i1,i2,i3 : integer32;
    f,g,h : Standard_Integer64_Vectors.Vector(A'range(1));
    ind : integer32;

  begin
   -- put("step 1, first edge inner normal :"); put(w); new_line;
   -- Write_Supported_Points(A,w);
    rnk := Rank_of_Supported_Points(A,w);
    if integer32(rnk) = A'last(1)-1 then
     -- put_line("Rank of supported points is full at edge, done!");
      v := w;
    else
     -- put("rank at supported points is "); put(rnk,1); new_line;
      Standard_Lattice_3d_Facets.Initial_Facet_Normal(A,i1,i2,i3,f);
     -- put("step 2, first 3D facet inner normal :"); put(f); new_line;
     -- Write_Supported_Points(A,f);
      rnk := Rank_of_Supported_Points(A,f);
      if integer32(rnk) = A'last(1)-1 then
       -- put_line("Rank of supported points is full at 3d facet, done!");
        v := f;
      else
       -- put("rank at supported points is "); put(rnk,1); new_line;
        for i in 1..A'last(1)-1 loop
          g := Normal(A,f);
         -- put("step "); put(i+2,1); put(", ");
         -- put("normal to previous face normal :"); put(g); new_line;
         -- if Is_Zero(g)
         --  then put_line("zero vector, bug !?"); exit;
         -- end if;
         -- put("inner product with previous face normal : "); 
         -- put(Standard_Lattice_Supports.Inner_Product(f,g),1); new_line;
          ind := Largest_Angle(A,f,g);
         -- put("index of largest angle : "); put(ind,1); new_line;
          h := Normal(A,f,ind);
         -- Write_Supported_Points(A,h);
          rnk := Rank_of_Supported_Points(A,h);
          if integer32(rnk) = A'last(1)-1 then
           -- put_line("Rank of supported points is full, done!");
            exit;
          end if;
         -- put("Continue with step "); put(i+3,1); put_line(" ...");
          f := h;
        end loop;
        v := h;
      end if;
    end if;
    rnk := Rank(A);
    if integer32(rnk) < A'last(1)-1 then
      put_line("support is not of full rank !???  BUG!"); 
      raise CONSTRAINT_ERROR;
    end if;
   -- put("Initial facet normal returns v = "); put(v); new_line;
   -- Write_Supported_Points(A,v);
  end Initial_Facet_Normal;

end Standard_Lattice_Polytopes;
