with Multprec_Integer_Matrices_io; use Multprec_Integer_Matrices_io;

with text_io;                            use text_io;
--with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Lattice_Supports;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
--with Multprec_Integer_Numbers_io;        use Multprec_Integer_Numbers_io;
--with Multprec_Integer_Vectors_io;        use Multprec_Integer_Vectors_io;
with Multprec_Integer_Linear_Solvers;    use Multprec_Integer_Linear_Solvers;
with Multprec_Lattice_Supports;
with Multprec_Lattice_3d_Facets;
--with Multprec_Lattice_3d_Facets_io;      use Multprec_Lattice_3d_Facets_io;
with Multprec_Integer_Orthogonals;

package body Multprec_Lattice_Polytopes is

  procedure Normalize ( v : in out Multprec_Integer_Vectors.Vector ) is

    zero : constant integer32 := 0;

  begin
    for i in v'range loop
      if Empty(v(i)) then
        v(i) := Multprec_Integer_Numbers.Create(zero);
      elsif Negative(v(i)) then
        if Equal(v(i),0) then
          Clear(v(i));
          v(i) := Multprec_Integer_Numbers.Create(zero);
        end if;
      end if;
    end loop;
  end Normalize;

  procedure Normalize ( A : in out Matrix ) is

    zero : constant integer32 := 0;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        if Empty(A(i,j)) then
          A(i,j) := Multprec_Integer_Numbers.Create(zero);
        elsif Negative(A(i,j)) then
          if Equal(A(i,j),0) then
            Clear(A(i,j));
            A(i,j) := Multprec_Integer_Numbers.Create(zero);
          end if;
        end if;
      end loop;
    end loop;
  end Normalize;

-- INITIAL FACET :

  function Shift ( A : Matrix; k : integer32 ) return Matrix is

  -- DESCRIPTION :
  --   The matrix on return has the k-th column removed and the coordinates
  --   of all other columns have been subtracted by the k-th column.

    res : Matrix(A'range(1),A'first(2)..A'last(2)-1);

  begin
   -- put_line("The matrix before the shift : "); put(A);
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
    Normalize(res);
   -- put_line("The matrix after the shift : "); put(res);
    return res;
  end Shift;

  function Rank_of_Upper ( A : Matrix ) return natural32 is

  -- DESCRIPTION :
  --   Returns the rank of an upper triangular matrix.

    res : natural32 := 0;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        if not Equal(A(i,j),0)
         then res := res + 1; exit;
        end if;
      end loop;
    end loop;
    return res;
 -- exception 
 --   when others => put_line("exception in Rank_of_Upper"); raise;
  end Rank_of_Upper;

  function Rank ( A : Matrix ) return natural32 is

    res : natural32;
    k : constant integer32 := Multprec_Lattice_3d_Facets.Lowest(A);
    V : Matrix(A'range(1),A'first(2)..A'last(2)-1) := Shift(A,k);

  begin
   -- put_line("in rank with V = "); put(V);
    Upper_Triangulate(V);
   -- put_line("V after Upper Triangulate :"); put(V);
    res := Rank_of_Upper(V);
    return res;
 -- exception
 --   when others => put_line("exception in Rank, V = "); put(V); raise;
  end Rank;

  function Normal ( A : Matrix; v : Multprec_Integer_Vectors.Vector )
                  return Multprec_Integer_Vectors.Vector is

    res : Multprec_Integer_Vectors.Vector(v'range);
    s : constant Standard_Integer_Vectors.Vector
      := Multprec_Lattice_Supports.Support(A,v);
    B : constant Matrix(A'range(1),s'range)
      := Multprec_Lattice_Supports.Support_Points(A,s);
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
   -- put_line("calling complement ...");
    res := Multprec_Integer_Orthogonals.Complement(C);
    return res;
 -- exception
 --   when others => put("exception occurred in Normal with v = ");
 --                  put(v); new_line; raise;
  end Normal;

  procedure Inner ( A : in Matrix; k : in integer32;
                    v : in out Multprec_Integer_Vectors.Vector ) is

    m : Integer_Number := Multprec_Lattice_Supports.Minimum(A,v);
    p : Integer_Number := Multprec_Lattice_Supports.Inner_Product(A,k,v);

  begin
    if p > m
     then Multprec_Integer_Vectors.Min(v);
    end if;
    Clear(m); Clear(p);
  end Inner;

  function Normal ( A : Matrix; v : Multprec_Integer_Vectors.Vector;
                    k : integer32 ) return Multprec_Integer_Vectors.Vector is

    res : Multprec_Integer_Vectors.Vector(v'range);
    s : constant Standard_Integer_Vectors.Vector
      := Multprec_Lattice_Supports.Support(A,v);
    B : constant Matrix(A'range(1),s'range)
      := Multprec_Lattice_Supports.Support_Points(A,s);
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
    res := Multprec_Integer_Orthogonals.Complement(C);
    Inner(A,k,res);
    return res;
  end Normal;

  function Largest_Angle
              ( A : Matrix; f,g : Multprec_Integer_Vectors.Vector )
              return integer32 is

    res,ind : integer32 := 0;
    s : constant Standard_Integer_Vectors.Vector
      := Multprec_Lattice_Supports.Support(A,f);
    w : Multprec_Integer_Vectors.Vector(A'range(1));
    x1,y1,x2,y2,p,q : Integer_Number;

  begin
    for k in A'range(2) loop
      if Standard_Lattice_Supports.Member(s,k) < s'first then
        w := Multprec_Lattice_3d_Facets.Shift(A,s(s'first),k);
        x1 := Multprec_Lattice_Supports.Inner_Product(w,f);
        y1 := Multprec_Lattice_Supports.Inner_Product(w,g);
        if y1 < 0 then Multprec_Integer_Numbers.Min(y1); end if;
        ind := k; res := k;
      end if;
      exit when (ind > 0);
    end loop;
    for k in ind+1..A'last(2) loop
      if Standard_Lattice_Supports.Member(s,k) < s'first then
        w := Multprec_Lattice_3d_Facets.Shift(A,s(s'first),k);
        x2 := Multprec_Lattice_Supports.Inner_Product(w,f);
        y2 := Multprec_Lattice_Supports.Inner_Product(w,g);
       -- put("  x1 = "); put(x1,1); put("  y1 = "); put(y1,1);
       -- put("  x2 = "); put(x2,1); put("  y2 = "); put(y2,1);
       -- put("  res = "); put(res,1); new_line;
        if y2 < 0 then Multprec_Integer_Numbers.Min(y2); end if;
        p := x1*y2; q := x2*y1;
        if p > q
         then x1 := x2; y1 := y2; res := k;
        end if;
        Clear(p); Clear(q);
      end if;
    end loop;
    return res;
  end Largest_Angle;

  function Is_Zero ( v : Multprec_Integer_Vectors.Vector ) return boolean is
  begin
    for i in v'range loop
      if not Equal(v(i),0)
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  function Rank_of_Supported_Points
              ( A : Matrix; v : Multprec_Integer_Vectors.Vector )
              return natural32 is

  -- DESCRIPTON :
  --   Returns the rank of affine plane spanned by those points of A
  --   supported in the direction of v.

    res : natural32;
    s : constant Standard_Integer_Vectors.Vector
      := Multprec_Lattice_Supports.Support(A,v);
    B : Matrix(A'range(1),s'range) 
      := Multprec_Lattice_Supports.Support_Points(A,s);

  begin
    res := Rank(B);
    Multprec_Integer_Matrices.Clear(B);
    return res;
 -- exception
 --   when others => put_line("exception in Rank_of_Supported_Points");
 --                  raise;
  end Rank_of_Supported_Points;

  procedure Initial_Facet_Normal
              ( A : in Matrix; rnk : out natural32;
                v : out Multprec_Integer_Vectors.Vector ) is

    k : constant integer32 := Multprec_Lattice_3d_Facets.Lowest(A);
    e : constant integer32 := Multprec_Lattice_3d_Facets.Initial_Edge(A,k);
    w : constant Multprec_Integer_Vectors.Vector(A'range(1))
      := Multprec_Lattice_3d_Facets.Edge_Normal(A,k,e);
    i1,i2,i3 : integer32;
    f,g,h : Multprec_Integer_Vectors.Vector(A'range(1));
    ind : integer32;

  begin
   -- put("step 1, first edge inner normal :"); put(w); new_line;
   -- Write_Supported_Points(A,w);
    rnk := Rank_of_Supported_Points(A,w);
    if integer32(rnk) = A'last(1)-1 then
     -- put_line("Rank of supported points at edge is full, done!");
      v := w;
    else
     -- put("rank at supported points is "); put(rnk,1); new_line;
      Multprec_Lattice_3d_Facets.Initial_Facet_Normal(A,i1,i2,i3,f);
     -- put("step 2, first 3D facet inner normal :"); put(f); new_line;
     -- Write_Supported_Points(A,f);
      rnk := Rank_of_Supported_Points(A,f);
      if integer32(rnk) = A'last(1)-1 then
       -- put_line("Rank of supported points at 3d facet is full, done!");
        v := f;
      else
       -- put("rank at supported points is "); put(rnk,1); new_line;
        for i in A'range loop
         -- put_line("Calling Normal ...");
          g := Normal(A,f);
         -- put("step "); put(i+2,1); put(", ");
         -- put("normal to previous face normal :"); put(g); new_line;
         -- if Is_Zero(g)
         --  then put_line("zero vector, bug !?"); exit;
         -- end if;
         -- put("inner product with previous face normal : "); 
         -- put(Multprec_Lattice_Supports.Inner_Product(f,g),1); new_line;
          ind := Largest_Angle(A,f,g);
         -- put("index of largest angle : "); put(ind,1); new_line;
          h := Normal(A,f,ind);
         -- put("New normal "); put(h); new_line;
         -- Write_Supported_Points(A,h);
          rnk := Rank_of_Supported_Points(A,h);
          if integer32(rnk) = A'last(1)-1 then
           -- put_line("Rank of supported points is full, done!");
            exit;
          end if;
         -- put("Continue with next step "); put(i,1); put_line(" ...");
          Multprec_Integer_Vectors.Copy(h,f);
        end loop;
        Multprec_Integer_Vectors.Copy(h,v);
      end if;
    end if;
    rnk := Rank(A);
 -- exception
 --   when others => put_line("exception occurred in Initial_Facet_Normal");
 --                  raise;
  end Initial_Facet_Normal;

end Multprec_Lattice_Polytopes;
