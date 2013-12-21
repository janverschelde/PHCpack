with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer64_Vectors_io;      use Standard_Integer64_Vectors_io;
with Standard_Integer64_Matrices_io;     use Standard_Integer64_Matrices_io;
with Standard64_Common_Divisors;         use Standard64_Common_Divisors;

package body Standard_Lattice_Polygons is

-- I. convert set into lexicographically decreasing matrix

  procedure Convert ( s : in List; A : out Matrix ) is

    cnt : integer32 := 0;
    t : List := s;
    v : Standard_Integer_Vectors.Link_to_Vector;

  begin
    while not Is_Null(t) loop
      v := Head_Of(t);
      cnt := cnt + 1;
      A(1,cnt) := integer64(v(1));
      A(2,cnt) := integer64(v(2));
      t := Tail_Of(t);
    end loop;
  end Convert;

  procedure Convert ( A : in Matrix; s : out List ) is

    last : List;
    v : Standard_Integer_Vectors.Vector(A'range(1));

  begin
    for j in A'range(2) loop
      for i in A'range(1) loop
        v(i) := integer32(A(i,j));
      end loop;
      Append(s,last,v);
    end loop;
  end Convert;

  function Is_Larger ( A : Matrix; i,j : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if column i of A is larger than column j of A,
  --   using the second coordinate as a tie breaker.

  begin
    if A(1,i) > A(1,j) then
      return true;
    elsif A(1,i) = A(1,j) then
      return (A(2,i) > A(2,j));
    else
      return false;
    end if;
  end Is_Larger;

  procedure Swap ( A : in out Matrix; i,j : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps columns i and j of the matrix A.

    tmp : integer64;

  begin
    tmp := A(1,i); A(1,i) := A(1,j); A(1,j) := tmp; 
    tmp := A(2,i); A(2,i) := A(2,j); A(2,j) := tmp; 
  end Swap;

  procedure Lexicographic_Decreasing_Sort ( A : in out Matrix ) is

    ind : integer32;

  begin
    for i in A'first(2)..A'last(2)-1 loop
      ind := i;
      for j in i+1..A'last(2) loop
        if not Is_Larger(A,ind,j)
         then ind := j;
        end if;
      end loop;
      if ind /= i
       then Swap(A,i,ind);
      end if;
    end loop;
  end Lexicographic_Decreasing_Sort;

-- II. compute vertex set

  function Position ( A : Matrix; i,j,k : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Determines the position of the j-th point with respect to i and k.

  -- REQUIRED : i,j,k belong to A'range(2).

  -- ON ENTRY :
  --   A        2-by-m matrix with in its columns x and y coordinates; 
  --   i        index to the column in A with coordinates of a point p;
  --   j        index to the column in A with coordinates of a point q;
  --   k        index to the column in A with coordinates of a point r.

  -- ON RETURN :
  --   -1       if the point q lies below the line spanned by p and r,
  --            or to the left of the vertical line spanned by p and r;
  --    0       if the point q lies on the line spanned by p and r;
  --   +1       if the point q lies above the line spanned by p and r,
  --            or to the right of the vertical line spanned by p and r.

    xij,xik,yij,yik : integer64;

  begin
    if A(1,j) = A(1,i) then
      if A(1,j) = A(1,k) then      -- collinear
        return 0;
      elsif A(2,j) > A(2,i) then   -- k-th point is irrelevant
        return +1;
      else
        return -1;
      end if;
    elsif A(1,i) = A(1,k) then     -- p and r span vertical line 
      if A(1,j) < A(1,i)           -- collinear case covered above
       then return -1;
       else return +1;
      end if;
    else
      xij := A(1,j) - A(1,i); xik := A(1,k) - A(1,i);
      yij := A(2,j) - A(2,i); yik := A(2,k) - A(2,i);
      if yij*xik = xij*yik then
        return 0;
      else
        if xij < 0 then yij := -yij; end if;  -- normalize slopes
        if xik < 0 then yik := -yik; end if;
        if yij*xik < xij*yik
         then return -1;
         else return +1;
        end if;
      end if;
    end if;
  end Position;

  function Vertices ( v : Standard_Integer_Vectors.Vector ) return natural32 is

  -- DESCRIPTION :
  --   Returns the number of vertices of the convex hull as indicated by
  --   the number of nonzero entries in v.

    res : natural32 := 0;

  begin
    for i in v'range loop
      if v(i) /= 0
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Vertices;

  function Hull ( A : Standard_Integer64_Matrices.Matrix;
                  v : Standard_Integer_Vectors.Vector )
                return Standard_Integer64_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the ordered matrix of points that span the convex hull,
  --   starting at the rightmost highest point running counterclockwise.

  -- ON ENTRY :
  --   A        columns of A contain x and y coordinates of points;
  --   v        v(k) = +1 if k-th point is on upper hull,
  --                 = -1 if k-th point is on lower hull,
  --                 = 0 if k-th point is an interior point.

    k : constant natural32 := Vertices(v);
    res : Matrix(1..2,1..integer32(k));
    ind : integer32 := 0;

  begin
    for i in v'range loop       -- first do upper hull
      if v(i) = +1 then
        ind := ind + 1;
        res(1,ind) := A(1,i);
        res(2,ind) := A(2,i);
      end if;
    end loop;
    for i in reverse v'range loop       -- then do lower hull
      if v(i) = -1 then
        ind := ind + 1;
        res(1,ind) := A(1,i);
        res(2,ind) := A(2,i);
      end if;
    end loop;
    return res;
  end Hull;

  function Next_Pos ( v : Standard_Integer_Vectors.Vector;
                      k : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the first element of v with index > k and positive entry.
  --   If there are no next posivite entries, v'last+1 is returned.

    res : integer32 := k+1;

  begin
    while res <= v'last loop
      if v(res) > 0
       then return res;
       else res := res + 1;
      end if;
    end loop;
    return res;
  end Next_Pos;

  function Next_Neg ( v : Standard_Integer_Vectors.Vector;
                      k : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the first element of v with index > k and negative entry.
  --   If there are no next negative entries, v'last is returned.

    res : integer32 := k+1;

  begin
    while res <= v'last loop
      if v(res) < 0
       then return res;
       else res := res + 1;
      end if;
    end loop;
    return v'last;
  end Next_Neg;

  function Prev_Pos ( v : Standard_Integer_Vectors.Vector;
                      k : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the first element of v with index < k and positive entry.
  --   If there are no next positive entries, then v'first-1 is returned.

    res : integer32 := k-1;

  begin
    while res >= v'first loop
      if v(res) > 0
       then return res;
       else res := res - 1;
      end if;
    end loop;
    return res;
  end Prev_Pos;

  function Prev_Neg ( v : Standard_Integer_Vectors.Vector;
                      k : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the first element of v with index < k and negative entry.
  --   If there are no next negative entries, then v'first is returned.

    res : integer32 := k-1;

  begin
    while res >= v'first loop
      if v(res) < 0
       then return res;
       else res := res - 1;
      end if;
    end loop;
    return v'first;
  end Prev_Neg;

  function Convex_Hull_2D ( A : Matrix ) return Matrix is

    p,q,r : integer32;
    v : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => +1);

  begin
    if A'first(2)+2 > A'last(2) then
      return A;
    else
      p := A'first(2); q := p+1; r := q+1;      -- compute upper hull
      loop
        v(q) := Position(A,p,q,r);
        if v(q) > 0 then
          p := q; q := r; r := Next_Pos(v,r);
        else
          if p = A'first(2) then
            q := r; r := Next_Pos(v,r);
          else
            q := p; p := Prev_Pos(v,p);
          end if;
        end if;
        exit when (r > A'last(2));
      end loop;
      p := A'last(2); q := Prev_Neg(v,p);       -- compute lower hull
      if q > v'first then
        r := Prev_Neg(v,q);
        if r >= A'first(2) then
          loop
            v(q) := Position(A,p,q,r);
            if v(q) < 0 then
              p := q; q := r; r := Prev_Neg(v,r);
            else
              v(q) := 0;              -- do not count q on upper hull
              if p = A'last(2) then
                q := r; r := Prev_Neg(v,r);
              else
                q := p; p := Next_Neg(v,p);
              end if;
            end if;
            exit when (v(q) >= 0);
          end loop;
        end if;
      end if;
      return Hull(A,v);
    end if;
  end Convex_Hull_2D;

  procedure Normal ( V : in Matrix; i : in integer32; N : in out Matrix ) is

  -- DESCRIPTION :
  --   Updates the i-th column of N with the coordinates of a vector
  --   perpendicular to the edge spanned by the i-th and (i+1)-th column
  --   of V, with (i+1) taken with respect to V'last(2).
  --   The components of the normals are relative prime.

  -- REQUIRED : N'range = V'range and i belongs to V'range(2).

    x,y,d : integer64;

  begin
    if i < V'last(2) then
      x := V(1,i+1) - V(1,i);
      y := V(2,i+1) - V(2,i);
    else
      x := V(1,1) - V(1,i);
      y := V(2,1) - V(2,i);
    end if;
    if x = 0 then
      N(1,i) := 1;
      N(2,i) := 0;
    elsif y = 0 then
      N(1,i) := 0;
      N(2,i) := 1;
    else
      d := gcd(x,y);
      N(1,i) := (-y)/d;
      N(2,i) := x/d;
    end if;
  end Normal;

  procedure Orient ( V : in Matrix; i : in integer32; N : in out Matrix ) is

  -- DESCRIPTION :
  --   Updates the orientation of the normal in the i-th column of N
  --   so that the normal is a proper inner normal.

  -- REQUIRED : V'last(2) > V'first(2) + 1.

    p : constant integer64 := V(1,i)*N(1,i) + V(2,i)*N(2,i);
    q : integer64;

  begin
    if i-1 >= V'first(1) then
      q := V(1,i-1)*N(1,i) + V(2,i-1)*N(2,i);
    else
      q := V(1,V'last(2))*N(1,i) + V(2,V'last(2))*N(2,i);
    end if;
    if q < p then
      N(1,i) := -N(1,i);
      N(2,i) := -N(2,i);
    end if;
  end Orient;

  function Inner_Normals ( V : Matrix ) return Matrix is

    res : Matrix(V'range(1),V'range(2));

  begin
    if V'last(2) = V'first(2) then
      res(1,res'first(2)) := 0;
      res(2,res'first(2)) := 0;
    elsif V'last(2) = V'first(2) + 1 then
      Normal(V,V'first(2),res);
      res(1,V'first(2)+1) := -res(1,V'first(2));
      res(2,V'first(2)+1) := -res(2,V'first(2));
    else
      for i in V'range(2) loop
        Normal(V,i,res);
        Orient(V,i,res);
      end loop;
    end if;
    return res;
  end Inner_Normals;

  function Inner_Normals ( s : List ) return List is

    res : List;
    A : Matrix(1..2,1..integer32(Length_Of(s)));

  begin
    Convert(s,A);
    Lexicographic_Decreasing_Sort(A);
    declare
      V : constant Matrix := Convex_Hull_2D(A);
      N : constant Matrix := Inner_Normals(V);
    begin
      Convert(N,res);
    end;
    return res;
  end Inner_Normals;

-- III. hyperplane representations

  function Rank ( A : Matrix; i : integer32; v : Vector ) return integer64 is
  begin
    return A(1,i)*v(1) + A(2,i)*v(2);
  end Rank;

  function Rank ( A : Matrix; v : Vector ) return Vector is

    res : Vector(A'range(2));

  begin
    for i in A'range(2) loop
      res(i) := Rank(A,i,v);
    end loop;
    return res;
  end Rank;

  function Minimum ( A : Matrix; v : Vector ) return integer64 is

     res : integer64 := Rank(A,A'first(2),v);
     val : integer64;

  begin
    for i in A'first(2)+1..A'last(2) loop
      val := Rank(A,i,v);
      if val < res
       then res := val;
      end if;
    end loop;
    return res;
  end Minimum;

  function Minima ( A,N : Matrix ) return Vector is
   
    res : Vector(N'range(2));
    v : Vector(1..2);

  begin
    for i in N'range(2) loop
      v(1) := N(1,i);
      v(2) := N(2,i);
      res(i) := Minimum(A,v);
    end loop;
    return res;
  end Minima;

  function Number_of_Minima ( v : Vector ) return natural32 is

    res : natural32 := 1;
    min : integer64 := v(v'first);

  begin
    for i in v'first+1..v'last loop
      if v(i) < min then
        min := v(i);
        res := 1;
      elsif v(i) = min then
        res := res + 1;
      end if;
    end loop;
    return res;
  end Number_of_Minima;

  function Rank ( A,N : Matrix ) return Matrix is

    res : Matrix(N'range(2),A'range(2));
    v : Vector(1..2);

  begin
    for i in N'range(2) loop
      v(1) := N(1,i);
      v(2) := N(2,i);
      for j in A'range(2) loop
        res(i,j) := Rank(A,j,v);
      end loop;
    end loop;
    return res;
  end Rank;

-- IV. sanity check

  procedure Check ( A,V,N : in Matrix; output : in boolean;
                    bug : out boolean ) is

    r : Vector(1..2);
    m : natural32;
    nA,nV : Vector(N'range(2));

  begin
    if output then
      put_line("The vertex set of "); put(A); put_line("is "); put(V);
      put_line("with inner normals :"); put(N);
      put_line("and ranking matrix :"); put(Rank(V,N));
    end if;
    bug := false;
    for i in V'range(2) loop
      r(1) := N(1,i); r(2) := N(2,i);
      m := Number_of_Minima(Rank(V,r));
      if output then
        put("normal"); put(r); put(" has ");
        put(m,1); put(" vertex minima");
      end if;
      if m < 2 then
        bug := true;
        if output 
         then put_line("  BUG!");
        end if;
      else
        m := Number_of_Minima(Rank(A,r));
        if output
         then put(" and "); put(m,1); put(" support minima ");
        end if;
        if m < 2
         then bug := true;
        end if;
        if output then
          if m < 2
           then put_line("  BUG!");
           else put_line("  okay");
          end if;
        end if;
      end if;
    end loop;
    nA := Minima(A,N); nV := Minima(V,N);
    if not Equal(nA,nV)
     then bug := true;
    end if;
    if output then
      put("Minima with vertices : "); put(nA); new_line;
      put("Minima with supports : "); put(nV);
      if not Equal(nA,nV)
       then put_line("  BUG!");
       else put_line("  okay");
      end if;
    end if;
  end Check;

end Standard_Lattice_Polygons;
