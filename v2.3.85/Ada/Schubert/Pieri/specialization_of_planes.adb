with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;

package body Specialization_of_Planes is

  function Random_Upper_Triangular
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(1..n,1..n);
    
  begin
    for j in 1..n loop                        -- assign values to jth column
      for i in 1..n-j loop
        res(i,j) := Random1;                  -- randoms above anti-diagonal
      end loop;
      res(n-j+1,j) := Create(1.0);            -- 1 = anti-diagonal element
      for i in n-j+2..n loop
        res(i,j) := Create(0.0);              -- zeros under anti-diagonal
      end loop;
    end loop;
    return res;
  end Random_Upper_Triangular;

  function Random_Lower_Triangular
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(1..n,1..n);

  begin
    for j in 1..n loop                         -- assign values to jth column
      for i in 1..(j-1) loop
        res(i,j) := Create(0.0);               -- zeros above diagonal
      end loop;
      res(j,j) := Create(1.0);                 -- 1 = diagonal element
      for i in (j+1)..n loop
        res(i,j) := Random1;                   -- randoms under diagonal
      end loop;
    end loop;
    return res;
  end Random_Lower_Triangular;

  function U_Matrix ( F : Standard_Complex_Matrices.Matrix; b : Bracket )
                    return Standard_Complex_Matrices.Matrix is

    m : constant integer32 := F'length(1) - b'length;
    res : Standard_Complex_Matrices.Matrix(F'range(1),1..m);
    rvf : Standard_Complex_Matrices.Matrix(F'range(1),F'range(2)) := F;
    rng : constant integer32 := F'length(2) - integer32(b(b'last));
    tmp : Complex_Number;
    ind : integer32 := 1;
    cnt : integer32 := 0;

  begin
    for j in 1..(rng/2) loop                   -- reverse last columns
      for i in F'range(1) loop
        tmp := rvf(i,F'last(2)-j+1);
        rvf(i,F'last(2)-j+1) := rvf(i,F'last(2)-rng+j);
        rvf(i,F'last(2)-rng+j) := tmp;
      end loop;
    end loop;
    for j in F'range(2) loop                   -- remove columns indexed by b
      if ((ind <= b'last) and then (natural32(j) = b(ind))) then
        ind := ind+1;
      else
        cnt := cnt+1;
        for i in F'range(1) loop
          res(i,cnt) := rvf(i,j);
        end loop;
      end if;
    end loop;
    return res;
  end U_Matrix;

  function Special_Plane ( m : integer32; b : Bracket )
                         return Standard_Complex_Matrices.Matrix is

    p : constant integer32 := b'length;
    n : constant integer32 := m+p;
    res : Standard_Complex_Matrices.Matrix(1..n,1..m);
    row,col : integer32;

  begin
    row := 1; col := 0;
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Create(0.0);
      end loop;
      if ((row <= p) and then (b(row) = natural32(i)))
       then row := row + 1;
       else col := col + 1; res(i,col) := Create(1.0);
      end if;
    end loop;
    return res;
  end Special_Plane;

  function Special_Bottom_Plane ( m : integer32; b : Bracket )
                                return Standard_Complex_Matrices.Matrix is

    p : constant integer32 := b'length;
    n : constant integer32 := m+p;
    res : Standard_Complex_Matrices.Matrix(1..n,1..m);
    row,col : integer32;

  begin
    row := 1; col := 0;
    for i in res'range(1) loop
      if ((row <= p) and then (b(row) = natural32(i))) then
        row := row + 1;
      else
        col := col + 1;
        for k in 1..i-1 loop             -- randoms above the diagonal
          res(k,col) := Random1;
        end loop;
        res(i,col) := Create(1.0);
        for k in i+1..n loop             -- zeros below the diagonal
          res(k,col) := Create(0.0);
        end loop;
      end if;
    end loop;
    return res;
  end Special_Bottom_Plane;

  function Special_Top_Plane ( m : integer32; b : Bracket )
                             return Standard_Complex_Matrices.Matrix is

    p : constant integer32 := b'length;
    n : constant integer32 := m+p;
    res : Standard_Complex_Matrices.Matrix(1..n,1..m);
    row,col : integer32;

  begin
    row := 1; col := 0;
    for i in res'range(1) loop
      if ((row <= p) and then (b(row) = natural32(i))) then
        row := row + 1;
      else
        col := col + 1;
        for k in 1..i-1 loop              -- zeros above the diagonal
          res(k,col) := Create(0.0);
        end loop;
        res(i,col) := Create(1.0);
        for k in i+1..n loop              -- randoms below the diagonal
          res(k,col) := Random1;
        end loop;
      end if;
    end loop;
    return res;
  end Special_Top_Plane;

  function Special_Plane
              ( n,m,k : integer32; b : Bracket;
                special : in Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(1..n,1..m+1-k);
    ran : Complex_Number;

  begin
    for i in res'range(1) loop                   -- initialize
      for j in res'range(2) loop
        res(i,j) := Create(0.0);
      end loop;
    end loop;
    for j in res'range(2) loop                   -- build j-th column
      for k in b'range loop
        ran := Random1;                          -- random for column k of mat
        for i in special'range(1) loop
          res(i,j) := res(i,j) + ran*special(i,integer32(b(k)));
        end loop;
      end loop;
    end loop;
    return res;
  end Special_Plane;

--  function Special_Bottom_Plane ( n,m,k : natural; b : Bracket )
--                                return Standard_Complex_Matrices.Matrix is
--
--    mat : Standard_Complex_Matrices.Matrix(1..n,b'range);
--
--  begin
--    for j in b'range loop               -- j-th column of matrix
--      for i in 1..b(j) loop
--        mat(i,j) := Random1;            -- random numbers above row b(j)
--      end loop;
--      for i in b(j)+1..n loop
--        mat(i,j) := Create(0.0);        -- zeros below row b(j)
--      end loop;
--    end loop;
--    return Special_Plane(n,m,k,b,mat);
--  end Special_Bottom_Plane;

--  function Special_Top_Plane ( n,m,k : natural; b : Bracket )
--                             return Standard_Complex_Matrices.Matrix is
--
--    mat : Standard_Complex_Matrices.Matrix(1..n,b'range);
--
--  begin
--    for j in b'range loop               -- j-th column of matrix
--      for i in 1..b(j)-1 loop
--        mat(i,j) := Create(0.0);        -- zeros below row b(j)
--      end loop;
--      for i in b(j)..n loop
--        mat(i,j) := Random1;            -- random numbers below row b(j)
--      end loop;
--    end loop;
--    return Special_Plane(n,m,k,b,mat);
--  end Special_Top_Plane;

  function Homotopy ( dim : integer32; start,target : Complex_Number )
                    return Poly is

  -- DESCRIPTION :
  --   Returns the polynomial start*(1-t) + target*t, where t is the
  --   is the last variable of index dim.
  --   This procedure is an auxiliary to building the moving U-matrices.

    res : Poly;
    t : Term;
    tdg : Degrees := new Standard_Natural_Vectors.Vector'(1..dim => 0);

  begin
    t.cf := start;
    t.dg := tdg;
    res := Create(t);               -- res = start
    tdg(tdg'last) := 1;             -- introduce t
    t.dg := tdg;
    Sub(res,t);                     -- res = (1-t)*start
    t.cf := target;
    Add(res,t);                     -- res = (1-t)*start + t*target
    Clear(tdg);                 
    return res;
  end Homotopy;

  function Constant_Poly ( dim : integer32; c : Complex_Number ) return Poly is

  -- DESCRIPTION :
  --   Returns the constant c represented as polynomial with as many
  --   variables as the given dimension.

    res : Poly;
    t : Term;
    tdg : constant Degrees 
        := new Standard_Natural_Vectors.Vector'(1..dim => 0);

  begin
    t.cf := c;
    t.dg := tdg;
    res := Create(t);
    return res;
  end Constant_Poly;

  function Moving_U_Matrix
             ( n : integer32; U,L : Standard_Complex_Matrices.Matrix ) 
             return Standard_Complex_Poly_Matrices.Matrix is

    res : Standard_Complex_Poly_Matrices.Matrix(L'range(1),L'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Homotopy(n,U(i,j),L(i,j));
      end loop;
    end loop;
    return res;
  end Moving_U_Matrix;

  function Moving_U_Matrix
             ( U : Standard_Complex_Matrices.Matrix;
               i,r : integer32; b : bracket ) 
             return Standard_Complex_Poly_Matrices.Matrix is

    p : constant integer32 := b'last;
    m : constant integer32 := U'length(1) - p;
    dim : constant integer32 := (m+p)*p+1;
    res : Standard_Complex_Poly_Matrices.Matrix(U'range(1),1..m+1-r);

  begin
    for j in res'range(2) loop
      if j+i-1 < integer32(b(b'last)) - b'last then
        for k in res'range(1) loop
          res(k,j) := Homotopy(dim,U(k,j+i),U(k,j+i-1));
        end loop;
      elsif j+i-1 = integer32(b(b'last)) - b'last then
        for k in res'range(1) loop
          res(k,j) := Homotopy(dim,U(k,m+1+i-r),U(k,j+i-1));
        end loop;
      else
        for k in res'range(1) loop
          res(k,j) := Constant_Poly(dim,U(k,j+i-1));
        end loop;
      end if;
    end loop;
    return res;
  end Moving_U_Matrix;

  function Slice ( M : Standard_Complex_Poly_Matrices.Matrix;
                   ind : integer32 )
                 return Standard_Complex_Poly_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the columns of M up to the given index.

    res : Standard_Complex_Poly_Matrices.Matrix(M'range(1),M'first(2)..ind);

  begin
    for j in res'range(2) loop
      for i in res'range(1) loop
        Copy(M(i,j),res(i,j));
      end loop;
    end loop;
    return res;
  end Slice;

  function Lower_Section
             ( M : Standard_Complex_Poly_Matrices.Matrix;
               row : integer32 )
             return Standard_Complex_Poly_Matrices.Matrix is

    res : Standard_Complex_Poly_Matrices.Matrix(M'range(1),M'range(2));
    cnt : integer32 := M'first(2)-1;
    add : boolean;

  begin
    for j in M'range(2) loop
      for i in (row+1)..M'last(1) loop
        add := (M(i,j) = Null_Poly);
        exit when not add;
      end loop;
      if add
       then cnt := cnt+1;
            for i in M'range(1) loop
              res(i,cnt) := M(i,j);
            end loop;
      end if;
    end loop;
    return Slice(res,cnt);
  end Lower_Section;

  function Upper_Section
             ( M : Standard_Complex_Poly_Matrices.Matrix;
               row : integer32 )
             return Standard_Complex_Poly_Matrices.Matrix is

    res : Standard_Complex_Poly_Matrices.Matrix(M'range(1),M'range(2));
    cnt : integer32 := M'first(2)-1;
    add : boolean;

  begin
    for j in M'range(2) loop
      for i in M'first(1)..(row-1) loop
        add := (M(i,j) = Null_Poly);
        exit when not add;
      end loop;
      if add
       then cnt := cnt+1;
            for i in M'range(1) loop
              res(i,cnt) := M(i,j);
            end loop;
      end if;
    end loop;
    return Slice(res,cnt);
  end Upper_Section;

end Specialization_of_Planes;
