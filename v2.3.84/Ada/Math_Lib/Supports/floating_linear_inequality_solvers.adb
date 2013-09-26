package body Floating_Linear_Inequality_Solvers is

-- AUXILIARIES :

  function Is_In ( v : Standard_Integer_Vectors.Vector;
                   i : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the number i appears in one of the elements of v.

  begin
    for j in v'range loop
      if v(j) = i
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Is_All_In ( v : Standard_Integer_Vectors.Vector; i : integer32 )
                     return boolean is

  -- DESCRIPTION :
  --   Returns true if v(k) = i, for k in v'range, false otherwise.

  begin
    for k in v'range loop
      if v(k) /= i
       then return false;
      end if;
    end loop;
    return true;
  end Is_All_In;

  procedure Copy ( m1 : in Matrix; m2 : out Matrix ) is

  -- REQUIRED : ranges must be equal.

  -- DESCRIPTION :
  --   Copies the first matrix m1 onto the second matrix m2.
  --   The statement `m2 := m1' does not produce the same result when
  --   compiled with the VADS compiler on IBM RS/6000.

  begin
    for i in m1'range(1) loop
      for j in m1'range(2) loop
        m2(i,j) := m1(i,j);
      end loop;
    end loop;
  end Copy;

  function Inner_Product ( m : Matrix; i,j : integer32 ) return double_float is

  -- DESCRIPTION :
  --   Returns the inner product of the normals to the ith and jth hyperplane.

    res : double_float := 0.0;

  begin
    for k in m'first(1)..m'last(1)-1 loop
      res := res + m(k,i)*m(k,j);
    end loop;
    return res;
  end Inner_Product;

-- AUXILIARIES FOR RECURSION ON THE DIMENSION :

  function Pivot ( m : Matrix; i : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the index of the coefficient with largest absolute value
  --   of the ith inequality.

    res : integer32 := m'first(1);
    max : double_float := abs(m(res,i));
    tmpabs : double_float;

  begin
    for j in m'first(1)+1..m'last(1)-1 loop
      tmpabs := abs(m(j,i));
      if tmpabs > max
       then max := tmpabs; res := j;
      end if;
    end loop;
    return res;
  end Pivot;

  procedure Eliminate ( m1 : in Matrix; i,j,piv : in integer32;
                        tol : in double_float; m2 : out Matrix ) is

  -- DESCRIPTION :
  --   Performs the elimination on the jth column of the original matrix m1
  --   of inequalities to the reduced matrix m2.

  -- REQUIRED : m(piv,i) /= 0 and dimensions of m2 must match.

    fac : double_float;

  begin
    if abs(m1(piv,j)) < tol then
      for k in m1'range(1) loop
        if k < piv then
          m2(k,j) := m1(k,j);
        elsif k > piv then
          m2(k-1,j) := m1(k,j);
        end if;
      end loop;
    else
      fac := m1(piv,j)/m1(piv,i);
      for k in m1'range(1) loop
        if k < piv then
          m2(k,j) := m1(k,j) - fac*m1(k,i);
        elsif k > piv then
          m2(k-1,j) := m1(k,j) - fac*m1(k,i);
        end if;
      end loop;
    end if;
  end Eliminate;

  function Eliminate ( m : Matrix; mlast2,i,piv : integer32;
                       tol : in double_float ) return Matrix is

  -- DESCRIPTION :
  --   Eliminates one unknown from the system of inequalities.

  -- REQUIRED : m(piv,i) /= 0 and m'last(1) > 2.

    res : Matrix(m'first(1)..m'last(1)-1,m'first(2)..mlast2);

  begin
    for j in res'range(2) loop
      Eliminate(m,i,j,piv,tol,res);
    end loop;
    return res;
  end Eliminate;

  function Eliminate ( m : Matrix; mlast2,i : integer32; tol : double_float )
                     return Matrix is

  -- DESCRIPTION :
  --   Computes the pivot for the ith inequality and eliminates one unknown.
  --   The matrix has as columns m'first(2)..mlast2.

  -- REQUIRED : Inner_Product(m,i,i) > 0 and m'last(1) > 2.

    piv : constant integer32 := Pivot(m,i);

  begin
    return Eliminate(m,mlast2,i,piv,tol);
  end Eliminate;

  function Back_Substitute ( m : Matrix; i,piv : integer32; x : Vector )
                           return Vector is

  -- DESCRIPTION :
  --   Returns the back substitution of the vector x which lies on the ith
  --   hyperplane of original inequality system before elimination.

  -- REQUIRED : abs(m(piv,i)) > tolerance.

    res : Vector(x'first..x'last+1);

  begin
    res(res'first..piv-1) := x(x'first..piv-1);
    res(piv+1..res'last) := x(piv..x'last);
    res(piv) := m(m'last(1),i);
    for k in m'first(1)..m'last(1)-1 loop
      if k < piv 
       then res(piv) := res(piv) - m(k,i)*x(k);
       elsif k > piv
           then res(piv) := res(piv) - m(k,i)*x(k-1);
      end if;
    end loop;
    res(piv) := res(piv)/m(piv,i);
    return res;
  end Back_Substitute;

-- SELECTORS :

  function Evaluate ( m : Matrix; i : integer32; x : Vector )
                    return double_float is

  -- DESCRIPTION :
  --   Evaluates the vector x in the ith inequality.

    res : double_float := 0.0;

  begin
    for j in x'range loop
      res := res + m(j,i)*x(j);
    end loop;
    return res;
  end Evaluate;

  function Satisfies ( m : Matrix; i : integer32; x : Vector;
                       tol : double_float ) return boolean is
  begin
    return (Evaluate(m,i,x) >= m(m'last(1),i) - tol);
  end Satisfies;

  function Satisfies ( m : Matrix; x : Vector; tol : double_float )
                     return boolean is
  begin
    for i in m'range(2) loop
      if not Satisfies(m,i,x,tol)
       then return false;
      end if;
    end loop;
    return true;
  end Satisfies;

  function Satisfies ( m : Matrix; fst,lst : integer32;
                       x : Vector; tol : double_float ) return boolean is
  begin
    for i in fst..lst loop
      if not Satisfies(m,i,x,tol)
       then return false;
      end if;
    end loop;
    return true;
  end Satisfies;

  function First_Violated ( m : Matrix; x : Vector; tol : double_float )
                          return integer32 is
  begin
    for i in m'range(2) loop
      if not Satisfies(m,i,x,tol)
       then return i;
      end if;
    end loop;
    return (m'last(2)+1);
  end First_Violated;

  function First_Violated ( m : Matrix; fst,lst : integer32;
                            x : Vector; tol : double_float )
                          return integer32 is
  begin
    for i in fst..lst loop
      if not Satisfies(m,i,x,tol)
       then return i;
      end if;
    end loop;
    return (lst+1);
  end First_Violated;

-- INCONSISTENCY CHECKS :

  function Inconsistent ( m : Matrix; i : integer32; tol : double_float )
                        return boolean is
  begin
    for j in m'first(1)..m'last(1)-1 loop
      if abs(m(j,i)) > tol
       then return false;
      end if;
    end loop;
    return (m(m'last(1),i) > tol);
  end Inconsistent;

  function Inconsistent ( m : Matrix; cols : Standard_Integer_Vectors.Vector;
                          x : Vector; tol : double_float ) return boolean is

  -- ALGORITHM : checks whether the inequality 0 >= c, with c > tol,
  --             can be derived from the combination of the inequalities.

    sum : double_float;

  begin
    for i in x'range loop            -- x should be a positive combination
      if x(i) < 0.0
       then return false;
      end if;
    end loop;
    for i in m'first(1)..m'last(1)-1 loop     -- check zero left hand side
      sum := 0.0;
      for j in x'range loop
        sum := sum + x(j)*m(i,cols(j));
      end loop;
      if abs(sum) > tol
       then return false;
      end if;
    end loop;
    sum := 0.0;                        -- right hand side must be positive
    for j in x'range loop
      sum := sum + x(j)*m(m'last(1),cols(j));
    end loop;
    return (sum > tol);
  end Inconsistent;

  function Inconsistent ( m : Matrix; i : integer32;
                          cols : Standard_Integer_Vectors.Vector; x : Vector;
                          tol : double_float ) return boolean is

    allcols : Standard_Integer_Vectors.Vector(cols'first..cols'last+1);
    allcoef : Vector(x'first..x'last+1);

  begin
    allcols(cols'range) := cols;      allcoef(x'range) := x;
    allcols(allcols'last) := i;
    if Is_In(cols,i)
     then allcoef(allcoef'last) := 0.0;
     else allcoef(allcoef'last) := 1.0;
    end if;
    return Inconsistent(m,allcols,allcoef,tol);
  end Inconsistent;

  function Inconsistent2D ( m : Matrix; cols : Standard_Integer_Vectors.Vector;
                            tol : double_float ) return Vector is

  -- REQUIRED : the normals in the corresponding inequality are opposite.

  -- DESCRIPTION :
  --   Returns the coefficients in the inconsistency proof.

    res : Vector(cols'range) := (cols'range => 1.0);
    f1,f2 : double_float;

  begin
    for i in m'first(1)..m'last(1)-1 loop
      f1 := abs(m(i,cols(1)));
      if f1 > tol then
        f2 := abs(m(i,cols(2)));
        if f1 > f2
         then res(2) := 1.0; res(1) := f2/f1;
         else res(1) := 1.0; res(2) := f1/f2;
        end if;
      end if;
      exit when (f1 > tol);
    end loop;
    return res;
  end Inconsistent2D;

  function Inconsistent2D ( m : Matrix; i : integer32;
                            cols : Standard_Integer_Vectors.Vector;
                            tol : double_float ) return Vector is

  -- REQUIRED : the inequalities of the columns define a nonsingular system.

  -- DESCRIPTION :
  --   Returns the coefficients in the inconsistency proof.

    res : Vector(cols'range) := (cols'range => 1.0);
    detm12 : constant double_float
           := m(1,cols(1))*m(2,cols(2)) - m(2,cols(1))*m(1,cols(2));

  begin
    if abs(detm12) > tol then
      res(1) := (m(2,i)*m(1,cols(2)) - m(1,i)*m(2,cols(2)))/detm12;
      res(2) := (m(1,i)*m(2,cols(1)) - m(2,i)*m(1,cols(1)))/detm12;
    end if;
    return res;
  end Inconsistent2D;

  function InconsistentnD
             ( m : Matrix; i,piv : integer32; x : Vector;
               cols : Standard_Integer_Vectors.Vector;
               k : integer32; tol : double_float ) return Vector is

  -- DESCRIPTION :
  --   Computes the coefficients of the inconsistency proof, based on the
  --   inconsistency proof of the eliminated problem.

    res : Vector(x'first..x'last+1);
    fac : double_float;

    procedure Compute_Factor is
    begin
      for j in x'range loop
        fac := fac + x(j)*m(piv,cols(j));
      end loop;
      fac := -fac/m(piv,i);
      if abs(fac) > tol then
        for j in x'range loop
          res(j) := x(j)/fac;
        end loop;
      else
        res(x'range) := x;
      end if;
    end Compute_Factor;

  begin
    if Is_In(cols,k) then
      res := (res'range => 0.0);
      if Is_All_In(cols,k)
       then res(res'first) := -m(piv,i)/m(piv,k);
       else fac := 0.0;
            Compute_Factor;
      end if;
    else
      fac := m(piv,k);
      Compute_Factor;
      if abs(fac) > tol
       then res(res'last) := 1.0/fac;
       else res(res'last) := 1.0;
      end if;
    end if;
    return res;
  end InconsistentnD;

-- CONSTRUCTORS :

  procedure Scale ( m : in out Matrix; i : in integer32;
                    tol : in double_float ) is

    ip : double_float := abs(m(Pivot(m,i),i)); --Inner_Product(m,i,i);

  begin
    if ip > tol then
      --ip := SQRT(ip);
      for j in m'range(1) loop
        m(j,i) := m(j,i)/ip;
      end loop;
    end if;
  end Scale;

  procedure Scale ( m : in out Matrix; tol : in double_float ) is
  begin
    for i in m'range(2) loop
      Scale(m,i,tol);
    end loop;
  end Scale;

  procedure Center ( m : in out Matrix; x : in Vector ) is
  begin
    for i in m'range(2) loop
      m(m'last(1),i) := m(m'last(1),i) - Evaluate(m,i,x);
    end loop;
   -- put_line("The centered inequality system : "); put(m,3,3,3);
  end Center;

  function Center ( m : Matrix; x : Vector ) return Matrix is

    mx : Matrix(m'range(1),m'range(2));   -- := m;  problems on RS/6000 ??

  begin
    Copy(m,mx);
    Center(mx,x);
    return mx;
  end Center;

  procedure Intersect2D ( m : in Matrix; i,j : in integer32;
                          tol : in double_float;
                          x : out Vector; sing : out boolean ) is

    detmij : constant double_float := m(1,i)*m(2,j) - m(2,i)*m(1,j);

  begin
    if abs(detmij) <= tol then
      sing := true;
    else
      sing := false;
      x(1) := (m(3,i)*m(2,j) - m(2,i)*m(3,j))/detmij;
      x(2) := (m(1,i)*m(3,j) - m(3,i)*m(1,j))/detmij;
    end if;
   -- put(" i : "); put(i,1); put("  j : "); put(j,1);
   -- put("  detmij : "); put(detmij,3,3,3); new_line;
  end Intersect2D;

-- SOLVE BY INTERSECTION :

  procedure Solve_Intersect2D
                ( m : in Matrix; i : in integer32; tol : in double_float;
                  x : in out Vector; fail : out boolean;
                  cols : out Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Enumerates all intersection points of the ith hyperplane with all the
  --   other previous ones.  The enumeration stops when either the current
  --   intersection point satisfies the system of inequalities or when the
  --   inequalities for the inconsistency proof are detected.

  -- REQUIRED : m'last(1) = 3 and the inequalities are centered.

    firstviol,j : integer32;
    ffail,sing : boolean := true;
    columns : Standard_Integer_Vectors.Vector(x'range) := (x'range => 0);

  begin
    j := m'first(2);
    loop
      Intersect2D(m,j,i,tol,x,sing);
      if not sing then
        firstviol := First_Violated(m,m'first(2),i-1,x,tol);
        if firstviol < j
         then columns(1) := firstviol; columns(2) := j; sing := true;
              x := Inconsistent2D(m,i,columns,tol);
        end if;
      else
        if Inner_Product(m,i,j) < -tol then
          sing := false; x := (x'range => 0.0);
          for k in m'first(1)..m'last(1)-1 loop
            if abs(m(k,i)) > tol then
              x(k) := m(3,i)/m(k,i);
              if m(k,i) > 0.0
               then sing := ((x(k) - m(3,j)/m(k,j)) > tol);
               else sing := ((m(3,j)/m(k,j) - x(k)) > tol);
              end if;
            end if;
            exit when (abs(m(k,i)) > tol);
          end loop;
          if sing then
            columns(1) := j; columns(2) := i;
            x := Inconsistent2D(m,columns,tol);
          end if;
        else
          sing := false;
        end if;
        if sing
         then firstviol := j;
         else firstviol := j+1;
        end if;
      end if;
      ffail := (firstviol < i);
      exit when not ffail or sing; 
      j := firstviol;
    end loop;
    fail := ffail;
    cols := columns;
  end Solve_Intersect2D;

  procedure Solve_IntersectnD
              ( m : in Matrix; i : in integer32; tol : in double_float;
                x : in out Vector; fail : out boolean;
                cols : out Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Eliminates one unknown by intersecting with the ith hyperplane.

  -- REQUIRED : m'last(1) > 3 and the inequalities are centered.

    m2 : Matrix(m'first(1)..m'last(1)-1,m'first(2)..i-1);
    piv : constant integer32 := Pivot(m,i);
    x2 : Vector(x'first..x'last-1);
    cols2 : Standard_Integer_Vectors.Vector(x2'range);
    k2 : integer32;

  begin
    if abs(m(piv,i)) > tol then
      m2 := Eliminate(m,i-1,i,piv,tol);
      x2 := (x2'range => 0.0);
      Solve(m2,tol,x2,k2,cols2);
      if k2 >= i then
        fail := false;
        x := Back_Substitute(m,i,piv,x2);
      else
        fail := true;
        cols(cols2'range) := cols2;
        cols(cols2'last+1) := k2;
        x := InconsistentnD(m,i,piv,x2,cols2,k2,tol);
      end if;
    end if;
  end Solve_IntersectnD;

  procedure Solve_Intersect
              ( m : in Matrix; i : in integer32; tol : in double_float;
                x : in out Vector; fail : out boolean;
                cols : out Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Intersects the inequalities with the ith hyperplane.

  -- REQUIRED : the inequalities are centered,
  --            the ith inequality is first one that is violated.

  begin
    if m'last(1) = 3 then
      Solve_Intersect2D(m,i,tol,x,fail,cols);
    elsif m'last(1) > 3 then
      Solve_IntersectnD(m,i,tol,x,fail,cols);
    end if;
  end Solve_Intersect;

  procedure Solve ( m : in Matrix; i : in integer32; tol : in double_float;
                    x : in out Vector; fail : out boolean;
                    cols : out Standard_Integer_Vectors.Vector ) is

    inx,xx : Vector(x'range) := (x'range => 0.0);
    wrk : Matrix(m'range(1),m'range(2));
    sing : boolean := true;
    fac : double_float;

  begin
    if abs(Inner_Product(m,i,i)) < tol then
      if m(m'last(1),i) > tol
       then fail := true; cols := (x'range => i); 
       else fail := false;
      end if;
    else
       if i = m'first(2) then
         x := (x'range => 0.0);
         for j in x'range loop
           if abs(m(j,i)) > tol then
             sing := false;
             x(j) := m(m'last(1),i)/m(j,i);
           end if;
           exit when not sing;
         end loop;
         fail := sing;
       else
         Copy(m,wrk);
         if x /= inx
          then Center(wrk,x);
         end if;
         fac := wrk(wrk'last(1),i)/Inner_Product(wrk,i,i);
         for j in inx'range loop
           inx(j) := fac*wrk(j,i);
         end loop;
         if First_Violated(wrk,inx,tol) > i then
           xx := x + inx;
           if Satisfies(m,xx,tol)
            then sing := false;
           end if;
         end if;
         if not sing then
           fail := false; x := xx;
         else
           Solve_Intersect(wrk,i,tol,inx,sing,cols);
           fail := sing;
           if not sing
            then Add(x,inx);
            else x := inx;
           end if;
         end if;
       end if;
    end if;
  end Solve;

  procedure Solve ( m : in Matrix; tol : in double_float; x : in out Vector;
                    k : out integer32;
                    cols : out Standard_Integer_Vectors.Vector ) is

    indviol : integer32 := First_Violated(m,x,tol);
    fail : boolean := false;

  begin
    while indviol <= m'last(2) loop
      Solve(m,indviol,tol,x,fail,cols);
      exit when fail;
      indviol := First_Violated(m,indviol+1,m'last(2),x,tol);
    end loop;
    if fail
     then k := indviol;
     else k := m'last(2) + 1;
    end if;
  end Solve;

  procedure Iterated_Solve
              ( m : in Matrix; tol : in double_float;
                x : in out Vector; k : out integer32;
                cols : out Standard_Integer_Vectors.Vector ) is

    indviol : integer32 := First_Violated(m,x,tol);
    fail : boolean := false;
    continue : boolean := true;

  begin
    while indviol <= m'last(2) loop
      Report(x,indviol-1,continue);
      exit when not continue;
      Solve(m,indviol,tol,x,fail,cols);
      exit when fail;
      indviol := First_Violated(m,indviol+1,m'last(2),x,tol);
    end loop;
    if fail
     then k := indviol;
     else k := m'last(2) + 1;
    end if;
  end Iterated_Solve;

end Floating_Linear_Inequality_Solvers;
