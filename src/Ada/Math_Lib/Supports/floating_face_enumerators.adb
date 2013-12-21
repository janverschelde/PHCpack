with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Givens_Rotations;                   use Givens_Rotations;
with Floating_Linear_Inequalities;       use Floating_Linear_Inequalities;

package body Floating_Face_Enumerators is

-- AUXILIARIES :

  function Is_In ( x : integer32; v : Standard_Integer_Vectors.Vector )
                 return boolean is

  -- DESCRIPTION :
  --   Returns true if there exists an entry k in v, with v(k) = x.

  begin
    for k in v'range loop
      if v(k) = x
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Is_Zero ( v : Standard_Floating_Vectors.Vector; tol : double_float )
                   return boolean is

  -- DESCRIPTION :
  --   Returns true if abs(v(i)) < tol, for all i in v'range.

  begin
    for i in v'range loop
      if abs(v(i)) > tol
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  function In_Edge ( x,a,b : Standard_Floating_Vectors.Vector;
                     tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the vector x lies between a and b.

    ba,xa : double_float;
    piv : integer32;

  begin
    for i in b'range loop     -- look for first i: b(i) - a(i) /= 0
      ba := b(i) - a(i);
      if abs(ba) > tol
       then piv := i; exit;
      end if;
    end loop;
    if abs(ba) < tol then
      return Equal(x,a);  -- in this case b = a, so test x = a
    else
      if ba < 0.0 then
        ba := -ba;
        xa := x(piv) - b(piv);
      else
        xa := x(piv) - a(piv);
      end if;
      if xa*ba >= 0.0 and then xa <= ba
       then return true;
       else return false;
      end if;
    end if;
  end In_Edge;

  function Convex_Decomposition 
             ( k : integer32; face : Standard_Integer_Vectors.Vector;
               x : Standard_Floating_Vectors.Vector;
               pts : Standard_Floating_VecVecs.VecVec;
               tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Decomposes the vector x in terms of the k+1 points of pts(face(*)).
  --   Returns true if x can be written as a convex combination.

    mat : Matrix(x'range,1..k);
    rhs : Standard_Floating_Vectors.Vector(x'range) := x;
    sol : Standard_Floating_Vectors.Vector(1..k) := (1..k => 0.0);
    ipvt : Standard_Integer_Vectors.Vector(1..k);
    sum : double_float := 0.0;

  begin
    for j in 1..k loop
      for i in mat'range(1) loop
        mat(i,j) := pts(face(j))(i) - pts(face(0))(i);
      end loop;
      ipvt(j) := j;
    end loop;
    Upper_Triangulate(mat,rhs,tol,ipvt);
    for j in k+1..rhs'last loop
      if abs(rhs(j)) > tol
       then return false;
      end if;
    end loop;
    Solve(mat,rhs,tol,sol);
    for j in sol'range loop
      if sol(j) < -tol
       then return false;
       else sum := sum + sol(j);
      end if;
    end loop;
    return (abs(sum-1.0) < tol);
  end Convex_Decomposition;

  function In_Face ( k : in integer32; face : Standard_Integer_Vectors.Vector;
                     x : Standard_Floating_Vectors.Vector; 
                     pts : Standard_Floating_VecVecs.VecVec;
                     tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Returns true if x lies in the given k-face, which contains entries
  --   to its elements in pts.  This means that we have to verify whether
  --   the vector x can be written as a convex combination of the points
  --   in the face.

  begin
    if k = 0 then
      return Equal(pts(face(face'first)).all,x);
    elsif k = 1 then
      return In_Edge(x,pts(face(face'first)).all,
                       pts(face(face'first+1)).all,tol);
    else
      return Convex_Decomposition(k,face,x,pts,tol);
    end if;
  end In_Face;

-- ELIMINATE TO MAKE UPPER TRIANGULAR :

  procedure Eliminate ( pivot : in integer32;
                        elim : in Standard_Floating_Vectors.Vector;
                        target : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Makes target(pivot) = 0.0 by means of making a linear
  --   combination of the vectors target and elim.
  --   Note that target is viewed as an inequality, so that we make sure
  --   to multiply it only with a positive number.

  -- REQUIRED : abs(target(pivot)) > tol and abs(elim(pivot)) > tol.

    a,b,tarfac,elifac : double_float;

  begin
    a := elim(pivot); b := target(pivot);
    if a > 0.0
     then tarfac := a;  elifac := b;
     else tarfac := -a; elifac := -b;
    end if;
    for j in target'range loop
      target(j) := tarfac*target(j) - elifac*elim(j);
    end loop;
  end Eliminate;

  procedure Eliminate ( l : in integer32;
                        pivots : in Standard_Integer_Vectors.Vector; 
                        elim : in Standard_Floating_VecVecs.VecVec;
                        tol : in double_float;
                        target : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Makes target(pivots(i)) = 0 by means of making a linear
  --   combination of the vectors target and elim(i), for i in 1..l.

  -- REQUIRED : elim(i)(pivots(i))) > tol.

  begin
    for i in 1..l loop
      if abs(target(pivots(i))) > tol
       then Eliminate(pivots(i),elim(i).all,target);
      end if;
    end loop;
  end Eliminate;

  function Pivot_after_Elimination
             ( l,k : in integer32; point : Standard_Floating_Vectors.Vector;
               face,pivots : Standard_Integer_Vectors.Vector;
               elim : Standard_Floating_VecVecs.VecVec; tol : double_float )
             return integer32 is

  -- DESCRIPTION :
  --   Returns the first nonzero element of the given point after elimination
  --   w.r.t. the entries in the face with lower index.

    work : Standard_Floating_Vectors.Vector(point'range)
         := point - elim(1-k).all;
    pivot : integer32;

  begin
    for i in (face'first+1)..face'last loop
      if (face(i) < l) and then (abs(work(pivots(i))) > tol)
       then Eliminate(pivots(i),elim(i).all,work);
      end if;
      exit when (face(i) > l);
    end loop;
    pivot := 0;
    for i in work'range loop
      if abs(work(pivots(i))) > tol
       then pivot := i;
      end if;
      exit when (pivot /= 0);
    end loop;
    return pivot;
  end Pivot_after_Elimination;

  procedure Update_Pivots
               ( point : in Standard_Floating_Vectors.Vector;
                 l : in integer32;
                 pivots : in out Standard_Integer_Vectors.Vector;
                 tol : in double_float; pivot : out integer32 ) is

  -- DESCRIPTION :
  --   Searches in the point(l..point'last) for the first nonzero entry
  --   and updates the pivoting information.

    temp,piv : integer32;

  begin
    piv := 0;
    for i in l..point'last loop
      if abs(point(pivots(i))) > tol
       then piv := i;
      end if;
      exit when (piv /= 0);
    end loop;
    if piv /= 0 and then (piv /= l) then
      temp := pivots(l);
      pivots(l) := pivots(piv);
      pivots(piv) := temp;
    end if;
    pivot := piv;
  end Update_Pivots;

  procedure Update_Eliminator
               ( elim : in out Standard_Floating_VecVecs.VecVec;
                 l : in integer32;
                 pivots : in out Standard_Integer_Vectors.Vector;
                 point : in Standard_Floating_Vectors.Vector;
                 tol : in double_float;
                 pivot : out integer32 ) is

  -- DESCRIPTION :
  --   Updates the vector of eliminators, by adding the lth elimination
  --   equation.  This procedure ensures the invariant condition on the
  --   eliminator, which has to be upper triangular.  If this cannot be
  --   achieved, degeneracy is indicated by pivot = 0.

  -- ON ENTRY :
  --   elim      vector of elimination equations: elim(i)(pivots(i)) /= 0
  --             and for j in 1..(i-1) : elim(i)(pivots(j)) = 0,
  --             for i in 1..(l-1), elim(0) contains the basis point;
  --   l         index of current elimination equation to be made;
  --   pivots    contains the pivoting information;
  --   tol       tolerance on the precision;
  --   point     new point to make the equation `point - elim(0) = 0'.

  -- ON RETURN :
  --   elim      if not degen, then elim(l)(pivots(l)) /= 0 and
  --             for i in 1..(l-1): elim(l)(pivots(i)) = 0;
  --   pivots    updated pivot information;
  --   pivot     the pivot that has been used for elim(l)(pivots(l)) /= 0;
  --             piv = 0 when the new elimination equation elim(l)
  --             became entirely zero after ensuring the invariant condition.

  begin
    elim(l) := new Standard_Floating_Vectors.Vector'(point - elim(0).all);
    Eliminate(l-1,pivots,elim,tol,elim(l).all);
    Update_Pivots(elim(l).all,l,pivots,tol,pivot);
  end Update_Eliminator;

  procedure Update_Eliminator_for_Sum
               ( elim : in out Standard_Floating_VecVecs.VecVec;
                 l : in integer32;
                 pivots : in out Standard_Integer_Vectors.Vector;
                 point : in Standard_Floating_Vectors.Vector;
                 index : in integer32;
                 tol : in double_float; pivot : out integer32 ) is

  -- DESCRIPTION :
  --   Updates the vector of eliminators, by adding the lth elimination
  --   equation.  This procedure ensures the invariant condition on the
  --   eliminator, which has to be upper triangular.  If this cannot be
  --   achieved, degeneracy is indicated by pivot = 0.

  -- ON ENTRY :
  --   elim      vector of elimination equations: elim(i)(pivots(i)) /= 0
  --             and for j in 1..(i-1) : elim(i)(pivots(j)) = 0,
  --             for i in 1..(l-1), elim(1-index) contains the basis point;
  --   l         index of current elimination equation to be made;
  --   pivots    contains the pivoting information;
  --   point     new point to make the equation `point - elim(1-index) = 0';
  --   index     indicates the current polytope;
  --   tol       tolerance for precision.

  -- ON RETURN :
  --   elim      if not degen, then elim(l)(pivots(l)) /= 0 and
  --             for i in 1..(l-1): elim(l)(pivots(i)) = 0;
  --   pivots    updated pivot information;
  --   piv       the pivot that has been used for elim(l)(pivots(l)) /= 0;
  --             piv = 0 when the new elimination equation elim(l)
  --             became entirely zero after ensuring the invariant condition.

  begin
    elim(l) := new Standard_Floating_Vectors.Vector'(point - elim(1-index).all);
    Eliminate(l-1,pivots,elim,tol,elim(l).all);
    Update_Pivots(elim(l).all,l,pivots,tol,pivot);
  end Update_Eliminator_for_Sum;

-- CREATE TABLEAU OF INEQUALITIES :

  procedure Create_Tableau_for_Vertices
               ( i,n : in integer32;
                 pts : in Standard_Floating_VecVecs.VecVec;
                 tab : out Matrix;
                 rhs : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Creates a system of linear inequalities that is feasible when
  --   the ith point is spanned by the other points in pts.

    column : integer32 := tab'first(2);

  begin
    for j in pts'range loop
      if j /= i then
        for k in pts(j)'range loop
          tab(k,column) := pts(j)(k);
        end loop;
        tab(tab'last(1),column) := 1.0;
        column := column + 1;
      end if;
    end loop;
    for k in pts(i)'range loop
      rhs(k) := pts(i)(k);
    end loop;
    rhs(n+1) := 1.0;
  end Create_Tableau_for_Vertices;

  procedure Create_Tableau_for_Edges
               ( i,j,n,pivot : in integer32;
                 elim : in Standard_Floating_Vectors.Vector;
                 pts : in Standard_Floating_VecVecs.VecVec; 
                 tol : in double_float; tab : out matrix; 
                 rhs : out Standard_Floating_Vectors.Vector;
                 lastcol : out integer32; degenerate : out boolean ) is

  -- DESCRIPTION :
  --   Creates a system of linear inequalities that express the fact
  --   that the points i and j span an edge of the convex hull.
  --   Degenerate means that there exists a point on the edge spanned
  --   by the points i and j that contains the tested edge.

    column : integer32 := 1;
    ineq : Standard_Floating_Vectors.Vector(1..n);

  begin
    degenerate := false;
    for k in pts'range loop
      if (k/=i) and then (k/=j) then
        ineq := pts(k).all - pts(i).all;        -- compute inequality
        if abs(ineq(pivot)) > tol               -- make ineq(pivot) = 0
         then Eliminate(pivot,elim,ineq);
        end if;
        if Is_Zero(ineq,tol) then                    -- check if degenerate
          if not In_Edge(pts(k).all,pts(i).all,pts(j).all,tol)
           then degenerate := true; return;
          end if;
        else
          for l in 1..n loop                -- fill in the column
            if l < pivot then
              tab(l,column) := ineq(l);
            elsif l > pivot then
              tab(l-1,column) := ineq(l);
            end if;
          end loop;
          tab(tab'last(1),column) := 1.0;
          column := column + 1;
        end if;
      end if;
    end loop;
    for k in tab'first(1)..(tab'last(1)-1) loop  -- right hand side
      rhs(k) := 0.0;
    end loop;
    rhs(tab'last(1)) := 1.0;
    lastcol := column-1;
  end Create_Tableau_for_Edges;

  procedure Create_Tableau_for_Faces
               ( k,n : in integer32;
                 face,pivots : in Standard_Integer_Vectors.Vector;
                 pts,elim : in Standard_Floating_VecVecs.VecVec;
                 tol : in double_float; tab : out Matrix;
                 rhs : out Standard_Floating_Vectors.Vector;
                 lastcol : out integer32; degenerate : out boolean ) is

  -- DESCRIPTION :
  --   Creates a system of linear inequalities that express the fact
  --   that the points pts(face(*)) span a face of the convex hull.
  --   Degenerate means that there exists a point on the face spanned
  --   by the points in the face that contains the tested face.

    column : integer32 := 1;
    ineq : Standard_Floating_Vectors.Vector(1..n);

  begin
    degenerate := false;
    for l in pts'range loop
      if not Is_In(l,face) then
        ineq := pts(l).all - elim(0).all;  -- new inequality
        Eliminate(k,pivots,elim,tol,ineq); -- ineq(pivots(i)) = 0, i=1,2,..k
        if Is_Zero(ineq,tol) then
          if not In_Face(k,face,pts(l).all,pts,tol)
            and then l < face(face'last)    -- lexicographic enumeration
            and then (Pivot_after_Elimination
                          (l,1,pts(l).all,face,pivots,elim,tol) /= 0)
           then degenerate := true; return;
          end if;
        else
          for ll in (k+1)..n loop              -- fill in the column
            tab(ll-k,column) := ineq(pivots(ll));
          end loop;
          tab(tab'last(1),column) := 1.0;
          column := column + 1;
        end if;
      end if;
    end loop;
    for l in tab'first(1)..(tab'last(1)-1) loop  -- right hand side
      rhs(l) := 0.0;
    end loop;
    rhs(tab'last(1)) := 1.0;
    lastcol := column-1;
  end Create_Tableau_for_Faces;

  procedure Create_Tableau_for_Faces_of_Sum
               ( k,n,i,r : in integer32;
                 ind,pivots : in Standard_Integer_Vectors.Vector;
                 pts,elim : Standard_Floating_VecVecs.VecVec;
                 tol : double_float;
                 face : in Standard_Integer_VecVecs.VecVec;
                 tab : out Matrix; rhs : out Standard_Floating_Vectors.Vector;
                 lastcol : out integer32; degenerate : out boolean ) is

  -- DESCRIPTION :
  --   Creates the table of inequalities for determining whether the given
  --   candidate face spans a face of the sum polytope.

  -- ON ENTRY :
  --   k         dimension of the face on the sum;
  --   n         dimension of the points;
  --   i         current polytope;
  --   r         number of polytopes;
  --   ind       indicate the beginning vector in pts of each polytope;
  --   etc...

    column : integer32 := 1;
    ineq : Standard_Floating_Vectors.Vector(1..n);
    last : integer32;

  begin
    degenerate := false;
    for l1 in face'first..i loop
      if l1 = r
       then last := pts'last;
       else last := ind(l1+1)-1;
      end if;
      for l2 in ind(l1)..last loop
        if not Is_In(l2,face(l1).all) then
          ineq := pts(l2).all - elim(1-l1).all;  -- new inequality
          Eliminate(k,pivots,elim,tol,ineq); -- ineq(pivots(i)) = 0, i=1,2,..k
          if Is_Zero(ineq,tol) then
            if not In_Face(face(l1)'length-1,face(l1).all,pts(l2).all,pts,tol)
              and then face(l1)(face(l1)'first) <= l2
              and then l2 < face(l1)(face(l1)'last)
                                            -- lexicographic enumeration
              and then (Pivot_after_Elimination
                          (l2,l1,pts(l2).all,face(l1).all,pivots,elim,tol) /= 0)
             then degenerate := true; return;
            end if;
          else
            for ll in (k+1)..n loop                 -- fill in the column
              tab(ll-k,column) := ineq(pivots(ll));
            end loop;
            tab(tab'last(1),column) := 1.0;
            column := column + 1;
          end if;
        end if;
      end loop;
    end loop;
    for l in tab'first(1)..(tab'last(1)-1) loop  -- right hand side
      rhs(l) := 0.0;
    end loop;
    rhs(tab'last(1)) := 1.0;
    lastcol := column-1;
  end Create_Tableau_for_Faces_of_Sum;

  procedure Create_Tableau_for_Lower_Edges
               ( i,j,n,pivot : in integer32;
                 elim : in Standard_Floating_Vectors.Vector;
                 pts : in Standard_Floating_VecVecs.VecVec;
                 tol : in double_float; tab : out Matrix;
                 rhs : out Standard_Floating_Vectors.Vector;
                 lastcol : out integer32; degenerate : out boolean ) is

  -- DESCRIPTION :
  --   Creates a system of linear inequalities that express the fact
  --   that the points i and j spann an edge of the lower hull.
  --   Degenerate means that there exists a point on the edge spanned
  --   by the points i and j that contains the tested edge.

    column : integer32 := 1;
    ineq : Standard_Floating_Vectors.Vector(1..n);

  begin
    degenerate := false;
    for k in pts'range loop
      if (k/=i) and then (k/=j) then
        ineq := pts(k).all - pts(i).all;         -- compute inequality
        if abs(ineq(pivot)) > tol                -- make ineq(pivot) = 0
         then Eliminate(pivot,elim,ineq);
        end if;
        if Is_Zero(ineq,tol) then                -- check if degenerate
          if not In_Edge(pts(k).all,pts(i).all,pts(j).all,tol)
           then degenerate := true; return;
          end if;
        else
          for l in 1..n loop                 -- fill in the column
            if l < pivot then
              tab(l,column) := ineq(l);
            elsif l > pivot then
              tab(l-1,column) := ineq(l);
            end if;
          end loop;
          column := column + 1;
        end if;
      end if;
    end loop;
    for k in tab'first(1)..(tab'last(1)-1) loop   -- right hand side
      rhs(k) := 0.0;
    end loop;
    rhs(tab'last(1)) := -1.0;
    lastcol := column-1;
  end Create_Tableau_for_Lower_Edges;

  procedure Create_Tableau_for_Lower_Faces
               ( k,n : in integer32;
                 face,pivots : in Standard_Integer_Vectors.Vector;
                 pts,elim : in Standard_Floating_VecVecs.VecVec;
                 tol : in double_float; tab : out Matrix;
                 rhs : out Standard_Floating_Vectors.Vector;
                 lastcol : out integer32; degenerate : out boolean ) is

    column : integer32 := 1;
    ineq : Standard_Floating_Vectors.Vector(1..n);

  begin
    degenerate := false;
    for l in pts'range loop
      if not Is_In(l,face) then
        ineq := pts(l).all - elim(0).all;   -- new inequality
        Eliminate(k,pivots,elim,tol,ineq); -- ineq(pivots(i)) = 0, i=1,2,..k
        if Is_Zero(ineq,tol) then
          if not In_Face(k,face,pts(l).all,pts,tol)
            and then l < face(face'last)     -- lexicographic enumeration
            and then (Pivot_after_Elimination
                        (l,1,pts(l).all,face,pivots,elim,tol) /= 0)
           then degenerate := true; return;
          end if;
        else 
          for ll in (k+1)..n loop             -- fill in the column
            tab(ll-k,column) := ineq(pivots(ll));
          end loop;
          column := column + 1;
        end if;
      end if;
    end loop;
    for l in tab'first(1)..(tab'last(1)-1) loop  -- right hand side
      rhs(l) := 0.0;
    end loop;
    rhs(tab'last(1)) := -1.0;
    lastcol := column-1;
  end Create_Tableau_for_Lower_Faces;

-- DETERMINE WHETHER THE CANDIDATE IS VERTEX, SPANS AN EDGE OR A K-FACE :

  function Is_Vertex ( i,m,n : integer32;
                       pts : Standard_Floating_VecVecs.VecVec;
                       tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   This function determines whether the ith point is a vertex of
  --   the m n-dimensional points in pts, by deciding on the feasibility
  --   of a linear-inequality system.

    tableau : Matrix(1..n+1,1..m);
    rhs : Standard_Floating_Vectors.Vector(tableau'range(1));
    sol : Standard_Floating_Vectors.Vector(tableau'range(1));
    col : Standard_Integer_Vectors.Vector(tableau'range(1));
    feasible : boolean;

  begin
    Create_Tableau_for_Vertices(i,n,pts,tableau,rhs);
    Complementary_Slackness(tableau,rhs,tol,sol,col,feasible);
    return not feasible;
  end Is_Vertex;

  function Is_Edge ( i,j,m,n : integer32;
                     pts : Standard_Floating_VecVecs.VecVec;
                     tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   This function determines whether the ith and jth point span an
  --   edge of the convex hull of m n-dimensional points in pts,
  --   by deciding on the feasibility of a linear-inequality system.

    elim : Standard_Floating_Vectors.Vector(1..n) := pts(i).all - pts(j).all;
    pivot : integer32;

  begin
    pivot := elim'first - 1;
    for k in elim'range loop
      if abs(elim(k)) > tol
       then pivot := k;
      end if;
      exit when pivot >= elim'first;
    end loop;
    if pivot < elim'first then
      return false;
    else
      declare
        tab : Matrix(1..n,1..m-1);
        rhs : Standard_Floating_Vectors.Vector(tab'range(1));
        sol : Standard_Floating_Vectors.Vector(tab'range(1));
        col : Standard_Integer_Vectors.Vector(tab'range(1));
        deg,feasible : boolean;
        lst : integer32;
      begin
        Create_Tableau_for_Edges(i,j,n,pivot,elim,pts,tol,tab,rhs,lst,deg);
        if deg then
          return false;
        elsif lst = 0 then
          return true;
        else
          Complementary_Slackness(tab,lst,rhs,tol,sol,col,feasible);
          return not feasible;
        end if;
      end;
    end if;
  end Is_Edge;

  function Is_Face ( k,n,m : integer32;
                     elim,pts : Standard_Floating_VecVecs.VecVec;
                     pivots,face : Standard_Integer_Vectors.Vector;
                     tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Applies Complementary Slackness to determine whether the given
  --   candidate face is a face of the polytope.

  -- ON ENTRY :
  --   k             dimension of the candidate face;
  --   n             dimension of the vector space;
  --   m             number of points which span the polytope;
  --   elim          elimination equations in upper triangular form;
  --   pts           the points which span the polytope;
  --   pivots        pivoting information for the elimination equations;
  --   face          entries of the points which span the candidate face;
  --   tol           tolerance on the precision.

    nn : constant integer32 := n-k+1;
    tab : Matrix(1..nn,1..(m-k));
    rhs : Standard_Floating_Vectors.Vector(tab'range(1));
    sol : Standard_Floating_Vectors.Vector(tab'range(1));
    col : Standard_Integer_Vectors.Vector(tab'range(1));
    deg,feasible : boolean;
    lst : integer32;

  begin
   -- put("The face being tested : "); put(face); new_line;
   -- put("The pivots : "); put(pivots); new_line;
   -- put_line("The elimination equations : ");
   -- for i in 1..elim'last loop
   --   put(elim(i)); put(" = 0 ");
   -- end loop;
   -- new_line;
    if m - k <= 1 then
      return true;
    else
      Create_Tableau_for_Faces(k,n,face,pivots,pts,elim,tol,tab,rhs,lst,deg);
      if deg then
        return false;
      elsif lst = 0 then
        return true;
      else 
        Complementary_Slackness(tab,lst,rhs,tol,sol,col,feasible);
        return not feasible;
      end if;
    end if;
  end Is_Face;

  function Is_Face_of_Sum
                ( k,n,m,i,r : integer32;
                  elim,pts : Standard_Floating_VecVecs.VecVec;
                  tol : double_float;
                  face : Standard_Integer_VecVecs.VecVec;
                  ind,pivots : Standard_Integer_Vectors.Vector )
                return boolean is

  -- DESCRIPTION :
  --   Applies Complementary Slackness to determine whether the given
  --   candidate face is a face of the polytope.

  -- ON ENTRY :
  --   k          dimension of the candidate face;
  --   n          dimension of the vector space;
  --   m          number of total points which span the polytope;
  --   i          current polytope;
  --   r          number of different polytopes;
  --   elim       elimination equations in upper triangular form;
  --   pts        the points which span the polytope;
  --   tol        tolerance on precision;
  --   face       entries of the points which span the candidate face;
  --   ind        indicates the starting vector in pts of each polytope;
  --   pivots     pivoting information for the elimination equations.

    nn : constant integer32 := n-k+1;
    tab : Matrix(1..nn,1..(m-k));
    rhs : Standard_Floating_Vectors.Vector(tab'range(1));
    sol : Standard_Floating_Vectors.Vector(tab'range(1));
    col : Standard_Integer_Vectors.Vector(tab'range(1));
    deg,feasible : boolean;
    lst : integer32;

  begin
    if m - k <= 1 then
      return true;
    else
      Create_Tableau_for_Faces_of_Sum
        (k,n,i,r,ind,pivots,pts,elim,tol,face,tab,rhs,lst,deg);
      if deg then
        return false;
      elsif lst = 0 then
        return true;
      else
        Complementary_Slackness(tab,lst,rhs,tol,sol,col,feasible);
        return not feasible;
      end if;
    end if;
  end Is_Face_of_Sum;

  function Is_Lower_Edge ( i,j,m,n : integer32;
                           pts : Standard_Floating_VecVecs.VecVec;
                           tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   This function determines whether the ith and jth point span an
  --   edge of the lower convex hull of m n-dimensional points in pts,
  --   by deciding on the feasibility of a linear-inequality system.

    elim : Standard_Floating_Vectors.Vector(1..n) := pts(i).all - pts(j).all;
    pivot : integer32;

  begin
    pivot := elim'first - 1;
    for k in elim'range loop
      if abs(elim(k)) > tol
       then pivot := k;
      end if;
      exit when pivot >= elim'first;
    end loop;
    if pivot < elim'first or else (pivot = elim'last) then
      return false;
    else
      declare
        tab : Matrix(1..n-1,1..m-1);
        rhs : Standard_Floating_Vectors.Vector(tab'range(1));
        sol : Standard_Floating_Vectors.Vector(tab'range(1));
        col : Standard_Integer_Vectors.Vector(tab'range(1));
        deg,feasible : boolean;
        lst : integer32;
      begin
        Create_Tableau_for_Lower_Edges
          (i,j,n,pivot,elim,pts,tol,tab,rhs,lst,deg);
        if deg then
          return false;
        elsif lst = 0 then
          return true;
        else
          Complementary_Slackness(tab,lst,rhs,tol,sol,col,feasible);
          return not feasible;
        end if;
      end;
    end if;
  end Is_Lower_Edge;

  function Is_Lower_Face 
                ( k,n,m : in integer32;
                  elim,pts : Standard_Floating_VecVecs.VecVec;
                  pivots,face : Standard_Integer_Vectors.Vector;
                  tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Applies Complementary Slackness to determine whether the given
  --   candidate face is a face of the lower hull of the polytope.

  -- ON ENTRY :
  --   k          dimension of the candidate face;
  --   n          dimension of the vector space;
  --   m          number of points which span the polytope;
  --   elim       elimination equations in upper triangular form;
  --   pts        the points which span the polytope;
  --   pivots     pivoting information for the elimination equations;
  --   face       entries of the points which span the candidate face.

    nn : constant integer32 := n-k;
    tab : Matrix(1..nn,1..(m-k));
    rhs : Standard_Floating_Vectors.Vector(tab'range(1));
    sol : Standard_Floating_Vectors.Vector(tab'range(1));
    col : Standard_Integer_Vectors.Vector(tab'range(1));
    deg,feasible : boolean;
    lst : integer32;

  begin
   -- put("The pivots : "); put(pivots); new_line;
   -- put_line("The elimination equations : ");
   -- for i in 1..elim'last loop
   --   put(elim(i)); put(" = 0 ");
   -- end loop;
   -- new_line;
    if m - k <= 1 then
      return true;
    else
      Create_Tableau_for_Lower_Faces
        (k,n,face,pivots,pts,elim,tol,tab,rhs,lst,deg);
      if deg then
        return false;
      elsif lst = 0 then
        return true;
      else
        Complementary_Slackness(tab,lst,rhs,tol,sol,col,feasible);
        return not feasible;
      end if;
    end if;
  end Is_Lower_Face;

-- TARGET ROUTINES :

  procedure Enumerate_Vertices ( pts : in Standard_Floating_VecVecs.VecVec;
                                 tol : in double_float ) is

    continue : boolean := true;
    n : constant integer32 := pts(pts'first).all'length;
    m : constant integer32 := pts'length - 1;

  begin
    for i in pts'range loop
      if Is_Vertex(i,m,n,pts,tol)
       then Process(i,continue);
      end if;
      exit when not continue;
    end loop;
  end Enumerate_Vertices;

  procedure Enumerate_Edges ( pts : in Standard_Floating_VecVecs.VecVec;
                              tol : in double_float ) is
  
    continue : boolean := true;
    n : constant integer32 := pts(pts'first).all'length;
    m : constant integer32 := pts'length;

    procedure Candidate_Edges ( i,n : in integer32 ) is
    begin
      for j in (i+1)..pts'last loop    -- enumerate all candidates
       -- put("Verifying "); put(pts(i).all); put(" and");
       -- put(pts(j).all); put(" :"); new_line;
        if Is_Edge(i,j,m,n,pts,tol)
         then Process(i,j,continue);
        end if;
        exit when not continue;
      end loop;
    end Candidate_Edges;

  begin
    for i in pts'first..(pts'last-1) loop
      Candidate_Edges(i,n);
    end loop;
  end Enumerate_Edges;

  procedure Enumerate_Lower_Edges 
              ( pts : in Standard_Floating_VecVecs.VecVec;
                tol : in double_float ) is

    continue : boolean := true;
    n : constant integer32  := pts(pts'first).all'length;
    m : constant integer32 := pts'length;

    procedure Candidate_Lower_Edges ( i,n : in integer32 ) is
    begin
      for j in (i+1)..pts'last loop    -- enumerate all candidates
       -- put("Verifying "); put(pts(i).all); put(" and");
       -- put(pts(j).all); put(" :"); new_line;
        if Is_Lower_Edge(i,j,m,n,pts,tol)
         then Process(i,j,continue);
        end if;
        exit when not continue;
      end loop;
    end Candidate_Lower_Edges;

  begin
    for i in pts'first..(pts'last-1) loop
      Candidate_Lower_Edges(i,n);
    end loop;
  end Enumerate_Lower_Edges;

  procedure Enumerate_Faces ( k : in integer32;
                              pts : in Standard_Floating_VecVecs.VecVec;
                              tol : in double_float ) is

    m : constant integer32 := pts'length;
    n : constant integer32 := pts(pts'first).all'length;
    candidate : Standard_Integer_Vectors.Vector(0..k);
    elim : Standard_Floating_VecVecs.VecVec(0..k);
    continue : boolean := true;

    procedure Candidate_Faces ( ipvt : in Standard_Integer_Vectors.Vector;
                                i,l : in integer32 ) is

    -- DESCRIPTION :
    --   Enumerates all candidate k-faces in lexicographic order.

      piv : integer32;
      pivots : Standard_Integer_Vectors.Vector(1..n);

    begin
      if l > k then
        if (k = m)
          or else Is_Face(k,n,m,elim,pts,ipvt,candidate,tol)
         then Process(candidate,continue);
        end if;
      else
        for j in (i+1)..pts'last loop
          candidate(l) := j;
          pivots := ipvt;
          Update_Eliminator(elim,l,pivots,pts(j).all,tol,piv);
          if piv /= 0
           then Candidate_Faces(pivots,j,l+1);
          end if;
          Clear(elim(l));
          exit when not continue;
        end loop;
      end if;
    end Candidate_Faces;

  begin
    if k <= m then
      declare
        ipvt : Standard_Integer_Vectors.Vector(1..n);
      begin
        for i in ipvt'range loop
          ipvt(i) := i;
        end loop;
        for i in pts'first..(pts'last-k) loop
          candidate(0) := i; 
          elim(0) := pts(i);  -- basis point
          Candidate_Faces(ipvt,i,1);
        end loop;
      end;
    end if;
  end Enumerate_Faces;

  procedure Enumerate_Lower_Faces 
              ( k : in integer32;
                pts : in Standard_Floating_VecVecs.VecVec;
                tol : in double_float ) is

    m : constant integer32 := pts'length;
    n : constant integer32 := pts(pts'first).all'length;
    candidate : Standard_Integer_Vectors.Vector(0..k);
    elim : Standard_Floating_VecVecs.VecVec(0..k);
    continue : boolean := true;

    procedure Candidate_Lower_Faces
                ( ipvt : in Standard_Integer_Vectors.Vector;
                  i,l : in integer32 ) is

    -- DESCRIPTION :
    --   Enumerates all candidate k-faces in lexicographic order.

      piv : integer32;
      pivots : Standard_Integer_Vectors.Vector(1..n);

    begin
      if l > k then
        --put_line("Testing the following candidate face :");
        --for ii in candidate'first..candidate'last-1 loop
        --  put(pts(candidate(ii))); put(" & ");
        --end loop;
        --put(pts(candidate(candidate'last))); new_line;
        if (k = m) or else Is_Lower_Face(k,n,m,elim,pts,ipvt,candidate,tol)
         then Process(candidate,continue);
        end if;
      else
        for j in (i+1)..pts'last loop
          candidate(l) := j;
          pivots := ipvt;
          -- put("Picking "); put(pts(j));
          -- put(" Pivots : "); put(pivots); new_line;
          Update_Eliminator(elim,l,pivots,pts(j).all,tol,piv);
          -- put(" update of eliminator piv = "); put(piv,1);
          -- put(" Pivots : "); put(pivots); new_line;
          if (piv /= 0) and (piv /= n)
           then Candidate_Lower_Faces(pivots,j,l+1);
          end if;
          Clear(elim(l));
          exit when not continue;
        end loop;
      end if;
    end Candidate_Lower_Faces;

  begin
    if k <= m then
      declare
        ipvt : Standard_Integer_Vectors.Vector(1..n);
      begin
        for i in ipvt'range loop
          ipvt(i) := i;
        end loop;
        for i in pts'first..(pts'last-k) loop
          candidate(0) := i;
          elim(0) := pts(i);  -- basis point
          Candidate_Lower_Faces(ipvt,i,1);
        end loop;
      end;
    end if;
  end Enumerate_Lower_Faces;

  procedure Enumerate_Faces_of_Sum
              ( ind,typ : in Standard_Integer_Vectors.Vector;
                k : in integer32;
                pts : in Standard_Floating_VecVecs.VecVec;
                tol : in double_float ) is

    m : constant integer32 := pts'length;             -- number of points
    n : constant integer32 := pts(pts'first)'length;  -- dimension of vectors
    r : constant integer32 := ind'length;             -- number of polytopes
    candidates : Standard_Integer_VecVecs.VecVec(1..r);
    elim : Standard_Floating_VecVecs.VecVec(1-r..k);
    pivots : Standard_Integer_Vectors.Vector(1..n);
    continue : boolean := true;
    lasti,sum : integer32;

    procedure Candidate_Faces_of_Sum
                ( ipvt : in Standard_Integer_Vectors.Vector;
                  i : in integer32 ) is

    -- DESCRIPTION :
    --   Enumerates all faces of the given type on the sum,
    --   i indicates the current polytope.

      procedure Candidate_Faces
                  ( ipvt : in Standard_Integer_Vectors.Vector;
                    start,l : in integer32 ) is
    
      -- DESCRIPTION :
      --   Enumerates all candidate k-faces, with k = typ(i).
      --   The parameter l indicates the current element of pts to be chosen.
      --   The previously chosen point is given by start.

        piv,last : integer32;

      begin
        if i = r
         then last := m;
         else last := ind(i+1)-1;
        end if;
        if l > typ(i) then
          if (typ(i) = last-ind(i)+1)
            or else Is_Face_of_Sum
                       (sum,n,last-i+1,i,r,elim,pts(pts'first..last),tol,
                        candidates,ind,ipvt)
           then
             if i = r
              then Process(candidates,continue);
              else Candidate_Faces_of_Sum(ipvt,i+1);
             end if;
          end if;
        else
          for j in (start+1)..(last-typ(i)+l) loop
            candidates(i)(l) := j;
            if l = 0 then
              Candidate_Faces(ipvt,j,l+1);
            else
              pivots := ipvt;
              Update_Eliminator_for_Sum
                (elim,sum-typ(i)+l,pivots,pts(j).all,i,tol,piv);
              if piv /= 0
               then Candidate_Faces(pivots,j,l+1);
              end if;
              Clear(elim(sum-typ(i)+l));
            end if;
            exit when not continue;
          end loop;
        end if;
      end Candidate_Faces;

    begin
      candidates(i) := new Standard_Integer_Vectors.Vector(0..typ(i));
      if i = r
       then lasti := pts'last;
       else lasti := ind(i+1)-1;
      end if;
      sum := sum + typ(i);
      for j in ind(i)..(lasti-typ(i)) loop
        candidates(i)(0) := j;
        elim(1-i) := pts(j);
        Candidate_Faces(ipvt,j,1);
      end loop;
      sum := sum - typ(i);
     -- for j in (sum+1)..(sum+typ(i)) loop
     --   Clear(elim(j));
     -- end loop;
      Clear(candidates(i));
    end Candidate_Faces_of_Sum;

  begin
    declare
      ipvt : Standard_Integer_Vectors.Vector(1..n);
    begin
      for i in ipvt'range loop
        ipvt(i) := i;
      end loop;
      sum := 0;
      Candidate_Faces_of_Sum(ipvt,1);
    end;
  end Enumerate_Faces_of_Sum;

end Floating_Face_Enumerators;
