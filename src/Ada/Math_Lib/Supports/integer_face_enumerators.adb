with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Standard_Common_Divisors;           use Standard_Common_Divisors;
with Standard_Integer_Linear_Solvers;    use Standard_Integer_Linear_Solvers;
with Face_Enumerators_Utilities;         use Face_Enumerators_Utilities;
with Integer_Linear_Inequalities;        use Integer_Linear_Inequalities;

package body Integer_Face_Enumerators is

-- ELIMINATE TO MAKE UPPER TRIANGULAR :

  procedure Eliminate ( pivot : in integer32; elim : in Vector;
                        target : in out Vector ) is

  -- DESCRIPTION :
  --   Makes target(pivot) = 0 by means of making a linear
  --   combination of the vectors target and elim.

  -- REQUIRED ON ENTRY :
  --   target(pivot) /= 0 and elim(pivot) /= 0

    a,b,lcmab,faca,facb : integer32;

  begin
    a := elim(pivot); b := target(pivot);
    lcmab := lcm(a,b);
    if lcmab < 0 then lcmab := -lcmab; end if;
    faca := lcmab/a;  facb := lcmab/b;
    if facb < 0
     then facb := -facb; faca := -faca;
    end if;
    for j in target'range loop
      target(j) := facb*target(j) - faca*elim(j);
    end loop;
    Scale(target);
  end Eliminate;

  procedure Eliminate ( l : in integer32; pivots : in Vector; 
                        elim : in VecVec; target : in out Vector ) is

  -- DESCRIPTION :
  --   Makes target(pivots(i)) = 0 by means of making a linear
  --   combination of the vectors target and elim(i), for i in 1..l.

  -- REQUIRED ON ENTRY :
  --   elim(i)(pivots(i)) /= 0

  begin
    for i in 1..l loop
      if target(pivots(i)) /= 0
       then Eliminate(pivots(i),elim(i).all,target);
      end if;
    end loop;
  end Eliminate;

  function Pivot_after_Elimination
             ( l,k : in integer32; point,face,pivots : in Vector;
               elim : in VecVec ) return integer32 is

  -- DESCRIPTION :
  --   Returns the first nonzero element of the given point after elimination
  --   w.r.t. the entries in the face with lower index.

    work : Vector(point'range) := point - elim(1-k).all;
    pivot : integer32;

  begin
    for i in (face'first+1)..face'last loop
      if (face(i) < l) and then (work(pivots(i)) /= 0)
       then Eliminate(pivots(i),elim(i).all,work);
      end if;
      exit when (face(i) > l);
    end loop;
    pivot := 0;
    for i in work'range loop
      if work(pivots(i)) /= 0
       then pivot := i;
      end if;
      exit when (pivot /= 0);
    end loop;
    return pivot;
  end Pivot_after_Elimination;

  procedure Update_Pivots ( point : in Vector; l : in integer32;
                            pivots : in out Vector; pivot : out integer32 ) is

  -- DESCRIPTION :
  --   Searches in the point(l..point'last) for the first nonzero entry
  --   and updates the pivoting information.

    temp,piv : integer32;

  begin
    piv := 0;
    for i in l..point'last loop
      if point(pivots(i)) /= 0
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

  procedure Update_Eliminator ( elim : in out VecVec; l : in integer32;
                                pivots : in out Vector;
                                point : in Vector; pivot : out integer32 ) is

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
  --   point     new point to make the equation `point - elim(0) = 0'.

  -- ON RETURN :
  --   elim      if not degen, then elim(l)(pivots(l)) /= 0 and
  --             for i in 1..(l-1): elim(l)(pivots(i)) = 0;
  --   pivots    updated pivot information;
  --   pivot     the pivot that has been used for elim(l)(pivots(l)) /= 0;
  --             piv = 0 when the new elimination equation elim(l)
  --             became entirely zero after ensuring the invariant condition.

  begin
    elim(l) := new Vector'(point - elim(0).all);
    Eliminate(l-1,pivots,elim,elim(l).all);
    Update_Pivots(elim(l).all,l,pivots,pivot);
  end Update_Eliminator;

  procedure Update_Eliminator_for_Sum
               ( elim : in out VecVec; l : in integer32;
                 pivots : in out Vector; point : in Vector;
                 index : in integer32; pivot : out integer32 ) is

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
  --   index     indicates the current polytope.

  -- ON RETURN :
  --   elim      if not degen, then elim(l)(pivots(l)) /= 0 and
  --             for i in 1..(l-1): elim(l)(pivots(i)) = 0;
  --   pivots    updated pivot information;
  --   piv       the pivot that has been used for elim(l)(pivots(l)) /= 0;
  --             piv = 0 when the new elimination equation elim(l)
  --             became entirely zero after ensuring the invariant condition.

  begin
    elim(l) := new Vector'(point - elim(1-index).all);
    Eliminate(l-1,pivots,elim,elim(l).all);
    Update_Pivots(elim(l).all,l,pivots,pivot);
  end Update_Eliminator_for_Sum;

-- CREATE TABLEAU OF INEQUALITIES :

  procedure Create_Tableau_for_Vertices
               ( i,n : in integer32; pts : in VecVec; tab : out matrix ) is

    column : integer32 := tab'first(2);

  begin
    for j in pts'range loop
      if j /= i then
        for k in pts(j)'range loop
          tab(k,column) := pts(j)(k);
        end loop;
        tab(tab'last(1),column) := 1;
        column := column + 1;
      end if;
    end loop;
    for k in pts(i)'range loop
      tab(k,tab'last(2)) := pts(i)(k);
    end loop;
    tab(tab'last(1),tab'last(2)) := 1;
  end Create_Tableau_for_Vertices;

  procedure Create_Tableau_for_Edges
               ( i,j,n,pivot : in integer32; elim : in Vector;
                 pts : in VecVec; tab : out matrix;
                 lastcol : out integer32; degenerate : out boolean ) is

    column : integer32 := 1;
    ineq : Vector(1..n);

  begin
    degenerate := false;
    for k in pts'range loop
      if (k/=i) and then (k/=j) then
        ineq := pts(k).all - pts(i).all;   -- compute inequality
        -- put("Inequality before eliminate : "); put(ineq); new_line;
        if ineq(pivot) /= 0                -- make ineq(pivot) = 0
         then Eliminate(pivot,elim,ineq);
        end if;
        -- put("Inequality after eliminate : "); put(ineq); new_line;
        if Is_Zero(ineq) then              -- check if degenerate
          if not In_Edge(pts(k).all,pts(i).all,pts(j).all)
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
          tab(tab'last(1),column) := 1;
          column := column + 1;
        end if;
      end if;
    end loop;
    for k in tab'first(1)..(tab'last(1)-1) loop  -- right hand side
      tab(k,tab'last(2)) := 0;
    end loop;
    tab(tab'last(1),tab'last(2)) := 1;
    lastcol := column-1;
  end Create_Tableau_for_Edges;

  procedure Create_Tableau_for_Faces
               ( k,n : in integer32; face,pivots : in Vector;
                 pts,elim : in VecVec; tab : out matrix;
                 lastcol : out integer32; degenerate : out boolean ) is

    column : integer32 := 1;
    ineq : Vector(1..n);

  begin
    degenerate := false;
    for l in pts'range loop
      if not Is_In(l,face) then
        ineq := pts(l).all - elim(0).all;   -- new inequality
        -- put("Inequality before : "); put(ineq); new_line;
        Eliminate(k,pivots,elim,ineq);      -- ineq(pivots(i)) = 0, i=1,2,..k
        -- put("and after elimination : "); put(ineq); new_line;
        if Is_Zero(ineq) then
          if --not In_Face(k,face,pts(l).all,pts)
             --and then 
             l < face(face'last)       -- lexicographic enumeration
            and then (Pivot_after_Elimination
                             (l,1,pts(l).all,face,pivots,elim) /= 0)
           then -- put("Degenerate for l = "); put(l,1); new_line;
                degenerate := true; return;
          end if;
        else
          for ll in (k+1)..n loop              -- fill in the column
            tab(ll-k,column) := ineq(pivots(ll));
          end loop;
          tab(tab'last(1),column) := 1;
          column := column + 1;
        end if;
      end if;
    end loop;
    for l in tab'first(1)..(tab'last(1)-1) loop  -- right hand side
      tab(l,tab'last(2)) := 0;
    end loop;
    tab(tab'last(1),tab'last(2)) := 1;
    lastcol := column-1;
  end Create_Tableau_for_Faces;

  procedure Create_Tableau_for_Faces_of_Sum
               ( k,n,i,r : in integer32; ind,pivots : in Vector;
                 pts,elim,face : in VecVec; tab : out matrix;
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
    ineq : Vector(1..n);
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
          Eliminate(k,pivots,elim,ineq);     -- ineq(pivots(i)) = 0, i=1,2,..k
          if Is_Zero(ineq) then
            if --not In_Face(face(l1)'length-1,face(l1).all,pts(l2).all,pts)
               --and then
                face(l1)(face(l1)'first) <= l2
              and then l2 < face(l1)(face(l1)'last)
                                            -- lexicographic enumeration
              and then (Pivot_after_Elimination
                           (l2,l1,pts(l2).all,face(l1).all,pivots,elim) /= 0)
             then degenerate := true; return;
            end if;
          else
            for ll in (k+1)..n loop                 -- fill in the column
              tab(ll-k,column) := ineq(pivots(ll));
            end loop;
            tab(tab'last(1),column) := 1;
            column := column + 1;
          end if;
        end if;
      end loop;
    end loop;
    for l in tab'first(1)..(tab'last(1)-1) loop  -- right hand side
      tab(l,tab'last(2)) := 0;
    end loop;
    tab(tab'last(1),tab'last(2)) := 1;
    lastcol := column-1;
  end Create_Tableau_for_Faces_of_Sum;

  procedure Create_Tableau_for_Lower_Edges
               ( i,j,n,pivot : in integer32; elim : in Vector;
                 pts : in VecVec; tab : out matrix;
                 lastcol : out integer32; degenerate : out boolean ) is

    column : integer32 := 1;
    ineq : Vector(1..n);

  begin
   -- put("The elimination vector : "); put(elim); new_line;
    degenerate := false;
    for k in pts'range loop
      if (k/=i) and then (k/=j) then
        ineq := pts(k).all - pts(i).all;   -- compute inequality
        -- put("Inequality before eliminate : "); put(ineq); new_line;
        if ineq(pivot) /= 0                -- make ineq(pivot) = 0
         then Eliminate(pivot,elim,ineq);
        end if;
        -- put("Inequality after eliminate : "); put(ineq); new_line;
        if Is_Zero(ineq) then                -- check if degenerate
          if not In_Edge(pts(k).all,pts(i).all,pts(j).all)
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
      tab(k,tab'last(2)) := 0;
    end loop;
    tab(tab'last(1),tab'last(2)) := -1;
    lastcol := column-1;
  end Create_Tableau_for_Lower_Edges;

  procedure Create_Tableau_for_Lower_Faces
               ( k,n : in integer32; face,pivots : in Vector;
                 pts,elim : in VecVec; tab : out matrix;
                 lastcol : out integer32; degenerate : out boolean ) is

    column : integer32 := 1;
    ineq : Vector(1..n);

  begin
    degenerate := false;
    for l in pts'range loop
      if not Is_In(l,face) then
        ineq := pts(l).all - elim(0).all;   -- new inequality
        -- put("Inequality before : "); put(ineq); new_line;
        Eliminate(k,pivots,elim,ineq);      -- ineq(pivots(i)) = 0, i=1,2,..k
        -- put("and after elimination : "); put(ineq); new_line;
        if Is_Zero(ineq) then
          if --not In_Face(k,face,pts(l).all,pts)
             --and then
               l < face(face'last)       -- lexicographic enumeration
            and then (Pivot_after_Elimination
                          (l,1,pts(l).all,face,pivots,elim) /= 0)
           then --put_line("Degenerate choice");
                degenerate := true; return;
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
      tab(l,tab'last(2)) := 0;
    end loop;
    tab(tab'last(1),tab'last(2)) := -1;
    lastcol := column-1;
  end Create_Tableau_for_Lower_Faces;

-- DETERMINE WHETHER THE CANDIDATE IS VERTEX, SPANS AN EDGE OR A K-FACE :

  function Is_Vertex ( i,m,n : integer32; pts : VecVec ) return boolean is

    tableau : matrix(1..n+1,1..m);
    feasible : boolean;

  begin
    Create_Tableau_for_Vertices(i,n,pts,tableau);
   -- put_line("The tableau :"); put(tableau);
    Integer_Complementary_Slackness(tableau,feasible);
    return not feasible;
  end Is_Vertex;

  function Is_Edge ( i,j,m,n : integer32; pts : VecVec ) return boolean is

    elim : Vector(1..n) := pts(i).all - pts(j).all;
    pivot : integer32;

  begin
    pivot := elim'first - 1;
    for k in elim'range loop
      if elim(k) /= 0
       then pivot := k;
      end if;
      exit when pivot >= elim'first;
    end loop;
    if pivot < elim'first then
      return false;
    else
      Scale(elim);
      declare
        tab : matrix(1..n,1..m-1);
        deg,feasible : boolean;
        lst : integer32;
      begin
        -- put("The elimination vector : "); put(elim); new_line;
        Create_Tableau_for_Edges(i,j,n,pivot,elim,pts,tab,lst,deg);
        -- put_line("The tableau :"); put(tab);
        if deg then
          return false;
        elsif lst = 0 then
          return true;
        else
          Integer_Complementary_Slackness(tab,lst,feasible);
          return not feasible;
        end if;
      end;
    end if;
  end Is_Edge;

  function Is_Face ( k,n,m : integer32; elim,pts : VecVec;
                     pivots,face : Vector ) return boolean is

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
  --   face          entries of the points which span the candidate face.

    nn : constant integer32 := n-k+1;
    tab : matrix(1..nn,1..(m-k));
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
      Create_Tableau_for_Faces(k,n,face,pivots,pts,elim,tab,lst,deg);
      -- put_line("The tableau of inequalities : "); put(tab);
      if deg then
        -- put_line("Tableau is degenerate: no solution :");
        return false;
      elsif lst = 0 then
        return true;
      else
        Integer_Complementary_Slackness(tab,lst,feasible);
        -- if feasible
        --  then put_line(" is feasible");
        --  else put_line(" is not feasible");
        -- end if;
        return not feasible;
      end if;
    end if;
  end Is_Face;

  function Is_Face_of_Sum
                ( k,n,m,i,r : integer32; elim,pts,face : VecVec;
                  ind,pivots : Vector ) return boolean is

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
  --   face       entries of the points which span the candidate face;
  --   ind        indicates the starting vector in pts of each polytope;
  --   pivots     pivoting information for the elimination equations.

    nn : constant integer32 := n-k+1;
    tab : matrix(1..nn,1..(m-k));
    deg,feasible : boolean;
    lst : integer32;

  begin
    if m - k <= 1 then
      return true;
    else
      Create_Tableau_for_Faces_of_Sum
        (k,n,i,r,ind,pivots,pts,elim,face,tab,lst,deg);
      -- put_line("The tableau of inequalities : "); put(tab);
      if deg then
        --put_line("Tableau is degenerate: no solution :");
        return false;
      elsif lst = 0 then
        return true;
      else
        Integer_Complementary_Slackness(tab,lst,feasible);
        --if feasible
        -- then put_line(" is feasible");
        -- else put_line(" is not feasible");
        --end if;
        return not feasible;
      end if;
    end if;
  end Is_Face_of_Sum;

  function Is_Lower_Edge ( i,j,m,n : integer32; pts : VecVec )
                         return boolean is

    elim : Vector(1..n) := pts(i).all - pts(j).all;
    pivot : integer32;

  begin
    pivot := elim'first - 1;
    for k in elim'range loop
      if elim(k) /= 0
       then pivot := k;
      end if;
      exit when pivot >= elim'first;
    end loop;
    if pivot < elim'first or else (pivot = elim'last) then
      return false;
    else
      Scale(elim);
      declare
        tab : matrix(1..n-1,1..m-1);
        deg,feasible : boolean;
        lst : integer32;
      begin
        Create_Tableau_for_Lower_Edges(i,j,n,pivot,elim,pts,tab,lst,deg);
        -- put_line("The tableau :"); put(tab);
        if deg then
          return false;
        elsif lst = 0 then
          return true;
        else
          Integer_Complementary_Slackness(tab,lst,feasible);
          return not feasible;
        end if;
      end;
    end if;
  end Is_Lower_Edge;

  function Is_Lower_Face 
                ( k,n,m : in integer32; elim,pts : VecVec;
                  pivots,face : Vector ) return boolean is

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
    tab : matrix(1..nn,1..(m-k));
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
      Create_Tableau_for_Lower_Faces(k,n,face,pivots,pts,elim,tab,lst,deg);
      -- put_line("The tableau of inequalities : "); put(tab);
      if deg then
        --put_line("Degenerate tableau");
        return false;
      elsif lst = 0 then
        --put_line("lst = 0");
        return true;
      else
        Integer_Complementary_Slackness(tab,lst,feasible);
        --if feasible
        -- then put_line(" is feasible");
        -- else put_line(" is not feasible");
        --end if;
        return not feasible;
      end if;
    end if;
  end Is_Lower_Face;

-- TARGET ROUTINES :

  procedure Enumerate_Vertices ( pts : in VecVec ) is

    continue : boolean := true;
    n : constant integer32 := pts(pts'first).all'length;
    m : constant integer32 := pts'length;

  begin
    for i in pts'range loop
      if Is_Vertex(i,m,n,pts)
       then Process(i,continue);
      end if;
      exit when not continue;
    end loop;
  end Enumerate_Vertices;

  procedure Enumerate_Edges ( pts : in VecVec ) is
  
    continue : boolean := true;
    n : constant integer32 := pts(pts'first).all'length;
    m : constant integer32 := pts'length;

    procedure Candidate_Edges ( i,n : in integer32 ) is
    begin
      for j in (i+1)..pts'last loop    -- enumerate all candidates
       -- put("Verifying "); put(pts(i).all); put(" and");
       -- put(pts(j).all); put(" :"); new_line;
        if Is_Edge(i,j,m,n,pts)
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

  procedure Enumerate_Lower_Edges ( pts : in VecVec ) is

    continue : boolean := true;
    n : constant integer32 := pts(pts'first).all'length;
    m : constant integer32 := pts'length;

    procedure Candidate_Lower_Edges ( i,n : in integer32 ) is
    begin
      for j in (i+1)..pts'last loop    -- enumerate all candidates
       -- put("Verifying "); put(pts(i).all); put(" and");
       -- put(pts(j).all); put(" :"); new_line;
        if Is_Lower_Edge(i,j,m,n,pts)
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

  procedure Enumerate_Faces ( k : in integer32; pts : in VecVec ) is

    m : constant integer32 := pts'length;
    n : constant integer32 := pts(pts'first).all'length;
    candidate : Vector(0..k);
    elim : VecVec(0..k);
    continue : boolean := true;

    procedure Candidate_Faces ( ipvt : in Vector; i,l : in integer32 ) is

    -- DESCRIPTION :
    --   Enumerates all candidate k-faces in lexicographic order.

      piv : integer32;
      pivots : Vector(1..n);

    begin
      if l > k then
        if (k = m)
           or else Is_Face(k,n,m,elim,pts,ipvt,candidate)
         then Process(candidate,continue);
        end if;
      else
        for j in (i+1)..pts'last loop
          candidate(l) := j;
          pivots := ipvt;
          Update_Eliminator(elim,l,pivots,pts(j).all,piv);
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
        ipvt : Vector(1..n);
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

  procedure Enumerate_Lower_Faces ( k : in integer32; pts : VecVec ) is

    m : constant integer32 := pts'length;
    n : constant integer32 := pts(pts'first).all'length;
    candidate : Vector(0..k);
    elim : VecVec(0..k);
    continue : boolean := true;

    procedure Candidate_Lower_Faces ( ipvt : in Vector; i,l : in integer32 ) is

    -- DESCRIPTION :
    --   Enumerates all candidate k-faces in lexicographic order.

      piv : integer32;
      pivots : Vector(1..n);

    begin
      if l > k then
        --put_line("Testing the following candidate face :");
        --for ii in candidate'first..candidate'last-1 loop
        --  put(pts(candidate(ii))); put(" & ");
        --end loop;
        --put(pts(candidate(candidate'last))); new_line;
        if (k = m) or else Is_Lower_Face(k,n,m,elim,pts,ipvt,candidate)
         then Process(candidate,continue);
        end if;
      else 
        for j in (i+1)..pts'last loop
          candidate(l) := j;
          pivots := ipvt;
          -- put("Picking "); put(pts(j));
          -- put(" Pivots : "); put(pivots); new_line;
          Update_Eliminator(elim,l,pivots,pts(j).all,piv);
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
        ipvt : Vector(1..n);
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
                  ( ind,typ : in Vector; k : in integer32; pts : in VecVec ) is

    m : constant integer32 := pts'length;             -- number of points
    n : constant integer32 := pts(pts'first)'length;  -- dimension of vectors
    r : constant integer32 := ind'length;             -- number of polytopes
    candidates : VecVec(1..r);
    elim : VecVec(1-r..k);
    pivots : Vector(1..n);
    continue : boolean := true;
    lasti,sum : integer32;

    procedure Candidate_Faces_of_Sum ( ipvt : in Vector; i : in integer32 ) is

    -- DESCRIPTION :
    --   Enumerates all faces of the given type on the sum,
    --   i indicates the current polytope.

      procedure Candidate_Faces ( ipvt : in Vector;
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
                         (sum,n,last-i+1,i,r,elim,pts(pts'first..last),
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
                     (elim,sum-typ(i)+l,pivots,pts(j).all,i,piv);
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
      candidates(i) := new Vector(0..typ(i));
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
      ipvt : Vector(1..n);
    begin
      for i in ipvt'range loop
        ipvt(i) := i;
      end loop;
      sum := 0;
      Candidate_Faces_of_Sum(ipvt,1);
    end;
  end Enumerate_Faces_of_Sum;

end Integer_Face_Enumerators;
