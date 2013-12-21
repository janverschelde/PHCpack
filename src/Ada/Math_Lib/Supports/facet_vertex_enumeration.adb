with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Norms_Equals;     use Standard_Floating_Norms_Equals;
with Linear_Minimization;                use Linear_Minimization;

package body Facet_Vertex_Enumeration is

-- AUXILIARIES :

  function Sort ( v : Standard_Integer_Vectors.Vector )
                return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Return the sorted vector v.

    res : Standard_Integer_Vectors.Vector(v'range) := v;
    min,ind : integer32;

  begin
    for i in res'first..res'last-1 loop
      ind := i;
      min := res(ind);
      for j in i+1..res'last loop
        if res(j) < min
         then ind := j; min := res(j);
        end if;
      end loop;
      if ind /= i
       then res(ind) := res(i); res(i) := min;
      end if;
    end loop;
    return res;
  end Sort;
 
  function Is_In ( L : Lists_of_Floating_Vectors.List;
                   v : Standard_Floating_Vectors.Vector;
                   tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Returns true if there is a vector in L that is equal to the given
  --   vector v within the given tolerance.

    use Standard_Floating_Vectors;
    use Lists_of_Floating_Vectors;
    tmp : List := L;
    lv : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lv := Head_Of(tmp);
      if Equal(lv.all,v,tol)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  procedure Append_Diff ( first,last : in out Lists_of_Floating_Vectors.List; 
                          v : in Standard_Floating_Vectors.Vector;
                          tol : in double_float ) is

  -- DESCRIPTION :
  --   Appends the vector v to the list first when all vectors in first
  --   are different from v up to the given tolerance.

    use Lists_of_Floating_Vectors;

  begin
    if not Is_In(first,v,tol)
     then Append(first,last,v);
    end if;
  end Append_Diff;

-- DATA MANIPULATION :

  function List_to_Matrix
             ( n : integer32; L : Lists_of_Floating_Vectors.List ) 
             return Standard_Floating_Matrices.Matrix is

    use Standard_Floating_Vectors;
    use Lists_of_Floating_Vectors;

    res : Standard_Floating_Matrices.Matrix(1..n,1..integer32(Length_Of(L)));
    tmp : List := L;
    lv : Link_to_Vector;

  begin
    for j in res'range(2) loop
      lv := Head_Of(tmp);
      for i in 1..n loop
        res(i,j) := lv(i);
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end List_to_Matrix;

  function List_to_Vector
             ( i : integer32; L : Lists_of_Floating_Vectors.List )
             return Standard_Floating_Vectors.Vector is

    use Standard_Floating_Vectors;
    use Lists_of_Floating_Vectors;

    res : Vector(1..integer32(Length_Of(L)));
    tmp : List := L;
    lv : Link_to_Vector;

  begin
    for j in res'range loop
      lv := Head_Of(tmp);
      res(j) := lv(i);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end List_to_Vector;

  function Matrix_to_List
             ( n : integer32; m : Standard_Floating_Matrices.Matrix )
             return Lists_of_Floating_Vectors.List is

    use Standard_Floating_Vectors;
    use Lists_of_Floating_Vectors;

    res,res_last : List;
    v : Vector(1..n);

  begin
    for j in m'range(2) loop
      for i in 1..n loop
        v(i) := m(i,j);
        Append(res,res_last,v);
      end loop;
    end loop;
    return res;
  end Matrix_to_List;

-- AUXILIARY TO Enumerate_Facets :

  function Barycenter ( pts : Standard_Floating_Matrices.Matrix )
                      return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the barycenter of the points in the columns of the matrix.

    use Standard_Floating_Vectors;

    res : Vector(pts'range(1)) := (pts'range(1) => 0.0);
    m : constant integer32 := pts'last(2);
   -- sum : double_float := 0.0;

  begin
    for i in res'range loop
      for j in 1..m loop
        res(i) := res(i) + pts(i,j);
      end loop;
      res(i) := res(i)/double_float(m);
    end loop;
    return res;
  end Barycenter;

  procedure Facet_Enumeration_Model
               ( pts : in Standard_Floating_Matrices.Matrix;
                 cff : out Standard_Floating_Matrices.Matrix;
                 cost,rhs : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   This procedure returns the linear minimization model to enumerate
  --   all facets of a given point configuration.
  --   The facets of the given point configuration are in one-to-one
  --   correspondence with the facets of a cone with apex the barycenter
  --   of the point configuration.  The original points are placed
  --   at height one and the barycenter is pulled down to level zero.
  --   The facets of the cone that are of interest are lower facets,
  --   with last component of the inner normal equal to one.
  --   So we can eliminate one unknown of this extended configuration
  --   so that the height one gets into the right-hand side vector.
  --   Also here the cost vector is a randomly chosen vector.

  -- ON ENTRY :
  --   pts       coordinates of the points are in the columns of the matrix.

  -- ON RETURN :
  --   cff       coefficients of the inequalities;
  --   cost      random cost vector;
  --   rhs       right-hand side vector, all entries equal to one.

    center : constant Standard_Floating_Vectors.Vector := Barycenter(pts);

  begin
    for i in pts'range(1) loop
      for j in pts'range(2) loop
        cff(i,j) := center(i)-pts(i,j);
      end loop;
    end loop;
    for i in cost'range loop
      cost(i) := abs(Random);
    end loop;
    rhs := (rhs'range => -1.0);
  end Facet_Enumeration_Model;

  procedure Lower_Facet_Enumeration_Model
               ( pts : in Standard_Floating_Matrices.Matrix;
                 cff : out Standard_Floating_Matrices.Matrix;
                 cost,rhs : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Returns the linear minimization model to enumerate lower facets.
  --   The points are embedded into projective space by adding one additional
  --   coordinate which is equal to one for each point.
  --   Because we are looking for lower facets, the inner normal must have
  --   last component positive, or without loss of generality, equal to one.
  --   By elimination the dimension equals again the original dimension.
  --   The cost function is here not chosen at random.
 
  -- ON ENTRY :
  --   pts       coordinates of the points are in the columns of the matrix.

  -- ON RETURN :
  --   cff       coefficients of the inequalities;
  --   cost      random cost vector;
  --   rhs       right-hand side vector, all entries equal to one.

    n : constant integer32 := pts'last(1);
    center : constant Standard_Floating_Vectors.Vector := Barycenter(pts);

  begin
    for i in pts'range(1) loop
      for j in pts'range(2) loop
        cff(i,j) := center(i)-pts(i,j);
      end loop;
    end loop;
   -- for i in cost'range loop
   --   cost(i) := abs(Random);
   -- end loop;
    rhs := (rhs'range => -1.0);
   -- for i in 1..n-1 loop
   --   for j in pts'range(2) loop
   --     cff(i,j) := pts(i,j);
   --   end loop;
   -- end loop;
   -- for j in pts'range(2) loop
   --   cff(n,j) := 1.0;
   -- end loop;
   -- for j in pts'range(2) loop
   --   rhs(j) := -abs(Random); -- -pts(n,j);  
   -- end loop;
    cost(1..n-1) := (1..n-1 => 0.0);
    cost(n) := -1.0;
  end Lower_Facet_Enumeration_Model;

-- ENUMERATORS OF VERTICES :

  procedure Enumerate_Vertices
               ( cff : in Standard_Floating_Matrices.Matrix;
                 rhs : in Standard_Floating_Vectors.Vector;
                 tol : in double_float;
                 points : out Lists_of_Floating_Vectors.List;
                 labels : out Lists_of_Integer_Vectors.List;
                 fail,infty : out boolean ) is

  -- NOTE ON THE OPTIMIZATION MODEL :
  --   The cost vector must be randomly chosen to ensure a unique sink.

    use Lists_of_Floating_Vectors;

    pts,pts_last : Lists_of_Floating_Vectors.List;
    lbl,lbl_last : Lists_of_Integer_Vectors.List;
    n : constant integer32 := cff'last(1);
    m : constant integer32 := cff'last(2);
    cost : Standard_Floating_Vectors.Vector(1..n);

    procedure Collect ( binv : in Standard_Floating_Matrices.Matrix;
                        active : in Standard_Integer_Vectors.Vector;
                        sol : in Standard_Floating_Vectors.Vector;
                        continue : out boolean ) is
    begin
      Append_Diff(pts,pts_last,sol,tol);
      Lists_of_Integer_Vectors.Append_Diff(lbl,lbl_last,Sort(active));
      continue := true;
    end Collect;
    procedure Collect_Vertices is new Enumerate_Feasible_Vertices(Collect);

  begin
    for i in 1..n loop
      cost(i) := abs(Random);
    end loop;
    Collect_Vertices(n,m,cff,cost,rhs,tol,fail,infty);
    points := pts;
    labels := lbl;
  end Enumerate_Vertices;

  function Enumerate_Vertex_Points
               ( cff : Standard_Floating_Matrices.Matrix;
                 rhs : Standard_Floating_Vectors.Vector; tol : double_float )
               return Lists_of_Floating_Vectors.List is

  -- NOTE ON THE OPTIMIZATION MODEL :
  --   The cost vector must be randomly chosen to ensure a unique sink.

    use Lists_of_Floating_Vectors;

    res,res_last : List;
    n : constant integer32 := cff'last(1);
    m : constant integer32 := cff'last(2);
    cost : Standard_Floating_Vectors.Vector(1..n);
    fail,infty : boolean;

    procedure Collect ( binv : in Standard_Floating_Matrices.Matrix;
                        active : in Standard_Integer_Vectors.Vector;
                        sol : in Standard_Floating_Vectors.Vector;
                        continue : out boolean ) is
    begin
      Append_Diff(res,res_last,sol,tol);
      continue := true;
    end Collect;
    procedure Collect_Vertices is new Enumerate_Feasible_Vertices(Collect);

  begin
    for i in 1..n loop
      cost(i) := abs(Random);
    end loop;
    Collect_Vertices(n,m,cff,cost,rhs,tol,fail,infty);
    return res;
  end Enumerate_Vertex_Points;

  function Enumerate_Vertex_Labels
               ( cff : Standard_Floating_Matrices.Matrix;
                 rhs : Standard_Floating_Vectors.Vector; tol : double_float )
               return Lists_of_Integer_Vectors.List is

  -- NOTE ON THE OPTIMIZATION MODEL :
  --   The cost vector must be randomly chosen to ensure a unique sink.

    use Lists_of_Integer_Vectors;

    res,res_last : List;
    n : constant integer32 := cff'last(1);
    m : constant integer32 := cff'last(2);
    cost : Standard_Floating_Vectors.Vector(1..n);
    fail,infty : boolean;

    procedure Collect ( binv : in Standard_Floating_Matrices.Matrix;
                        active : in Standard_Integer_Vectors.Vector;
                        sol : in Standard_Floating_Vectors.Vector;
                        continue : out boolean ) is
    begin
      Append_Diff(res,res_last,Sort(active));
      continue := true;
    end Collect;
    procedure Collect_Vertices is new Enumerate_Feasible_Vertices(Collect);

  begin
    for i in 1..n loop
      cost(i) := abs(Random);
    end loop;
    Collect_Vertices(n,m,cff,cost,rhs,tol,fail,infty);
    return res;
  end Enumerate_Vertex_Labels;

-- ENUMERATORS OF FACETS :

  procedure Enumerate_Facets
               ( pts : in Standard_Floating_Matrices.Matrix;
                 tol : in double_float;
                 facets : out Lists_of_Floating_Vectors.List;
                 labels : out Lists_of_Integer_Vectors.List;
                 fail,infty : out boolean ) is

    fcs,fcs_last : Lists_of_Floating_Vectors.List;
    lbl,lbl_last : Lists_of_Integer_Vectors.List;
    cff : Standard_Floating_Matrices.Matrix(pts'range(1),pts'range(2));
    cost : Standard_Floating_Vectors.Vector(pts'range(1));
    rhs : Standard_Floating_Vectors.Vector(pts'range(2));

    procedure Collect ( binv : in Standard_Floating_Matrices.Matrix;
                        active : in Standard_Integer_Vectors.Vector;
                        sol : in Standard_Floating_Vectors.Vector;
                        continue : out boolean ) is

      extsol : Standard_Floating_Vectors.Vector(sol'first..sol'last+1);

    begin
      extsol(sol'range) := sol;
      extsol(sol'last+1) := Eval(sol'last,active(1),pts,sol);
      Append_Diff(fcs,fcs_last,extsol,tol);
      Lists_of_Integer_Vectors.Append_Diff(lbl,lbl_last,Sort(active));
      continue := true;
    end Collect;
    procedure Collect_Facets is new Enumerate_Feasible_Vertices(Collect);

  begin
    Facet_Enumeration_Model(pts,cff,cost,rhs);
    Collect_Facets(pts'last(1),pts'last(2),cff,cost,rhs,tol,fail,infty);
    facets := fcs;
    labels := lbl;
  end Enumerate_Facets;

  function Enumerate_Facet_Inequalities
               ( pts : Standard_Floating_Matrices.Matrix; tol : double_float )
               return Lists_of_Floating_Vectors.List is

    res,res_last : Lists_of_Floating_Vectors.List;
    cff : Standard_Floating_Matrices.Matrix(pts'range(1),pts'range(2));
    cost : Standard_Floating_Vectors.Vector(pts'range(1));
    rhs : Standard_Floating_Vectors.Vector(pts'range(2));
    fail,infty : boolean;

    procedure Collect ( binv : in Standard_Floating_Matrices.Matrix;
                        active : in Standard_Integer_Vectors.Vector;
                        sol : in Standard_Floating_Vectors.Vector;
                        continue : out boolean ) is

      extsol : Standard_Floating_Vectors.Vector(sol'first..sol'last+1);

    begin
      extsol(sol'range) := sol;
      extsol(sol'last+1) := Eval(sol'last,active(1),pts,sol);
      Append_Diff(res,res_last,extsol,tol);
      continue := true;
    end Collect;
    procedure Collect_Facets is new Enumerate_Feasible_Vertices(Collect);

  begin
    Facet_Enumeration_Model(pts,cff,cost,rhs);
    Collect_Facets(pts'last(1),pts'last(2),cff,cost,rhs,tol,fail,infty);
    return res;
  end Enumerate_Facet_Inequalities;

  function Enumerate_Facet_Labels
               ( pts : Standard_Floating_Matrices.Matrix; tol : double_float )
               return Lists_of_Integer_Vectors.List is

    res,res_last : Lists_of_Integer_Vectors.List;
    cff : Standard_Floating_Matrices.Matrix(pts'range(1),pts'range(2));
    cost : Standard_Floating_Vectors.Vector(pts'range(1));
    rhs : Standard_Floating_Vectors.Vector(pts'range(2));
    fail,infty : boolean;

    procedure Collect ( binv : in Standard_Floating_Matrices.Matrix;
                        active : in Standard_Integer_Vectors.Vector;
                        sol : in Standard_Floating_Vectors.Vector;
                        continue : out boolean ) is
    begin
      Lists_of_Integer_Vectors.Append_Diff(res,res_last,Sort(active));
      continue := true;
    end Collect;
    procedure Collect_Facets is new Enumerate_Feasible_Vertices(Collect);

  begin
    Facet_Enumeration_Model(pts,cff,cost,rhs);
    Collect_Facets(pts'last(1),pts'last(2),cff,cost,rhs,tol,fail,infty);
    return res;
  end Enumerate_Facet_Labels;

-- ENUMERATORS OF LOWER FACETS :

  procedure Enumerate_Lower_Facets
               ( pts : in Standard_Floating_Matrices.Matrix;
                 tol : in double_float;
                 facets : out Lists_of_Floating_Vectors.List;
                 labels : out Lists_of_Integer_Vectors.List;
                 fail,infty : out boolean ) is

    fcs,fcs_last : Lists_of_Floating_Vectors.List;
    lbl,lbl_last : Lists_of_Integer_Vectors.List;
    n : constant integer32 := pts'last(1);
    m : constant integer32 := pts'last(2);
    cff : Standard_Floating_Matrices.Matrix(1..n,1..m);
    cost : Standard_Floating_Vectors.Vector(1..n);
    rhs : Standard_Floating_Vectors.Vector(1..m);

    procedure Collect ( binv : in Standard_Floating_Matrices.Matrix;
                        active : in Standard_Integer_Vectors.Vector;
                        sol : in Standard_Floating_Vectors.Vector;
                        continue : out boolean ) is

      extsol : Standard_Floating_Vectors.Vector(1..n+1);

    begin
     -- extsol(1..n-1) := sol;
     -- extsol(n) := 1.0;
     -- extsol(n+1) := Eval(n,active(1),pts,extsol(1..n));
     -- Append_Diff(fcs,fcs_last,extsol,tol);
     -- Lists_of_Integer_Vectors.Append_Diff(lbl,lbl_last,Sort(active));
     -- continue := true;
      if sol(sol'last) > tol then
        extsol(sol'range) := sol;
        extsol(n+1) := Eval(n,active(1),pts,sol);
        Append_Diff(fcs,fcs_last,extsol,tol);
        Lists_of_Integer_Vectors.Append_Diff(lbl,lbl_last,Sort(active));
        continue := true;
      else
        continue := false;
      end if;
    end Collect;
    procedure Collect_Lower_Facets is new Enumerate_Feasible_Vertices(Collect);

  begin
    Lower_Facet_Enumeration_Model(pts,cff,cost,rhs);
    Collect_Lower_Facets(n,m,cff,cost,rhs,tol,fail,infty);
    facets := fcs;
    labels := lbl;
  end Enumerate_Lower_Facets;

  function Enumerate_Lower_Facet_Inequalities
               ( pts : Standard_Floating_Matrices.Matrix; tol : double_float )
               return Lists_of_Floating_Vectors.List is

    res,res_last : Lists_of_Floating_Vectors.List;
    n : constant integer32 := pts'last(1);
    m : constant integer32 := pts'last(2);
    cff : Standard_Floating_Matrices.Matrix(1..n,1..m);
    cost : Standard_Floating_Vectors.Vector(1..n);
    rhs : Standard_Floating_Vectors.Vector(1..m);
    fail,infty : boolean;

    procedure Collect ( binv : in Standard_Floating_Matrices.Matrix;
                        active : in Standard_Integer_Vectors.Vector;
                        sol : in Standard_Floating_Vectors.Vector;
                        continue : out boolean ) is

      extsol : Standard_Floating_Vectors.Vector(1..n+1);

    begin
     -- extsol(1..n-1) := sol;
     -- extsol(n) := 1.0;
     -- extsol(n+1) := Eval(n,active(1),pts,extsol(1..n));
     -- Append_Diff(res,res_last,extsol,tol);
     -- continue := true;
      if sol(sol'last) > tol
       then extsol(sol'range) := sol;
            extsol(n+1) := Eval(n,active(1),pts,sol);
            Append_Diff(res,res_last,extsol,tol);
            continue := true;
       else continue := false;
      end if;
    end Collect;
    procedure Collect_Lower_Facets is new Enumerate_Feasible_Vertices(Collect);

  begin
    Lower_Facet_Enumeration_Model(pts,cff,cost,rhs);
    Collect_Lower_Facets(n,m,cff,cost,rhs,tol,fail,infty);
    return res;
  end Enumerate_Lower_Facet_Inequalities;

  function Enumerate_Lower_Facet_Labels
               ( pts : Standard_Floating_Matrices.Matrix; tol : double_float )
               return Lists_of_Integer_Vectors.List is

    res,res_last : Lists_of_Integer_Vectors.List;
    n : constant integer32 := pts'last(1);
    m : constant integer32 := pts'last(2);
    cff : Standard_Floating_Matrices.Matrix(1..n,1..m);
    cost : Standard_Floating_Vectors.Vector(1..n);
    rhs : Standard_Floating_Vectors.Vector(1..m);
    fail,infty : boolean;

    procedure Collect ( binv : in Standard_Floating_Matrices.Matrix;
                        active : in Standard_Integer_Vectors.Vector;
                        sol : in Standard_Floating_Vectors.Vector;
                        continue : out boolean ) is
    begin
     -- continue := true;
      if sol(sol'last) > tol
       then Lists_of_Integer_Vectors.Append_Diff(res,res_last,Sort(active));
            continue := true;
       else continue := false;
      end if;
    end Collect;
    procedure Collect_Lower_Facets is new Enumerate_Feasible_Vertices(Collect);

  begin
    Lower_Facet_Enumeration_Model(pts,cff,cost,rhs);
    Collect_Lower_Facets(n,m,cff,cost,rhs,tol,fail,infty);
    return res;
  end Enumerate_Lower_Facet_Labels;

end Facet_Vertex_Enumeration;
