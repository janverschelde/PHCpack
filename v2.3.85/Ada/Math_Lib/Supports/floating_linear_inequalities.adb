with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Givens_Rotations;                   use Givens_Rotations;
with Linear_Minimization;

package body Floating_Linear_Inequalities is

-- USAGE OF LP :

  procedure Is_In_Cone ( tab : in out Matrix; lastcol : in integer32;
                         rhs : in out Standard_Floating_Vectors.Vector;
                         tol : in double_float;
                         solution : out Standard_Floating_Vectors.Vector;
                         columns : out Standard_Integer_Vectors.Vector;
                         feasible : out boolean ) is

  -- ALGORITHM:
  --   The following linear program will be solved:
  --
  --    min u_1 + .. + u_n
  --
  --        l_1*p_1i + l_2*p_2i + .. + l_m*p_mi + u_i*q_i = q_i   i=1,2,..,n
  --
  --  to determine whether q belongs to the cone spanned by the
  --  vectors p_1,p_2,..,p_m.
  --  If all u_i are zero and all constraints are satisfied,
  --  then q belongs to the cone spanned by the vectors.

  -- CONSTANTS :

    n  : constant integer32 := rhs'last;
    m  : constant integer32 := 2*n;             -- number of constraints
    nb : constant integer32 := lastcol+n;       -- number of unknowns

  -- VARIABLES :

    mat : matrix(0..m,0..nb);
    sol : Standard_Floating_Vectors.Vector(1..nb) := (1..nb => 0.0);
   -- bas : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);
   -- slk : Standard_Integer_Vectors.Vector(1..nb) := (1..nb => 0);

    cnt : integer32;
    s : double_float;

  begin

   -- INITIALIZATION OF THE TARGET :

    for i in 0..(nb-n) loop
      mat(0,i) := 0.0;              -- sum of the lambda's
    end loop;
    for i in (nb-n+1)..nb loop
      mat(0,i) := -1.0;             -- sum of the mu's
    end loop;

   -- INITIALIZATION OF THE CONSTRAINTS :

    for i in 1..n loop
      for j in (nb-n+1)..nb loop
        if i /= j-nb+n
         then mat(i,j)   := 0.0;
              mat(i+n,j) := 0.0;
        end if;
      end loop;
    end loop;

    for j in tab'first(2)..lastcol loop
      for i in tab'range(1) loop
        mat(i,j)   := -tab(i,j);
        mat(i+n,j) :=  tab(i,j);
      end loop;
    end loop;

    for i in rhs'range loop
      mat(i,0)   := -rhs(i);
      mat(i+n,0) :=  rhs(i);
      mat(i,i+nb-n)   := -rhs(i);
      mat(i+n,i+nb-n) :=  rhs(i);
    end loop;

  -- SOLVE THE LINEAR PROGRAM :

   -- New_Tableau.Init_Basis(nb,bas,slk);
   -- Dual_Simplex(mat,sol,bas,slk,tol,false);

  -- RETURN AND CHECK THE SOLUTION :

    cnt := 0;
    for i in tab'first(2)..lastcol loop
      if abs(sol(i)) > tol then
        cnt := cnt + 1;
        solution(cnt) := sol(i);
        columns(cnt) := i;
      end if;
    end loop;

    s := 0.0;
    for i in (nb-n+1)..nb loop
      s := s + sol(i);
    end loop;
    if abs(s) > tol
     then feasible := false; return;
    end if;
    feasible := true;

  end Is_In_Cone;

-- AUXILIARIES :

  function Index_of_Maximum ( v : Standard_Floating_Vectors.Vector;
                              lst : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the index in the vector v(v'first..lst) for which the maximum
  --   is obtained.

    res : integer32 := v'first;
    max : double_float := v(res);

  begin
    for i in v'first+1..lst loop
      if v(i) > max
       then max := v(i); res := i;
      end if;
    end loop;
    return res;
  end Index_of_Maximum;

  function Norm ( v : Standard_Floating_Vectors.Vector ) return double_float is

  -- DESCRIPTION :
  --   Returns the 2-norm of the vector.

    res : double_float := 0.0;

  begin
    for j in v'range loop
      res := res + v(j)*v(j);
    end loop;
    if res + 1.0 = 1.0
     then return 0.0;
     else return SQRT(res);
    end if;
  end Norm;

  function Norm ( v : Standard_Floating_Vectors.vector;
                  frst : integer32 ) return double_float is

  -- DESCRIPTION :
  --   Returns the 2-norm of the vector v(frst..v'last).

    res : double_float := 0.0;
 
  begin
    for j in frst..v'last loop
      res := res + v(j)*v(j);
    end loop;
    if res + 1.0 = 1.0
     then return 0.0;
     else return SQRT(res);
    end if;
  end Norm;

  procedure Scale ( v : in out Standard_Floating_Vectors.vector ) is

  -- DESCRIPTION :
  --   Returns the normed vector, i.e. v/Norm(v).

    nrm : double_float := Norm(v);

  begin
    if nrm + 1.0 /= 1.0 then
      for i in v'range loop
        v(i) := v(i)/nrm;
      end loop;
    end if;
  end Scale;

  function Norm ( mat : matrix; i : integer32 ) return double_float is

  -- DESCRIPTION :
  --   Returns the 2-norm of the ith column in the matrix.

    res : double_float := 0.0;

  begin
    for j in mat'range(1) loop
      res := res + mat(j,i)*mat(j,i);
    end loop;
    if res + 1.0 = 1.0
     then return 0.0;
     else return SQRT(res);
    end if;
  end Norm;

  procedure Scale ( mat : in out matrix; lastcolumn : in integer32 ) is

  -- DESCRIPTION :
  --   Scales the columns in the matrix, by dividing by their norm.

    nrm : double_float;

  begin
    for i in mat'first(2)..lastcolumn loop
      nrm := Norm(mat,i);
      if nrm + 1.0 /= 1.0 then
        for j in mat'range(1) loop
          mat(j,i) := mat(j,i)/nrm;
        end loop;
      end if;
    end loop;
  end Scale;

  function Inner_Product ( mat : Matrix; i : integer32;
                           v : Standard_Floating_Vectors.Vector )
                         return double_float is
  -- DESCRIPTION :
  --   Computes the inner product of the given vector v and the vector
  --   in the ith column of the matrix.

    res : double_float := 0.0;

  begin
    for k in mat'range(1) loop
      res := res + mat(k,i)*v(k);
    end loop;
    return res;
  end Inner_Product;

  function Maximal_Product ( mat : matrix; i : integer32;
                             v : Standard_Floating_Vectors.Vector;
                             frst : integer32 ) return integer32 is
  -- DESCRIPTION :
  --   Returns the index in v(frst..v'last) for which mat(index,i)*v(index)
  --   becomes maximal.

    res : integer32 := frst;
    max : double_float := mat(res,i)*v(res);
    tmp : double_float;

  begin
    for j in frst+1..v'last loop
      tmp := mat(j,i)*v(j);
      if tmp > max
       then max := tmp; res := j;
      end if;
    end loop;
    return res;
  end Maximal_Product;

  function Inner_Product
             ( mat : Matrix; i,j : integer32 ) return double_float is

  -- DESCRIPTION :
  --   Computes the inner product of the vectors in the columns i and j
  --   of the matrix.

    res : double_float := 0.0;

  begin
    for k in mat'range(1) loop
      res := res + mat(k,i)*mat(k,j);
    end loop;
    return res;
  end Inner_Product;

  function Inner_Products ( mat : matrix; lastcol : integer32;
                            v : Standard_Floating_Vectors.Vector )
                          return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector with all inner products of the given vector
  --   and the column vectors of the matrix.

    res : Standard_Floating_Vectors.Vector(mat'first(2)..lastcol);

  begin
    for i in res'range loop
      res(i) := Inner_Product(mat,i,v);
    end loop;
    return res;
  end Inner_Products;

  function Positive ( v : Standard_Floating_Vectors.Vector;
                      lst : integer32; tol : double_float )
                    return boolean is
  -- DESCRIPTION :
  --   Returns true if all elements in v(v'first..lst) are >= 0.

  begin
    for i in v'first..lst loop
      if abs(v(i)) > tol and then (v(i) < 0.0)
       then return false;
      end if;
    end loop;
    return true;
  end Positive;

  function Positive
               ( mat : matrix; rhs : Standard_Floating_Vectors.Vector;
                 ind,col : integer32;
                 solind : double_float; tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   The given matrix is upper triangular up to the column indicated
  --   by ind.  This function returns true when by replacing column ind
  --   by column col, a positive solution would be obtained.
  --   The number solind is sol(ind) in the solution vector.

    sol : Standard_Floating_Vectors.Vector(1..ind) := (1..ind => 0.0);

  begin
    sol(ind) := solind;
    for i in reverse mat'first(1)..(ind-1) loop
      for j in i+1..(ind-1) loop
        sol(i) := sol(i) + mat(i,j)*sol(j);
      end loop;
      sol(i) := sol(i) + mat(i,col)*sol(ind);
      sol(i) := (rhs(i) - sol(i))/mat(i,i);
      if (abs(sol(i)) > tol) and then (sol(i) < 0.0)
       then --put_line("The solution before returning false : ");
            --put(sol,3,3,3); new_line;
            return false;
      end if;
    end loop;
    --put_line("The solution before returning true : ");
    --put(sol,3,3,3); new_line;
    return true;
  end Positive; 

  function Positive_Diagonal
               ( mat : Matrix; rhs : Standard_Floating_Vectors.Vector;
                 ind,col : integer32;
                 solind : double_float; tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   The given matrix is diagonal up to the column indicated
  --   by ind.  This function returns true when by replacing column ind
  --   by column col, a positive solution would be obtained.
  --   The number solind is sol(ind) in the solution vector.

    sol : Standard_Floating_Vectors.Vector(1..ind) := (1..ind => 0.0);

  begin
    sol(ind) := solind;
    for i in reverse mat'first(1)..(ind-1) loop
      sol(i) := (rhs(i) - mat(i,col)*sol(ind))/mat(i,i);
      if (abs(sol(i)) > tol) and then (sol(i) < 0.0)
       then --put_line("The solution before returning false : ");
            --put(sol,3,3,3); new_line;
            return false;
      end if;
    end loop;
    --put_line("The solution before returning true : ");
    --put(sol,3,3,3); new_line;
    return true;
  end Positive_Diagonal;

-- PIVOTING PROCEDURES :

  procedure Interchange ( v : in out Standard_Floating_Vectors.Vector;
                          i,j : in integer32 ) is

  -- DESCRIPTION :
  --   The entries v(i) and v(j) are interchanged.

    temp : double_float;

  begin
    temp := v(i); v(i) := v(j); v(j) := temp;
  end Interchange;

  procedure Interchange_Columns ( mat : in out matrix; i,j : in integer32 ) is

  -- DESCRIPTION :
  --   The columns i and j of the given matrix are interchanged.

    temp : double_float;

  begin
    for k in mat'range(1) loop
      temp := mat(k,i); mat(k,i) := mat(k,j); mat(k,j) := temp;
    end loop;
  end Interchange_Columns;

  procedure Interchange_Rows
                 ( mat : in out Matrix; lastcol : in integer32;
                   v : in out Standard_Floating_Vectors.Vector;
                   i,j : in integer32 ) is

  -- DESCRIPTION :
  --   The rows i and j of the given matrix and vector are interchanged.

    temp : double_float;

  begin
    for k in mat'first(2)..lastcol loop
      temp := mat(i,k); mat(i,k) := mat(j,k); mat(j,k) := temp;
    end loop;
    temp := v(i); v(i) := v(j); v(j) := temp;
  end Interchange_Rows;

  procedure Rotate ( mat : in out Matrix; lastcol,i : in integer32;
                     v : in out Standard_Floating_Vectors.Vector;
                     tol : in double_float ) is

  -- DESCRIPTION :
  --   Performs Givens rotations on the matrix and vector such that 
  --   mat(j,i) = 0, for all j > i.

  -- NOTE :
  --   A diagonalization procedure constructs an orthonormal basis
  --   and does not preserve the angles between the vectors.

    cosi,sine : double_float;

  begin
    for j in i+1..mat'last(1) loop
      if abs(mat(j,i)) > tol then
        Givens_Factors(mat,i,j,cosi,sine);
        Givens_Rotation(mat,lastcol,i,j,cosi,sine);
        Givens_Rotation(v,i,j,cosi,sine);
      end if;
    end loop;
  end Rotate;

  procedure Upper_Triangulate ( mat : in out Matrix; lastcol,i : in integer32;
                                v : in out Standard_Floating_Vectors.Vector;
                                tol : in double_float ) is

  -- DESCRIPTION :
  --   Makes the upper diagonal elements of the ith colume zero:
  --   mat(j,i) = 0, for all j /= i.

    fac : double_float;

  begin
    for j in mat'first(1)..(i-1) loop
      if (abs(mat(j,i)) > tol) then
        fac := mat(j,i)/mat(i,i);
        for k in i..lastcol loop
          mat(j,k) := mat(j,k) - fac*mat(i,k);
        end loop;
      end if;
    end loop;
  end Upper_Triangulate;

-- INITIALIZE : COMPUTE CLOSEST VECTOR TO RHS

  procedure Initialize
              ( tableau : in out Matrix; lastcol : in integer32;
                rhs,sol : in out Standard_Floating_Vectors.Vector;
                tol : in double_float;
                columns : in out Standard_Integer_Vectors.Vector;
                cosines : out Standard_Floating_Vectors.Vector;
                feasible : out boolean ) is

  -- DESCRIPTION :
  --   Scales the tableau and right hand side by dividing by its norm. 
  --   Then computes the vector of cosines.
  --   If no vector with positive cosine with the right hand side vector
  --   can be found, then the problem is not feasible, otherwise,
  --   the vector with largest positive cosine will be the first column.
  --   At last, this first column will be transformed to the unit vector,
  --   by means of Givens rotations.
  --   Then the first solution component equals the first component of
  --   the transformed right hand side vector.

    ips : Standard_Floating_Vectors.Vector(tableau'first(2)..lastcol);
    index : integer32;

  begin
    Scale(tableau,lastcol); Scale(rhs);
   -- put_line("The tableau after scaling : "); put(tableau,3,3,3);
   -- put_line(" with scaled right hand side : "); put(rhs,3,3,3); new_line;
    ips := Inner_Products(tableau,lastcol,rhs);
    index := Index_of_Maximum(ips,lastcol);
    if (abs(ips(index)) < tol) or else (ips(index) < 0.0) then
      if index <= rhs'last
       then feasible := (abs(rhs(index)) < tol);
       else feasible := false;
      end if;
    else
      feasible := true;
      columns(columns'first) := index;
      columns(index) := columns'first;
      Interchange_Columns(tableau,tableau'first(2),index);
      Interchange(ips,tableau'first(2),index);
      index := Maximal_Product(tableau,tableau'first(2),rhs,rhs'first);
      if index /= tableau'first(2)
       then Interchange_Rows(tableau,lastcol,rhs,tableau'first(2),index);
      end if;
      Rotate(tableau,lastcol,tableau'first(2),rhs,tol);
      -- Upper_Triangulate(tableau,lastcol,tableau'first(2),rhs,tol);
      sol(sol'first) := rhs(rhs'first);
    end if;
    cosines := ips;
   -- put_line("At the end of Initialize : ");
   -- put(" cosines : "); put(ips,3,3,3); new_line;
   -- put_line("The tableau : "); put(tableau,3,3,3);
   -- put_line(" with right hand side : "); put(rhs,3,3,3); new_line;
   -- put(" and last column "); put(lastcol,1); new_line;
  end Initialize;

  procedure Initialize
              ( tableau : in out Matrix; lastcol : in integer32;
                rhs,sol : in out Standard_Floating_Vectors.Vector;
                tol : in double_float;
                columns : in out Standard_Integer_Vectors.Vector;
                feasible : out boolean ) is

  -- DESCRIPTION :
  --   Scales the tableau and right hand side by dividing by its norm. 
  --   Then computes the vector of cosines.
  --   If no vector with positive cosine with the right hand side vector
  --   can be found, then the problem is not feasible, otherwise,
  --   the vector with largest positive cosine will be the first column.
  --   At last, this first column will be transformed to the unit vector,
  --   by means of Givens rotations.
  --   Then the first solution component equals the first component of
  --   the transformed right hand side vector.

    ips : Standard_Floating_Vectors.Vector(tableau'first(2)..lastcol);
    index : integer32;

  begin
    Scale(tableau,lastcol); Scale(rhs);
   -- put_line("The tableau after scaling : "); put(tableau,3,3,3);
   -- put_line(" with scaled right hand side : "); put(rhs,3,3,3); new_line;
    ips := Inner_Products(tableau,lastcol,rhs);
    index := Index_of_Maximum(ips,lastcol);
    if (abs(ips(index)) < tol) or else (ips(index) < 0.0) then
      if index <= rhs'last
       then feasible := (abs(rhs(index)) < tol);
       else feasible := false;
      end if;
    else
      feasible := true;
      columns(columns'first) := index;
      columns(index) := columns'first;
      Interchange_Columns(tableau,tableau'first(2),index);
      index := Maximal_Product(tableau,tableau'first(2),rhs,rhs'first);
      if index /= tableau'first(2)
       then Interchange_Rows(tableau,lastcol,rhs,tableau'first(2),index);
      end if;
      Rotate(tableau,lastcol,tableau'first(2),rhs,tol);
      -- Upper_Triangulate(tableau,lastcol,tableau'first(2),rhs,tol);
      sol(sol'first) := rhs(rhs'first);
    end if;
   -- put_line("At the end of Initialize : ");
   -- put_line("The tableau : "); put(tableau,3,3,3);
   -- put_line(" with right hand side : "); put(rhs,3,3,3); new_line;
   -- put(" and last column "); put(lastcol,1); new_line;
  end Initialize;

-- PERFORM THE NEXT STEPS :

  procedure Select_Column 
                  ( mat : in Matrix; lastcol : in integer32; 
                    rhs,cosines : in Standard_Floating_Vectors.Vector;
                    index : in integer32;
                    tol : in double_float; row,col : out integer32 ) is

  -- DESCRIPTION :
  --   Selects a column col with in matrix for which
  --    1) the cosine is maximal and
  --    2) row = Maximal_Product(mat,col,rhs,index), with mat(row,col) > 0.0
  --   If the resulting column is smaller than index, then no such column
  --   has been found.

    cl : integer32 := index;
    maxcos,prod : double_float;
    rw : integer32 := Maximal_Product(mat,cl,rhs,index);
    mp : integer32;

  begin
   -- put_line("The tableau when Select_Column :"); put(mat,3,3,3);
   -- put_line(" with right hand side : "); put(rhs,3,3,3); new_line;
   -- put_line("and vector of cosines : "); put(cosines,3,3,3); new_line;
   -- put("The index : "); put(index,1); new_line;
    prod := mat(rw,cl)*rhs(rw);
    if prod > tol
      and then Positive(mat,rhs,index,cl,rhs(rw)/mat(rw,cl),tol)
              -- Positive_Diagonal(mat,rhs,index,cl,rhs(rw)/mat(rw,cl),tol)
     then maxcos := cosines(cl);
     else maxcos := -2.0;
    end if;
    for i in index+1..lastcol loop
      mp := Maximal_Product(mat,i,rhs,index);
     -- put("mp : "); put(mp,1); put(" when testing "); put(i,1);
     -- put_line("th column");
      prod := mat(mp,i)*rhs(mp);
      if prod > tol and then cosines(i) > maxcos then
        if Positive(mat,rhs,index,i,rhs(mp)/mat(mp,i),tol)
           -- Positive_Diagonal(mat,rhs,index,i,rhs(mp)/mat(mp,i),tol)
         then cl := i; rw := mp; maxcos := cosines(i);
        end if;
      end if;
    end loop;
   -- put("maximal cosine : "); put(maxcos,3,3,3); new_line;
    if maxcos >= -1.0
     then col := cl; row := rw;
         -- put("column : "); put(cl,1);
         -- put(" and row : "); put(rw,1); new_line;
     else col := index-1; row := index-1;
         -- put("column : "); put(index-1,1);
         -- put(" and row : "); put(index-1,1); new_line;
    end if;
  end Select_Column;

  procedure Select_Column 
                  ( mat : in Matrix; lastcol : in integer32; 
                    rhs : in Standard_Floating_Vectors.Vector;
                    index,column : in integer32;
                    tol : in double_float; row,col : out integer32 ) is

  -- DESCRIPTION :
  --   Selects first column col with in matrix for which
  --    row = Maximal_Product(mat,col,rhs,index), with mat(row,col) > 0.0.
  --   If the resulting column is smaller than index, then no such column
  --   has been found.

    cl : integer32 := index-1;
    prod : double_float;
    rw,mp : integer32;

  begin
   -- put_line("The tableau when Select_Column :"); put(mat,3,3,3);
   -- put_line(" with right hand side : "); put(rhs,3,3,3); new_line;
   -- put("The index : "); put(index,1); new_line;
    for i in column..lastcol loop
      mp := Maximal_Product(mat,i,rhs,index);
     -- put("mp : "); put(mp,1); put(" when testing "); put(i,1);
     -- put_line("th column");
      prod := mat(mp,i)*rhs(mp);
      if prod > tol
        and then Positive(mat,rhs,index,i,rhs(mp)/mat(mp,i),tol)
                -- Positive_Diagonal(mat,rhs,index,i,rhs(mp)/mat(mp,i),tol)
       then cl := i; rw := mp;
      end if;
      exit when (cl >= index);
    end loop;
    if cl >= index
     then col := cl; row := rw;
         -- put("column : "); put(cl,1);
         -- put(" and row : "); put(rw,1); new_line;
     else col := index-1; row := index-1;
         -- put("column : "); put(index-1,1);
         -- put(" and row : "); put(index-1,1); new_line;
    end if;
  end Select_Column;

  procedure Next ( tableau : in out Matrix; lastcol : in integer32;
                   rhs,sol,cosines : in out Standard_Floating_Vectors.Vector;
                   pivots : in out Standard_Integer_Vectors.Vector;
                   tol : in double_float;
                   index : in integer32; feasi : in out boolean ) is

  -- DESCRIPTION :
  --   Selects one new column out of the tableau.

  -- ON ENTRY :
  --   tableau     upper triangular up to the row and column index-1;
  --   lastcol     index of the last column of interest in the tableau;
  --   rhs         right hand side vector;
  --   sol         solution vector for the components 1..index;
  --   cosines     vector of cosines;
  --   pivots      vector of column interchangements;
  --   tol         tolerance to decide whether a number equals zero;
  --   index       indicator of the current step;
  --   feasi       must be true on entry.

  -- ON RETURN :
  --   tableau     upper triangular up to the row and column index,
  --               under the condition that feasi remains true;
  --   rhs         transformed right hand side vector;
  --   sol         new solution, with additional meaningful component
  --               sol(index);
  --   cosines     eventually permuted vector of cosines:
  --               cosines(index) is the cosine of the (new) column,
  --               indicated by index, and the right hand side vector;
  --   pivots      updated vector of column interchangements;
  --   feasi       if true, then a new column which yields a positive
  --               contribution of the right hand side vector has been 
  --               found, if this was not the case, then feasi is false;

    col,row,tmp : integer32;

  begin
    Select_Column(tableau,lastcol,rhs,cosines,index,tol,row,col);
    if col < index then
      feasi := false;
    else
      feasi := true;
      if col /= index then
        tmp := pivots(col); pivots(col) := pivots(index);
        pivots(index) := tmp;
        Interchange_Columns(tableau,col,index);
        Interchange(cosines,col,index);
      end if;
      if row /= index
       then Interchange_Rows(tableau,lastcol,rhs,row,index);
      end if;
      Rotate(tableau,lastcol,index,rhs,tol);
      -- Upper_Triangulate(tableau,lastcol,index,rhs,tol);
      -- Solve(tableau,rhs,tol,index,sol);  -- only needed at end
      -- feasi := Positive(sol,index,tol);  -- double check...
    end if;
  end Next;

  procedure Next ( tableau : in out Matrix; lastcol : in integer32;
                   rhs,sol : in out Standard_Floating_Vectors.Vector;
                   pivots : in out Standard_Integer_Vectors.Vector;
                   tol : in double_float;
                   index : in integer32; rank : out integer32;
                   feasi : in out boolean ) is

  -- DESCRIPTION :
  --   Lexicographic enumeration of the candidates, stops when a feasible
  --   solution is encountered.

  -- ON ENTRY :
  --   tableau     upper triangular up to the row and column index-1;
  --   lastcol     index of the last column of interest in the tableau;
  --   rhs         right hand side vector;
  --   sol         solution vector for the components 1..index;
  --   pivots      vector of column interchangements;
  --   tol         tolerance to decide whether a number equals zero;
  --   index       indicator of the current step;
  --   feasi       must be true on entry.

  -- ON RETURN :
  --   tableau     upper triangular up to the row and column index,
  --               under the condition that feasi remains true;
  --   rhs         transformed right hand side vector;
  --   sol         new solution, with additional meaningful component
  --               sol(index);
  --   pivots      updated vector of column interchangements;
  --   rank        rank of the current tableau;
  --   feasi       if true, then a new column which yields a positive
  --               contribution of the right hand side vector has been 
  --               found, if this was not the case, then feasi is false;

    column,col,row,tmp : integer32;
    stop : boolean := false;

  begin
    column := index;
    loop
      Select_Column(tableau,lastcol,rhs,index,column,tol,row,col);
      if col < column then
        feasi := false; stop := true;
      else
        feasi := true;
        if col /= index then
          tmp := pivots(col); pivots(col) := pivots(index);
          pivots(index) := tmp;
          Interchange_Columns(tableau,col,index);
        end if;
        if row /= index
         then Interchange_Rows(tableau,lastcol,rhs,row,index);
        end if;
        Rotate(tableau,lastcol,index,rhs,tol);
       -- Upper_Triangulate(tableau,lastcol,index,rhs,tol);
        if (index = lastcol) or else (index = tableau'last(1))
                             or else (Norm(rhs,index+1) < tol) then
          stop := true; rank := index;
        else
          Next(tableau,lastcol,rhs,sol,pivots,tol,index+1,rank,feasi);
          if feasi
           then stop := true;
           else column := col+1;
          end if;
         end if;
      end if;
      exit when stop;
    end loop;
  end Next;

  procedure Complementary_Slackness
                 ( tableau : in out matrix;
                   rhs : in out Standard_Floating_Vectors.Vector;
                   tol : in double_float;
                   solution : out Standard_Floating_Vectors.Vector;
                   columns : out Standard_Integer_Vectors.Vector;
                   feasible : out boolean ) is
  begin
    Complementary_Slackness
       (tableau,tableau'last(2),rhs,tol,solution,columns,feasible);
  end Complementary_Slackness;

  procedure Complementary_Slackness1
                 ( tableau : in out Matrix; lastcol : in integer32;
                   rhs : in out Standard_Floating_Vectors.Vector;
                   tol : in double_float; 
                   solution : out Standard_Floating_Vectors.Vector;
                   columns : out Standard_Integer_Vectors.Vector;
                   feasible : out boolean ) is

  -- NOTE :
  --   This version is a straigth forward search and finds not always
  --   the feasible solution.

    feasi : boolean;
    pivots : Standard_Integer_Vectors.Vector(tableau'first(2)..lastcol);
    cosis : Standard_Floating_Vectors.Vector(tableau'first(2)..lastcol);
    sol : Standard_Floating_Vectors.Vector(rhs'range) := (rhs'range => 0.0);
    rank : integer32;

  begin
    for i in pivots'range loop  -- initialize vector of pivots
      pivots(i) := i;
    end loop;
    Initialize(tableau,lastcol,rhs,sol,tol,pivots,cosis,feasi);
    if feasi then
      if Norm(rhs,rhs'first+1) > tol then -- check whether residual > 0
        if tableau'first(1)+1 > lastcol then
          feasi := false;
        else
          rank := 1;
          for i in tableau'first(1)+1..tableau'last(1) loop
            Next(tableau,lastcol,rhs,sol,cosis,pivots,tol,i,feasi);
            exit when not feasi;
            exit when (i = lastcol);
            rank := rank + 1;
            exit when (Norm(rhs,i+1) < tol);
          end loop;
          if feasi
           then Solve(tableau,rhs,tol,rank,sol);
          end if;
          end if;
       end if;
       if feasi and then (tableau'last(1) > lastcol)
        then feasi := (Norm(rhs,lastcol+1) < tol);
       end if;
       if feasi and then (Norm(sol) < tol)
        then feasi := (Norm(rhs) < tol);
       end if;
       solution := sol;
       if tableau'last(1) <= lastcol then
         for i in tableau'range(1) loop
           columns(i) := pivots(i);
         end loop;
       else
         for i in tableau'first(2)..lastcol loop
           columns(i) := pivots(i);
         end loop;
       end if;
    end if;
    feasible := feasi;
  end Complementary_Slackness1;

 -- procedure Complementary_Slackness_the_Original
  procedure Complementary_Slackness
                 ( tableau : in out Matrix; lastcol : in integer32;
                   rhs : in out Standard_Floating_Vectors.Vector;
                   tol : in double_float;
                   solution : out Standard_Floating_Vectors.Vector;
                   columns : out Standard_Integer_Vectors.Vector;
                   feasible : out boolean ) is

  -- NOTE :
  --   This version enumerates in a lexicographic order all candidates for
  --   entering the basis and stops when a feasible solution has been found.

    feasi : boolean;
    pivots : Standard_Integer_Vectors.Vector(tableau'first(2)..lastcol);
    sol : Standard_Floating_Vectors.Vector(rhs'range) := (rhs'range => 0.0);
    ind,rank : integer32;

  begin
    for i in pivots'range loop  -- initialize vector of pivots
      pivots(i) := i;
    end loop;
    Initialize(tableau,lastcol,rhs,sol,tol,pivots,feasi);
    if feasi then
      if Norm(rhs,rhs'first+1) > tol then  -- check whether residual > 0
        if tableau'first(1)+1 > lastcol then
          feasi := false;
        else
          ind := tableau'first(1)+1;
          Next(tableau,lastcol,rhs,sol,pivots,tol,ind,rank,feasi);
          if feasi
           then Solve(tableau,rhs,tol,rank,sol);
          end if;
        end if;
      end if;
      if feasi and then (tableau'last(1) > lastcol)
       then feasi := (Norm(rhs,lastcol+1) < tol);
      end if;
      if feasi and then (Norm(sol) < tol)
       then feasi := (Norm(rhs) < tol);
      end if;
      solution := sol;
      if tableau'last(1) <= lastcol then
        for i in tableau'range(1) loop
          columns(i) := pivots(i);
        end loop;
      else
        for i in tableau'first(2)..lastcol loop
          columns(i) := pivots(i);
        end loop;
      end if;
    end if;
    feasible := feasi;
  end Complementary_Slackness;
 -- end Complementary_Slackness_the_Original;

  procedure Complementary_Slackness2
                 ( tableau : in out Matrix; lastcol : in integer32;
                   rhs : in out Standard_Floating_Vectors.Vector;
                   tol : in double_float;
                   solution : out Standard_Floating_Vectors.Vector;
                   columns : out Standard_Integer_Vectors.Vector;
                   feasible : out boolean ) is

  -- NOTE :
  --   This version explicitly uses linear programming.

  begin
    Is_In_Cone(tableau,lastcol,rhs,tol,solution,columns,feasible);
  end Complementary_Slackness2;

  procedure Complementary_Slackness_the_New
 -- procedure Complementary_Slackness
                 ( tableau : in out Matrix; lastcol : in integer32;
                   rhs : in out Standard_Floating_Vectors.Vector;
                   tol : in double_float;
                   solution : out Standard_Floating_Vectors.Vector;
                   columns : out Standard_Integer_Vectors.Vector;
                   feasible : out boolean ) is
  
  begin
    if rhs(rhs'last) > 0.0 then
      declare
        n : constant integer32 := tableau'last(1)-1;
        m : constant integer32 := lastcol;
        lincff : Standard_Floating_Matrices.Matrix(1..n,1..m);
        difrhs : Standard_Floating_Vectors.Vector(1..m);
        sol : Standard_Floating_Vectors.Vector(1..n+1);
        fail : boolean;
      begin
        columns := (columns'range => 0);
        solution := (solution'range => 0.0);
        for j in 1..m loop
          difrhs(j) := tableau(n+1,j)+abs(Random); --*10.0**(-8);
          for i in 1..n loop
            lincff(i,j) := tableau(i,j) - rhs(i);
          end loop;
        end loop;
        Linear_Minimization.Feasible(n,m,lincff,difrhs,tol,sol,fail);
        feasible := fail;
      end;
    else
      -- Complementary_Slackness_the_Original
      Complementary_Slackness
        (tableau,lastcol,rhs,tol,solution,columns,feasible);
--       declare
--         n : constant natural := tableau'last(1);
--         m : constant natural := lastcol+1;
--         lincff : Standard_Floating_Matrices.Matrix(1..n,1..m);
--         difrhs : Standard_Floating_Vectors.Vector(1..m);
--         sol : Standard_Floating_Vectors.Vector(1..n+1);
--         fail : boolean;
--       begin
--         columns := (columns'range => 0);
--         solution := (solution'range => 0.0);
--         for j in 1..m-1 loop
--           difrhs(j) := abs(Random); --*10.0**(-8);
--           for i in 1..n loop
--             lincff(i,j) := tableau(i,j);
--           end loop;
--         end loop;
--         difrhs(m) := abs(Random);
--         for i in 1..n-1 loop
--           lincff(i,m) := 0.0;
--         end loop;
--         lincff(n,m) := 1.0;
--         Linear_Minimization.Feasible(n,m,lincff,difrhs,tol,sol,fail);
--         feasible := fail;
--       end;
    end if;
 -- end Complementary_Slackness;
  end Complementary_Slackness_the_New;

end Floating_Linear_Inequalities;
