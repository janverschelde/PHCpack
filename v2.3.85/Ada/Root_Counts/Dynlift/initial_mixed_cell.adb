with Standard_Integer_VecVecs;           use Standard_Integer_VecVecs;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Standard_Integer_Linear_Solvers;    use Standard_Integer_Linear_Solvers;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Transforming_Integer32_Vector_Lists;
 use Transforming_Integer32_Vector_Lists;
with Standard_Integer32_Vertices;        use Standard_Integer32_Vertices;
with Frequency_Graph;                    use Frequency_Graph;

procedure Initial_Mixed_Cell
                ( n : in integer32; mix : in Vector; pts : in Array_of_Lists;
                  mic : out Mixed_Cell; rest : in out Array_of_Lists ) is

  res : Mixed_Cell;
  acc,tmp,rest_last : List;
  pt : Link_to_Vector;
  expts : Array_of_Lists(pts'range);
  perm : Vector(pts'range);
  mind,minl,dims,lent,index : integer32;
  nullvec : constant Vector := (1..n => 0);
  shiftvecs : VecVec(pts'range);
  grap : Matrix(1..n,pts'range);
  fail : boolean := false;
  rowcnt,cnt : integer32 := 0;
  mat : Matrix(1..n,1..n+1);
  ipvt : Vector(1..n+1);

-- AUXILIAIRIES

  procedure Add ( pt,shift : in Vector; L : in out List ) is

    newpt : constant Link_to_Vector := new Vector'(pt);

  begin
    Add(newpt.all,shift);
    Construct(newpt,L);
  end Add;

  procedure Add ( pt : in Vector; L : in out List ) is

    newpt : constant Link_to_Vector := new Vector'(pt);

  begin
    Construct(newpt,L);
  end Add;

  procedure Fill ( L : in List; rowcnt : in out integer32;
                   m : in out Matrix ) is

  -- DESCRIPTION :
  --   Fills the matrix with the elements in l.

  -- ON ENTRY :
  --   L        list of point vectors;
  --   rowcnt   m(m'first(1)..rowcnt) is already defined;
  --   m        matrix with m'range(2) = range of points in l.

  --   rowcnt   updated row counter;
  --   m        matrix with the nonzero points of L

    tmp : List := L;
    pt : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      declare
        nullvec : constant Vector(pt'range) := (pt'range => 0);
      begin
        if pt.all /= nullvec then
          rowcnt := rowcnt + 1;
          for k in pt'range loop
            m(rowcnt,k) := pt(k);
          end loop;
        end if;
        tmp := Tail_Of(tmp);
      end;
    end loop;
  end Fill;

  procedure Initialize ( L : in List; rowcnt : out integer32;
                         m : in out Matrix; ipvt : in out Vector ) is

  -- DESCRIPTION :
  --   Given a list of linearly independent points, eventually with
  --   0 in l, a matrix will be constructed which contains these points
  --   in upper triangular form.

    cnt : integer32 := 0;
    piv : integer32; 

  begin
    Fill(L,cnt,m);
    for i in 1..cnt loop
      m(i,m'last(2)) := i;
    end loop;
    for k in 1..cnt loop
      Triangulate(k,m,ipvt,piv);
    end loop;
    rowcnt := cnt;
  end Initialize;

  procedure Linearly_Independent
                   ( m : in out Matrix; rowcnt : in integer32;
                     ipvt : in Vector;
                     L : in List; len : out integer32;
                     indep,indep_last : out List ) is

  -- DESCRIPTION :
  --   Computes the list of points which are linearly independent w.r.t.
  --   the matrix m.

  -- ON ENTRY :
  --   m             m(1..rowcnt,m'range(2)) is upper triangular,
  --                 can be used as work space;
  --   rowcnt        counter for the number of rows in m:
  --   ipvt          vector with the pivoting information;
  --   l             list of points to consider.

  -- ON RETURN :
  --   len           length of the list indep;
  --   indep         the list of linearly independent points;
  --   indep_last    pointer to the last element of indep.

    tmp,res,res_last : List;
    pt : Link_to_Vector;
    wipvt : Vector(ipvt'range) := ipvt;
    piv,cnt : integer32;

  begin
    if rowcnt < m'last(2) then
      tmp := l; cnt := 0;
      while not Is_Null(tmp) loop
        pt := Head_Of(tmp);
        for i in pt'range loop
          m(rowcnt+1,i) := pt(i);
        end loop;
        m(rowcnt+1,m'last(2)) := rowcnt+1;
        Triangulate(rowcnt+1,m,wipvt,piv);
        if m(rowcnt+1,rowcnt+1) /= 0 then
          cnt := cnt + 1;
          Append(res,res_last,pt.all);
        end if;
        wipvt := ipvt;
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    len := cnt;
    indep := res; indep_last := res_last;
  end Linearly_Independent;

  procedure Construct_Basis ( L : in List; rowcnt : in out integer32;
                              m : in out Matrix; ipvt : in out Vector ) is

  -- DESCRIPTION :
  --   Constructs a triangular basis for the vectors in L.

  -- REQUIRED :  dimensions of the vectors must match.

  -- ON ENTRY :
  --   L         list of vector;
  --   m         triangular basis, stored in the rows of the matrix;
  --   rowcnt    index to the last significant row in m, could be 0,
  --             when there is no initial basis to take into account;
  --   ipvt      initial pivoting vector, must be equal to the identity
  --             permutation vector, when rowcnt = 0.

  -- ON RETURN :
  --   m         triangular basis, stored in the rows of the matrix;
  --   rowcnt    index to the last significant row in m;
  --   ipvt      contains pivoting information.

    tmp : List := L;
    wipvt : Vector(ipvt'range);
    pt : Link_to_Vector;
    piv : integer32;

  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      for i in pt'range loop
        m(rowcnt+1,i) := pt(i);
      end loop;
      wipvt := ipvt;
      m(rowcnt+1,m'last(2)) := rowcnt+1;
      Triangulate(rowcnt+1,m,wipvt,piv);
      if m(rowcnt+1,rowcnt+1) /= 0
       then rowcnt := rowcnt + 1;
            ipvt := wipvt;
      end if;
      exit when (rowcnt >= m'last(2));
      tmp := Tail_Of(tmp);
    end loop;
  end Construct_Basis;

--  function Linearly_Independent ( l : List; x : Vector ) return boolean is
--
--  -- DESCRIPTION :
--  --   Returns true if the given vector x is linearly independent w.r.t.
--  --   the vectors in l.
--
--    basis : Matrix(1..Length_Of(l)+1,x'first..x'last+1);
--    ipvt : Vector(basis'range(2));
--    rowcnt : natural := 0;
--    piv : natural;
--
--  begin
--    for i in ipvt'range loop
--      ipvt(i) := i;
--      for j in basis'range(1) loop
--        basis(j,i) := 0;
--      end loop;
--    end loop;
--    Construct_Basis(l,rowcnt,basis,ipvt);
--    for i in x'range loop
--      basis(rowcnt+1,i) := x(i);
--    end loop;
--    basis(rowcnt+1,basis'last(2)) := rowcnt+1;
--    Triangulate(rowcnt+1,basis,ipvt,piv);
--    return (basis(rowcnt+1,rowcnt+1) /= 0);
--  end Linearly_Independent;

  procedure Linearly_Independent
                ( L,xl : in List; indep,indep_last : in out List ) is

  -- DESCRIPTION :
  --   Returns those points in xl that are linearly independent from the
  --   points in L.

  begin
    if not Is_Null(xl) then
      declare
        x : constant Link_to_Vector := Head_Of(xl);
        basis : Matrix(1..integer32(Length_Of(L))+1,x'first..x'last+1);
        ipvt : Vector(basis'range(2));
        rowcnt,len : integer32 := 0;
      begin
        for i in ipvt'range loop
          ipvt(i) := i;
          for j in basis'range(1) loop
            basis(j,i) := 0;
          end loop;
        end loop;
        Construct_Basis(l,rowcnt,basis,ipvt);
        Linearly_Independent(basis,rowcnt,ipvt,xl,len,indep,indep_last);
      end;
    end if;
  end Linearly_Independent;

  function Construct_Rest ( L : Array_of_Lists; start : integer32 )
                          return List is

  -- DESCRIPTION :
  --   Collects all remaining points of L(i) into one single list,
  --   for i in start..L'last.

    tmp,res,res_last : List;
    pt : Link_to_Vector;

  begin
    for i in start..L'last loop
      tmp := L(i);
      while not Is_Null(tmp) loop
        pt := Head_Of(tmp);
        if not Is_In(res,pt)
         then Append(res,res_last,pt.all);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Construct_Rest;

  procedure Sort_Co_Linearly_Independent
               ( L : in out List; al : in Array_of_Lists;
                 start : in integer32 ) is

  -- DESCRIPTION :
  --   Sorts the points in l, by putting those points that are linearly
  --   independent w.r.t. the rest of the poins in al, in front of the list.

    rest : List := Construct_Rest(al,start);
    tmp,indep,indep_last : List;
    pt : Link_to_Vector;

  begin
   -- put_line("The set of points in al just after Construct_Rest :"); put(al);
   -- put_line("The rest of the points : "); put(rest);
    Linearly_Independent(rest,l,indep,indep_last);
    if Is_Null(indep) then
      null;                              -- nothing to put in front of l
    else
      tmp := L;                 -- append the other points of l to indep
      -- put_line("Linearly co-independent points : "); put(indep);
      -- put_line(" w.r.t. the rest : "); put(rest);
      while not Is_Null(tmp) loop
        pt := Head_Of(tmp);
        if not Is_In(indep,pt)
         then Append(indep,indep_last,pt.all);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      Copy(indep,l); Deep_Clear(indep);
    end if;
    Deep_Clear(rest);
   -- put_line("The list of points in al after Sort_Co : "); put(al);
  end Sort_Co_Linearly_Independent;

  procedure Incremental_Dimension
               ( m : in out Matrix; rowcnt : in integer32; ipvt : in Vector;
                 l : in List; dim : out integer32 ) is

  -- DESCRIPTION :
  --   Computes the number of points which are linearly independent w.r.t.
  --   the matrix m.

  -- ON ENTRY :
  --   m         m(1..rowcnt,m'range(2)) is upper triangular,
  --             can be used as work space;
  --   rowcnt    counter for the number of rows in m:
  --   ipvt      vector with the pivoting information;
  --   l         list of points to consider.

  -- ON RETURN :
  --   dim       the number of linearly independent points.

   tmp : List;
   pt : Link_to_Vector;
   cnt,piv : integer32;
   newipvt,wipvt : Vector(ipvt'range);

  begin
    wipvt := ipvt; newipvt := ipvt;
    tmp := l; cnt := 0;
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      cnt := cnt + 1;
      for i in pt'range loop
        m(rowcnt+cnt,i) := pt(i);
      end loop;
      m(rowcnt+cnt,m'last(2)) := rowcnt+cnt;
      Triangulate(rowcnt+cnt,m,newipvt,piv);
      if m(rowcnt+cnt,rowcnt+cnt) /= 0
       then wipvt := newipvt;
       else newipvt := wipvt;
            cnt := cnt - 1;
      end if;
      tmp := Tail_Of(tmp);
      exit when ((rowcnt + cnt) >= m'last(1));
    end loop;
    dim := cnt;
  end Incremental_Dimension;

  procedure Incremental_Dimension
               ( m : in Matrix; rowcnt : in integer32; ipvt : in Vector;
                 L : in out List; dim,len : out integer32 ) is

  -- DESCRIPTION :
  --   Computes the increase in dimension by considering the points
  --   in the list l, w.r.t. the matrix m.

  -- ON ENTRY :
  --   m         m(1..rowcnt,m'range(2)) is upper triangular,
  --   rowcnt    counter for the number of rows in m:
  --   ipvt      vector with the pivoting information;
  --   L         list of points to consider.
 
  -- ON RETURN :
  --   L         list of linearly independent points w.r.t. m;
  --   dim       increase in dimension;
  --   len       length of the list of linearly independent points
  --             w.r.t. the rows in the matrix m.

  -- ALGORITHM :
  --   works in two stages:
  --   1. determination of all linearly independent points;
  --   2. determination of the dimension.

    workm : Matrix(m'range(1),m'range(2));
    indep,indep_last : List;

  begin
   -- Initialize workm :
   -- put_line("The list to investigate : "); put(l);
    for i in 1..rowcnt loop
      for j in m'range(2) loop
        workm(i,j) := m(i,j);
      end loop;
    end loop;
   -- Determine linearly independent points :
    Linearly_Independent(workm,rowcnt,ipvt,L,len,indep,indep_last);
   -- Determine the incremental dimension :
    Incremental_Dimension(workm,rowcnt,ipvt,indep,dim);
    Copy(indep,L); Deep_Clear(indep);
  end Incremental_Dimension;

  procedure Next_Point ( acc,L : in List; pt : out Link_to_Vector;
                         fail : out boolean; rowcnt : in out integer32;
                         m : in out Matrix; ipvt : in out Vector ) is

  -- DESCRIPTION :
  --   A new point out of L, not in the list acc will be chosen.
  --   The point pt has to be linearly independent from the other,
  --   already chosen points.  Therefore, the upper triangular matrix m
  --   will be used and properly updated.

    res : Link_to_Vector;
    tmp : List := L;
    newipvt : Vector(ipvt'range) := ipvt;
    done : boolean := false;
    piv : integer32;

  begin
    while not Is_Null(tmp) loop
      res := Head_Of(tmp);
      if not Is_In(acc,res) then
        --put("Checking point "); put(res); new_line;
        rowcnt := rowcnt + 1;
        for i in res'range loop
          m(rowcnt,i) := res(i);
        end loop;
        m(rowcnt,m'last(2)) := rowcnt;
        Triangulate(rowcnt,m,newipvt,piv);
        if m(rowcnt,rowcnt) = 0 then
          rowcnt := rowcnt - 1;
          newipvt := ipvt;
        else
          ipvt := newipvt;
          pt := res;
          done := true;
        end if;
      end if;
      exit when done;
      tmp := Tail_Of(tmp);
    end loop;
    fail := not done;
  end Next_Point;

begin
 -- INITIALIZE : compute for each list the basis points
  res.nor := new Vector'(1..n+1 => 0); res.nor(n+1) := 1;
  res.pts := new Array_of_Lists(mix'range);
  for i in pts'range loop
    expts(i) := Max_Extremal_Points(n,pts(i));  perm(i) := i;
   -- put("Extremal points for component "); put(i,1); put_line(" :");
   -- put(expts(i));
  end loop;
 -- INITIALIZE : order the lists according to occurencies
  grap := Graph(n,expts);
 -- put_line("The graph matrix of the extremal points : "); put(grap);
  for i in expts'range loop
    Sort(expts(i),grap);
   -- put("Ordered extremal points for component "); put(i,1); put_line(" :");
   -- put(expts(i));
  end loop;
 -- INITIALIZE : choose one anchor point for each component,
 --              shift when necessary :
  for i in expts'range loop
    if Is_In(pts(i),nullvec) then
      Add(nullvec,res.pts(i));
      shiftvecs(i) := null;
    else -- ADD LAST POINT = LEAST IMPORTANT
      tmp := expts(i);
      if not Is_Null(tmp) then
        while not Is_Null(Tail_Of(tmp)) loop
          tmp := Tail_Of(tmp);
        end loop;
        pt := Head_Of(tmp);
        Add(pt.all,res.pts(i));
        shiftvecs(i) := new Vector'(pt.all);
        -- put("The shift vector : "); put(shiftvecs(i)); new_line;
        Shift(expts(i),shiftvecs(i));
        -- put_line("The shifted list of points : "); put(expts(i));
      end if;
    end if;
  end loop;
 -- INITIALIZE : intermediate data structures mat and ipvt :
  for i in ipvt'range loop
    ipvt(i) := i;
  end loop;
  for i in mat'range(1) loop
    for j in mat'range(2) loop
      mat(i,j) := 0;
    end loop;
  end loop;
  rowcnt := 0;
 -- COMPUTE CELL : based on lists with extremal points
  for i in expts'range loop   -- MAIN LOOP
   -- LOOK FOR SMALLEST SET :
    index := i;
    if i = 1 then
      mind := integer32(Length_Of(expts(i)))-1;
      for j in i+1..expts'last loop
        dims := integer32(Length_Of(expts(j)))-1;
        if dims < mind
         then index := j; mind := dims;
        end if;
      end loop;
    else
      --put("Calling Incremental_Dimension for i = "); put(i,1); new_line;
      Incremental_Dimension(mat,rowcnt,ipvt,expts(i),mind,minl);
      --put(" dimension : "); put(mind,1); 
      --put(" length : " ); put(minl,1); new_line;
      for j in i+1..expts'last loop
        --put("Calling Incremental_Dimension for j = "); put(j,1); new_line;
        Incremental_Dimension(mat,rowcnt,ipvt,expts(j),dims,lent);
        --put(" dimension : "); put(dims,1); 
        --put(" length : " ); put(lent,1); new_line;
        if (dims < mind) or else ((dims = mind) and then (lent < minl))
         then index := j; mind := dims; minl := lent;
        end if;
      end loop;
    end if;
   -- put("The index : "); put(index,1); new_line;
   -- put(" with minimal dimension : "); put(mind,1); new_line;
   -- put("perm before permute : "); put(perm); new_line;
    if index /= i then
      tmp := expts(index); expts(index) := expts(i); expts(i) := tmp;
      mind := perm(i);  perm(i) := perm(index);  perm(index) := mind;
    end if;
   -- SORT ACCORDING TO OCCURRENCIES AND TO CO-LINEARLY :
    Sort(expts(i),grap,i-1,perm);
   -- put_line("After first ordering : "); put(expts(i));
    Sort_Co_Linearly_Independent(expts(i),expts,i+1);
   -- put_line("The ordered list : "); put(expts(i));
   -- put("perm after permute : "); put(perm); new_line;
   -- CHOOSE THE POINTS :
    if i = 1 then
      Add(nullvec,acc);
      tmp := expts(i);
      for j in 1..mix(perm(i)) loop
        pt := Head_Of(tmp);
        if pt.all = nullvec then  -- skip the null vector
          tmp := Tail_Of(tmp);
          if not Is_Null(tmp)
           then pt := Head_Of(tmp);
          end if;
        end if;
        exit when Is_Null(tmp);
        Add(pt.all,acc); 
        -- put_line("frequency graph before ignore :"); put(grap);
        Ignore(grap,pt.all);
        -- put_line("frequency graph after ignore  :"); put(grap);
        if shiftvecs(perm(i)) /= null
         then Add(pt.all,shiftvecs(perm(i)).all,res.pts(perm(i)));
         else Add(pt.all,res.pts(perm(i)));
        end if;
        tmp := Tail_Of(tmp);
        exit when Is_Null(tmp);
      end loop;
      cnt := integer32(Length_Of(acc));
      Initialize(acc,rowcnt,mat,ipvt);
      -- put_line("The list acc : "); put(acc);
      -- put_line("The matrix mat : "); put(mat,1,rowcnt);
      -- put_line("The vector ipvt : "); put(ipvt); new_line;
      fail := (cnt <= mix(perm(i)));
    else
      for j in 1..mix(perm(i)) loop
        Next_Point(acc,expts(i),pt,fail,rowcnt,mat,ipvt);
        if not fail then
          if shiftvecs(perm(i)) /= null
           then Add(pt.all,shiftvecs(perm(i)).all,res.pts(perm(i)));
           else Add(pt.all,res.pts(perm(i)));
          end if;
          Add(pt.all,acc); cnt := cnt + 1; 
          -- put_line("frequency graph before ignore :"); put(grap);
          Ignore(grap,pt.all);
          -- put_line("frequency graph after ignore  :"); put(grap);
        end if;
        exit when fail;
      end loop;
      -- put_line("The list acc : "); put(acc);
      -- put_line("The matrix mat : "); put(mat,1,rowcnt);
      -- put_line("The vector ipvt : "); put(ipvt); new_line;
      fail := (integer32(Length_Of(res.pts(perm(i)))) <= mix(perm(i)));
    end if;
    exit when fail;
  end loop;  -- END OF MAIN LOOP
  Deep_Clear(acc);
 -- COMPUTE THE REST OF THE POINT LISTS :
  if not fail then
    for i in pts'range loop
      tmp := pts(i);
      rest_last := rest(i);
      while not Is_Null(tmp) loop
        pt := Head_Of(tmp);
        if not Is_In(res.pts(i),pt)
         then Append(rest(i),rest_last,pt);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end if;
 -- GIVE THE POINTS IN THE INITIAL SIMPLEX LIFTING VALUE 0 :
  for i in res.pts'range loop
    tmp := res.pts(i);
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      declare
        lpt : Link_to_Vector;
      begin
        lpt := new Vector(1..n+1);
        lpt(pt'range) := pt.all; lpt(n+1) := 0;
        Set_Head(tmp,lpt);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end loop;
  mic := res;
end Initial_Mixed_Cell;
