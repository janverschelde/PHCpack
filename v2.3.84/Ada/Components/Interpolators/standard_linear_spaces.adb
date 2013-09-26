with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;

package body Standard_Linear_Spaces is

-- AUXILIARY :

  function Is_In ( v : Standard_Integer_Vectors.Vector; k : integer32 )
                 return boolean is

  -- DESCRIPTION :
  --   Returns true if there is an index i in v such that v(i) = k.

  begin
    for i in v'range loop
      if v(i) = k
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

-- TARGET ROUTINES :

  procedure Rank ( nbvecs,n : in integer32;
                   pts : in Matrix; tol : in double_float;
                   trivec : out Matrix; rnk : out natural32 ) is

    col : integer32 := 1;

  begin
    for i in 1..nbvecs loop
      for j in 1..n loop
        trivec(i,j) := pts(i,j) - pts(nbvecs+1,j);
      end loop;
    end loop;
    Triangulate(trivec,tol,nbvecs,n);
    rnk := 0;
    for i in 1..nbvecs loop
      while (col <= n) and then (AbsVal(trivec(i,col)) < tol) loop
        col := col + 1;
      end loop;
      exit when (col > n);
      rnk := rnk + 1;
    end loop;
  end Rank;

  function Pivots ( nbvecs,n : integer32; 
                    trivec : Matrix; rnk : natural32; tol : double_float )
                  return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..integer32(rnk));
    col : integer32 := 1;

  begin
    for i in 1..nbvecs loop
      while (col <= n)
          and then (AbsVal(trivec(i,col)) < tol) loop
        col := col+1;
      end loop;
      exit when (col > n);
      res(i) := col;
    end loop;
    return res;
  end Pivots;

  function Kernel ( pts,trivec : Matrix; rnk : natural32;
                    pivots : Standard_Integer_Vectors.Vector;
                    tol : double_float ) return VecVec is

    n : constant integer32 := pts'last(2);
    res : VecVec(1..n-integer32(rnk));
    mat : Matrix(1..integer32(rnk),1..integer32(rnk));
    rhs : Standard_Complex_Vectors.Vector(1..integer32(rnk));
    ind : integer32 := 0;
    info : integer32;
    ipvt : Standard_Integer_Vectors.Vector(1..integer32(rnk));
    eva : Complex_Number;

  begin
    for i in res'range loop                      -- initialize
      res(i) := new Standard_Complex_Vectors.Vector'(0..n => Create(0.0));
    end loop;
    for i in mat'range(1) loop                   -- set up linear system
      for j in mat'range(2) loop
        mat(i,j) := trivec(i,pivots(j));
      end loop;
    end loop;
    lufac(mat,integer32(rnk),ipvt,info);
    for j in 1..n loop                           -- compute kernel vectors
      if not Is_In(pivots,j) then
        ind := ind+1;
        res(ind)(j) := Create(1.0);              -- assign ones at nonpivots
        for i in 1..integer32(rnk) loop
          rhs(i) := -trivec(i,j);
        end loop;
        lusolve(mat,integer32(rnk),ipvt,rhs);    -- solve for pivot entries
        for i in 1..integer32(rnk) loop
          res(ind)(pivots(i)) := rhs(i);         -- assign pivot entries
        end loop;
        eva := Create(0.0);                      -- evaluate for constant term
        for k in pts'range(2) loop
          eva := eva + pts(pts'first,k)*res(ind)(k);
        end loop;
        res(ind)(0) := -eva;
      end if;
    end loop;
    return res;
  end Kernel;

end Standard_Linear_Spaces;
