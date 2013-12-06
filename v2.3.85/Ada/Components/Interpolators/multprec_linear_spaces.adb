with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_Linear_Solvers;    use Multprec_Complex_Linear_Solvers;

--with integer_io;                         use integer_io;
--with Multprec_Complex_Number_Tools;      use Multprec_Complex_Number_Tools;
--with Standard_Complex_Numbers;
--with Standard_Complex_Matrices;
--with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
--with Standard_Complex_Linear_Solvers;
--with Multprec_Complex_Matrices_io;       use Multprec_Complex_Matrices_io;

package body Multprec_Linear_Spaces is

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

  procedure Normalize_Rows
              ( mat : in out Matrix; n,m : in integer32; 
                size : in natural32; tol : in double_float ) is

  -- DESCRIPTION :
  --   Divides each row by its first nonzero element.

    temp : Floating_Number;
    quot : Complex_Number;
    done : boolean;

  begin
    for i in 1..n loop
      done := false;
      for j in 1..m loop
        temp := AbsVal(mat(i,j));
        if temp > tol then
          for k in j+1..m loop
            Div(mat(i,k),mat(i,j));
          end loop;
          -- Clear(mat(i,j));
          -- mat(i,j) := Create(1);
          -- Set_Size(mat(i,j),size);
          quot := mat(i,j)/mat(i,j);
          Copy(quot,mat(i,j));
          Clear(quot);
          done := true;
        end if;
        Clear(temp);
        exit when done;
      end loop;
    end loop;
  end Normalize_Rows;

  procedure Rank ( nbvecs,n : in integer32; size : in natural32;
                   pts : in Matrix; tol : in double_float;
                   trivec : in out Matrix; rnk : out natural32 ) is

    col : integer32;
    continue : boolean;
    temp : Floating_Number;

   -- use Standard_Complex_Numbers;

   -- stamat : Standard_Complex_Matrices.Matrix(1..nbvecs+1,1..n);
   -- stavec : Standard_Complex_Matrices.Matrix(1..nbvecs,1..n);

  begin
   -- for i in 1..nbvecs+1 loop
   --   for j in 1..n loop
   --     stamat(i,j) := Round(pts(i,j));
   --   end loop;
   -- end loop;
    for i in 1..nbvecs loop
      for j in 1..n loop
        trivec(i,j) := pts(i,j) - pts(nbvecs+1,j);
       -- stavec(i,j) := Round(trivec(i,j)); --stamat(i,j) - stamat(nbvecs+1,j);
      end loop;
    end loop;
   -- Normalize_Rows(trivec,nbvecs,n,size,tol);
   -- put_line(file,"The matrix before triangulation : "); put(file,trivec);
   -- put_line(file,"the standard matrix : "); put(file,stavec);
   -- Triangulate(file,trivec,tol,size,nbvecs,n);
    Triangulate(trivec,tol,size,nbvecs,n);
   -- put_line(file,"The matrix after triangulation : "); put(file,trivec);
   -- Standard_Complex_Linear_Solvers.Triangulate(stavec,tol,nbvecs,n);
   -- put_line(file,"The standard triangulation : "); put(file,stavec);
   -- put_line(file,"Structure of the standard triangulation : ");
   -- for i in 1..nbvecs loop
   --   for j in 1..n loop
   --     if Standard_Complex_Numbers.AbsVal(stavec(i,j)) > tol
   --      then put(file,"*");
   --      else put(file,"0");
   --     end if;
   --   end loop;
   --   new_line(file);
   -- end loop;
    rnk := 0;
    col := 1;
    for i in 1..nbvecs loop
      continue := true;
      while (col <= n) and continue loop
        temp := AbsVal(trivec(i,col));
        if temp < tol
         then col := col + 1;
         else continue := false;
        end if;
        Clear(temp);
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
    continue : boolean;
    temp : Floating_Number;

  begin
    for i in 1..nbvecs loop
      continue := true;
      while (col <= n) and continue loop
        temp := AbsVal(trivec(i,col));
        if temp < tol
         then col := col + 1;
         else continue := false;
        end if;
        Clear(temp);
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
    rhs : Multprec_Complex_Vectors.Vector(1..integer32(rnk));
    ind : integer32 := 0;
    info : integer32;
    ipvt : Standard_Integer_Vectors.Vector(1..integer32(rnk));
    acc,eva : Complex_Number;

  begin
    for i in res'range loop                     -- initialize
      res(i) := new Multprec_Complex_Vectors.Vector'
                      (0..n => Create(integer(0)));
    end loop;
    for i in mat'range(1) loop                  -- set up linear system
      for j in mat'range(2) loop
        Copy(trivec(i,pivots(j)),mat(i,j));
      end loop;
    end loop;
    lufac(mat,integer32(rnk),ipvt,info);
    for j in 1..n loop                          -- compute kernel vectors
      if not Is_In(pivots,j) then
        ind := ind+1;
        res(ind)(j) := Create(integer(1));      -- assign ones at nonpivots
        for i in 1..integer32(rnk) loop
          Copy(trivec(i,j),rhs(i));
          Min(rhs(i));
        end loop;
        lusolve(mat,integer32(rnk),ipvt,rhs);   -- solve for pivot entries
        for i in 1..integer32(rnk) loop
          Copy(rhs(i),res(ind)(pivots(i)));     -- assign pivot entries
        end loop;
        eva := Create(integer(0));              -- evaluate for constant term
        for k in pts'range(2) loop
          acc := pts(pts'first,k)*res(ind)(k);
          Add(eva,acc);
          Clear(acc);
        end loop;
        Min(eva);
        Copy(eva,res(ind)(0));
        Clear(eva);
      end if;
    end loop;
    Clear(mat);
    return res;
  end Kernel;

end Multprec_Linear_Spaces;
