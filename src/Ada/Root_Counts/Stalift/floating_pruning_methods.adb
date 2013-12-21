with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;
with Floating_Linear_Inequalities;       use Floating_Linear_Inequalities;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;

package body Floating_Pruning_Methods is

-- GENERAL AUXILIARIES :

  procedure Normalize ( v : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION : Divides every entry by v(v'last).

  begin
    for i in v'range loop
      v(i) := v(i)/v(v'last);
    end loop;
  end Normalize;

  function Convert ( fa : Face_Array ) return Array_of_Lists is

  -- DESCRIPTION :
  --   Converts the array of faces into an array of lists by
  --   converting the first element of each list of faces.

    res : Array_of_Lists(fa'range);

  begin
    for k in fa'range loop
      res(k) := Shallow_Create(fa(k).all);
    end loop;
    return res;
  end Convert;

-- AUXILIARIES FOR THE PRUNING ALGORITHMS :

  procedure Initialize ( n : in integer32; mat : out Matrix; 
                         ipvt : out Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Initializes an n*(n+1) matrix with zeroes and the pivoting vector
  --   with the entries 1..n.

  begin
    for i in 1..n loop
      for j in 1..n+1 loop
        mat(i,j) := 0.0;
      end loop;
    end loop;
    for i in 1..n loop
      ipvt(i) := i;
    end loop;
  end Initialize;

  function Number_of_Inequalities
             ( mix : Vector; lifted : Array_of_Lists ) return natural32 is

  -- DESCRIPTION :
  --   Returns the maximal number of inequalities for pruning.

    res : natural32 := 0;

  begin
    for k in lifted'range loop
      res := res + Length_Of(lifted(k)) - natural32(mix(k)) - 1;
    end loop;
    return res;
  end Number_of_Inequalities;

  procedure Ordered_Inequalities ( k : in integer32; mat : in out Matrix ) is

  -- DESCRIPTION :
  --   Defines k inequalities mat(k,k) - mat(k+1,k) >= 0.

  begin
    for i in mat'first(1)..k loop
      for j in mat'range(2) loop
        mat(i,j) := 0.0;
      end loop;
      mat(i,i) := 1.0; mat(i,i+1) := -1.0;
    end loop;
  end Ordered_Inequalities;

  procedure Check_and_Update
              ( mic : in Face_Array; lifted : in Array_of_Lists;
                m : in Matrix; ipvt : in Standard_Integer_Vectors.Vector;
                tol : in double_float;
                mixsub,mixsub_last : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Computes the normal to the points in pts, by solving the 
  --   linear system defined by m and ipvt.
  --   If the computed normal is an inner normal w.r.t. the lifted points,
  --   then the mixed subdivision will be updated with a new cell.

    use Standard_Floating_Vectors;
    v : Standard_Floating_Vectors.Vector(m'range(2)) := Solve(m,tol,ipvt);
    pts : Array_of_Lists(mic'range);

  begin
    if abs(v(v'last)) > tol then
      Normalize(v);
      if v(v'last) < 0.0
       then Min(v);
      end if;
     -- pts := Convert(mic);
     -- Update(pts,v,mixsub,mixsub_last);
      declare
        cell : Mixed_Cell;
      begin
        cell.nor := new Standard_Floating_Vectors.Vector'(v);
        cell.pts := new Array_of_Lists'(Convert(mic));
        cell.sub := null;
        Append(mixsub,mixsub_last,cell);
      end;
    end if;
  end Check_and_Update;

  procedure Create_Equalities
              ( n : in integer32; f : in Face; mat,ineq : in Matrix;
                tol : in double_float;
                ipvt : in Standard_Integer_Vectors.Vector;
                row,rowineq : in integer32;
                newmat,newineq : in out Matrix;
                newipvt : in out Standard_Integer_Vectors.Vector;
                newrow,newrowineq : in out integer32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Creates new equalities and uses them to eliminate unknowns in the
  --   matrices of equalities and inequalities.  Failure is reported when
  --   a zero row is encountered.  On entry, all new* variables must be
  --   initialized with the corresponding *-ones.

    shi : Standard_Floating_Vectors.Vector(1..n+1) := f(f'first).all;
    fl : boolean := false;
    pivot : integer32;

  begin
    for i in f'first+1..f'last loop
      newrow := newrow + 1;
      for j in f(i)'range loop
        newmat(newrow,j) := f(i)(j) - shi(j);
      end loop;
      Switch(newipvt,newrow,newmat);
      Upper_Triangulate(newrow,newmat,tol,newipvt,pivot);
      fl := (pivot = 0);
      exit when fl;
      Switch(newrow,pivot,ineq'first,rowineq,newineq);
    end loop;
    fail := fl;
  end Create_Equalities;

  procedure Complementary_Slackness
              ( tableau : in out matrix;
                tol : in double_float; feasible : out boolean ) is

    lastcol : constant integer32 := tableau'last(2)-1;
    rhs,sol : Standard_Floating_Vectors.Vector(tableau'range(1));
    columns : Vector(sol'range);

  begin
    for i in rhs'range loop
      rhs(i) := double_float(tableau(i,tableau'last(2)));
    end loop;
    Complementary_Slackness(tableau,lastcol,rhs,tol,sol,columns,feasible);
  end Complementary_Slackness;

  function Check_Feasibility
             ( k,m,n : integer32; ineq : Matrix;
               tol : double_float ) return boolean is
 
  -- DESCRIPTION :
  --   Returns true if -v(n+1) > 0 can be derived, false otherwise.

  -- ON ENTRY :
  --   k       current unknown that has been eliminated,
  --           for all i <= k : ineq(l,i) = 0, for l in ineq'first..m;
  --   m       number of inequalities;
  --   n       dimension;
  --   ineq    matrix of inequalities.

    tableau : Matrix(1..n-k+1,1..m+1);
    feasi : boolean;

  begin
    if m = 0 or n-k+1 < 2  then -- second clause in or is patched up
      feasi := false;
    else
      for i in k+1..n+1 loop
        for j in 1..m loop
          tableau(i-k,j) := ineq(j,i);
        end loop;
        tableau(i-k,m+1) := 0.0;
      end loop;
      tableau(n-k+1,m+1) := -1.0;
      Complementary_Slackness(tableau,tol,feasi);
    end if;
    return feasi;
  end Check_Feasibility;

  procedure Update_Inequalities
              ( k,rowmat1,rowmat2,n : in integer32;
                mat : in Matrix; ipvt : in Standard_Integer_Vectors.Vector;
                tol : in double_float;
                rowineq : in out integer32; ineq : in out Matrix;
                lifted : in Array_of_Lists; mic : in out Face_Array ) is

  -- DESCRIPTION :
  --   The inequalities will be updated w.r.t. the equality
  --   constraints on the inner normal.

    tmp : List;
    pt,shi : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in ineq'first..rowineq loop      -- update the old inequalities
      for j in rowmat1..rowmat2 loop
        Upper_Triangulate(j,mat,tol,i,ineq);
      end loop;
    end loop;
    shi := mic(k)(mic(k)'first);
    tmp := lifted(k);                      -- make new inequalities
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if not Is_In(mic(k),pt.all) then
        rowineq := rowineq + 1;
        for j in pt'range loop
          ineq(rowineq,j) := pt(j) - shi(j);
        end loop;
        Switch(ipvt,rowineq,ineq);
        for i in 1..rowmat2 loop
          Upper_Triangulate(i,mat,tol,rowineq,ineq);
        end loop;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Update_Inequalities;

-- CONSTRUCTION WITH PRUNING :

  procedure Gen1_Create
              ( n : in integer32; mix : in Vector; fa : in Array_of_Faces; 
                lifted : in Array_of_Lists; tol : in double_float;
                nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
                mixsub : out Mixed_Subdivision ) is

    res,res_last : Mixed_Subdivision;
    accu : Face_Array(fa'range);
    ma : Matrix(1..n,1..n+1);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    ineqrows : integer32;

    procedure Compute_Mixed_Cells
                 ( k,row : in integer32; mat : in Matrix;
                   ipvt : in Standard_Integer_Vectors.Vector;
                   rowineq : in integer32; ineq : in Matrix;
                   mic : in out Face_Array; continue : out boolean );

    -- DESCRIPTION :
    --   Backtrack recursive procedure to enumerate face-face combinations.

    -- ON ENTRY :
    --   k         index for current point set;
    --   row       number of rows already in the matrix;
    --   mat       matrix which determines the inner normal;
    --   ipvt      contains the pivoting information;
    --   rowineq   number of inequality constraints already in ineq;
    --   ineq      matrix for the inequality constraints on the
    --             inner normal v, of type <.,v> >= 0;
    --   mic       contains the current selected faces, up to k-1.

    -- ON RETURN :
    --   mic       updated selected faces.
    --   continue  indicates whether to continue the creation or not.

    procedure Process_Inequalities
                 ( k,rowmat1,rowmat2 : in integer32;
                   mat : in Matrix; ipvt : in Standard_Integer_Vectors.Vector;
                   rowineq : in out integer32; ineq : in out Matrix;
                   mic : in out Face_Array; cont : out boolean ) is

    -- DESCRIPTION :
    --   Updates inequalities and checks feasibility before proceeding.

    begin
      Update_Inequalities
        (k,rowmat1,rowmat2,n,mat,ipvt,tol,rowineq,ineq,lifted,mic);
      if Check_Feasibility(rowmat2,rowineq,n,ineq,tol)
       then nbfail(k) := nbfail(k) + 1.0;
            cont := true;
       else nbsucc(k) := nbsucc(k) + 1.0;
            Compute_Mixed_Cells(k+1,rowmat2,mat,ipvt,rowineq,ineq,mic,cont);
      end if;
    end Process_Inequalities;

    procedure Compute_Mixed_Cells
                 ( k,row : in integer32; mat : in Matrix;
                   ipvt : in Standard_Integer_Vectors.Vector;
                   rowineq : in integer32; ineq : in Matrix;
                   mic : in out Face_Array; continue : out boolean ) is

      old : constant Mixed_Subdivision := res_last;
      cont : boolean := true;
      tmpfa : Faces;

    begin
      if k > mic'last then
        Check_and_Update(mic,lifted,mat,ipvt,tol,res,res_last);
        if old /= res_last
         then Process(Head_Of(res_last),continue);
         else continue := true;
        end if;
      else
        tmpfa := fa(k);
        while not Is_Null(tmpfa) loop  -- enumerate faces of kth polytope
          mic(k) := Head_Of(tmpfa);
          declare                                 -- update the matrices
            fl : boolean;
            newipvt : Standard_Integer_Vectors.Vector(ipvt'range) := ipvt;
            newmat : Matrix(mat'range(1),mat'range(2)) := mat;
            newineq : Matrix(ineq'range(1),ineq'range(2)) := ineq;
            newrow : integer32 := row;
            newrowineq : integer32 := rowineq;
          begin
            Create_Equalities
              (n,mic(k),mat,ineq,tol,ipvt,row,rowineq,newmat,newineq,
               newipvt,newrow,newrowineq,fl);
            if fl
             then nbfail(k) := nbfail(k) + 1.0;
             else Process_Inequalities(k,row+1,newrow,newmat,newipvt,
                                       newrowineq,newineq,mic,cont);
            end if;
          end;
          tmpfa := Tail_Of(tmpfa);
          exit when not cont;
        end loop;
        continue := cont;
      end if;
    end Compute_Mixed_Cells;

  begin
    Initialize(n,ma,ipvt);
    ineqrows := integer32(Number_of_Inequalities(mix,lifted));
    declare
      ineq : matrix(1..ineqrows,1..n+1);
      cont : boolean;
    begin
      ineq(1,1) := 0.0;
      Compute_Mixed_Cells(accu'first,0,ma,ipvt,0,ineq,accu,cont);
    end;
    mixsub := res;
  end Gen1_Create;

  procedure Create
              ( n : in integer32; mix : in Vector; fa : in Array_of_Faces; 
                lifted : in Array_of_Lists; tol : in double_float;
                nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
		mixsub : out Mixed_Subdivision ) is

    res,res_last : Mixed_Subdivision;
    accu : Face_Array(fa'range);
    ma : Matrix(1..n,1..n+1);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    ineqrows : integer32;

    procedure Compute_Mixed_Cells
                 ( k,row : in integer32; mat : in Matrix;
                   ipvt : in Standard_Integer_Vectors.Vector;
                   rowineq : in integer32; ineq : in Matrix;
                   mic : in out Face_Array );

    -- DESCRIPTION :
    --   Backtrack recursive procedure to enumerate face-face combinations.

    -- ON ENTRY :
    --   k         index for current point set;
    --   row       number of rows already in the matrix;
    --   mat       matrix which determines the inner normal;
    --   ipvt      contains the pivoting information;
    --   rowineq   number of inequality constraints already in ineq;
    --   ineq      matrix for the inequality constraints on the
    --             inner normal v, of type <.,v> >= 0;
    --   mic       contains the current selected faces, up to k-1.

    -- ON RETURN :
    --   mic       updated selected faces.

    procedure Process_Inequalities
                 ( k,rowmat1,rowmat2 : in integer32;
                   mat : in Matrix; ipvt : in Standard_Integer_Vectors.Vector;
                   rowineq : in out integer32; ineq : in out matrix;
                   mic : in out Face_Array) is

    -- DESCRIPTION :
    --   Updates inequalities and checks feasibility before proceeding.

    begin
      Update_Inequalities
        (k,rowmat1,rowmat2,n,mat,ipvt,tol,rowineq,ineq,lifted,mic);
      if Check_Feasibility(rowmat2,rowineq,n,ineq,tol)
       then nbfail(k) := nbfail(k) + 1.0;
       else nbsucc(k) := nbsucc(k) + 1.0;
            Compute_Mixed_Cells(k+1,rowmat2,mat,ipvt,rowineq,ineq,mic);
      end if;
    end Process_Inequalities;

    procedure Compute_Mixed_Cells 
                 ( k,row : in integer32; mat : in Matrix;
                   ipvt : in Standard_Integer_Vectors.Vector;
                   rowineq : in integer32; ineq : in Matrix;
                   mic : in out Face_Array ) is

      tmpfa : Faces;

    begin
      if k > mic'last then
        Check_and_Update(mic,lifted,mat,ipvt,tol,res,res_last);
      else
        tmpfa := fa(k);
        while not Is_Null(tmpfa) loop  -- enumerate faces of kth polytope
          mic(k) := Head_Of(tmpfa);
          declare                                      -- update matrices
            fl : boolean;
            newipvt : Standard_Integer_Vectors.Vector(ipvt'range) := ipvt;
            newmat : Matrix(mat'range(1),mat'range(2)) := mat;
            newineq : Matrix(ineq'range(1),ineq'range(2)) := ineq;
            newrow : integer32 := row;
            newrowineq : integer32 := rowineq;
          begin
            Create_Equalities
                (n,mic(k),mat,ineq,tol,ipvt,row,rowineq,newmat,newineq,
                 newipvt,newrow,newrowineq,fl);
            if fl 
             then nbfail(k) := nbfail(k) + 1.0;
             else Process_Inequalities
                    (k,row+1,newrow,newmat,newipvt,newrowineq,newineq,mic);
            end if;
          end;
          tmpfa := Tail_Of(tmpfa);
        end loop;
      end if;
    end Compute_Mixed_Cells;

  begin
    Initialize(n,ma,ipvt);
    ineqrows := integer32(Number_of_Inequalities(mix,lifted));
    declare
      ineq : Matrix(1..ineqrows,1..n+1);
    begin
      ineq(1,1) := 0.0;
      Compute_Mixed_Cells(accu'first,0,ma,ipvt,0,ineq,accu);
    end;
    mixsub := res;
  end Create;

end Floating_Pruning_Methods;
