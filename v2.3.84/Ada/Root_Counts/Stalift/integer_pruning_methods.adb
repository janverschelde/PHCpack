with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Matrices;
with Standard_Integer_Norms;             use Standard_Integer_Norms;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Standard_Integer_Linear_Solvers;    use Standard_Integer_Linear_Solvers;
with Standard_Integer_Linear_Equalities; use Standard_Integer_Linear_Equalities;
with Integer_Linear_Inequalities;        use Integer_Linear_Inequalities;
with Floating_Linear_Inequality_Solvers; use Floating_Linear_Inequality_Solvers;

package body Integer_Pruning_Methods is

-- GENERAL AUXILIARIES :

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
                         ipvt : out Vector ) is

  -- DESCRIPTION :
  --   Initializes the matrix of the equality constraints on the normals
  --   and sets the pivoting vector ipvt to 1..n+1.

  begin
    for i in 1..n+1 loop
      ipvt(i) := i;
      for j in 1..n+1 loop
        mat(i,j) := 0;
      end loop;
    end loop;
  end Initialize;

  function Number_of_Inequalities
             ( mix : Vector; lifted : Array_of_Lists ) return natural32 is

  -- DESCRIPTION :
  --   Returns the maximal number of inequalities on the inner normal.

    res : natural32 := 0;

  begin
    for k in lifted'range loop
      res := res + Length_Of(lifted(k)) - natural32(mix(k)) - 1;
    end loop;
    return res;
  end Number_of_Inequalities;

  procedure Ordered_Inequalities ( k : in integer32; mat : in out matrix ) is

  -- DESCRIPTION :
  --   Defines k inequalities mat(k,k) - mat(k+1,k) >= 0.

  begin
    for i in mat'first(1)..k loop
      for j in mat'range(2) loop
        mat(i,j) := 0;
      end loop;
      mat(i,i) := 1; mat(i,i+1) := -1;
    end loop;
  end Ordered_Inequalities;

  procedure Check_and_Update
                ( mic : in Face_Array; lifted : in Array_of_Lists;
                  m : in matrix; ipvt : in Vector;
                  mixsub,mixsub_last : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Computes the normal to the points in mic, by solving the 
  --   linear system defined by m and ipvt.  Updates the mixed subdivision
  --   if the computed normal is an inner normal.

    v : Vector(m'range(2)) := Solve(m,ipvt);
    pts : Array_of_Lists(mic'range);

  begin
    if v(v'last) /= 0 then
      Normalize(v);
      if v(v'last) < 0
       then Min(v);
      end if;
      pts := Convert(mic);
      Update(pts,v,mixsub,mixsub_last);
    end if;
  end Check_and_Update;

  procedure Create_Equalities
                ( n : in integer32; f : in Face; mat,ineq : in matrix;
                  ipvt : in Vector; row,rowineq : in integer32;
                  newmat,newineq : in out matrix; newipvt : in out Vector;
                  newrow,newrowineq : in out integer32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Creates new equalities and uses them to eliminate unknowns in the
  --   matrices of equalities and inequalities.  Failure is reported when
  --   a zero row is encountered.  On entry, all new* variables must be
  --   initialized with the corresponding *-ones.

    shi : constant Vector(1..n+1) := f(f'first).all;
    fl : boolean := false;
    pivot : integer32;

  begin
    for i in f'first+1..f'last loop
      newrow := newrow + 1;
      for j in f(i)'range loop
        newmat(newrow,j) := f(i)(j) - shi(j);
      end loop;
      Triangulate(newrow,newmat,newipvt,pivot);
      fl := (pivot = 0);
      exit when fl;
      Switch(newrow,pivot,ineq'first,rowineq,newineq);
      --Scale(newrow,newmat);
    end loop;
    fail := fl;
  end Create_Equalities;

  function Check_Feasibility
              ( k,m,n : integer32; ineq : matrix ) return boolean is
 
  -- DESCRIPTION :
  --   Returns true if -v(n+1) > 0 can be derived, false otherwise.

  -- ON ENTRY :
  --   k      current unknown that has been eliminated,
  --           for all i <= k : ineq(l,i) = 0, for l in ineq'first..m;
  --   m      number of inequalities;
  --   n      dimension;
  --   ineq   matrix of inequalities.

    tableau : matrix(1..n-k+1,1..m+1);
    feasi : boolean;

  begin
    if m = 0 then
      feasi := false;
    else
      for i in k+1..n+1 loop
        for j in 1..m loop
          tableau(i-k,j) := ineq(j,i);
        end loop;
        tableau(i-k,m+1) := 0;
      end loop;
      tableau(n-k+1,m+1) := -1;
      Integer_Complementary_Slackness(tableau,feasi);
    end if;
    return feasi;
  end Check_Feasibility;

  function New_Check_Feasibility
              ( k,m,n : integer32; ineq : matrix; tol : double_float;
                solu : Standard_Floating_Vectors.Vector ) return boolean is

  -- DESCRIPTION :
  --   Returns true if -v(n+1) > 0 can be derived, false otherwise.

  -- ON ENTRY :
  --   k      current unknown that has been eliminated,
  --           for all i <= k : ineq(l,i) = 0, for l in ineq'first..m;
  --   m      number of inequalities;
  --   n      dimension;
  --   ineq   matrix of inequalities.

    tableau : Standard_Floating_Matrices.matrix(1..n-k+1,1..m);
    tabsolu : Standard_Floating_Vectors.Vector(1..n-k);
    cols : Vector(tabsolu'range);
    kfail : integer32;
    max,tmpabs : double_float;
         -- used for scaling of the last component of the inner normal

  begin
    for i in k+1..n loop
      for j in 1..m loop
        tableau(i-k,j) := double_float(ineq(j,i));
      end loop;
      tabsolu(i-k) := solu(i);
    end loop;
    max := 0.0;
    for j in 1..m loop
      tableau(n-k+1,j) := -double_float(ineq(j,n+1));
      tmpabs := abs(tableau(n-k+1,j));
      if tmpabs > max
       then max := tmpabs;
      end if;
    end loop;
    if max > 1.0 then
      for j in 1..m loop
        tableau(n-k+1,j) := tableau(n-k+1,j)/max;
      end loop;
    end if;
    Scale(tableau,tol);
    Solve(tableau,tol,tabsolu,kfail,cols);
    return (kfail <= m);
  end New_Check_Feasibility;

  procedure Update_Inequalities
               ( k,rowmat1,rowmat2,n : in integer32;
                 mat : in matrix; ipvt : in vector; 
                 rowineq : in out integer32; ineq : in out matrix;
                 lifted : in Array_of_Lists; mic : in out Face_Array ) is

  -- DESCRIPTION :
  --   The inequalities will be updated w.r.t. the equality
  --   constraints on the inner normal.

    tmp : List;
    pt,shi : Link_to_Vector;

  begin
    for i in ineq'first..rowineq loop   -- triangulate old inequalities
      for j in rowmat1..rowmat2 loop
        Triangulate(j,mat,i,ineq);
      end loop;
      --Scale(i,ineq);
    end loop;
    tmp := lifted(k);                         -- create and triangulate 
    shi := mic(k)(mic(k)'first);                    -- new inequalities
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if not Is_In(mic(k),pt.all)
       then rowineq := rowineq + 1;
            for j in pt'range loop
              ineq(rowineq,j) := pt(j) - shi(j);
            end loop;
            Switch(ipvt,rowineq,ineq);
            for i in 1..rowmat2 loop
              Triangulate(i,mat,rowineq,ineq);
              --Scale(rowineq,ineq);
            end loop;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Update_Inequalities;

  procedure Eliminate
              ( k : in integer32; m : in Matrix; tol : in double_float;
                x : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Eliminates the kth component of x by using the matrix.

   fac : double_float;

  begin
    if abs(x(k)) > tol then
      fac := x(k)/double_float(m(k,k));
      for j in k..x'last loop
        x(j) := x(j) - fac*double_float(m(k,j));
      end loop;
    end if;
  end Eliminate;

  procedure Switch ( l,pivot : in integer32;
                     x : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Applies the pivotation information to the vector x.

    tmp : double_float;

  begin
    if l /= pivot
     then tmp := x(l); x(l) := x(pivot); x(pivot) := tmp;
    end if;
  end Switch;

  procedure Eliminate ( rowmat1,rowmat2 : in integer32; mat : in Matrix;
                        ipvt : in Vector; tol : in double_float;
                        x : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Returns an eliminated solution vector.

  begin
    for k in rowmat1..rowmat2 loop
      Switch(k,ipvt(k),x);
      Eliminate(k,mat,tol,x);
    end loop;
  end Eliminate;

-- GENERAL CONSTRUCTOR in several versions :

  -- procedure Compute_Mixed_Cells ( );  -- general specification.

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
  --   testsolu  test solution to check current set of inequalilities;
  --   mic       contains the current selected faces, up to k-1.

  -- ON RETURN :
  --   mic       updated selected faces;
  --   issupp    true if current tuple of faces is supported;
  --   continue  indicates whether to continue the creation or not.

-- CONSTRUCTION WITH PRUNING :

  procedure Gen1_Create_CS
               ( n : in integer32; mix : in Vector; fa : in Array_of_Faces; 
                 lifted : in Array_of_Lists;
                 nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
                 mixsub : out Mixed_Subdivision ) is

    res,res_last : Mixed_Subdivision;
    accu : Face_Array(fa'range);
    ma : matrix(1..n+1,1..n+1);
    ipvt : vector(1..n+1);
    ineqrows : integer32;

    procedure Compute_Mixed_Cells
                 ( k,row : in integer32; mat : in matrix; ipvt : in vector;
                   rowineq : in integer32; ineq : in matrix;
                   mic : in out Face_Array; continue : out boolean );

    -- DESCRIPTION :
    --   Backtrack recursive procedure to enumerate face-face combinations.

    procedure Process_Inequalities
                 ( k,rowmat1,rowmat2 : in integer32;
                   mat : in matrix; ipvt : in vector;
                   rowineq : in out integer32; ineq : in out matrix;
                   mic : in out Face_Array; cont : out boolean ) is

    -- DESCRIPTION :
    --   Updates inequalities and checks feasibility before proceeding.

    begin
      Update_Inequalities(k,rowmat1,rowmat2,n,mat,ipvt,rowineq,ineq,lifted,mic);
      if Check_Feasibility(rowmat2,rowineq,n,ineq)
       then nbfail(k) := nbfail(k) + 1.0;
            cont := true;
       else nbsucc(k) := nbsucc(k) + 1.0;
            Compute_Mixed_Cells(k+1,rowmat2,mat,ipvt,rowineq,ineq,mic,cont);
      end if;
    end Process_Inequalities;

    procedure Compute_Mixed_Cells
                 ( k,row : in integer32; mat : in matrix; ipvt : in vector;
                   rowineq : in integer32; ineq : in matrix;
                   mic : in out Face_Array; continue : out boolean ) is

      old : Mixed_Subdivision := res_last;
      cont : boolean := true;
      tmpfa : Faces;

    begin
      if k > mic'last then
        Check_and_Update(mic,lifted,mat,ipvt,res,res_last);
        if old /= res_last
         then Process(Head_Of(res_last),continue);
         else continue := true;
        end if;
      else
        tmpfa := fa(k);
        while not Is_Null(tmpfa) loop  -- enumerate faces of kth polytope
          mic(k) := Head_Of(tmpfa);
          declare                                    -- update matrices
            fl : boolean;
            newipvt : vector(ipvt'range) := ipvt;
            newmat : matrix(mat'range(1),mat'range(2)) := mat;
            newineq : matrix(ineq'range(1),ineq'range(2)) := ineq;
            newrow : integer32 := row;
            newrowineq : integer32 := rowineq;
          begin
            Create_Equalities(n,mic(k),mat,ineq,ipvt,row,rowineq,
                              newmat,newineq,newipvt,newrow,newrowineq,fl);
            if fl
             then nbfail(k) := nbfail(k) + 1.0;
             else Process_Inequalities(k,row+1,newrow,newmat,newipvt,
                                       newrowineq,newineq,mic,cont);
            end if;
          end;
          exit when not cont;
          tmpfa := Tail_Of(tmpfa);
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
      ineq(1,1) := 0;
      Compute_Mixed_Cells(accu'first,0,ma,ipvt,0,ineq,accu,cont);
    end;
    mixsub := res;
  end Gen1_Create_CS;

  procedure Create_CS
              ( n : in integer32; mix : in Vector; fa : in Array_of_Faces; 
                lifted : in Array_of_Lists;
                nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
                mixsub : out Mixed_Subdivision ) is

    res,res_last : Mixed_Subdivision;
    accu : Face_Array(fa'range);
    ma : matrix(1..n+1,1..n+1);
    ipvt : vector(1..n+1);
    ineqrows : integer32;

    procedure Compute_Mixed_Cells
                 ( k,row : in integer32; mat : in Matrix; ipvt : in Vector;
                   rowineq : in integer32; ineq : in Matrix;
                   mic : in out Face_Array );

    -- DESCRIPTION :
    --   Backtrack recursive procedure to enumerate face-face combinations.

    procedure Process_Inequalities
                 ( k,rowmat1,rowmat2 : in integer32;
                   mat : in matrix; ipvt : in vector;
                   rowineq : in out integer32; ineq : in out matrix;
                   mic : in out Face_Array ) is

    -- DESCRIPTION :
    --   Updates inequalities and checks feasibility before proceeding.

    begin
      Update_Inequalities(k,rowmat1,rowmat2,n,mat,ipvt,rowineq,ineq,lifted,mic);
      if Check_Feasibility(rowmat2,rowineq,n,ineq)
       then nbfail(k) := nbfail(k) + 1.0;
       else nbsucc(k) := nbsucc(k) + 1.0;
            Compute_Mixed_Cells(k+1,rowmat2,mat,ipvt,rowineq,ineq,mic);
      end if;
    end Process_Inequalities;

    procedure Compute_Mixed_Cells 
                 ( k,row : in integer32; mat : in Matrix; ipvt : in Vector;
                   rowineq : in integer32; ineq : in Matrix;
                   mic : in out Face_Array ) is

      tmpfa : Faces;

    begin
      if k > mic'last then
        Check_and_Update(mic,lifted,mat,ipvt,res,res_last);
      else
        tmpfa := fa(k);
        while not Is_Null(tmpfa) loop   -- enumerate faces of kth polytype
          mic(k) := Head_Of(tmpfa);
          declare                                     -- update matrices
            fl : boolean;
            newipvt : vector(ipvt'range) := ipvt;
            newmat : matrix(mat'range(1),mat'range(2)) := mat;
            newineq : matrix(ineq'range(1),ineq'range(2)) := ineq;
            newrow : integer32 := row;
            newrowineq : integer32 := rowineq;
          begin
            Create_Equalities(n,mic(k),mat,ineq,ipvt,row,rowineq,newmat,
                              newineq,newipvt,newrow,newrowineq,fl);
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
      ineq : matrix(1..ineqrows,1..n+1);
    begin
      ineq(1,1) := 0;
      Compute_Mixed_Cells(accu'first,0,ma,ipvt,0,ineq,accu);
    end;
    mixsub := res;
  end Create_CS;

  procedure New_Create_CS
              ( n : in integer32; mix : in Vector; fa : in Array_of_Faces; 
                lifted : in Array_of_Lists;
                nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
		mixsub : out Mixed_Subdivision ) is

    res,res_last : Mixed_Subdivision;
    accu : Face_Array(fa'range);
    ma : matrix(1..n+1,1..n+1);
    ipvt : vector(1..n+1);
    tol : constant double_float := double_float(n)*10.0**(-11);
    ineqrows : integer32;

    procedure Compute_Mixed_Cells
                 ( k,row : in integer32; mat : in matrix; ipvt : in vector;
                   rowineq : in integer32; ineq : in matrix;
                   testsolu : in Standard_Floating_Vectors.Vector;
                   mic : in out Face_Array );

    -- DESCRIPTION :
    --   Backtrack recursive procedure to enumerate face-face combinations.

    procedure Process_Inequalities
                 ( k,rowmat1,rowmat2 : in integer32;
                   mat : in matrix; ipvt : in vector;
                   rowineq : in out integer32; ineq : in out matrix;
                   testsolu : in Standard_Floating_Vectors.Vector;
                   mic : in out Face_Array ) is

    -- DESCRIPTION :
    --   Updates the inequalities and checks feasibility before proceeding.

      newsolu : Standard_Floating_Vectors.Vector(testsolu'range) := testsolu;

    begin
      Update_Inequalities
        (k,rowmat1,rowmat2,n,mat,ipvt,rowineq,ineq,lifted,mic);
      Eliminate(rowmat1,rowmat2,mat,ipvt,tol,newsolu);
      if New_Check_Feasibility(rowmat2,rowineq,n,ineq,tol,newsolu)
       then nbfail(k) := nbfail(k) + 1.0;
       else nbsucc(k) := nbsucc(k) + 1.0;
            Compute_Mixed_Cells
              (k+1,rowmat2,mat,ipvt,rowineq,ineq,testsolu,mic);
      end if;
    end Process_Inequalities;

    procedure Compute_Mixed_Cells 
                 ( k,row : in integer32; mat : in matrix; ipvt : in vector;
                   rowineq : in integer32; ineq : in matrix;
                   testsolu : in Standard_Floating_Vectors.Vector;
                   mic : in out Face_Array ) is

    -- DESCRIPTION :
    --   Backtrack recursive procedure to enumerate face-face combinations.

      tmpfa : Faces;

    begin
      if k > mic'last then
        Check_and_Update(mic,lifted,mat,ipvt,res,res_last);
      else
        tmpfa := fa(k); 
        while not Is_Null(tmpfa) loop -- enumerate faces of kth polytope
          mic(k) := Head_Of(tmpfa);
          declare                                     -- update matrices
            fl : boolean;
            newipvt : vector(ipvt'range) := ipvt;
            newmat : matrix(mat'range(1),mat'range(2)) := mat;
            newineq : matrix(ineq'range(1),ineq'range(2)) := ineq;
            newrow : integer32 := row;
            newrowineq : integer32 := rowineq;
          begin
            Create_Equalities(n,mic(k),mat,ineq,ipvt,row,rowineq,newmat,
                              newineq,newipvt,newrow,newrowineq,fl);
            if fl 
             then nbfail(k) := nbfail(k) + 1.0;
             else Process_Inequalities(k,row+1,newrow,newmat,newipvt,
                                       newrowineq,newineq,testsolu,mic);
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
      ineq : matrix(1..ineqrows,1..n+1);
      solu : Standard_Floating_Vectors.Vector(1..n+1);
    begin
      ineq(1,1) := 0;
      solu := (solu'range => 0.0);
      Compute_Mixed_Cells(accu'first,0,ma,ipvt,0,ineq,solu,accu);
    end;
    mixsub := res;
  end New_Create_CS;

-- AUXILIARIES FOR THE CONSTRAINT PRUNING :

  function Is_Supported ( f : Face; normal : Vector; supp : integer32 )
                        return boolean is

    ip : integer32;

  begin
    for i in f'range loop
      ip := f(i).all*normal;
      if ip /= supp
       then return false;
      end if;
    end loop;
    return true;
  end Is_Supported;

  function Is_Supported ( f : Face; k : integer32; normals,suppvals : List )
                        return boolean is

  -- DESCRIPTION :
  --   Returns true if there exists a normal for which the inner product
  --   with each point in the face equals the corresponding supporting value
  --   of the kth component.

    tmpnor : List := normals;
    tmpsup : List := suppvals;
    support : integer32;

  begin
    while not Is_Null(tmpnor) loop
      support := Head_Of(tmpsup)(k);
      if Is_Supported(f,Head_Of(tmpnor).all,support)
       then return true;
       else tmpnor := Tail_Of(tmpnor);
            tmpsup := Tail_Of(tmpsup);
      end if;
    end loop;
    return false;
  end Is_Supported;
 
  procedure Update ( mic : in Mixed_Cell; normals,suppvals : in out List ) is

  -- DESCRIPTION :
  --   Updates the list of normals and supporting values with the information
  --   of a mixed cell.

    normal : Link_to_Vector := new Vector'(mic.nor.all);
    suppval : Link_to_Vector := new Vector(mic.pts'range);

  begin
    Construct(normal,normals);
    for i in suppval'range loop
      suppval(i) := normal*Head_Of(mic.pts(i));
    end loop;
    Construct(suppval,suppvals);
  end Update;

  procedure Create_CCS
                 ( n : in integer32; mix : in Vector; fa : in Array_of_Faces; 
                   lifted : in Array_of_Lists;
                   nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
                   normals,suppvals : in out List;
                   mixsub : out Mixed_Subdivision ) is

    res,res_last : Mixed_Subdivision;
    accu : Face_Array(fa'range);
    ma : matrix(1..n+1,1..n+1);
    ipvt : vector(1..n+1);
    ineqrows : integer32;

    procedure Compute_Mixed_Cells
                 ( k,row : in integer32; mat : in matrix; ipvt : in vector;
                   rowineq : in integer32; ineq : in matrix;
                   mic : in out Face_Array; issupp : in out boolean );

    -- DESCRIPTION :
    --   Backtrack recursive procedure to enumerate face-face combinations.

    procedure Process_Inequalities
                 ( k,rowmat1,rowmat2 : in integer32;
                   mat : in matrix; ipvt : in vector;
                   rowineq : in out integer32; ineq : in out matrix;
                   mic : in out Face_Array; issupp : in out boolean ) is

    -- DESCRIPTION :
    --   Updates inequalities and checks feasibility before proceeding.

    begin
      Update_Inequalities(k,rowmat1,rowmat2,n,mat,ipvt,rowineq,ineq,lifted,mic);
      if issupp
       then issupp := Is_Supported(mic(k),k,normals,suppvals);
      end if;
      if not issupp and then Check_Feasibility(rowmat2,rowineq,n,ineq)
       then nbfail(k) := nbfail(k) + 1.0;
       else nbsucc(k) := nbsucc(k) + 1.0;
            Compute_Mixed_Cells(k+1,rowmat2,mat,ipvt,rowineq,ineq,mic,issupp);
      end if;
    end Process_Inequalities;

    procedure Compute_Mixed_Cells 
                 ( k,row : in integer32; mat : in matrix; ipvt : in vector;
                   rowineq : in integer32; ineq : in matrix;
                   mic : in out Face_Array; issupp : in out boolean ) is

    -- DESCRIPTION :
    --   Backtrack recursive procedure to enumerate face-face combinations.

      old : Mixed_Subdivision;
      tmpfa : Faces;

    begin
      if k > mic'last then
        old := res_last;
        Check_and_Update(mic,lifted,mat,ipvt,res,res_last);
        if old /= res_last
         then Update(Head_Of(res_last),normals,suppvals);
        end if;
      else
        tmpfa := fa(k);
        while not Is_Null(tmpfa) loop  -- enumerate faces of kth polytope
          mic(k) := Head_Of(tmpfa);
          declare                                      -- update matrices
            fl : boolean;
            newipvt : vector(ipvt'range) := ipvt;
            newmat : matrix(mat'range(1),mat'range(2)) := mat;
            newineq : matrix(ineq'range(1),ineq'range(2)) := ineq;
            newrow : integer32 := row;
            newrowineq : integer32 := rowineq;
          begin
            Create_Equalities(n,mic(k),mat,ineq,ipvt,row,rowineq,newmat,
                              newineq,newipvt,newrow,newrowineq,fl);
            if fl 
             then nbfail(k) := nbfail(k) + 1.0;
             else Process_Inequalities(k,row+1,newrow,newmat,newipvt,
                                       newrowineq,newineq,mic,issupp);
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
      ineq : matrix(1..ineqrows,1..n+1);
      supported : boolean := true;
    begin
      ineq(1,1) := 0;
      Compute_Mixed_Cells(accu'first,0,ma,ipvt,0,ineq,accu,supported);
    end;
    mixsub := res;
  end Create_CCS;

end Integer_Pruning_Methods;
