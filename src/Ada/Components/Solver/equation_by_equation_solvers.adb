with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Solution_Splitters;        use Standard_Solution_Splitters;
with Extrinsic_Diagonal_Continuation;    use Extrinsic_Diagonal_Continuation;
with Intrinsic_Diagonal_Continuation;    use Intrinsic_Diagonal_Continuation;
with Hypersurfaces_and_Filters;          use Hypersurfaces_and_Filters;
with Flow_Tables;                        use Flow_Tables;
with Intrinsic_Witness_Sets_io;          use Intrinsic_Witness_Sets_io;

package body Equation_by_Equation_Solvers is 

-- AUXILIARIES :

  function Create_Plane ( b,v : Vector ) return Matrix is

  -- DESCRIPTION :
  --   Returns a matrix with two columns at 0: b, and at 1: v.

    res : Matrix(b'range,0..1);

  begin
    for i in b'range loop
      res(i,0) := b(i);
      res(i,1) := v(i);
    end loop;
    return res;
  end Create_Plane;

  generic

    with function f ( k : integer32; x : Vector ) return Complex_Number;

  procedure G_Hypersurface_Roots
               ( file : in file_type; n,i,d : in integer32;
                 b,v : in Vector; s : out Solution_List );

  -- DESCRIPTION :
  --   Generic version to compute witness set for i-th hypersurface
  --   of degree d in n variables, cut with b + t*v.

  procedure G_Hypersurface_Roots
               ( file : in file_type; n,i,d : in integer32;
                 b,v : in Vector; s : out Solution_List ) is

    nrm : double_float;

    function Eval ( x : Vector ) return Complex_Number is
    begin
      return f(i,x);
    end Eval;
   -- procedure R_Find_Roots is new RG_Hypersurface_Witness_Set(Eval);
    procedure S_Find_Roots is new SG_Hypersurface_Witness_Set(Eval);
              
  begin
    put(file,"Computing witness set for equation "); put(file,i,1);
    put_line(file," ...");
   -- R_Find_Roots(file,n,d,b,v,s,nrm);
    S_Find_Roots(natural32(n),natural32(d),b,v,s,nrm);
    put(file," ... found "); put(file,Length_Of(s),1);
    put(file," roots with max residual norm :");
    put(file,nrm,3); put_line(file,".");
   -- put(file,Length_Of(s(i)),1,s(i));
  end G_Hypersurface_Roots;

  procedure P_Hypersurface_Roots
               ( file : in file_type; n,i : in integer32;
                 p : in Poly; ep : in Eval_Poly; b,v : in Vector;
                 s : out Solution_List ) is

  -- DESCRIPTON :
  --   Calls the routine to compute a witness set for the i-th hypersurface.
  --   Writes some extra diagnostics on screen.

    nrm : double_float;
              
  begin
    put(file,"Computing witness set for equation "); put(file,i,1);
    put_line(file," ...");
    SP_Hypersurface_Witness_Set(natural32(n),natural32(Degree(p)),ep,b,v,s,nrm);
    put(file," ... found "); put(file,Length_Of(s),1);
    put(file," roots with max residual norm :");
    put(file,nrm,3); put_line(file,".");
   -- put(file,Length_Of(s(i)),1,s(i));
  end P_Hypersurface_Roots;

  generic

    with function f ( k : integer32; x : Vector ) return Complex_Number;

  procedure G_Filter
               ( file : in file_type; i : in integer32;
                 s : in out Solution_List; b,v : in Vector;
                 tol : in double_float );

  -- DESCRIPTION :
  --   Filters a witness set for one hypersurface using previous
  --   equations, discarding those points which satisfy any of
  --   those previous equations.

  -- ON ENTRY :
  --   file      to write diagnostics on;
  --   i         number of the current polynomial equation;
  --   s         points in the current witness set;
  --   b         offset vector of a random line;
  --   v         direction vector of a random line;
  --   ep        equations used for testing witness set;
  --   tol       tolerance to decide whether number is zero.

  -- ON RETURN :
  --   s         filtered witness set.

  procedure G_Filter
               ( file : in file_type; i : in integer32;
                 s : in out Solution_List; b,v : in Vector;
                 tol : in double_float ) is

   -- procedure Filter is new RG_Filter(f);
    procedure Filter is new SG_Filter(f);

  begin
    put(file,"Filtering witness set "); put(file,i,1);
    put_line(file," with previous equations ...");
   -- Filter(file,i-1,s,b,v,tol);
    Filter(natural32(i-1),s,b,v,tol);
    put(file,"  ... retained "); put(file,Length_Of(s),1);
    put_line(file," points.");
  end G_Filter;

  procedure P_Filter
               ( file : in file_type; i : in integer32;
                 s : in out Solution_List; b,v : in Vector;
                 ep : in Eval_Poly_Sys; tol : in double_float ) is

  -- DESCRIPTION :
  --   Filters a witness set for one hypersurface using previous
  --   equations, discarding those points which satisfy any of
  --   those previous equations.

  -- ON ENTRY :
  --   file      to write diagnostics on;
  --   i         number of the current polynomial equation;
  --   s         points in the current witness set;
  --   b         offset vector of a random line;
  --   v         direction vector of a random line;
  --   ep        equations used for testing witness set;
  --   tol       tolerance to decide whether number is zero.

  -- ON RETURN :
  --   s         filtered witness set.

  begin
    put(file,"Filtering witness set "); put(file,i,1);
    put_line(file," with previous equations ...");
    SP_Filter(s,ep,b,v,tol);
   -- RP_Filter(file,s,ep,b,v,tol); -- reporting version
    put(file,"  ... retained "); put(file,Length_Of(s),1);
    put_line(file," points.");
--  exception
--    when others => put_line("exception happened in P_Filter ..."); raise;
  end P_Filter;

  generic

    with function f ( k : integer32; x : Vector ) return Complex_Number;

  procedure G_Split_Filter
               ( file : in file_type; i,j : in integer32;
                 s : in out Array_of_Solution_Lists; p : in Matrix;
                 tol : in double_float; qws : out Solution_List );

  -- DESCRIPTION :
  --   Filters the j-th witness set with the i-th hypersurface.

  procedure G_Split_Filter
               ( file : in file_type; i,j : in integer32;
                 s : in out Array_of_Solution_Lists; p : in Matrix;
                 tol : in double_float; qws : out Solution_List ) is

  -- DESCRIPTION :
  --   Filters the j-th witness set with the i-th hypersurface.

    function Eval ( x : Vector ) return Complex_Number is
    begin
      return f(i,x);
    end Eval;
   -- procedure R_Split_Filter is new RG_SPlit_Filter(Eval);
    procedure S_Split_Filter is new SG_SPlit_Filter(Eval);

    len : constant natural32 := Length_Of(s(j));

  begin
    put(file,"Filtering witness set "); put(file,j,1);
    put(file," using current equation "); put(file,i,1);
    put_line(file," ...");
    S_Split_Filter(tol,s(j),qws,p);
   -- R_Split_Filter(file,tol,s(j),qws,p);
    put(file,"  ... filtered "); put(file,len,1); put(file," points : ");
    put(file,Length_Of(s(j)),1); put(file," remained and ");
    put(file,Length_Of(qws),1);  put_line(file," new.");
  end G_Split_Filter;

  procedure P_Split_Filter
               ( file : in file_type; i,j : in integer32; ep : in Eval_Poly;
                 s : in out Array_of_Solution_Lists; p : in Matrix;
                 tol : in double_float; qws : out Solution_List ) is

  -- DESCRIPTION :
  --   Filters the j-th witness set with the i-th hypersurface.

    len : constant natural32 := Length_Of(s(j));

  begin
    put(file,"Filtering witness set "); put(file,j,1);
    put(file," using current equation "); put(file,i,1);
    put_line(file," ...");
    SP_Split_Filter(ep,tol,s(j),qws,p);
   -- RP_Split_Filter(file,ep,tol,s(j),qws,p);
    put(file,"  ... filtered "); put(file,len,1); put(file," points : ");
    put(file,Length_Of(s(j)),1); put(file," remained and ");
    put(file,Length_Of(qws),1);  put_line(file," new.");
  end P_Split_Filter;

  procedure Singular_Filter
              ( file : in file_type;
                sols : in out Solution_List; tol : in double_float ) is

    len : constant natural32 := Length_Of(sols);
    sing,regu : Solution_List;

  begin
    Reporting_Singular_Filter(file,sols,tol,sing,regu);
   -- Silent_Singular_Filter(sols,tol,sing,regu);
    put(file,"Tested "); put(file,len,1); put(file," solutions, found ");
    put(file,Length_Of(regu),1); put(file," regular and ");
    put(file,Length_Of(sing),1); put_line(file," singular.");
    Clear(sols);
    sols := regu;
  end Singular_Filter;

  generic

    with function f ( k : integer32; x : Vector ) return Complex_Number;
    with function jf ( k : integer32; x : Vector ) return Vector;

  procedure G_Diagonal_Homotopy
              ( file : in file_type; n,i,j : in integer32;
                qws : in Solution_List;
                s : in out Array_of_Solution_Lists;
                p : in out VecMat; line : Matrix );

  -- DESCRIPTION :
  --   Executes the intrinsic diagonal homotopy to intersect the i-th
  --   witness set of a hypersurface with the j-th witness set.

  -- ON ENTRY :
  --   file     for diagnostics and intermediate results;
  --   n        dimension of the ambient space;
  --   i        number of the hypersurface;
  --   j        second component in the intersection, j<i;
  --   qws      witness points from j-th witness set;
  --   s        witness sets for all components;
  --   p        affine linear spaces defining all witness sets;
  --   line     random line defines witness set of hypersurface.

  -- ON RETURN :
  --   s        updated i-th component with witness set;
  --   p        new plane to define i-th witness set.

  procedure G_Diagonal_Homotopy
              ( file : in file_type; n,i,j : in integer32;
                qws : in Solution_List;
                s : in out Array_of_Solution_Lists;
                p : in out VecMat; line : Matrix ) is

    tol : constant double_float := 1.0E-10;
    sols : Solution_List;
    m : constant integer32 := j + 1;
    plane : Matrix(1..n,0..m);
    len : natural32;

    function fA ( x : Vector ) return Vector is

      res : Vector(1..1);

    begin
      res(1) := f(i,x);
      return res;
    end fA;

    function jfA ( x : Vector ) return Matrix is

      res : Matrix(1..1,x'range);
      val : constant Vector(x'range) := jf(i,x);

    begin
      for i1 in x'range loop
        res(1,i1) := val(i1);
      end loop;
      return res;
    end jfA;

    function fB ( x : Vector ) return Vector is

      res : Vector(1..j);

    begin
      for i1 in 1..j loop
        res(i1) := f(i1,x);
      end loop;
      return res;
    end fB;

    function jfB ( x : Vector ) return Matrix is

      res : Matrix(1..j,x'range);
      val : Vector(x'range);

    begin
      for i1 in 1..j loop
        val := jf(i1,x);
        for i2 in x'range loop
          res(i1,i2) := val(i2);
        end loop;
      end loop;
      return res;
    end jfB;
    
    procedure Diagonal_Homotopy is
      new Generic_Diagonal_Homotopy(fA,jfA,fB,jfB);

  begin
   -- put_line("intersecting in G_Diagonal_Homotopy ...");
    put(file,"Intersecting witness set "); put(file,i,1);
    put(file," with witness set "); put(file,j,1);
    put_line(file," ...");
    Diagonal_Homotopy(file,1,j,n,n-1,n-j,s(i),qws,line,p(j).all,sols,plane);
    len := Length_Of(sols);
    put(file,"  ... computed "); put(file,len,1);
    put_line(file," candidate intersection points.");
    Clear(s(i));
    Singular_Filter(file,sols,tol);
    put(file,"  computing singular filter : ");
    put(file,Length_Of(sols),1);
    put_line(file," regular solutions remained.");
    s(i) := sols;
    p(i) := new Matrix'(plane);
  end G_Diagonal_Homotopy;

  generic

    with function f ( k : integer32; x : Vector ) return Complex_Number;
    with function jf ( k : integer32; x : Vector ) return Vector;
    with function Q ( x : Vector ) return Vector;

  procedure GQ_Diagonal_Homotopy
              ( file : in file_type; n,i,j : in integer32;
                qws : in Solution_List;
                s : in out Array_of_Solution_Lists;
                p : in out VecMat; line : Matrix );

  -- DESCRIPTION :
  --   Executes the intrinsic diagonal homotopy to intersect the i-th
  --   witness set of a hypersurface with the j-th witness set.
  --   Filters solutions using the discriminating Q equations.

  -- ON ENTRY :
  --   file     for diagnostics and intermediate results;
  --   n        dimension of the ambient space;
  --   i        number of the hypersurface;
  --   j        second component in the intersection, j<i;
  --   qws      witness points from j-th witness set;
  --   s        witness sets for all components;
  --   p        affine linear spaces defining all witness sets;
  --   line     random line defines witness set of hypersurface.

  -- ON RETURN :
  --   s        updated i-th component with witness set;
  --   p        new plane to define i-th witness set.

  procedure GQ_Diagonal_Homotopy
              ( file : in file_type; n,i,j : in integer32;
                qws : in Solution_List;
                s : in out Array_of_Solution_Lists;
                p : in out VecMat; line : Matrix ) is

    tol : constant double_float := 1.0E-10;
    sols : Solution_List;
    m : constant integer32 := j + 1;
    plane : Matrix(1..n,0..m);
    len : natural32;

    function fA ( x : Vector ) return Vector is

      res : Vector(1..1);

    begin
      res(1) := f(i,x);
      return res;
    end fA;

    function jfA ( x : Vector ) return Matrix is

      res : Matrix(1..1,x'range);
      val : constant Vector(x'range) := jf(i,x);

    begin
      for i1 in x'range loop
        res(1,i1) := val(i1);
      end loop;
      return res;
    end jfA;

    function fB ( x : Vector ) return Vector is

      res : Vector(1..j);

    begin
      for i1 in 1..j loop
        res(i1) := f(i1,x);
      end loop;
      return res;
    end fB;

    function jfB ( x : Vector ) return Matrix is

      res : Matrix(1..j,x'range);
      val : Vector(x'range);

    begin
      for i1 in 1..j loop
        val := jf(i1,x);
        for i2 in x'range loop
          res(i1,i2) := val(i2);
        end loop;
      end loop;
      return res;
    end jfB;
    
    procedure Diagonal_Homotopy is
      new Generic_Diagonal_Homotopy(fA,jfA,fB,jfB);

   -- procedure Q_Filter is new QRG_Filter(Q);
    procedure Q_Filter is new QSG_Filter(Q);

  begin
   -- put_line("intersecting in GQ_Diagonal_Homotopy ...");
    put(file,"Intersecting witness set "); put(file,i,1);
    put(file," with witness set "); put(file,j,1);
    put_line(file," ...");
    Diagonal_Homotopy(file,1,j,n,n-1,n-j,s(i),qws,line,p(j).all,sols,plane);
    len := Length_Of(sols);
    put(file,"  ... computed "); put(file,len,1);
    put_line(file," candidate intersection points.");
    Clear(s(i));
    Q_Filter(sols,plane,tol); -- Q_Filter(file,sols,plane,tol);
    put(file,"  computing Q-filter : ");
    put(file,Length_Of(sols),1); put_line(file," remained.");
    Singular_Filter(file,sols,tol);
    put(file,"  computing singular filter : ");
    put(file,Length_Of(sols),1);
    put_line(file," regular solutions remained.");
    s(i) := sols;
    p(i) := new Matrix'(plane);
  end GQ_Diagonal_Homotopy;

  procedure P_Diagonal_Homotopy
              ( file : in file_type; report : in boolean;
                n,i,j : in integer32;
                f : in Poly_Sys; qws : in Solution_List;
                s : in out Array_of_Solution_Lists;
                p : in out VecMat; line : Matrix ) is

  -- DESCRIPTION :
  --   Executes the intrinsic diagonal homotopy to intersect the i-th
  --   witness set of a hypersurface with the j-th witness set.

  -- ON ENTRY :
  --   file     for diagnostics and intermediate results;
  --   report   indicates if path trackers report diagnostics,
  --            if false, then path trackers are silent;
  --   n        dimension of the ambient space;
  --   i        number of the hypersurface;
  --   j        second component in the intersection, j<i;
  --   f        polynomial system;
  --   qws      witness points from j-th witness set;
  --   s        witness sets for all components;
  --   p        affine linear spaces defining all witness sets;
  --   line     random line defines witness set of hypersurface.

  -- ON RETURN :
  --   s        updated i-th component with witness set;
  --   p        new plane to define i-th witness set.

    tol : double_float := 1.0E-8;
    sols : Solution_List;
    m : constant integer32 := j + 1;
    plane : Matrix(1..n,0..m);
    len : natural32;

  begin
   -- put_line("intersecting in P_Diagonal_Homotopy ...");
    put(file,"Intersecting witness set "); put(file,i,1);
    put(file," with witness set "); put(file,j,1);
    put_line(file," ...");
   -- put_line("calling Intrinsic_Diagonal_Homotopy ...");
    Intrinsic_Diagonal_Homotopy
      (file,report,natural32(n-1),natural32(n-j),
       f(i..i),f(1..j),s(i),qws,line,p(j).all,sols,plane);
   -- put_line("returning from Intrinsic_Diagonal_Homotopy");
    len := Length_Of(sols);
    put(file,"  ... computed "); put(file,len,1);
    put_line(file," candidate intersection points.");
    Clear(s(i));
    Singular_Filter(file,sols,tol);
    s(i) := sols;
    p(i) := new Matrix'(plane);
  end P_Diagonal_Homotopy;

  procedure P_Solve_Equation_by_Equation
              ( file : in file_type; filename : in string; 
                report,step : in boolean;
                ne,nv,k : in integer32; p : in Poly_Sys;
                witset : in out Array_of_Solution_Lists;
                planes : in out VecMat ) is

    timer : Timing_Widget;
    tol : constant double_float := 1.0E-8;
    ep : constant Eval_Poly_Sys(p'range) := Create(p);
    b : constant Vector(1..nv) := Random_Vector(1,nv);
    v : constant Vector(1..nv) := Random_Vector(1,nv);
    line : constant Matrix(1..nv,0..1) := Create_Plane(b,v);
    qws : Solution_List;
    superfluous : boolean;
    flow : Flow_Table := Create(ne,nv);
    start_equ : integer32;

  begin
   -- put_line("inside P_Solve_Equation_by_Equation ...");
    tstart(timer);
    if k = 0 then
      P_Hypersurface_Roots(file,nv,1,p(1),ep(1),b,v,witset(1));
      planes(1) := new Matrix'(line);
      Store_Hypersurface_Degree(flow,1,integer32(Length_Of(witset(1))));
      start_equ := 2;
    else
      start_equ := k+1;
    end if;
    for i in start_equ..ne loop
     -- put("inside main solver loop with i = "); put(i,1); put_line(" ...");
      P_Hypersurface_Roots(file,nv,i,p(i),ep(i),b,v,witset(i));
      Store_Hypersurface_Degree(flow,i,integer32(Length_Of(witset(i))));
      if i = 2 then Update(flow,1,witset); end if;
     -- put_line("calling P_filter ...");
      P_Filter(file,i,witset(i),b,v,ep(1..i-1),tol);
     -- put_line(" ... returned from P_filter");
      superfluous := true;
      for j in 1..i-1 loop
        if not Is_Null(witset(j)) then
          P_Split_Filter(file,i,j,ep(i),witset,planes(j).all,tol,qws);
          if not Is_Null(qws) then
            superfluous := false;
            P_Diagonal_Homotopy(file,report,nv,i,j,p,qws,witset,planes,line);
            if step and not Is_Null(witset(i)) then
              Write_Witness_Stone
                (file,filename,natural32(nv),natural32(nv-i),
                 p(1..i),witset(i),planes(i).all);
            end if;
            Clear(qws);
          end if;
        end if;
      end loop;
      if superfluous then
        put(file,"Equation "); put(file,i,1);
        put_line(file," is superfluous.");
        Clear(witset(i));
      end if;
      Update(flow,i,witset);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Equation-by-Equation Solver");
    new_line(file);
    put_line(file,"The flow table :"); Write(file,flow);
    Clear(flow);
    Write_Witness_Sets(file,filename,natural32(nv),p,witset,planes);
  end P_Solve_Equation_by_Equation;

  procedure G_Solve_Equation_by_Equation
              ( file : in file_type; filename : in string;
                report,step : in boolean;
                ne,nv,k : in integer32; d : in Standard_Natural_Vectors.Vector;
                witset : out Array_of_Solution_Lists; planes : out VecMat ) is

    timer : Timing_Widget;
    tol : constant double_float := 1.0E-8;
    b : constant Vector(1..nv) := Random_Vector(1,nv);
    v : constant Vector(1..nv) := Random_Vector(1,nv);
    line : constant Matrix(1..nv,0..1) := Create_Plane(b,v);
    qws : Solution_List;
    superfluous : boolean;
    flow : Flow_Table := Create(ne,nv);
    start_equ : integer32;

    procedure Roots is new G_Hypersurface_Roots(f);
    procedure Filter is new G_Filter(f);
    procedure Split_Filter is new G_Split_Filter(f);
    procedure Diagonal_Homotopy is new G_Diagonal_Homotopy(f,jf);

  begin
    tstart(timer);
    if k = 0 then
      Roots(file,nv,1,integer32(d(1)),b,v,witset(1));
      planes(1) := new Matrix'(line);
      Store_Hypersurface_Degree(flow,1,integer32(Length_Of(witset(1))));
      start_equ := 2;
    else
      start_equ := k+1;
    end if;
    for i in start_equ..ne loop
      Roots(file,nv,i,integer32(d(i)),b,v,witset(i));
      Store_Hypersurface_Degree(flow,i,integer32(Length_Of(witset(i))));
      if i = 2 then Update(flow,1,witset); end if;
      Filter(file,i,witset(i),b,v,tol);
      superfluous := true;
      for j in 1..i-1 loop
        if not Is_Null(witset(j)) then
          Split_Filter(file,i,j,witset,planes(j).all,tol,qws);
          if not Is_Null(qws) then
            superfluous := false;
            Diagonal_Homotopy(file,nv,i,j,qws,witset,planes,line);
            if step and not Is_Null(witset(i)) then
              Write_Witness_Stone
                (file,filename,natural32(nv),natural32(nv-i),
                 witset(i),planes(i).all);
            end if;
            Clear(qws);
          end if;
        end if;
      end loop;
      if superfluous then
        put(file,"Equation "); put(file,i,1);
        put_line(file," is superfluous.");
        Clear(witset(i));
      end if;
      Update(flow,i,witset);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Equation-by-Equation Solver");
    new_line(file);
    put_line(file,"The flow table :"); Write(file,flow);
    Clear(flow);
    Write_Witness_Sets(file,filename,natural32(nv),witset,planes);
  end G_Solve_Equation_by_Equation;

  procedure GQ_Solve_Equation_by_Equation
              ( file : in file_type; filename : in string;
                report,step : in boolean;
                ne,nv,k : in integer32; d : in Standard_Natural_Vectors.Vector;
                witset : out Array_of_Solution_Lists; planes : out VecMat ) is

    timer : Timing_Widget;
    tol : constant double_float := 1.0E-8;
    b : constant Vector(1..nv) := Random_Vector(1,nv);
    v : constant Vector(1..nv) := Random_Vector(1,nv);
    line : constant Matrix(1..nv,0..1) := Create_Plane(b,v);
    qws : Solution_List;
    superfluous : boolean;
    flow : Flow_Table := Create(ne,nv);
    start_equ : integer32;

    procedure Roots is new G_Hypersurface_Roots(f);
    procedure Filter is new G_Filter(f);
    procedure Split_Filter is new G_Split_Filter(f);
    procedure Diagonal_Homotopy is new GQ_Diagonal_Homotopy(f,jf,Q);

  begin
    tstart(timer);
    if k = 0 then
      Roots(file,nv,1,integer32(d(1)),b,v,witset(1));
      planes(1) := new Matrix'(line);
      Store_Hypersurface_Degree(flow,1,integer32(Length_Of(witset(1))));
      start_equ := 2;
    else
      start_equ := k+1;
    end if;
    for i in start_equ..ne loop
      Roots(file,nv,i,integer32(d(i)),b,v,witset(i));
      Store_Hypersurface_Degree(flow,i,integer32(Length_Of(witset(i))));
      if i = 2 then Update(flow,1,witset); end if;
      Filter(file,i,witset(i),b,v,tol);
      superfluous := true;
      for j in 1..i-1 loop
        if not Is_Null(witset(j)) then
          Split_Filter(file,i,j,witset,planes(j).all,tol,qws);
          if not Is_Null(qws) then
            superfluous := false;
            Diagonal_Homotopy(file,nv,i,j,qws,witset,planes,line);
            if step and not Is_Null(witset(i)) then
              Write_Witness_Stone
                (file,filename,natural32(nv),natural32(nv-i),
                 witset(i),planes(i).all);
            end if;
            Clear(qws);
          end if;
        end if;
      end loop;
      if superfluous then
        put(file,"Equation "); put(file,i,1);
        put_line(file," is superfluous.");
        Clear(witset(i));
      end if;
      Update(flow,i,witset);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Equation-by-Equation Solver");
    new_line(file);
    put_line(file,"The flow table :"); Write(file,flow);
    Clear(flow);
    Write_Witness_Sets(file,filename,natural32(nv),witset,planes);
  end GQ_Solve_Equation_by_Equation;

end Equation_by_Equation_Solvers;
