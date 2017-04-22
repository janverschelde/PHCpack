with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Natural_Matrices;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with Standard_Complex_VecMats;          use Standard_Complex_VecMats;
with Standard_Complex_VecMats_io;       use Standard_Complex_VecMats_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Root_Refiners;            use Standard_Root_Refiners;
with Brackets;                          use Brackets;
with Matrix_Indeterminates;
with Determinantal_Systems;
with Localization_Posets;               use Localization_Posets;
with Localization_Poset_Strings;        use Localization_Poset_Strings;
with Deformation_Posets;                use Deformation_Posets;
with Symbolic_Minor_Equations;
with Curves_into_Grassmannian;
with Curves_into_Grassmannian_io;
with Drivers_for_Input_Planes;          use Drivers_for_Input_Planes;
with Drivers_for_Pieri_Homotopies;      use Drivers_for_Pieri_Homotopies; 
with Pieri_Root_Count;
with Pieri_Homotopy;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_PolySys_Container;
with Standard_Solutions_Container;

function use_c2pieri ( job : integer32;
                       a : C_intarrs.Pointer;
		       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer ) return integer32 is

  procedure Write_Menu is
  begin
    new_line;
    put_line("General MENU to use Pieri homotopies from C :");
    put_line("  0. display this menu;");
    put_line("  1. initialize dimensions (m,p,q);");
    put_line("  2. initialize m*p + q*(m+p) input m-planes;");
    put_line("  3. initialize m*p + q*(m+p) interpolation points;");
    put_line("  4. store pivots of pattern at start solution curve;");
    put_line("  5. store pivots of pattern at target solution curve;");
    put_line("  6. store coefficients of start solution curve;");
    put_line("  7. retrieve coefficients of target solution curve;");
    put_line("  8. track solution path without intermediate output;");
    put_line("  9. track solution path with output diagnostics;");
    put_line(" 10. verify intersection conditions without output;");
    put_line(" 11. verify intersection conditions with extra output;");
    put_line(" 12. destroy the state machine;");
    put_line(" 13. compute the combinatorial Pieri root count;");
    put_line(" 14. return localization poset as string in b;");
    put_line(" 15. run the Pieri homotopies on random input data;");
    put_line(" 16. generate real planes osculating a rational curve;");
    put_line(" 17. put Schubert polynomial system in container.");
  end Write_Menu;

  procedure Get_Dimensions
              ( a : C_intarrs.Pointer; m,p,q : out integer32 ) is

  -- DESCRIPTION :
  --   Extracts the dimension (m,p,q) from the array a.

    v : constant C_Integer_Array(0..2)
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));

  begin
    m := integer32(v(0));
    p := integer32(v(1));
    q := integer32(v(2));
   -- put("Ada: m = "); put(m,1);
   -- put("  p = "); put(p,1);
   -- put("  q = "); put(q,1); new_line;
  end Get_Dimensions;

  procedure Get_All_Dimensions
              ( a,b : C_intarrs.Pointer; m,p,q,n : out integer32 ) is

  -- DESCRIPTION : 
  --   Calls Get_Dimensions on the array a and extracts n from b.

    v : constant C_Integer_Array(0..0) := C_intarrs.Value(b);

  begin
    Get_Dimensions(a,m,p,q);
    n := integer32(v(0));
   -- put("Ada: n = "); put(n,1); new_line;
  end Get_All_Dimensions;

  procedure Get_and_Initialize_Input_Planes
              ( c : C_dblarrs.Pointer; m,p,q,n : in integer32 ) is

  -- DESCRIPTION :
  --   Extracts the coefficients of the input planes from c
  --   and initializes the Pieri homotopy state machine with it.

    nbcff : constant integer32 := 2*(m+p)*m*n;
    d : constant Interfaces.C.size_t := Interfaces.C.size_t(nbcff);
    v : C_Double_Array(0..d-1)
      := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(d));
    ind : Interfaces.C.size_t := 0;
    planes : VecMat(1..n);

  begin
   -- put_line("The coefficients of the input planes : ");
   -- for i in v'range loop
   --   put(double_float(v(i))); new_line;
   -- end loop;
    for k in planes'range loop
      declare
        plane : Standard_Complex_Matrices.Matrix(1..m+p,1..m);
      begin
        for i in 1..m+p loop
          for j in 1..m loop
            plane(i,j) := Create(double_float(v(ind)),
                                 double_float(v(ind+1)));
            ind := ind + 2;
          end loop;
        end loop;
        planes(k) := new Standard_Complex_Matrices.Matrix'(plane);
       -- put("input plane "); put(k,1); put_line(" :"); put(plane);
      end;
    end loop;
    Pieri_Homotopy.Initialize_Input_Planes(planes);
  end Get_and_Initialize_Input_Planes;

  procedure Get_and_Initialize_Interpolation_Points
              ( c : C_dblarrs.Pointer; m,p,q,n : in integer32 ) is

  -- DESCRIPTION :
  --   Extracts the coefficients of the interpolation points from c
  --   and initializes the Pieri homotopy state machine with it.

    nbcff : constant integer32 := 2*n;
    d : constant Interfaces.C.size_t := Interfaces.C.size_t(nbcff);
    v : constant C_Double_Array(0..d-1)
      := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(d));
    ind : Interfaces.C.size_t := 0;
    points : Standard_Complex_Vectors.Vector(1..n);
    
  begin
   -- put_line("The coefficients of the interpolation points : ");
   -- for i in v'range loop
   --   put(double_float(v(i))); new_line;
   -- end loop;
    for k in points'range loop
      points(k) := Create(double_float(v(ind)),double_float(v(ind+1)));
      ind := ind + 2;
    end loop;
   -- put_line("The interpolation points : "); put_line(points);
    Pieri_Homotopy.Initialize_Interpolation_Points(points);
  end Get_and_Initialize_Interpolation_Points;

  function Get_Complex_Vector
             ( c : C_dblarrs.Pointer; n : in integer32 )
             return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Extracts from c the coefficients of a complex n-vector.

    res : Standard_Complex_Vectors.Vector(1..n);
    nbcff : constant integer32 := 2*n;
    d : constant Interfaces.C.size_t := Interfaces.C.size_t(nbcff);
    v : C_Double_Array(0..d-1)
      := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(d));
    ind : Interfaces.C.size_t := 0;
    
  begin
   -- put_line("The coefficients of the start plane : ");
   -- for i in v'range loop
   --   put(double_float(v(i))); new_line;
   -- end loop;
    for k in res'range loop
      res(k) := Create(double_float(v(ind)),double_float(v(ind+1)));
      ind := ind + 2;
    end loop;
    return res;
  end Get_Complex_Vector;

  function Get_Number ( a : C_intarrs.Pointer ) return integer32 is

  -- DESCRIPTION :
  --   Extracts a number from the array a and returns its value.

    v : constant C_Integer_Array(0..0) := C_intarrs.Value(a);
    p : constant integer32 := integer32(v(0));

  begin
    return p;
  end Get_Number;

  procedure Get_Pivots
              ( a : in C_intarrs.Pointer; p : in integer32;
                top,bottom : out Standard_Natural_Vectors.Vector ) is

    d : constant Interfaces.C.size_t := Interfaces.C.size_t(2*p);
    vp : constant C_Integer_Array(0..d-1)
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(d));
    ind : Interfaces.C.size_t := 0;

  begin
    for i in 1..p loop
      top(i) := natural32(vp(ind));
      ind := ind + 1;
    end loop;
    for i in 1..p loop
      bottom(i) := natural32(vp(ind));
      ind := ind + 1;
    end loop;
  end Get_Pivots;

  procedure Write_Pivots
              ( top,bottom : in Standard_Natural_Vectors.Vector ) is
  begin
    put("Ada: [");
    put(top);
    put(" ],[");
    put(bottom);
    put_line(" ]");
  end Write_Pivots;

  function Job1 return integer32 is -- initialize dimensions

    m,p,q : integer32;

  begin
    Get_Dimensions(a,m,p,q);
    Pieri_Homotopy.Initialize_Dimensions
      (natural32(m),natural32(p),natural32(q));
    return 0;
  end Job1;

  function Job2 return integer32 is -- initialize input planes

    m,p,q,n : integer32;

  begin
    Get_All_Dimensions(a,b,m,p,q,n);
    Get_and_Initialize_Input_Planes(c,m,p,q,n);
    return 0;
  end Job2;

  function Job3 return integer32 is -- initialize interpolation points

    m,p,q,n : integer32;

  begin
    Get_All_Dimensions(a,b,m,p,q,n);
    Get_and_Initialize_Interpolation_Points(c,m,p,q,n);
    return 0;
  end Job3;

  function Job4 return integer32 is -- store pivots at start

    p : constant integer32 := Get_Number(a);
    top,bottom : Standard_Natural_Vectors.Vector(1..p);

  begin
    Get_Pivots(b,p,top,bottom);
   -- Write_Pivots(top,bottom);
    Pieri_Homotopy.Store_Start_Pivots(top,bottom);
    return 0;
  end Job4;

  function Job5 return integer32 is -- store pivots at target

    p : constant integer32 := Get_Number(a);
    top,bottom : Standard_Natural_Vectors.Vector(1..p);

  begin
    Get_Pivots(b,p,top,bottom);
   -- Write_Pivots(top,bottom);
    Pieri_Homotopy.Store_Target_Pivots(top,bottom);
    return 0;
  end Job5;

  function Job6 return integer32 is -- store coefficients at start

    n : constant integer32 := Get_Number(a);            
    x : Standard_Complex_Vectors.Vector(1..n);

  begin
   -- put("The number of variables : "); put(n,1); new_line;
    if n > 0
     then x := Get_Complex_Vector(c,n);
    end if;
    Pieri_Homotopy.Store_Start(x);
    return 0;
  end Job6;

  function Job7 return integer32 is -- retrieve target solution

    n : constant integer32 := Get_Number(a);
    x : Standard_Complex_Vectors.Vector(1..n);

  begin
    Pieri_Homotopy.Retrieve_Target(x);
    Assign(x,c);
    return 0;
  end Job7;

  function Job13 return integer32 is -- compute Pieri root count

    m,p,q : integer32;
    r : natural32;

  begin
    Get_Dimensions(a,m,p,q);
    r := Pieri_Root_Count(natural32(m),natural32(p),natural32(q));
    Assign(integer32(r),b);
    return 0;
  end Job13;

  function Job14 return integer32 is -- localization poset string

    m,p,q : integer32;

  begin
    Get_Dimensions(a,m,p,q);
    declare
      root : constant Node(p)
           := Trivial_Root(natural32(m),natural32(p),natural32(q));
      lnkroot : Link_to_Node := new Node'(root);
      nq : constant integer32 := m*p + q*(m+p);
      level_poset : Array_of_Nodes(0..nq);
    begin
      Q_Top_Bottom_Create(lnkroot,root.bottom(p),natural32(m+p));
      level_poset := Create_Leveled_Poset(lnkroot);
      Count_Roots(level_poset);
      declare
        s : constant string := Poset_to_String(level_poset);
        sv : constant Standard_Integer_Vectors.Vector
           := String_to_Integer_Vector(s);
      begin
        Assign(integer32(s'last),a);
        Assign(sv,b);
      end;
    end;
    return 0;
  end Job14;

  procedure Show_Input_Planes_and_Points
              ( n,m,p,q : in integer32;
                c : in Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Writes the coefficients c of the n input m-planes to screen
  --   for testing purposes.
  --   If q > 0, then the interpolation points are shown as well.

    ind : integer32 := 0;

  begin
    for i in 1..n loop
      put("input plane "); put(i,1); put_line(" :");
      for j in 1..(m+p) loop
        for k in 1..m loop
          put("(");
          ind := ind + 1; put(c(ind)); put(" ");
          ind := ind + 1; put(c(ind)); put(")");
         -- if ind mod 2 = 0
         --  then new_line;
         -- end if;
        end loop;
        new_line;
      end loop;
    end loop;
    if q > 0 then
      put_line("the interpolation points :");
      for i in 1..n loop
        ind := ind + 1; put(c(ind)); put(" ");
        ind := ind + 1; put(c(ind)); new_line;
      end loop;
    end if;
  end Show_Input_Planes_and_Points;

  function Construct_Input_Planes
             ( n,m,p : in integer32;
               c : in Standard_Floating_Vectors.Vector ) return VecMat is

  -- DESCRIPTION :
  --   Returns the n input m-planes in (m+p)-space as
  --   a proper vector of matrices.

    res : VecMat(1..n);
    dim : constant integer32 := m+p;
    ind : integer32 := 0;

  begin
    for i in 1..n loop
      declare
        A : Standard_Complex_Matrices.Matrix(1..dim,1..m);
      begin
        for j in 1..(m+p) loop
          for k in 1..m loop
            A(j,k) := Create(c(ind+1),c(ind+2));
            ind := ind + 2;
          end loop;
        end loop;
        res(i) := new Standard_Complex_Matrices.Matrix'(A);
      end;
    end loop;
    return res;
  end Construct_Input_Planes;

  procedure Refine_Roots ( p : in Poly_Sys; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Applies Newton's method to all solution of p
  --   to enscure there is data for (err,rco,res).

    epsxa : constant double_float := 1.0E-8;
    epsfa : constant double_float := 1.0E-8;
    tolsing : constant double_float := 1.0E-8;
    max : constant natural32 := 3;
    numit : natural32 := 0;
    deflate : boolean := false;

  begin
    Silent_Root_Refiner(p,sols,epsxa,epsfa,tolsing,numit,max,deflate);
  end Refine_Roots;

  procedure Extract_System_and_Solutions
              ( m,p : in integer32; planes : in VecMat;
                index_poset : in Array_of_Array_of_Nodes;
                deform_poset : in Array_of_Array_of_VecMats ) is

  -- DESCRIPTION :
  --   Extracts the solutions and the polynomial system from the input.
  --   The solution container will contain the solutions
  --   and the systems container the system.

    dim : constant integer32 := m*p;
    target_level : constant integer32 := dim;
    xpm : constant Standard_Complex_Poly_Matrices.Matrix(1..m+p,1..p)
        := Symbolic_Minor_Equations.Localization_Pattern
             (natural32(m+p),index_poset(dim)(1).top,
                             index_poset(dim)(1).bottom);
    solplanes : constant VecMat := deform_poset(target_level)(1).all;
    locmap : constant Standard_Natural_Matrices.Matrix(1..m+p,1..p)
           := Determinantal_Systems.Standard_Coordinate_Frame
                (xpm,solplanes(1).all);
    locsys : constant Poly_Sys(planes'range)
           := Create_Hypersurface_System(natural32(dim),locmap,xpm,planes);
    sols : Solution_List := Solution_Planes(locmap,solplanes);

  begin
    Matrix_Indeterminates.Reduce_Symbols(locmap);
    Refine_Roots(locsys,sols);
   -- put_line("THE GENERIC SYSTEM : ");
   -- put_line(locsys);
    Standard_PolySys_Container.Initialize(locsys);
   -- put_line("THE SOLUTIONS :");
   -- put(standard_output,Length_Of(sols),Head_Of(sols).n,sols);
    Standard_Solutions_Container.Initialize(sols);
  end Extract_System_and_Solutions;

  procedure Solve_Hypersurface_Schubert_Problem
              ( n,m,p,tdim : in integer32; rc : out natural32 ) is

  -- DESCRIPTION :
  --   Extracts the input data from c and solves the
  --   corresponding hypersurface Schubert problem.
                
    cff : Standard_Floating_Vectors.Vector(1..tdim);
    planes : VecMat(1..n);
    root : constant Node(p) := Trivial_Root(natural32(m),natural32(p));
    lnkroot : Link_to_Node := new Node'(root);
    level_poset : Array_of_Nodes(0..m*p);
    index_poset : Array_of_Array_of_Nodes(0..m*p);
    deform_poset : Array_of_Array_of_VecMats(index_poset'range);
    target_level : constant integer32 := m*p;
    npaths : Standard_Natural_Vectors.Vector(1..target_level)
           := (1..target_level => 0);
    timings : Duration_Array(1..target_level) := (1..target_level => 0.0);

  begin
    Assign(natural32(tdim),c,cff);
   -- Show_Input_Planes_and_Points(n,m,p,0,cff);
    planes := Construct_Input_Planes(n,m,p,cff);
   -- put_line("The input planes :"); put(planes);
    Top_Bottom_Create(lnkroot,natural32(m+p));
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    rc := natural32(level_poset(m*p).roco);
   -- put("number of roots : "); put(rc,1); new_line;
    deform_poset := Create(index_poset);
    Matrix_Indeterminates.Initialize_Symbols(natural32(m+p),natural32(p));
    for i in index_poset(target_level)'range loop
      declare
        root : constant Node := index_poset(target_level)(i).all;
      begin
       -- Solve(standard_output,m+p,deform_poset,root,planes(1..target_level),
       --       true,true,npaths,timings);
        Solve(natural32(m+p),deform_poset,root,
              planes(1..target_level),npaths,timings);
        Extract_System_and_Solutions(m,p,planes,index_poset,deform_poset);
      end;
    end loop;
  end Solve_Hypersurface_Schubert_Problem;

  procedure Extract_System_and_Solutions
              ( n,m,p,q : in integer32; planes : in VecMat;
                points : in Standard_Complex_Vectors.Vector;
                index_poset : in Array_of_Array_of_Nodes;
                deform_poset : in Array_of_Array_of_VecMats ) is

  -- DESCRIPTION :
  --   Extracts the solutions and the polynomial system from the parameters.

    dim : constant integer32 := m*p+q*(m+p);
    top : constant Bracket := index_poset(dim)(1).top;
    bot : constant Bracket := index_poset(dim)(1).bottom; 
    xpm : constant Standard_Complex_Poly_Matrices.Matrix(1..m+p,1..p)
        := Curves_into_Grassmannian.Symbolic_Create
             (natural32(m),natural32(p),natural32(q),top,bot);
    solplanes : constant VecMat := deform_poset(dim)(1).all;
    locmap : constant Standard_Natural_Matrices.Matrix(1..(m+p)*(q+1),1..p)
           := Curves_into_Grassmannian.Standard_Coordinate_Frame
                (natural32(m),natural32(p),natural32(q),
                 top,bot,solplanes(1).all);
    sols : Solution_List
         := Determinantal_Systems.Solution_Planes
              (top,bot,locmap,solplanes);
    locsys : constant Poly_Sys(planes'range)
           := Determinantal_Systems.Create_Polynomial_System
                (top,bot,locmap,xpm,points,planes);

  begin
    Curves_into_Grassmannian_io.Two_Set_up_Symbol_Table
      (natural32(m),natural32(p),natural32(q),top,bot);
   -- put_line("THE GENERIC SYSTEM :");
   -- put_line(locsys);
    Standard_PolySys_Container.Initialize(locsys);
   -- put_line("THE SOLUTIONS :");
   -- put(standard_output,Length_Of(sols),Head_Of(sols).n,sols);
    Standard_Solutions_Container.Initialize(sols);
  end Extract_System_and_Solutions;

  procedure Solve_Quantum_Schubert_Problem
              ( n,m,p,q,tdim : in integer32; rc : out natural32 ) is

  -- DESCRIPTION :
  --   Extracts the input data from c and solves the
  --   corresponding quantum Schubert problem.

    cff : Standard_Floating_Vectors.Vector(1..tdim+2*n);
    planes : VecMat(1..n);
    pts : Standard_Complex_Vectors.Vector(1..n);
    ind : integer32 := tdim;
    level_poset : Array_of_Nodes(0..n);
    index_poset : Array_of_Array_of_Nodes(0..n);
    deform_poset : Array_of_Array_of_VecMats(index_poset'range);
    the_root : constant Node(p)
             := Trivial_Root(natural32(m),natural32(p),natural32(q));
    lnkroot : Link_to_Node := new Node'(the_root);
    npaths : Standard_Natural_Vectors.Vector(1..n) := (1..n => 0);
    timings : Duration_Array(1..n) := (1..n => 0.0);

  begin
    Assign(natural32(tdim+2*n),c,cff);
   -- Show_Input_Planes_and_Points(n,m,p,q,cff);
    planes := Construct_Input_Planes(n,m,p,cff);
   -- put_line("The input planes :"); put(planes);
    for i in 1..n loop
      pts(i) := Create(cff(ind+1),cff(ind+2));
      ind := ind + 2;
    end loop;
   -- put_line("the interpolation points :"); put_line(pts);
    Q_Top_Bottom_Create(lnkroot,the_root.bottom(p),natural32(m+p));
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    rc := natural32(level_poset(n).roco);
   -- put("number of roots : "); put(rc,1); new_line;
    deform_poset := Create(index_poset);
    declare
      root : constant Node := index_poset(n)(1).all;
    begin
     -- Solve(standard_output,m+p,q,rc,deform_poset,root,planes,
     --       pts,true,true,npaths,timings);
      Solve(natural32(m+p),natural32(q),rc,deform_poset,root,
            planes,pts,npaths,timings);
      Extract_System_and_Solutions
        (n,m,p,q,planes,pts,index_poset,deform_poset);
    end;
  end Solve_Quantum_Schubert_Problem;

  function Job15 return integer32 is -- run Pieri homotopies

    m,p,q,n,mdim,tdim : integer32;
    rc : natural32;

  begin
   -- put_line("running the Pieri homotopies ...");
    Get_Dimensions(a,m,p,q);
    n := m*p + q*(m+p);
    mdim := m*(m+p);
    tdim := 2*n*mdim;
    if q = 0
     then Solve_Hypersurface_Schubert_Problem(n,m,p,tdim,rc);
     else Solve_Quantum_Schubert_Problem(n,m,p,q,tdim,rc);
    end if;
    Assign(integer32(rc),b);
    return 0;
  end Job15;

  function Flatten ( n : integer32; v : VecMat )
                   return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns one vector of range 1..n which contains all coefficients
  --   of the (real) planes in v.

    res : Standard_Floating_Vectors.Vector(1..n);
    mat : Standard_Complex_Matrices.Link_to_Matrix;
    ind : integer32 := 0;

  begin
    for k in v'range loop
      mat := v(k);
      for i in mat'range(1) loop
        for j in mat'range(2) loop
          ind := ind + 1;
          res(ind) := REAL_PART(mat(i,j));
        end loop;
      end loop;
    end loop;
    return res;
  end Flatten;

  function Job16 return integer32 is -- generate real osculating planes

  -- DESCRIPTION :
  --   Generates n real m-planes in d-space osculating a rational curve
  --   at n given points.

    n,m,p,q,d : integer32;

  begin
    Get_Dimensions(a,m,p,q);
    d := m+p;
    n := m*p + q*d;
   -- put("c2pieri expects "); put(n,1); put(" points ...");
    declare
      s : Standard_Floating_Vectors.Vector(1..n);
      n1 : Interfaces.C.size_t := Interfaces.C.size_t(n-1);
      i1 : Interfaces.C.size_t;
      v : C_Double_Array(0..n1)
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(n));
      planes : VecMat(1..n);
      coeffs : Standard_Floating_Vectors.Vector(1..n*m*d);
    begin
      for i in 1..n loop
        i1 := Interfaces.C.size_t(i-1);
        s(i) := double_float(v(i1));
      end loop;
     -- put_line("the interpolation points :"); put_line(s);
      planes := Osculating_Input_Planes
                  (natural32(m),natural32(p),natural32(q),s);
      coeffs := Flatten(coeffs'last,planes);
     -- put_line("the coefficients of the planes :"); put_line(coeffs);
      Assign(coeffs,c);
    end;
    return 0;
  end Job16;

  function Standard_Localization_Map
             ( m,p : integer32; top,bottom : Bracket )
             return Standard_Natural_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the standard localization map (with ones defined
  --   by the top pivots) for an m-plane in (m+p)-space.

    res : Standard_Natural_Matrices.Matrix(1..m+p,1..p);

  begin
    for j in res'range(2) loop
      for i in res'range(1) loop
        if i < integer32(top(j)) then
          res(i,j) := 0;
        elsif i = integer32(top(j)) then
          res(i,j) := 1;
        elsif i <= integer32(bottom(j)) then
          res(i,j) := 2;
        else
          res(i,j) := 0;
        end if;
      end loop;
    end loop;
    return res;
  end Standard_Localization_Map;

  procedure Hypersurface_Target_System ( m,p : in integer32 ) is

  -- DESCRIPTION :
  --   Makes the target system for the hypersurface Pieri case.

    d : constant integer32 := m+p;         -- ambient dimension
    n : constant integer32 := m*p;         -- number of input planes
    mdim : constant integer32 := m*d;      -- number of coefficients in m-plane
    tdim : constant integer32 := 2*n*mdim; -- number of real numbers
    cff : Standard_Floating_Vectors.Vector(1..tdim);
    planes : VecMat(1..n);
    root : constant Node(p) := Trivial_Root(natural32(m),natural32(p));
    xpm : constant Standard_Complex_Poly_Matrices.Matrix(1..d,1..p)
        := Symbolic_Minor_Equations.Localization_Pattern
             (natural32(d),root.top,root.bottom);
    locmap : Standard_Natural_Matrices.Matrix(1..m+p,1..p)
           := Standard_Localization_Map(m,p,root.top,root.bottom);
    locsys : Poly_Sys(planes'range);
   
  begin
    Assign(natural32(tdim),c,cff);
    planes := Construct_Input_Planes(n,m,p,cff);
   -- put_line("The input planes :"); put(planes);
    Matrix_Indeterminates.Initialize_Symbols(natural32(m+p),natural32(p));
    Matrix_Indeterminates.Reduce_Symbols(locmap);
    locsys := Create_Hypersurface_System(natural32(n),locmap,xpm,planes);
    Standard_PolySys_Container.Initialize(locsys);
  end Hypersurface_Target_System;

  function Job17 return integer32 is -- make Schubert polynomial system

  -- DESCRIPTION :
  --   Gets the dimensions and coefficients of the m-planes
  --   and puts the target system in the container.

    m,p,q : integer32;

  begin
    Get_Dimensions(a,m,p,q);
    if q = 0
     then Hypersurface_Target_System(m,p);
    end if;
    return 0;
  end Job17;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => Write_Menu;
      when 1 => return Job1; -- initialize dimensions
      when 2 => return Job2; -- initialize input planes
      when 3 => return Job3; -- initialize interpolation points
      when 4 => return Job4; -- store pivots at start
      when 5 => return Job5; -- store pivots at target
      when 6 => return Job6; -- store coefficients of start solution
      when 7 => return Job7; -- retrieve target solution
      when 8 => Pieri_Homotopy.Track_Path;
      when 9 => Pieri_Homotopy.Track_Path(Standard_Output);
      when 10 =>
        declare
          res : double_float := Pieri_Homotopy.Verify_Determinants;
        begin
          Assign(res,c);
        end;
      when 11 =>
        declare
          res : double_float
              := Pieri_Homotopy.Verify_Determinants(Standard_Output);
        begin
          Assign(res,c);
        end;
      when 12 => Pieri_Homotopy.Clear;
      when 13 => return Job13; -- compute Pieri root count
      when 14 => return Job14; -- return localization poset string
      when 15 => return Job15; -- run Pieri homotopies
      when 16 => return Job16; -- generate real osculating planes
      when 17 => return Job17; -- make Schubert polynomial system
      when others => put_line("  Sorry.  Invalid operation.");
                     return 1;
    end case;
    return 0;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_c2pieri;
