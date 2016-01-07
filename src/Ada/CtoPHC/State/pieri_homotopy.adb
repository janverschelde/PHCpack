with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Poly_Matrices_io;  use Standard_Complex_Poly_Matrices_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Homotopy;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Continuation_Data;         use Standard_Continuation_Data;
with Process_io;                         use Process_io;
with Continuation_Parameters;
with Standard_Path_Trackers;             use Standard_Path_Trackers;
with Brackets;                           use Brackets;
with Bracket_Monomials;                  use Bracket_Monomials;
with Standard_Bracket_Polynomials;       use Standard_Bracket_Polynomials;
with Standard_Bracket_Polynomials_io;    use Standard_Bracket_Polynomials_io;
with Standard_Bracket_Systems;           use Standard_Bracket_Systems;
with Specialization_of_Planes;           use Specialization_of_Planes;
with Symbolic_Minor_Equations;           use Symbolic_Minor_Equations;
with Numeric_Minor_Equations;            use Numeric_Minor_Equations;
with Curves_into_Grassmannian;           use Curves_into_Grassmannian;
with Curves_into_Grassmannian_io;        use Curves_into_Grassmannian_io;
with Verification_with_Determinants;     use Verification_with_Determinants;
with Pieri_Homotopies;                   use Pieri_homotopies;

package body Pieri_Homotopy is

-- INTERNAL DATA :

  m_dim,p_dim,q_deg : integer32;
  input_planes : Link_to_VecMat;
  interpolation_points : Standard_Complex_Vectors.Link_to_Vector;
  start_top,start_bottom : Standard_Natural_Vectors.Link_to_Vector;
  target_top,target_bottom : Standard_Natural_Vectors.Link_to_Vector;
  start,target : Standard_Complex_Vectors.Link_to_Vector;

-- TARGET OPERATIONS :

  procedure Initialize_Dimensions ( m,p,q : in natural32 ) is
  begin
    m_dim := integer32(m);
    p_dim := integer32(p);
    q_deg := integer32(q);
  end Initialize_Dimensions;

  procedure Initialize_Input_Planes ( planes : in VecMat ) is
  begin
    input_planes := new VecMat'(planes);
  end Initialize_Input_Planes;

  procedure Initialize_Interpolation_Points
              ( points : in Standard_Complex_Vectors.Vector ) is
  begin
    interpolation_points := new Standard_Complex_Vectors.Vector'(points);
  end Initialize_Interpolation_Points;

  procedure Initialize ( m,p,q : in natural32; planes : in VecMat;
                         points : in Standard_Complex_Vectors.Vector ) is
  begin
    Initialize_Dimensions(m,p,q);
    Initialize_Input_Planes(planes);
    Initialize_Interpolation_Points(points);
  end Initialize;

  procedure Store_Start_Pivots
              ( top,bottom : in Standard_Natural_Vectors.Vector ) is
  begin
    Standard_Natural_Vectors.Clear(start_top);
    Standard_Natural_Vectors.Clear(start_bottom);
    start_top := new Standard_Natural_Vectors.Vector'(top);
    start_bottom := new Standard_Natural_Vectors.Vector'(bottom);
  end Store_Start_Pivots;

  procedure Store_Target_Pivots
              ( top,bottom : in Standard_Natural_Vectors.Vector ) is
  begin
    Standard_Natural_Vectors.Clear(target_top);
    Standard_Natural_Vectors.Clear(target_bottom);
    target_top := new Standard_Natural_Vectors.Vector'(top);
    target_bottom := new Standard_Natural_Vectors.Vector'(bottom);
  end Store_Target_Pivots;

  procedure Store_Start ( x : in Standard_Complex_Vectors.Vector ) is
  begin
    Standard_Complex_Vectors.Clear(start);
    start := new Standard_Complex_Vectors.Vector'(x);
  end Store_Start;

  procedure Retrieve_Target ( x : out Standard_Complex_Vectors.Vector ) is
  begin
    x := target.all;
  end Retrieve_Target;

  function Degree_of_Freedom
             ( top,bottom : Standard_Natural_Vectors.Vector )
             return natural32 is

  -- DESCRIPTION :
  --   Returns the degree of freedom in the localization pattern
  --   defined by top and bottom pivots.

    res : natural32 := 0;

  begin
    for i in top'range loop
      res := res + bottom(i) - top(i);
    end loop;
    return res;
  end Degree_of_Freedom;

  function Sum_of_Difference
             ( a,b : Standard_Natural_Vectors.Vector ) return natural32 is

  -- DESCRIPTION :
  --   Returns the sum of the entries in the difference a - b.

  -- REQUIRED : a-b must be a vector of natural, nonnegative numbers!

    use Standard_Natural_Vectors;
    d : constant Vector(a'range) := a - b;

  begin
    return Sum(d);
  end Sum_of_Difference;

  function Start_Solution_Plane return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Assigns the solution vector stored in the vector start
  --   to a matrix representing a solution plane.
  --   The assignment is done columnwise, following the convention 
  --   that the first nonzero element in each column equals one.

    mpp : constant integer32 := m_dim+p_dim;
    res : Standard_Complex_Matrices.Matrix(1..mpp*(q_deg+1),1..p_dim);
    cnt : integer32 := start'first - 1;

  begin
    for j in 1..p_dim loop
      for i in res'range(1) loop
        if i < integer32(start_top(j)) or i > integer32(start_bottom(j)) then
          res(i,j) := Create(0.0);
        elsif i = integer32(start_top(j)) then
          res(i,j) := Create(1.0);
        else
          cnt := cnt + 1;
          res(i,j) := start(cnt);
        end if;
      end loop;
    end loop;
    return res;
  end Start_Solution_Plane;

  function Target_Solution_Plane return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Assigns the solution vector stored in the vector target
  --   to a matrix representing a solution plane.
  --   The assignment is done columnwise, following the convention 
  --   that the first nonzero element in each column equals one.

    mpp : constant integer32 := m_dim+p_dim;
    res : Standard_Complex_Matrices.Matrix(1..mpp*(q_deg+1),1..p_dim);
    cnt : integer32 := target'first - 1;

  begin
    for j in 1..p_dim loop
      for i in res'range(1) loop
        if i < integer32(target_top(j)) or i > integer32(target_bottom(j)) then
          res(i,j) := Create(0.0);
        -- elsif i = target_top(j)
        elsif i = integer32(start_top(j)) then
          res(i,j) := Create(1.0);
        else
          cnt := cnt + 1;
          res(i,j) := target(cnt);
        end if;
      end loop;
    end loop;
    return res;
  end Target_Solution_Plane;

  procedure Call_Path_Tracker
              ( x : in out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Calls the path tracker with the homotopy in Standard_Homotopy
  --   and one solution x.

    procedure Sil_Cont1 is
      new Linear_Single_Normal_Silent_Continue
            (Max_Norm,Standard_Homotopy.Eval,Standard_Homotopy.Diff,
             Standard_Homotopy.Diff);
    procedure Sil_Cont2 is
      new Linear_Single_Conditioned_Silent_Continue
            (Max_Norm,Standard_Homotopy.Eval,Standard_Homotopy.Diff,
             Standard_Homotopy.Diff);

    sol : Solution(x'last-1);
    s : Solu_Info; 
    tol : constant double_float := 1.0E-12;
    err : double_float := 1.0;
    w : integer32 := 1;
    dum : Standard_Floating_Vectors.Link_to_Vector;
    pp1,pp2 : Continuation_Parameters.Pred_Pars;
    cp1,cp2 : Continuation_Parameters.Corr_Pars;

  begin
    sol.t := x(x'last);
    sol.m := 1;
    sol.v := x(x'first..x'last-1);
    sol.err := 0.0;
    sol.rco := 1.0;
    sol.res := 0.0;
    s := Deep_Create(sol);
    Continuation_Parameters.Tune(2);
    pp1 := Continuation_Parameters.Create_for_Path;
    pp2 := Continuation_Parameters.Create_End_Game;
    cp1 := Continuation_Parameters.Create_for_Path;
    cp2 := Continuation_Parameters.Create_End_Game;
    Sil_Cont1(s,Create(1.0),tol,false,pp1,cp1);
    Sil_Cont2(s,Create(1.0),tol,false,0,w,dum,err,pp2,cp2);
    x(x'first..x'last-1) := s.sol.v;
    x(x'last) := s.sol.t;
    Clear(s);
  end Call_Path_Tracker;

  procedure Call_Path_Tracker
              ( file : in file_type;
                x : in out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Calls the path tracker with the homotopy in Standard_Homotopy
  --   and one solution x.

    procedure Rep_Cont1 is
      new Linear_Single_Normal_Reporting_Continue
            (Max_Norm,Standard_Homotopy.Eval,Standard_Homotopy.Diff,
             Standard_Homotopy.Diff);
    procedure Rep_Cont2 is
      new Linear_Single_Conditioned_Reporting_Continue
            (Max_Norm,Standard_Homotopy.Eval,Standard_Homotopy.Diff,
             Standard_Homotopy.Diff);

    sol : Solution(x'last-1);
    s : Solu_Info; 
    tol : constant double_float := 1.0E-12;
    w : integer32 := 1;
    err : double_float := 1.0;
    dum : Standard_Floating_Vectors.Link_to_Vector;
    pp1,pp2 : Continuation_Parameters.Pred_Pars;
    cp1,cp2 : Continuation_Parameters.Corr_Pars;

  begin
    sol.t := x(x'last);
    sol.m := 1;
    sol.v := x(x'first..x'last-1);
    sol.err := 0.0;
    sol.rco := 1.0;
    sol.res := 0.0;
    s := Deep_Create(sol);
    Continuation_Parameters.Tune(2);
    pp1 := Continuation_Parameters.Create_for_Path;
    pp2 := Continuation_Parameters.Create_End_Game;
    cp1 := Continuation_Parameters.Create_for_Path;
    cp2 := Continuation_Parameters.Create_End_Game;
    Process_io.Set_Output_Code(spc);
    Rep_Cont1(file,s,Create(1.0),tol,false,pp1,cp1);
    Rep_Cont2(file,s,Create(1.0),tol,false,0,w,dum,err,pp2,cp2);
    x(x'first..x'last-1) := s.sol.v;
    x(x'last) := s.sol.t;
    Clear(s);
  end Call_Path_Tracker;

  procedure Display_Homotopy_Settings
              ( file : in file_type; nv,level : in integer32;
                top,bottom : in Bracket; mixed : in boolean;
                xpm : in Standard_Complex_Poly_Matrices.Matrix;
                locmap : in Standard_Natural_Matrices.Matrix;
                bp : in Bracket_Polynomial ) is

  -- DESCRIPTION :
  --   Displays all settings for the Pieri homotopy as diagnostic.

  begin
    new_line(file);
    put(file,"Tracking a path with");
    put(file,"  m = "); put(file,m_dim,1);
    put(file,"  p = "); put(file,p_dim,1);
    put(file,"  q = "); put(file,q_deg,1); new_line(file);
    put(file,"Pivots at start  : ");
    put(file,"["); put(file,start_top); put(file," ]");
    put(file,"["); put(file,start_bottom); put(file," ]"); new_line(file);
    put(file,"Pivots at target : ");
    put(file,"["); put(file,target_top); put(file," ]");
    put(file,"["); put(file,target_bottom); put(file," ]"); new_line(file);
    put(file,"Degree of freedom at start  : ");
    put(file,Degree_of_Freedom(start_top.all,start_bottom.all),1);
    new_line(file);
    put(file,"Degree of freedom at target : ");
    put(file,Degree_of_Freedom(target_top.all,target_bottom.all),1);
    new_line(file);
    put(file,"The number of variables : "); put(file,nv,1); new_line(file);
    if mixed
     then Two_Set_up_Symbol_Table(natural32(m_dim),natural32(p_dim),
                                  natural32(q_deg),top,bottom);
     else One_Set_up_Symbol_Table(natural32(m_dim),natural32(p_dim),
                                  natural32(q_deg),top,bottom);
    end if;
    put_line(file,"The symbolic localization pattern : "); put(file,xpm);
    put_line(file,"The standard coordinate frame : "); put(file,locmap);
    put_line(file,"The general form of the intersection conditions : ");
    put(file,bp);
    put(file,"The level of the homotopy : ");
    put(file,level,1); new_line(file);
    put_line(file,"Moving towards the input plane :");
    put(file,input_planes(level).all,3);
  end Display_Homotopy_Settings;

  procedure Track_Top_Path is

  -- IMPORTANT NOTE : only works for top pivots

    mpp : constant integer32 := m_dim+p_dim;
    top : constant Bracket(target_top'range) := Bracket(target_top.all); 
    bottom : constant Bracket(target_bottom'range)
           := Bracket(target_bottom.all); 
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..mpp,1..p_dim)
        := Symbolic_Create(natural32(m_dim),natural32(p_dim),
                           natural32(q_deg),top,bottom);
    solpla : constant Standard_Complex_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
           := Start_Solution_Plane;
    locmap : constant Standard_Natural_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
           := Standard_Coordinate_Frame
                (natural32(m_dim),natural32(p_dim),natural32(q_deg),
                 top,bottom,solpla);
    bm : constant Bracket_Monomial
       := Maximal_Minors(natural32(mpp),natural32(mpp));
    bs : constant Bracket_System(0..integer32(Number_of_Brackets(bm)))
       := Minor_Equations(natural32(mpp),natural32(m_dim),bm);
    bp : constant Bracket_Polynomial := bs(1);
    nv : constant integer32 := integer32(Number_of_Variables(top,bottom));
    nv2 : constant integer32 := nv + 2;
    special : constant Standard_Complex_Matrices.Matrix(1..mpp,1..m_dim)
            := Special_Plane(m_dim,Modulo(top,natural32(mpp)));
    level : constant integer32
          := integer32(Degree_of_Freedom(target_top.all,target_bottom.all));
    start_equ,target_equ : Poly;
    homsys,lochom : Poly_Sys(1..level+1);
    solloc : constant Standard_Complex_Vectors.Vector
           := Column_Vector_Rep(locmap,solpla);
    locsol : Standard_Complex_Vectors.Vector(1..solloc'last+2);
 
  begin
   -- set up the homotopy
    for i in 1..level-1 loop
      declare
        eva : Standard_Complex_Poly_Matrices.Matrix(xpm'range(1),xpm'range(2))
            := Eval(xpm,interpolation_points(i),Create(1.0));
      begin
        homsys(i) := Expanded_Minors(input_planes(i).all,eva,bp);
        Standard_Complex_Poly_Matrices.Clear(eva);
      end;
    end loop;
    Swap(xpm,nv2-1,nv2);
    start_equ := Expanded_Minors(special,xpm,bp);
    target_equ := Expanded_Minors(input_planes(level).all,xpm,bp);
    Standard_Complex_Poly_Matrices.Clear(xpm);
    homsys(level) := Linear_Interpolation(target_equ,start_equ,nv2);
    Clear(target_equ); Clear(start_equ);
    homsys(level+1)
      := Moving_Parameter(nv2,nv2-1,nv2,Create(1.0),
                          (Create(1.0)/interpolation_points(level)));
    Divide_Common_Factor(homsys(level),nv2);
    lochom := Column_Localize(top,bottom,locmap,homsys);
   -- do the path tracking
    locsol(solloc'range) := solloc;
    locsol(locsol'last-1) := Create(1.0); -- start value for s
    locsol(locsol'last) := Create(0.0);   -- start value for t
    Standard_Homotopy.Create(lochom,locsol'last);
    Call_Path_Tracker(locsol);
    Standard_Homotopy.Clear;
    Clear(homsys); Clear(lochom);
   -- store target vector
    Standard_Complex_Vectors.Clear(target);
    target := new Standard_Complex_Vectors.Vector'(locsol(1..level));
  end Track_Top_Path;

  procedure Track_Top_Path ( file : in file_type ) is

  -- IMPORTANT NOTE : only works for top pivots

    mpp : constant integer32 := m_dim+p_dim;
    top : constant Bracket(target_top'range) := Bracket(target_top.all); 
    bottom : constant Bracket(target_bottom'range)
           := Bracket(target_bottom.all); 
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..mpp,1..p_dim)
        := Symbolic_Create(natural32(m_dim),natural32(p_dim),
                           natural32(q_deg),top,bottom);
    solpla : constant Standard_Complex_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
           := Start_Solution_Plane;
    locmap : constant Standard_Natural_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
           := Standard_Coordinate_Frame
                (natural32(m_dim),natural32(p_dim),natural32(q_deg),
                 top,bottom,solpla);
    bm : constant Bracket_Monomial
       := Maximal_Minors(natural32(mpp),natural32(mpp));
    bs : constant Bracket_System(0..integer32(Number_of_Brackets(bm)))
       := Minor_Equations(natural32(mpp),natural32(m_dim),bm);
    bp : constant Bracket_Polynomial := bs(1);
    nv : constant integer32 := integer32(Number_of_Variables(top,bottom));
    nv2 : constant integer32 := nv + 2;
    special : constant Standard_Complex_Matrices.Matrix(1..mpp,1..m_dim)
            := Special_Plane(m_dim,Modulo(top,natural32(mpp)));
    level : constant integer32
          := integer32(Degree_of_Freedom(target_top.all,target_bottom.all));
    start_equ,target_equ : Poly;
    homsys,lochom : Poly_Sys(1..level+1);
    solloc : constant Standard_Complex_Vectors.Vector
           := Column_Vector_Rep(locmap,solpla);
    locsol : Standard_Complex_Vectors.Vector(1..solloc'last+2);
 
  begin
    Display_Homotopy_Settings(file,nv,level,top,bottom,false,xpm,locmap,bp);
   -- set up the homotopy
    for i in 1..level-1 loop
      declare
        eva : Standard_Complex_Poly_Matrices.Matrix(xpm'range(1),xpm'range(2))
            := Eval(xpm,interpolation_points(i),Create(1.0));
      begin
        homsys(i) := Expanded_Minors(input_planes(i).all,eva,bp);
        Standard_Complex_Poly_Matrices.Clear(eva);
      end;
    end loop;
    Swap(xpm,nv2-1,nv2);
    start_equ := Expanded_Minors(special,xpm,bp);
    target_equ := Expanded_Minors(input_planes(level).all,xpm,bp);
    Standard_Complex_Poly_Matrices.Clear(xpm);
    homsys(level) := Linear_Interpolation(target_equ,start_equ,nv2);
    Clear(target_equ); Clear(start_equ);
    homsys(level+1)
      := Moving_Parameter(nv2,nv2-1,nv2,Create(1.0),
                          (Create(1.0)/interpolation_points(level)));
    Divide_Common_Factor(homsys(level),nv2);
    put_line(file,"The homotopy : "); put(file,homsys);
    lochom := Column_Localize(top,bottom,locmap,homsys);
    Reduce_Symbols(top,bottom,locmap);
    put_line("The homotopy in local coordinates : "); put(file,lochom);
   -- do the path tracking
    put_line(file,"The solution plane at the start : ");
    put(file,solpla,3);
    locsol(solloc'range) := solloc;
    locsol(locsol'last-1) := Create(1.0); -- start value for s
    locsol(locsol'last) := Create(0.0);   -- start value for t
    put_line(file,"The solution vector at the start : ");
    put_line(file,locsol);
    put_line(file,"Evaluation of the start solution in the local homotopy : ");
    put_line(file,Eval(lochom,locsol));
    Standard_Homotopy.Create(lochom,locsol'last);
    Call_Path_Tracker(file,locsol);
    Standard_Homotopy.Clear;
    put_line(file,"The solution vector at the end : ");
    put_line(file,locsol);
    put_line(file,"Evaluation of the end solution in the local homotopy : ");
    put_line(file,Eval(lochom,locsol));
    Clear(homsys); Clear(lochom);
   -- store target vector
    Standard_Complex_Vectors.Clear(target);
    target := new Standard_Complex_Vectors.Vector'(locsol(1..level));
  end Track_Top_Path;

  procedure Track_Bottom_Path is

  -- IMPORTANT NOTE : only works for bottom pivots

    mpp : constant integer32 := m_dim+p_dim;
    top : constant Bracket(target_top'range) := Bracket(target_top.all); 
    bottom : constant Bracket(target_bottom'range)
           := Bracket(target_bottom.all); 
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..mpp,1..p_dim)
        := Symbolic_Create(natural32(m_dim),natural32(p_dim),
                           natural32(q_deg),top,bottom);
    solpla : constant Standard_Complex_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
           := Start_Solution_Plane;
    locmap : constant Standard_Natural_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
           := Standard_Coordinate_Frame
                (natural32(m_dim),natural32(p_dim),natural32(q_deg),
                 top,bottom,solpla);
    bm : constant Bracket_Monomial
       := Maximal_Minors(natural32(mpp),natural32(mpp));
    bs : constant Bracket_System(0..integer32(Number_of_Brackets(bm)))
       := Minor_Equations(natural32(mpp),natural32(m_dim),bm);
    bp : constant Bracket_Polynomial := bs(1);
    nv : constant integer32 := integer32(Number_of_Variables(top,bottom));
    nv2 : constant integer32 := nv + 2;
    special : constant Standard_Complex_Matrices.Matrix(1..mpp,1..m_dim)
            := Special_Plane(m_dim,Modulo(bottom,natural32(mpp)));
    level : constant integer32
          := integer32(Degree_of_Freedom(target_top.all,target_bottom.all));
    start_equ : Poly := Expanded_Minors(special,xpm,bp);
    target_equ : Poly := Expanded_Minors(input_planes(level).all,xpm,bp);
    homsys,lochom : Poly_Sys(1..level+1);
    solloc : constant Standard_Complex_Vectors.Vector
           := Column_Vector_Rep(locmap,solpla);
    locsol : Standard_Complex_Vectors.Vector(1..solloc'last+2);
 
  begin
   -- set up the homotopy
    for i in 1..level-1 loop
      declare
        eva : Standard_Complex_Poly_Matrices.Matrix(xpm'range(1),xpm'range(2))
            := Eval(xpm,interpolation_points(i),Create(1.0));
      begin
        homsys(i) := Expanded_Minors(input_planes(i).all,eva,bp);
        Standard_Complex_Poly_Matrices.Clear(eva);
      end;
    end loop;
    Standard_Complex_Poly_Matrices.Clear(xpm);
    homsys(level) := Linear_Interpolation(target_equ,start_equ,nv2);
    Clear(target_equ); Clear(start_equ);
    homsys(level+1) := Moving_Parameter(nv2,nv2-1,nv2,Create(1.0),
                                        interpolation_points(level));
    lochom := Column_Localize(top,bottom,locmap,homsys);
   -- do the path tracking
    locsol(solloc'range) := solloc;
    locsol(locsol'last-1) := Create(1.0); -- start value for s
    locsol(locsol'last) := Create(0.0);   -- start value for t
    Standard_Homotopy.Create(lochom,locsol'last);
    Call_Path_Tracker(locsol);
    Standard_Homotopy.Clear;
    Clear(homsys); Clear(lochom);
   -- store target vector
    Standard_Complex_Vectors.Clear(target);
    target := new Standard_Complex_Vectors.Vector'(locsol(1..level));
  end Track_Bottom_Path;

  procedure Track_Bottom_Path ( file : in file_type ) is

  -- IMPORTANT NOTE : only works for bottom pivots

    mpp : constant integer32 := m_dim+p_dim;
    top : constant Bracket(target_top'range) := Bracket(target_top.all); 
    bottom : constant Bracket(target_bottom'range)
           := Bracket(target_bottom.all); 
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..mpp,1..p_dim)
        := Symbolic_Create(natural32(m_dim),natural32(p_dim),
                           natural32(q_deg),top,bottom);
    solpla : constant Standard_Complex_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
           := Start_Solution_Plane;
    locmap : constant Standard_Natural_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
           := Standard_Coordinate_Frame
                (natural32(m_dim),natural32(p_dim),natural32(q_deg),
                 top,bottom,solpla);
    bm : constant Bracket_Monomial
       := Maximal_Minors(natural32(mpp),natural32(mpp));
    bs : constant Bracket_System(0..integer32(Number_of_Brackets(bm)))
       := Minor_Equations(natural32(mpp),natural32(m_dim),bm);
    bp : constant Bracket_Polynomial := bs(1);
    nv : constant integer32 := integer32(Number_of_Variables(top,bottom));
    nv2 : constant integer32 := nv + 2;
    special : constant Standard_Complex_Matrices.Matrix(1..mpp,1..m_dim)
            := Special_Plane(m_dim,Modulo(bottom,natural32(mpp)));
    level : constant integer32
          := integer32(Degree_of_Freedom(target_top.all,target_bottom.all));
    start_equ : Poly := Expanded_Minors(special,xpm,bp);
    target_equ : Poly := Expanded_Minors(input_planes(level).all,xpm,bp);
    homsys,lochom : Poly_Sys(1..level+1);
    solloc : constant Standard_Complex_Vectors.Vector
           := Column_Vector_Rep(locmap,solpla);
    locsol : Standard_Complex_Vectors.Vector(1..solloc'last+2);
 
  begin
    Display_Homotopy_Settings(file,nv,level,top,bottom,false,xpm,locmap,bp);
   -- set up the homotopy
    for i in 1..level-1 loop
      declare
        eva : Standard_Complex_Poly_Matrices.Matrix(xpm'range(1),xpm'range(2))
            := Eval(xpm,interpolation_points(i),Create(1.0));
      begin
        homsys(i) := Expanded_Minors(input_planes(i).all,eva,bp);
        Standard_Complex_Poly_Matrices.Clear(eva);
      end;
    end loop;
    Standard_Complex_Poly_Matrices.Clear(xpm);
    homsys(level) := Linear_Interpolation(target_equ,start_equ,nv2);
    Clear(target_equ); Clear(start_equ);
    homsys(level+1) := Moving_Parameter(nv2,nv2-1,nv2,Create(1.0),
                                        interpolation_points(level));
    put_line(file,"The homotopy : "); put(file,homsys);
    lochom := Column_Localize(top,bottom,locmap,homsys);
    Reduce_Symbols(top,bottom,locmap);
    put_line("The homotopy in local coordinates : "); put(file,lochom);
   -- do the path tracking
    put_line(file,"The solution plane at the start : ");
    put(file,solpla,3);
    locsol(solloc'range) := solloc;
    locsol(locsol'last-1) := Create(1.0); -- start value for s
    locsol(locsol'last) := Create(0.0);   -- start value for t
    put_line(file,"The solution vector at the start : ");
    put_line(file,locsol);
    put_line(file,"Evaluation of the start solution in the local homotopy : ");
    put_line(file,Eval(lochom,locsol));
    Standard_Homotopy.Create(lochom,locsol'last);
    Call_Path_Tracker(file,locsol);
    Standard_Homotopy.Clear;
    put_line(file,"The solution vector at the end : ");
    put_line(file,locsol);
    put_line(file,"Evaluation of the end solution in the local homotopy : ");
    put_line(file,Eval(lochom,locsol));
    Clear(homsys); Clear(lochom);
   -- store target vector
    Standard_Complex_Vectors.Clear(target);
    target := new Standard_Complex_Vectors.Vector'(locsol(1..level));
  end Track_Bottom_Path;

  procedure Track_Mixed_Path is

  -- IMPORTANT NOTE : only works for mixed top+bottom pivots

    mpp : constant integer32 := m_dim+p_dim;
    top : constant Bracket(target_top'range) := Bracket(target_top.all); 
    bottom : constant Bracket(target_bottom'range)
           := Bracket(target_bottom.all); 
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..mpp,1..p_dim)
        := Symbolic_Create(natural32(m_dim),natural32(p_dim),
                           natural32(q_deg),top,bottom);
    solpla : constant Standard_Complex_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
           := Start_Solution_Plane;
    locmap : constant Standard_Natural_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
           := Standard_Coordinate_Frame
                (natural32(m_dim),natural32(p_dim),natural32(q_deg),
                 top,bottom,solpla);
    bm : constant Bracket_Monomial
       := Maximal_Minors(natural32(mpp),natural32(mpp));
    bs : constant Bracket_System(0..integer32(Number_of_Brackets(bm)))
       := Minor_Equations(natural32(mpp),natural32(m_dim),bm);
    bp : constant Bracket_Polynomial := bs(1);
    nv : constant integer32 := integer32(Number_of_Variables(top,bottom));
    nv2 : constant integer32 := nv + 2;
    top_special : constant Standard_Complex_Matrices.Matrix(1..mpp,1..m_dim)
                := Special_Plane(m_dim,Modulo(top,natural32(mpp)));
    bot_special : constant Standard_Complex_Matrices.Matrix(1..mpp,1..m_dim)
                := Special_Plane(m_dim,Modulo(bottom,natural32(mpp)));
    level : constant integer32 
          := integer32(Degree_of_Freedom(target_top.all,target_bottom.all));
    top_start,bot_start,top_target,bot_target : Poly;
    homsys,lochom : Poly_Sys(1..level+2);
    solloc : constant Standard_Complex_Vectors.Vector
           := Column_Vector_Rep(locmap,solpla);
    locsol : Standard_Complex_Vectors.Vector(1..solloc'last+3);
    map : Standard_Complex_Poly_Matrices.Matrix(xpm'range(1),xpm'range(2))
        := Insert(xpm,nv2);
 
  begin
   -- set up the homotopy
    bot_target := Expanded_Minors(input_planes(level-1).all,map,bp);
    bot_start := Expanded_Minors(bot_special,map,bp);
    homsys(level-1) := Linear_Interpolation(bot_target,bot_start,nv2+1);
    Divide_Common_Factor(homsys(level-1),nv2+1);
    homsys(level+1)
      := Moving_Parameter(nv2+1,nv2-1,nv2+1,Create(1.0),
                          interpolation_points(level-1));
    Clear(bot_target); Clear(bot_start);
    Swap(map,nv2-1,nv2);
    for i in 1..level-2 loop
      declare
        eva : Standard_Complex_Poly_Matrices.Matrix(map'range(1),map'range(2))
            := Eval(map,interpolation_points(i),Create(1.0));
      begin
        homsys(i) := Expanded_Minors(input_planes(i).all,eva,bp);
        Standard_Complex_Poly_Matrices.Clear(eva);
      end;
    end loop;
    Swap(map,nv2,nv2+1);
    top_target := Expanded_Minors(input_planes(level).all,map,bp);
    top_start := Expanded_Minors(top_special,map,bp);
    homsys(level) := Linear_Interpolation(top_target,top_start,nv2+1);
    homsys(level+2)
      := Moving_Parameter(nv2+1,nv2,nv2+1,Create(1.0),
                          (Create(1.0)/interpolation_points(level)));
    Divide_Common_Factor(homsys(level),nv2+1);
    Clear(top_target); Clear(top_start);
    Standard_Complex_Poly_Matrices.Clear(xpm);
    Standard_Complex_Poly_Matrices.Clear(map);
    lochom := Column_Localize(top,bottom,locmap,homsys);
   -- do the path tracking
    locsol(solloc'range) := solloc;
    locsol(locsol'last-2) := Create(1.0); -- start value for s1
    locsol(locsol'last-1) := Create(1.0); -- start value for s2
    locsol(locsol'last) := Create(0.0);   -- start value for t
    Standard_Homotopy.Create(lochom,locsol'last);
    Call_Path_Tracker(locsol);
    Standard_Homotopy.Clear;
    Clear(homsys); Clear(lochom);
   -- store target vector
    Standard_Complex_Vectors.Clear(target);
    target := new Standard_Complex_Vectors.Vector'(locsol(1..level));
  end Track_Mixed_Path;

  procedure Track_Mixed_Path ( file : in file_type ) is

  -- IMPORTANT NOTE : only works for mixed top+bottom pivots

    mpp : constant integer32 := m_dim+p_dim;
    top : constant Bracket(target_top'range) := Bracket(target_top.all); 
    bottom : constant Bracket(target_bottom'range)
           := Bracket(target_bottom.all); 
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..mpp,1..p_dim)
        := Symbolic_Create(natural32(m_dim),natural32(p_dim),
                           natural32(q_deg),top,bottom);
    solpla : constant Standard_Complex_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
           := Start_Solution_Plane;
    locmap : constant Standard_Natural_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
           := Standard_Coordinate_Frame
                (natural32(m_dim),natural32(p_dim),
                 natural32(q_deg),top,bottom,solpla);
    bm : constant Bracket_Monomial
       := Maximal_Minors(natural32(mpp),natural32(mpp));
    bs : constant Bracket_System(0..integer32(Number_of_Brackets(bm)))
       := Minor_Equations(natural32(mpp),natural32(m_dim),bm);
    bp : constant Bracket_Polynomial := bs(1);
    nv : constant integer32 := integer32(Number_of_Variables(top,bottom));
    nv2 : constant integer32 := nv + 2;
    top_special : constant Standard_Complex_Matrices.Matrix(1..mpp,1..m_dim)
                := Special_Plane(m_dim,Modulo(top,natural32(mpp)));
    bot_special : constant Standard_Complex_Matrices.Matrix(1..mpp,1..m_dim)
                := Special_Plane(m_dim,Modulo(bottom,natural32(mpp)));
    level : constant integer32
          := integer32(Degree_of_Freedom(target_top.all,target_bottom.all));
    top_start,bot_start,top_target,bot_target : Poly;
    homsys,lochom : Poly_Sys(1..level+2);
    solloc : constant Standard_Complex_Vectors.Vector
           := Column_Vector_Rep(locmap,solpla);
    locsol : Standard_Complex_Vectors.Vector(1..solloc'last+3);
    map : Standard_Complex_Poly_Matrices.Matrix(xpm'range(1),xpm'range(2))
        := Insert(xpm,nv2);
 
  begin
    Display_Homotopy_Settings(file,nv,level,top,bottom,true,xpm,locmap,bp);
   -- set up the homotopy
    bot_target := Expanded_Minors(input_planes(level-1).all,map,bp);
    bot_start := Expanded_Minors(bot_special,map,bp);
    homsys(level-1) := Linear_Interpolation(bot_target,bot_start,nv2+1);
    Divide_Common_Factor(homsys(level-1),nv2+1);
    homsys(level+1)
      := Moving_Parameter(nv2+1,nv2-1,nv2+1,Create(1.0),
                          interpolation_points(level-1));
    Clear(bot_target); Clear(bot_start);
    Swap(map,nv2-1,nv2);
    for i in 1..level-2 loop
      declare
        eva : Standard_Complex_Poly_Matrices.Matrix(map'range(1),map'range(2))
            := Eval(map,interpolation_points(i),Create(1.0));
      begin
        homsys(i) := Expanded_Minors(input_planes(i).all,eva,bp);
        Standard_Complex_Poly_Matrices.Clear(eva);
      end;
    end loop;
    Swap(map,nv2,nv2+1);
    top_target := Expanded_Minors(input_planes(level).all,map,bp);
    top_start := Expanded_Minors(top_special,map,bp);
    homsys(level) := Linear_Interpolation(top_target,top_start,nv2+1);
    homsys(level+2)
      := Moving_Parameter(nv2+1,nv2,nv2+1,Create(1.0),
                          (Create(1.0)/interpolation_points(level)));
    Divide_Common_Factor(homsys(level),nv2+1);
    Clear(top_target); Clear(top_start);
    Standard_Complex_Poly_Matrices.Clear(xpm);
    Standard_Complex_Poly_Matrices.Clear(map);
    put_line(file,"The homotopy : "); put(file,homsys);
    lochom := Column_Localize(top,bottom,locmap,homsys);
    Reduce_Symbols(top,bottom,locmap);
    put_line("The homotopy in local coordinates : "); put(file,lochom);
   -- do the path tracking
    put_line(file,"The solution plane at the start : ");
    put(file,solpla,3);
    locsol(solloc'range) := solloc;
    locsol(locsol'last-2) := Create(1.0); -- start value for s1
    locsol(locsol'last-1) := Create(1.0); -- start value for s2
    locsol(locsol'last) := Create(0.0);   -- start value for t
    put_line(file,"The solution vector at the start : ");
    put_line(file,locsol);
    put_line(file,"Evaluation of the start solution in the local homotopy : ");
    put_line(file,Eval(lochom,locsol));
    Standard_Homotopy.Create(lochom,locsol'last);
    Call_Path_Tracker(file,locsol);
    Standard_Homotopy.Clear;
    Clear(homsys); Clear(lochom);
    put_line(file,"The solution vector at the end : ");
    put_line(file,locsol);
    put_line(file,"Evaluation of the end solution in the local homotopy : ");
    put_line(file,Eval(lochom,locsol));
   -- store target vector
    Standard_Complex_Vectors.Clear(target);
    target := new Standard_Complex_Vectors.Vector'(locsol(1..level));
  end Track_Mixed_Path;

  procedure Track_Path is

    d_top : constant natural32
          :=  Sum_of_Difference(start_top.all,target_top.all);
    d_bot : constant natural32
          :=  Sum_of_Difference(target_bottom.all,start_bottom.all);

  begin
    if d_top = 0 then
      Track_Bottom_Path;
    elsif d_bot = 0 then
      Track_Top_Path;
    else
      Track_Mixed_Path;
    end if;
  end Track_Path;

  procedure Track_Path ( file : in file_type ) is

    d_top : constant natural32
          :=  Sum_of_Difference(start_top.all,target_top.all);
    d_bot : constant natural32
          :=  Sum_of_Difference(target_bottom.all,start_bottom.all);

  begin
    if d_top = 0 then
      Track_Bottom_Path(file);
    elsif d_bot = 0 then
      Track_Top_Path(file);
    else
      Track_Mixed_Path(file);
    end if;
  end Track_Path;

  function Verify_Determinants return double_float is

    res : double_float := 0.0;
    mpp : constant integer32 := m_dim+p_dim;
    top : constant Bracket(target_top'range) := Bracket(target_top.all);
    bottom : constant Bracket(target_bottom'range)
           := Bracket(target_bottom.all);
    xpm : constant Standard_Complex_Poly_Matrices.Matrix(1..mpp,1..p_dim)
        := Symbolic_Create(natural32(m_dim),natural32(p_dim),
                           natural32(q_deg),top,bottom);
    level : constant integer32
          := integer32(Degree_of_Freedom(target_top.all,target_bottom.all));
    stapla : constant Standard_Complex_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
           := Start_Solution_Plane;
   -- solpla : constant Standard_Complex_Matrices.Matrix
   --                     (1..mpp*(q_deg+1),1..p_dim)
   --        := Target_Solution_Plane;
    locmap : constant Standard_Natural_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
        -- := Standard_Coordinate_Frame(m_dim,p_dim,q_deg,top,bottom,solpla);
           := Standard_Coordinate_Frame
                (natural32(m_dim),natural32(p_dim),natural32(q_deg),
                 top,bottom,stapla);
    lxpm : constant Standard_Complex_Poly_Matrices.Matrix
                      (xpm'range(1),xpm'range(2))
         := Column_Localize(top,bottom,locmap,xpm);
    xs : constant Standard_Complex_Poly_Matrices.Matrix
                    (xpm'range(1),xpm'range(2))
       := Substitute(lxpm,target(1..level));
    eva : Standard_Complex_Poly_Matrices.Matrix(xpm'range(1),xpm'range(2));
    mat : Standard_Complex_Matrices.Matrix(eva'range(1),eva'range(2));
    det : Complex_Number;

  begin
    for i in 1..level loop
      eva := Eval(xs,interpolation_points(i),Create(1.0));
      mat := Convert(eva);
      det := Determinant(mat,input_planes(i).all);
      res := res + AbsVal(det);
    end loop;
    return res;
  end Verify_Determinants;

  function Verify_Determinants ( file : in file_type ) return double_float is

    res : double_float := 0.0;
    mpp : constant integer32 := m_dim+p_dim;
    top : constant Bracket(target_top'range) := Bracket(target_top.all);
    bottom : constant Bracket(target_bottom'range)
           := Bracket(target_bottom.all);
    xpm : constant Standard_Complex_Poly_Matrices.Matrix(1..mpp,1..p_dim)
        := Symbolic_Create(natural32(m_dim),natural32(p_dim),
                           natural32(q_deg),top,bottom);
    level : constant integer32
          := integer32(Degree_of_Freedom(target_top.all,target_bottom.all));
    stapla : constant Standard_Complex_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
           := Start_Solution_Plane;
    solpla : constant Standard_Complex_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
           := Target_Solution_Plane;
    locmap : constant Standard_Natural_Matrices.Matrix
                        (1..mpp*(q_deg+1),1..p_dim)
        -- := Standard_Coordinate_Frame(m_dim,p_dim,q_deg,top,bottom,solpla);
           := Standard_Coordinate_Frame
                (natural32(m_dim),natural32(p_dim),natural32(q_deg),
                 top,bottom,stapla);
    lxpm : constant Standard_Complex_Poly_Matrices.Matrix
                      (xpm'range(1),xpm'range(2))
         := Column_Localize(top,bottom,locmap,xpm);
    xs : constant Standard_Complex_Poly_Matrices.Matrix
                    (xpm'range(1),xpm'range(2))
       := Substitute(lxpm,target(1..level));
    eva : Standard_Complex_Poly_Matrices.Matrix(xpm'range(1),xpm'range(2));
    mat : Standard_Complex_Matrices.Matrix(eva'range(1),eva'range(2));
    det : Complex_Number;

  begin
    One_Set_up_Symbol_Table(natural32(m_dim),natural32(p_dim),
                            natural32(q_deg),top,bottom);
    put_line(file,"The symbolic localization pattern : "); put(file,xpm);
    put_line(file,"The localization map : "); put(file,locmap);
    put_line(file,"The solution plane : "); put(file,solpla,3);
    put_line("The map after substitution : "); put(file,xs);
    put(file,"The current level : "); put(file,level,1); new_line(file);
    for i in 1..level loop
      eva := Eval(xs,interpolation_points(i),Create(1.0));
      mat := Convert(eva);
      det := Determinant(mat,input_planes(i).all);
      put(file,"determinant "); put(file,i,1); put(file," = ");
      put(file,det); new_line(file);
      res := res + AbsVal(det);
    end loop;
    return res;
  end Verify_Determinants;

  procedure Clear is
  begin
    Deep_Clear(input_planes);
    Standard_Complex_Vectors.Clear(interpolation_points);
    Standard_Natural_Vectors.Clear(start_top);
    Standard_Natural_Vectors.Clear(start_bottom);
    Standard_Natural_Vectors.Clear(target_top);
    Standard_Natural_Vectors.Clear(target_bottom);
    Standard_Complex_Vectors.Clear(start);
    Standard_Complex_Vectors.Clear(target);
  end Clear;

end Pieri_Homotopy;
