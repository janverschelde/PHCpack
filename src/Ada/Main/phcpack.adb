with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_Randomizers;  use Standard_Complex_Poly_Randomizers;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_Scaling;                   use Standard_Scaling;
with Reduction_of_Polynomial_Systems;    use Reduction_of_Polynomial_Systems;  
with Standard_Homotopy;
with Total_Degree_Start_Systems;         use Total_Degree_Start_Systems;
with Partitions_of_Sets_Of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with m_Homogeneous_Bezout_Numbers;       use m_Homogeneous_Bezout_Numbers;
with m_Homogeneous_Start_Systems;        use m_Homogeneous_Start_Systems;
with Set_Structure;
with Degree_Sets_Tables;                 use Degree_Sets_Tables;
with Standard_Linear_Product_System;
with Random_Product_Start_Systems;       use Random_Product_Start_Systems;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Black_Mixed_Volume_Computations;    use Black_Mixed_Volume_Computations;
with Black_Polyhedral_Continuations;     use Black_Polyhedral_Continuations;
with BKK_Bound_Computations;             use BKK_Bound_Computations;
with Continuation_Parameters;
with Standard_IncFix_Continuation;       use Standard_IncFix_Continuation; 
with Standard_Root_Refiners;             use Standard_Root_Refiners;

package body PHCpack is

-- GENERAL UTILITIES :

  function Minimum ( a,b : natural32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the minimum of a and b.

  begin
    if a <= b
     then return a;
     else return b;
    end if;
  end Minimum;

  function Minimum ( a,b,c : natural32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the minimum of the three numbers a, b, and c.

     minbc : constant natural32 := Minimum(b,c);

  begin
    if a <= minbc
     then return a;
     else return minbc;
    end if;
  end Minimum;

  function Minimum ( a,b,c,d : natural32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the minimum of the four numbers.

    minbcd : constant natural32 := Minimum(b,c,d);

  begin
    if a <= minbcd
     then return a;
     else return minbcd;
    end if;
  end Minimum;

-- 1. PRE-PROCESSING : SCALING AND REDUCTION

  procedure Equation_Scaling
                ( file : in file_type; p : in Poly_Sys; s : out Poly_Sys ) is

    res : Poly_Sys(p'range);

  begin
    Copy(p,res);
    Scale(res);
    put(file,res);
    s := res;
  end Equation_Scaling;

  procedure Linear_Reduction
              ( file : in file_type; p : in Poly_Sys; r : out Poly_Sys ) is

    res : Poly_Sys(p'range);
    success,inconsistent,infinite : boolean := false;

  begin
    Copy(p,res);
    reduce(res,success,inconsistent,infinite);
    if success then
      if inconsistent
       then put_line(file,"system is inconsistent");
      end if;
      if infinite
       then put_line(file,"system has infinite number of solutions");
      end if;
    end if;
    put(file,res);
    r := res;
  end Linear_Reduction;

-- 2. ROOT COUNTING AND CONSTRUCTION OF START SYSTEM

  procedure Total_Degree
              ( file : in file_type; p : in Poly_Sys; d : out natural32 ) is
  begin
    d := m_Homogeneous_Bezout_Numbers.Total_Degree(p);
  end Total_Degree;

  procedure Total_Degree
              ( file : in file_type; p : in Poly_Sys; d : out natural32;
                q : out Poly_Sys; qsols : out Solution_List ) is

    qq : Poly_Sys(p'range);
    qqsols : Solution_List;

  begin
    d := m_Homogeneous_Bezout_Numbers.Total_Degree(p);
    Start_System(p,qq,qqsols);
    q := qq; qsols := qqsols;
  end Total_Degree;

  function Set_Structure_Bound ( p : Poly_Sys ) return natural32 is
  begin
    Random_Product_Start_Systems.Build_Set_Structure(p);
    return natural32(Permanent(Degree_Sets_Tables.Create));
  end Set_Structure_Bound;

  procedure Implicit_Lifting
              ( file : in file_type; p : in Poly_Sys; mv : out natural32 ) is
  begin
    mv := BKK_by_Implicit_Lifting(p);
  end Implicit_Lifting;

  procedure Implicit_Lifting
              ( file : in file_type; p : in Poly_Sys; mv : out natural32;
                q : out Poly_Sys; qsols : out Solution_List ) is

    qq : constant Poly_Sys(p'range) := Complex_Randomize1(p);
    qqsols : Solution_List := Solve_by_Implicit_Lifting(file,qq);

  begin
    mv := Length_Of(qqsols);
    Set_Continuation_Parameter(qqsols,Create(0.0));
    q := qq; qsols := qqsols;
  end Implicit_Lifting;

  procedure Static_Lifting
              ( file : in file_type; p : in Poly_Sys; mv : out natural32 ) is
  begin
    mv := BKK_by_Static_Lifting(file,p);
  end Static_Lifting;

  procedure Static_Lifting
              ( file : in file_type; p : in Poly_Sys; mv : out natural32;
                q : out Poly_Sys; qsols : out Solution_List ) is

    qq : constant Poly_Sys(p'range) := Complex_Randomize1(p);
    qqsols : Solution_List := Solve_by_Static_Lifting(file,qq);

  begin
    mv := Length_Of(qqsols);
    Set_Continuation_Parameter(qqsols,Create(0.0));
    q := qq; qsols := qqsols;
  end Static_Lifting;

  procedure Count_Roots
              ( p : in out Poly_Sys; rc : out natural32;
                q : out Poly_Sys; qsols : out Solution_List ) is

    n : constant natural32 := natural32(p'length);
    d,m,bz,bs,mv,nl : natural32;
    bz64 : natural64;
    z : Partition(1..n);
    lifsup : Link_to_Array_of_Lists;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    mixsub : Mixed_Subdivision;

  begin
    d := m_Homogeneous_Bezout_Numbers.Total_Degree(p);
    PB(p,bz64,m,z);
    bz := natural32(bz64);  -- patch ...
    bs := Set_Structure_Bound(p);
    if m = 1
     then bz := d;
    end if;
    if n <= 10 then
      declare -- problems with systems with one monomial equations
      begin
        Black_Box_Mixed_Volume_Computation(p,mix,lifsup,mixsub,mv);
      exception
        when others => mv := 0;
      end;
    end if;
    if mv = 0
     then rc := Minimum(d,bz,bs);
     else rc := Minimum(d,bz,bs,mv);
    end if;
    if rc = d then
      Start_System(p,q,qsols);
    elsif rc = bz then
      m_Homogeneous_Start_System(p,z,q,qsols);
    elsif rc = bs then
      Standard_Linear_Product_System.Init(n);
      Build_Random_Product_System(n);
      q := Standard_Linear_Product_System.Polynomial_System;
      Standard_Linear_Product_System.Solve(qsols,nl);
      Set_Structure.Clear;
      Standard_Linear_Product_System.Clear; 
    else
      Black_Box_Polyhedral_Continuation(p,mix.all,lifsup.all,mixsub,q,qsols);
    end if;
  end Count_Roots;

-- 3. POLYNOMIAL CONTINUATION

  procedure Artificial_Parameter_Continuation
              ( p,q : in Poly_Sys; sols : in out Solution_List;
                k : in natural32 := 2;
                a : in Complex_Number := Create(1.0); 
                target : in Complex_Number := Create(1.0) ) is

    procedure Cont is 
      new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                          Standard_Homotopy.Eval,Standard_Homotopy.Diff);

  begin
    Standard_Homotopy.Create(p,q,k,a);
    Continuation_Parameters.Tune(0);
    Cont(sols,false,0,target);
    Standard_Homotopy.Clear;
  end Artificial_Parameter_Continuation;

  procedure Artificial_Parameter_Continuation
              ( file : in file_type; p,q : in Poly_Sys;
                sols : in out Solution_List;
                k : in natural32 := 2;
                a : in Complex_Number := Create(1.0); 
                target : in Complex_Number := Create(1.0) ) is

    procedure Cont is 
      new Reporting_Continue(Max_Norm,Standard_Homotopy.Eval,
                             Standard_Homotopy.Eval,Standard_Homotopy.Diff);

  begin
    Standard_Homotopy.Create(p,q,k,a);
    Continuation_Parameters.Tune(0);
    Cont(file,sols,false,0,target);
    Standard_Homotopy.Clear;
  end Artificial_Parameter_Continuation;

  procedure Natural_Parameter_Continuation
              ( file : in file_type; h : in Poly_Sys; k : in natural32;
                t0,t1 : in Complex_Number; sols : in out Solution_List ) is
  begin
    null;
  end Natural_Parameter_Continuation;

-- 4. POST-PROCESSING : VALIDATION

  procedure Refine_Roots
               ( p : in Poly_Sys; sols : in out Solution_List ) is

    epsxa,epsfa : constant double_float := 1.0E-8;   -- defaults
    tolsing : constant double_float := 1.0E-8;
    maxit : constant natural32 := 3;
    numit : natural32 := 0;
    deflate : boolean := false;

  begin
    Silent_Root_Refiner(p,sols,epsxa,epsfa,tolsing,numit,maxit,deflate);
  end Refine_Roots;

  procedure Refine_Roots
              ( file : in file_type; p : in Poly_Sys;
                sols : in out Solution_List ) is

    epsxa,epsfa : constant double_float := 1.0E-8;   -- defaults
    tolsing : constant double_float := 1.0E-8;
    maxit : constant natural32 := 3;
    numit : natural32 := 0;
    deflate : boolean := false;

  begin
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,false);
  end Refine_Roots;

end PHCpack;
