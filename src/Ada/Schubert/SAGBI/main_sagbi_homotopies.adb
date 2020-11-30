with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Numbers_io;                         use Numbers_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Complex_Matrices;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Symbol_Table;                       use Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Functions;    use Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_Jacomats;
with Exponent_Vectors;                   use Exponent_Vectors;
with Standard_Poly_Laur_Convertors;      use Standard_Poly_Laur_Convertors;
with Matrix_Indeterminates;
with Standard_Homotopy;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_IncFix_Continuation;       use Standard_IncFix_Continuation;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Main_Poly_Continuation;             use Main_Poly_Continuation;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Integer_Polyhedral_Continuation;    use Integer_Polyhedral_Continuation;
with Standard_Integer32_Triangulations;  use Standard_Integer32_Triangulations;
with Standard_Dynamic32_Triangulations;  use Standard_Dynamic32_Triangulations;
with Triangulations_and_Subdivisions;    use Triangulations_and_Subdivisions;
with Bracket_Expansions;                 use Bracket_Expansions;
with SAGBI_Homotopies;                   use SAGBI_Homotopies;
with Matrix_Homotopies;
with Matrix_Homotopies_io;
with Osculating_Planes;                  use Osculating_Planes;

package body Main_SAGBI_Homotopies is

  procedure Main ( n,d : in natural32 ) is

    file : file_type;

  begin
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Main(file,n,d);
  end Main;

  function Convert ( mat : Standard_Floating_Matrices.Matrix )
                   return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(mat'range(1),mat'range(2));

  begin
    for i in mat'range(1) loop
      for j in mat'range(2) loop
        res(i,j) := Create(mat(i,j));
      end loop;
    end loop;
    return res;
  end Convert;

  function Random_Osculating_SAGBI_Homotopy
             ( file : file_type; n,d,dim : natural32;
               locmap : Standard_Natural_Matrices.Matrix ) return Poly_Sys is

  -- DESCRIPTION :
  --   Generates planes that are osculating a rational normal curve.
  --   The system on return is the SAGBI homotopy.

    p : Poly;
    l : Poly_Sys(1..integer32(dim));
    s : double_float;
    inc : constant double_float := 2.0/(double_float(dim));
    mat : Standard_Floating_Matrices.Matrix(1..integer32(n),1..integer32(n-d));

  begin
   -- p := Lifted_Localized_Laplace_Expansion(n,d);
    p := Lifted_Localized_Laplace_Expansion(locmap);
    new_line(file);
    put_line(file,"The selected s-values : ");
    s := Random;
    for i in l'range loop
      s := s + inc;
      if s >= 1.0
       then s := s - 2.0;
      end if;      -- s lies in [-1.0,+1.0]
      put(file,s); new_line(file);
      mat := Orthogonal_Basis(n,n-d,s);
      Matrix_Homotopies.Add_Target(natural32(i),Convert(mat));
      l(i) := Intersection_Condition(mat,p);
    end loop;
    return l;
  end Random_Osculating_SAGBI_Homotopy;

  function Given_Osculating_SAGBI_Homotopy
             ( file : file_type; n,d,dim : natural32;
               locmap : Standard_Natural_Matrices.Matrix ) return Poly_Sys is

  -- DESCRIPTION :
  --   Reads s-values and generates planes that are osculating a 
  --   rational normal curve.
  --   The system on return is the SAGBI homotopy.

    p : Poly;
    l : Poly_Sys(1..integer32(dim));
    s : Standard_Floating_Vectors.Vector(1..integer32(dim));
    mat : Standard_Floating_Matrices.Matrix(1..integer32(n),1..integer32(n-d));

  begin
   -- p := Lifted_Localized_Laplace_Expansion(n,d);
    p := Lifted_Localized_Laplace_Expansion(locmap);
    new_line;
    put("Give "); put(dim,1); put_line(" distinct real values for s : ");
    for i in s'range loop
      Read_Double_Float(s(i));
    end loop;
    new_line(file);
    put_line(file,"The selected s-values : ");
    put_line(file,s);
    for i in l'range loop
      mat := Orthogonal_Basis(n,n-d,s(i));
      Matrix_Homotopies.Add_Target(natural32(i),Convert(mat));
      l(i) := Intersection_Condition(mat,p);
    end loop;
    return l;
  end Given_Osculating_SAGBI_Homotopy;

  function Read_Input_SAGBI_Homotopy
             ( file : file_type; n,d,dim : integer32;
               locmap : Standard_Natural_Matrices.Matrix ) return Poly_Sys is

  -- DESCRIPTION :
  --   Reads the input planes from file.
  --   The system on return is the SAGBI homotopy.

    planesfile : file_type;
    p : Poly;
    L : Poly_Sys(1..dim);
    mat : Standard_Floating_Matrices.Matrix(1..n,1..n-d);

  begin
   -- p := Lifted_Localized_Laplace_Expansion(n,d);
    p := Lifted_Localized_Laplace_Expansion(locmap);
    new_line;
    put_line("Reading the name of the file with the input planes.");
    Read_Name_and_Open_File(planesfile);
    new_line(file);
    put_line(file,"The input planes : ");
    for i in L'range loop
      get(planesfile,mat);
      new_line(file); put(file,mat);
      mat := Orthogonalize(mat);
      Matrix_Homotopies.Add_Target(natural32(i),Convert(mat));
      L(i) := Intersection_Condition(mat,p);
    end loop;
    Close(planesfile);
    return L;
  end Read_Input_SAGBI_Homotopy;

  function Complex_Random_SAGBI_Homotopy
             ( n,d,dim : integer32;
               locmap : Standard_Natural_Matrices.Matrix ) return Poly_Sys is

  -- DESCRIPTION :
  --   Generates a SAGBI homotopy by taking random complex (n-d)-planes.
  --   This system will be the start system in the Cheater's homotopy.

    p : Poly;
    L : Poly_Sys(1..dim);
    mat : Standard_Complex_Matrices.Matrix(1..n,1..n-d);

  begin
   -- p := Lifted_Localized_Laplace_Expansion(n,d);
    p := Lifted_Localized_Laplace_Expansion(locmap);
    for i in L'range loop
      mat := Random_Orthogonal_Matrix(natural32(n),natural32(n-d));
      Matrix_Homotopies.Add_Start(natural32(i),mat);
      L(i) := Intersection_Condition(mat,p);
    end loop;
    return l;
  end Complex_Random_SAGBI_Homotopy;

  function Start_System ( L : Poly_Sys ) return Poly_Sys is

  -- DESCRIPTION :
  --   Returns the system L after evaluation for t = 0.
  --   This will be the start system in the SAGBI homotopy, flat deformation.

    r : Poly_Sys(L'range);
    m : constant integer32 := integer32(Number_of_Unknowns(L(L'first)));

  begin
    for i in L'range loop
      r(i) := Eval(L(i),Create(0.0),m);
    end loop;
    return r;
  end Start_System; 

  function Target_System ( L : Poly_Sys ) return Poly_Sys is

  -- DESCRIPTION :
  --   Returns the system L after evaluation for t = 1.
  --   This will be the target system in the SAGBI homotopy, flat deformation.

    r : Poly_Sys(L'range);
    m : constant integer32 := integer32(Number_of_Unknowns(L(L'first)));

  begin
    for i in L'range loop
      r(i) := Eval(L(i),Create(1.0),m);
    end loop;
    return r;
  end Target_System;

  procedure Polyhedral_Homotopy_Continuation
               ( file : in file_type;
                 p : in Poly_Sys; sols : in out Solution_List;
                 t : in Triangulation; lifted : in List; rep : in boolean ) is

  -- DESCRIPTION :
  --   This is the first continuation stage, the resolution of the start system
  --   in the SAGBI homotopy by means of polyhedral homotopy continuation.

    use Standard_Complex_Laur_JacoMats;

    timer : Timing_Widget;
    lifted_lq,lq : Laur_Sys(p'range);
    mix : constant Standard_Integer_Vectors.Vector(1..1) := (1..1 => p'last);
    lif : constant Array_of_Lists(1..1) := (1..1 => lifted);
    h : Eval_Coeff_Laur_Sys(p'range);
    c : Standard_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    m : Mult_Factors(j'range(1),j'range(2));
    mixsub : Mixed_Subdivision;

  begin
    tstart(timer);
    mixsub := Shallow_Create(p'last,t);
    lq := Polynomial_to_Laurent_System(p);
    lifted_lq := Perform_Lifting(p'last,mix,lif,lq);
    h := Create(lq);
    for i in c'range loop
      declare
        coeff_lq : constant Standard_Complex_Vectors.Vector := Coeff(lq(i));
      begin
        c(i) := new Standard_Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          c(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
    e := Create(lq);
    Create(lq,j,m);
    if rep
     then Mixed_Solve(file,lifted_lq,lif,h,c,e,j,m,mix,mixsub,sols);
     else Mixed_Solve(lifted_lq,lif,h,c,e,j,m,mix,mixsub,sols);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Polyhedral Homotopy Continuation");
  end Polyhedral_Homotopy_Continuation;

  procedure Solve_Start_System
               ( file : in file_type;
                 p : in Poly_Sys; sols : in out Solution_List;
                 report : in boolean ) is

  -- DESCRIPTION :
  --   Computes the volume of the support of p.  This is the
  --   combinatorial-geometric set up for the first continuation stage.

    timer : Timing_Widget;
    support : constant List := Create(p(p'first));
    lifted,lifted_last : List;
    t : Triangulation;
    vol : natural32;

  begin
    new_line(file);
    put_line(file,"The support of the start system : ");
    new_line(file);
    put(file,support);
    tstart(timer);
    Dynamic_Lifting(support,false,false,0,lifted,lifted_last,t);
    new_line(file);
    put_line(file,"The lifted support : ");
    new_line(file);
    put(file,lifted);
    vol := Volume(t);
    tstop(timer);
    new_line(file);
    put(file,"The volume : "); put(file,vol,1); new_line(file);
    new_line(file);
    print_times(file,timer,"Triangulation and Volume Computation");
    Polyhedral_Homotopy_Continuation(file,p,sols,t,lifted,report);
    Clear(t);
  end Solve_Start_System;

  procedure Flat_Deformation
               ( file : in file_type; sh : in Poly_Sys;
                 sols : in out Solution_List; report : in boolean ) is

  -- DESCRIPTION :
  --   Performs flat deformation as defined by the SAGBI homotopy.
  --   This is the second stage in the continuation.

    use Standard_Complex_Jaco_Matrices;

    timer : Timing_Widget;
    sh_eval : Eval_Poly_Sys(sh'range);
    jac_mat : Jaco_Mat(sh'range,sh'first..sh'last+1);
    eva_jac : Eval_Jaco_Mat(jac_mat'range(1),jac_mat'range(2));

    function Eval ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                  return Standard_Complex_Vectors.Vector is

      xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      return Eval(sh_eval,xt);
    end Eval;

    function Diff ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                  return Standard_Complex_Matrices.Matrix is

      res : Standard_Complex_Matrices.Matrix(x'range,x'range);
      xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      for i in res'range(1) loop
        for j in res'range(2) loop
          res(i,j) := Eval(eva_jac(i,j),xt);
        end loop;
      end loop;
      return res;
    end Diff;

    function Diff ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                  return Standard_Complex_Vectors.Vector is

      res : Standard_Complex_Vectors.Vector(x'range);
      xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      for i in res'range loop
        res(i) := Eval(eva_jac(i,eva_jac'last(2)),xt);
      end loop;
      return res;
    end Diff;

    procedure Sil_Cont is new Silent_Continue(Max_Norm,Eval,Diff,Diff);
    procedure Rep_Cont is new Reporting_Continue(Max_Norm,Eval,Diff,Diff);

  begin
    tstart(timer);
    sh_eval := Create(sh);
    jac_mat := Create(sh);
    eva_jac := Create(jac_mat);
    Set_Continuation_Parameter(sols,Create(0.0));
    if report
     then Rep_Cont(file,sols,false,target=>Create(1.0));
     else Sil_Cont(sols,false,target=>Create(1.0));
    end if;
    tstop(timer);
    new_line(file); print_times(file,timer,"Flat Deformation");
    Clear(sh_eval);
    Clear(jac_mat);
    Clear(eva_jac);
  end Flat_Deformation;

  procedure Solve_Target_System
               ( file : in file_type; start,target : in Poly_Sys;
                 sols : in out Solution_List; report : in boolean ) is

  -- DESCRIPTION :
  --   This is the third and last continuation stage: Cheater's homotopy.
  --   It is implemented in the space of polynomials.

    timer : Timing_Widget;
    n : constant integer32 := target'last;
    a : constant Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    b : constant Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);

  begin
    tstart(timer);
    Standard_Homotopy.Create(target,start,2,a,b,true);       -- linear cheater
    Set_Continuation_Parameter(sols,Create(0.0));
    declare
      procedure Sil_Cont is
        new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                            Standard_Homotopy.Diff,Standard_Homotopy.Diff);
      procedure Rep_Cont is
        new Reporting_Continue(Max_Norm,Standard_Homotopy.Eval,
                               Standard_Homotopy.Diff,Standard_Homotopy.Diff);
    begin
      if report
       then Rep_Cont(file,sols,false,target=>Create(1.0));
       else Sil_Cont(sols,false,target=>Create(1.0));
      end if;
    end;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Cheater's homotopy to target system");
  end Solve_Target_System;

--  procedure Solve_Target_System
--               ( file : in file_type; dim : in natural;
--                 symtarget : in Poly; sols : in out Solution_List; 
--                 report : in boolean ) is
--
--  -- DESCRIPTION :
--  --   This is the third and last continuation stage: Cheater's homotopy.
--  --   It uses a coefficient-matrix homotopy and takes as input the target
--  --   system with coefficients as brackets.
--  --   This is far more expensive than the linear Cheater's homotopy.
--
--    use Standard_Complex_Jaco_Matrices;
--
--    timer : Timing_Widget;
--    h : constant Standard_Complex_Poly_Functions.Eval_Coeff_Poly
--      := Create(symtarget);
--    homsys : Poly_Sys(1..dim);
--    jacmat : Eval_Coeff_Jaco_Mat(1..dim,1..dim);
--    mulfac : Mult_Factors(1..dim,1..dim);
--    brkcff : Standard_Complex_Vectors.Vector := Coeff(symtarget);
--    syscff : Standard_Complex_VecVecs.VecVec(1..dim);
--
--  begin
--    tstart(timer);
--    Set_Continuation_Parameter(sols,Create(0.0));
--    for i in homsys'range loop
--      Copy(symtarget,homsys(i));
--    end loop;
--    Create(homsys,jacmat,mulfac);
--    for i in syscff'range loop
--      syscff(i) := new Standard_Complex_Vectors.Vector(brkcff'range);
--    end loop;
--    declare
--
--      function Eval ( x : Standard_Complex_Vectors.Vector;
--                      t : Complex_Number )
--                    return Standard_Complex_Vectors.Vector is
--
--        y : Standard_Complex_Vectors.Vector(x'range);
--        c : Standard_Complex_Vectors.Vector(brkcff'range);
--
--      begin
--        for i in y'range loop
--          c := Intersection_Coefficients(Matrix_Homotopies.Eval(i,t),brkcff);
--          y(i) := Eval(h,c,x);
--        end loop;
--        return y;
--      end Eval;
--
--      function dHt ( x : Standard_Complex_Vectors.Vector;
--                     t : Complex_Number ) 
--                   return Standard_Complex_Vectors.Vector is 
--
--        y : Standard_Complex_Vectors.Vector(x'range) := x;
--
--      begin
--        return y;
--      end dHt;
--
--      function dHx ( x : Standard_Complex_Vectors.Vector;
--                     t : Complex_Number ) 
--                   return Standard_Complex_Matrices.Matrix is 
--      begin
--        for i in syscff'range loop
--          syscff(i).all
--            := Intersection_Coefficients(Matrix_Homotopies.Eval(i,t),brkcff);
--        end loop;
--        return Eval(jacmat,mulfac,syscff,x);
--      end dHx;
--
--      procedure Sil_Cont is new Silent_Continue(Max_Norm,Eval,dHt,dHx);
--      procedure Rep_Cont is new Reporting_Continue(Max_Norm,Eval,dHt,dHx);
--
--    begin
--      if report
--       then Rep_Cont(file,sols,false,Create(1.0));
--       else Sil_Cont(sols,false,Create(1.0));
--      end if;
--    end;
--    Standard_Complex_VecVecs.Clear(syscff);
--    tstop(timer);
--    new_line(file);
--    print_times(file,timer,"Cheater's homotopy to target system");
--  end Solve_Target_System;

  procedure Refine_Roots ( file : in file_type;
                           p : in Poly_Sys; sols : in out Solution_List ) is

    epsxa : constant double_float := 1.0E-8;
    epsfa : constant double_float := 1.0E-8;
    tolsing : constant double_float := 1.0E-8;
    max : constant natural32 := 3;
    numit : natural32 := 0;
    deflate : boolean := false;

  begin
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,max,deflate,false);
  end Refine_Roots;

  procedure Set_Parameters ( file : in file_type; report : out boolean ) is

  -- DESCRIPTION :
  --   Interactive determination of the continuation and output parameters.

    oc : natural32;

  begin
    new_line;
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    report := not (oc = 0);
    new_line;
    put_line("See the output file for results.");
    new_line;
  end Set_Parameters;

  function t_Symbol return Symbol is

  -- DESCRIPTION :
  --   Returns the symbol to represent the continuation parameter t.

    res : Symbol;
  
  begin
    res(1) := 't';
    for i in 2..res'last loop
      res(i) := ' ';
    end loop;
    return res;
  end t_Symbol;

  function Lowest_Degree_Localization_Map
             ( n,d : integer32) return Standard_Natural_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the lowest-degree localization map.

    res : Standard_Natural_Matrices.Matrix(1..n,1..d);

  begin
    for i in 1..n-d loop
      for j in 1..d loop
        res(i,j) := 2;
      end loop;
    end loop;
    for i in n-d+1..n loop
      for j in 1..d loop
        if i = j+n-d
         then res(i,j) := 1;
         else res(i,j) := 0;
        end if;
      end loop;
    end loop;
    return res;
  end Lowest_Degree_Localization_Map;

  function Read_Localization_Map
             ( n,d : integer32 ) return Standard_Natural_Matrices.Matrix is

  -- DESCRIPTION :
  --   Reads a localization map from standard input.

    res : Standard_Natural_Matrices.Matrix(1..n,1..d);

  begin
    put_line("Reading localization map: 0,1 for I and 2 = free element.");
    put_line("The lowest-degree localization map looks like :");
    put(Lowest_Degree_Localization_Map(n,d));
    put_line("and the general localization map is :");
    put(Localization_Map(natural32(n),natural32(d)));
    put("Give a "); put(n,1); put("-by-"); put(d,1);
    put_line(" matrix to represent your own map :");
    get(res); skip_line;
    return res;
  end Read_Localization_Map;

  procedure Main ( file : in file_type; n,d : in natural32 ) is

  -- DESCRIPTION :
  --   There are four parts in the elaboration of the SAGBI homotopy :
  --     0. Set up of the polynomials in the SAGBI homotopy;
  --     1. Solve the start system by polyhedral continuation;
  --     2. Apply the flat deformations to solve complex instance;
  --     3. Cheater's homotopy from complex to real instance.

    dim : constant natural32 := (n-d)*d;
    timer,totaltimer : Timing_Widget;
    realsagbih,compsagbih,realtarget,comptarget,starts
      : Poly_Sys(1..integer32(dim));
    sols : Solution_List;
    report : boolean;
    inputchoice,localchoice : character;
    locmap : Standard_Natural_Matrices.Matrix(1..integer32(n),1..integer32(d));

  begin
    put(file,"SAGBI Homotopies for m = "); put(file,n-d,1);
    put(file," and p = "); put(file,d,1); new_line(file);
    new_line;
    put_line("MENU for a localization for the output planes : ");
    put_line("  1. Lowest-degree map, with I in lower-right corner.");
    put_line("  2. General map, able to represent any intersection.");
    put_line("  3. Give in your own localization map.");
    put("Type 1, 2, or 3 to select : "); Ask_Alternative(localchoice,"123");
    case localchoice is
      when '1' =>
        locmap := Lowest_Degree_Localization_Map(integer32(n),integer32(d));
      when '2' => locmap := Localization_Map(n,d);
      when '3' => locmap := Read_Localization_Map(integer32(n),integer32(d));
      when others => null;
    end case;
    new_line;
    put_line("MENU for constructing target system based on input planes : ");
    put_line("  1. Generate input planes osculating a rational normal curve.");
    put_line("  2. Interactively give s-values for the "
                               & "osculating input planes.");
    put_line("  3. Give the name of the file with input planes.");
    put("Type 1, 2, or 3 to select : "); Ask_Alternative(inputchoice,"123");
    tstart(totaltimer);
    tstart(timer);
    Matrix_Indeterminates.Initialize_Symbols(n,d);
    Matrix_Indeterminates.Reduce_Symbols(locmap);
    Matrix_Homotopies.Init(dim);
    case inputchoice is
      when '1' => 
        realsagbih := Random_Osculating_SAGBI_Homotopy(file,n,d,dim,locmap);
      when '2' => 
        realsagbih := Given_Osculating_SAGBI_Homotopy(file,n,d,dim,locmap);
      when '3' =>
        realsagbih := Read_Input_SAGBI_Homotopy
                        (file,integer32(n),integer32(d),integer32(dim),locmap);
      when others => null;
    end case;
    realtarget := Target_System(realsagbih);
    compsagbih := Complex_Random_SAGBI_Homotopy
                    (integer32(n),integer32(d),integer32(dim),locmap);
    comptarget := Target_System(compsagbih);
    starts := Start_System(compsagbih);
    Symbol_Table.Add(t_Symbol);
    new_line(file);
    put_line(file,"The target polynomial system with real coefficients : ");
    new_line(file);
    put(file,natural32(realtarget'last),realtarget);
    new_line(file);
    put_line(file,"The localization map : "); put(file,locmap);
    new_line(file);
    Matrix_Homotopies_io.Write(file);
    put_line(file,"The target polynomial system with complex coefficients : ");
    put_line(file,comptarget);
    new_line(file);
    put_line(file,"The SAGBI Homotopy as complex polynomial system : ");
    new_line(file);
    put_line(file,compsagbih);
    new_line(file);
    put_line(file,"The SAGBI Homotopy at t=0, the start system : ");
    new_line(file);
    put_line(file,starts);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Setting up the polynomials in the SAGBI homotopy");
    Set_Parameters(file,report);
    Solve_Start_System(file,starts,sols,report);
    Flat_Deformation(file,compsagbih,sols,report);
    Solve_Target_System(file,comptarget,realtarget,sols,report);
   -- Solve_Target_System(file,dim,symtarget,sols,report); 
   --  this is the determinantal cheater's homotopy, very expensive
    Refine_Roots(file,realtarget,sols);
    Matrix_Indeterminates.Clear_Symbols;
    tstop(totaltimer);
    new_line(file);
    print_times(file,totaltimer,"Total time for Solving the System");
  end Main;

end Main_SAGBI_Homotopies;
