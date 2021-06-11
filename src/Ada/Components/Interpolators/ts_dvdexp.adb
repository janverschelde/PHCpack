with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Multprec_Complex_Numbers;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_VecVecs;
with Multprec_Complex_VecVecs;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Symbol_Table;                       use Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Total_Degree_Start_Systems;         use Total_Degree_Start_Systems;
with Standard_Homotopy;
with Standard_Complex_Solutions;
with Continuation_Parameters;
with Standard_IncFix_Continuation;       use Standard_IncFix_Continuation;
with Standard_Polynomial_Interpolators;
with Multprec_Polynomial_Interpolators;
with Witness_Sets;                       use Witness_Sets;
with Sampling_Machine;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with Rectangular_Sample_Grids;           use Rectangular_Sample_Grids;
with Standard_Divided_Differences;
with Multprec_Divided_Differences;

procedure ts_dvdexp is

-- DESCRIPTION :
--   This procedure implements a test experiment on divided differences:
--   interpolate a dense polynomial of degree d with coefficients equal
--   to one with the direct method and with divided differences.

  procedure Add_Symbols is

  -- DESCRIPTION :
  --   Initializes the symbol table and adds "x" and "y" to it.

    xsb,ysb,zsb : Symbol;

  begin
    Symbol_Table.Init(3);
    xsb := (xsb'range => ' ');
    xsb(1) := 'x';
    Symbol_Table.Add(xsb);
    ysb := (ysb'range => ' ');
    ysb(1) := 'y';
    Symbol_Table.Add(ysb);
    zsb := (zsb'range => ' ');
    zsb(1) := 'z';
    Symbol_Table.Add(zsb);
  end Add_Symbols;

  function Is_Zero ( v : Standard_Natural_Vectors.Vector ) return boolean is
  begin
    for i in v'range loop
      if v(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  procedure Set_Constant_to_One 
              ( p : in out Standard_Complex_Polynomials.Poly ) is

    use Standard_Complex_Numbers,Standard_Complex_Polynomials;

    procedure Change_Coefficient ( t : in out Term; continue : out boolean ) is

    begin
      if Is_Zero(t.dg.all)
       then t.cf := Create(1.0);
            continue := false;
       else continue := true;
      end if;
    end Change_Coefficient;
    procedure Change_Constant is new Changing_Iterator(Change_Coefficient);

  begin
    Change_Constant(p);
  end Set_Constant_to_One;

  procedure Generate_Test_Standard_Polynomial
              ( file : in file_type; d,cff : in natural32;
                p : out Standard_Complex_Polynomials.Poly;
                m : out natural32 ) is

  -- DESCRIPTION :
  --   Generates a test polynomial and writes this on file.

  -- ON ENTRY :
  --   file     for results and diagnostics;
  --   d        degree of the test polynomial.

  -- ON RETURN :
  --   p        dense polynomial with all coefficients equal to one;
  --   m        number of monomials in p.

  begin
    put(file,"Generating dense polynomial of degree "); put(file,d,1);
    p := Standard_Polynomial_Interpolators.Create(d,2,cff);
    if cff /= 1
     then Set_Constant_to_One(p);
    end if;
   -- new_line(file);
   -- put_line(file,"The test polynomial : "); put(file,p);
   -- new_line(file);
    m := Standard_Complex_Polynomials.Number_of_Terms(p);
    put(file," with #terms : "); put(file,m,1); new_line(file);
  end Generate_Test_Standard_Polynomial;

  procedure Add_Slice
                ( p : in Standard_Complex_Polynomials.Poly;
                  sys : out Standard_Complex_Poly_Systems.Poly_Sys;
                  hypsli : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Generates a random hyperplane and adds this as slice in the
  --   system on return.  The coefficients of the slice are in hypsli.

  -- ON ENTRY :
  --   p          polynomial in two variables;

  -- ON RETURN :
  --   sys        polynomial system, sys(1) = p, sys(2) = hypsli;
  --   hypsli     random hyperplane through the origin.

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    hyp : Poly;
    slice : Standard_Complex_Vectors.Vector(0..2) := Random_Vector(0,2);
    t : Term;

  begin
    slice(0) := Create(0.0);
    t.dg := new Standard_Natural_Vectors.Vector'(1..2 => 0);
    t.cf := slice(0);
    Add(hyp,t);
    t.dg(1) := 1; t.cf := slice(1);
    Add(hyp,t);
    t.dg(1) := 0;
    t.dg(2) := 1; t.cf := slice(2);
    Add(hyp,t);
    Clear(t);
    sys(1) := p;
    sys(2) := hyp;
    hypsli := slice;
  end Add_Slice;

  procedure Write_Summary
              -- ( file : in file_type;
               ( sols : in Standard_Complex_Solutions.Solution_List;
                 maxerr,minrco,maxres : out double_float ) is

  -- DESCRIPTION :
  --   Writes the (err,rcond,res)-banner for every solution in the list,
  --   returns maximal errors, minimal rcond and maximal residual.

    use Standard_Complex_Solutions;
    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Link_to_Solution;

  begin
    maxerr := 0.0; minrco := 1.0E+10; maxres := 0.0;
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
     -- put(file,i,3); put(file," : ");
     -- put(file,"  err : "); put(file,ls.err,3);
     -- put(file,"  rco : "); put(file,ls.rco,3);
     -- put(file,"  res : "); put(file,ls.res,3);
     -- new_line(file);
      if ls.err > maxerr then maxerr := ls.err; end if;
      if ls.rco < minrco then minrco := ls.rco; end if;
      if ls.res > maxres then maxres := ls.res; end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Summary;

  procedure Generic_Points
               ( file : in file_type;
                 p : in Standard_Complex_Polynomials.Poly;
                 sys : out Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : out Standard_Complex_Solutions.Solution_List;
                 hyp : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes the generic points on the polynomial p.

  -- ON ENTRY :
  --   file      for results and diagnostics;
  --   p         polynomial in two variables;

  -- ON RETURN :
  --   sys       polynomial system, sys(1) = p, sys(2) = hyp;
  --   sols      solutions to the system;
  --   hyp       random hyperplane through the origin.

    use Standard_Complex_Numbers;

    q : Standard_Complex_Poly_Systems.Poly_Sys(1..2);
    epsxa : constant double_float := 1.0E-13;
    epsfa : constant double_float := 1.0E-13;
    tolsing : constant double_float := 1.0E-8;
    numit : natural32 := 0;
    max : constant natural32 := 5;
    maxerr,minrco,maxres : double_float;
    deflate : boolean := false;

    procedure Continue is
      new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                          Standard_Homotopy.Diff,Standard_Homotopy.Diff);
  begin
    Add_Slice(p,sys,hyp);
    Start_System(sys,q,sols);
    Standard_Homotopy.Create(sys,q,2,Random1);
    Continuation_Parameters.Tune(2);
    Continue(sols,false); -- ,Create(1.0));
    Silent_Root_Refiner(sys,sols,epsxa,epsfa,tolsing,numit,max,deflate);
   -- Write_Summary(file,sols,maxerr,minrco,maxres);
    Write_Summary(sols,maxerr,minrco,maxres);
    put(file,"GP Max err :"); put(file,maxerr,3);
    put(file,"  Min rco :"); put(file,minrco,3);
    put(file,"  Max res :"); put(file,maxres,3); new_line(file);
    Standard_Complex_Poly_Systems.Clear(q);
    Standard_Homotopy.Clear;
  end Generic_Points;

  procedure Create_Standard_Sample_Grid
               ( file : in file_type;
                 sys : Standard_Complex_Poly_Systems.Poly_Sys;
                 hyp : in Standard_Complex_Vectors.Vector;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 grid : out Array_of_Standard_Sample_Lists ) is

  -- DESCRIPTION :
  --   Creates a grid of samples.

  -- ON ENTRY :
  --   file      for results and diagnostics;
  --   sys       polynomial system, sys(1) = p, sys(2) = hyp;
  --   sols      solutions to the system;
  --   hyp       random hyperplane through origin.
 
  -- ON RETURN :
  --   grid      grid of samples.

    lhyp : constant Standard_Complex_Vectors.Link_to_Vector
         := new Standard_Complex_Vectors.Vector'(hyp);
    sli : constant Standard_Complex_VecVecs.VecVec(1..1) := (1..1 => lhyp);
    sps : constant Standard_Sample_List := Create(sols,sli);

  begin
    Sampling_Machine.Initialize(sys);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    grid := Create1(sps,natural32(grid'last));
    Sampling_Machine.Clear;
  end Create_Standard_Sample_Grid;

  procedure Create_Multprec_Sample_Grid
               ( file : in file_type; size : in natural32;
                 sys : Standard_Complex_Poly_Systems.Poly_Sys;
                 embsys : Standard_Complex_Poly_Systems.Poly_Sys;
                 hyp : in Standard_Complex_Vectors.Vector;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 grid : out Array_of_Multprec_Sample_Lists ) is

  -- DESCRIPTION :
  --   Creates a grid of samples.

  -- ON ENTRY :
  --   file      for results and diagnostics;
  --   sys       polynomial system, sys(1) = p, sys(2) = hyp;
  --   sols      solutions to the system;
  --   hyp       random hyperplane through origin.
 
  -- ON RETURN :
  --   grid      grid of samples.

    sli : constant Standard_Complex_VecVecs.VecVec := Slices(embsys,1);
    sps : constant Standard_Sample_List := Create(sols,sli);
    mpsys : Multprec_Complex_Poly_Systems.Poly_Sys(1..2);

  begin
    mpsys(1) := Convert(sys(1));
    Sampling_Machine.Initialize(embsys,mpsys,1,size);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    grid := Create1(sps,natural32(grid'last),size);
    Sampling_Machine.Clear;
  end Create_Multprec_Sample_Grid;

  procedure Standard_Newton_Interpolate
               ( file : in file_type;
                 grid : in Array_of_Standard_Sample_Lists;
                 exp : out Standard_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Computes the Newton from of the interpolator through p. 

  -- ON ENTRY :
  --   file      for results and diagnostics;
  --   grid      grid of sample points.

  -- ON RETRUN :
  --   exp       the expanded Newton form of the interpolator.

    use Standard_Complex_Numbers;
    use Standard_Divided_Differences;

    q : Newton_Interpolator1;
    eq : Newton_Form_Evaluator1;
    maxres : double_float;

  begin
    q := Create(grid,Create(0.0));
    eq := Create(q);
    exp := Expand(eq);
    maxres := Maximal_Error(q,grid);
    put(file,"Maximal residual of Newton form on grid :");
    put(file,maxres,3);
    put(file,"  |LOG10| : "); put(file,AbsVal(LOG10(maxres)),2,3,0);
    new_line(file);
   -- put_line(file,"The expanded Newton form :"); put_line(file,exp);
    Clear(q); Clear(eq);
  end Standard_Newton_Interpolate;

  procedure Multprec_Newton_Interpolate
               ( file : in file_type;
                 grid : in Array_of_Multprec_Sample_Lists;
                 exp : out Multprec_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Computes the Newton from of the interpolator through p.

  -- ON ENTRY :
  --   file      for results and diagnostics;
  --   grid      grid of sample points.

  -- ON RETRUN :
  --   exp       the expanded Newton form of the interpolator.

    use Multprec_Complex_Numbers;
    use Multprec_Divided_Differences;

    q : Newton_Interpolator1;
    eq : Newton_Form_Evaluator1;
    tmp : Floating_Number;
    maxres : double_float;

  begin
    q := Create(grid,Create(integer32(0)));
    eq := Create(q);
    exp := Expand(eq);
    tmp := Maximal_Error(q,grid); maxres := Round(tmp); Clear(tmp);
    put(file,"Maximal residual of Newton form on grid :");
    put(file,maxres,3);
    put(file,"  |LOG10| : "); put(file,AbsVal(LOG10(maxres)),2,3,0);
    new_line(file);
   -- put_line(file,"The expanded Newton form :"); put_line(file,exp);
    Clear(q); Clear(eq);
  end Multprec_Newton_Interpolate;

  procedure Standard_Direct_Interpolate
               ( file : in file_type; 
                 p : in Standard_Complex_Polynomials.Poly;
                 grid : in Array_of_Standard_Sample_Lists;
                 ip : out Standard_Complex_Polynomials.Poly;
                 rcond : out double_float ) is

  -- DESCRIPTION :
  --   Computes the interpolating polynomial with the direct approach. 

  -- ON ENTRY :
  --   file      for results and diagnostics;
  --   grid      grid of sample points.

  -- ON RETRUN :
  --   ip        interpolating polynomial.

    use Standard_Polynomial_Interpolators;
    samples : constant Standard_Complex_VecVecs.VecVec
            := Extract_Samples(grid);

  begin
    Interpolate(p,samples,ip,rcond);
  end Standard_Direct_Interpolate;

  procedure Multprec_Direct_Interpolate
               ( file : in file_type; 
                 p : in Multprec_Complex_Polynomials.Poly;
                 grid : in Array_of_Multprec_Sample_Lists;
                 ip : out Multprec_Complex_Polynomials.Poly;
                 rcond : out double_float ) is

  -- DESCRIPTION :
  --   Computes the interpolating polynomial with the direct approach. 

  -- ON ENTRY :
  --   file      for results and diagnostics;
  --   grid      grid of sample points.

  -- ON RETRUN :
  --   ip        interpolating polynomial.

    use Multprec_Polynomial_Interpolators;
    samples : constant Multprec_Complex_VecVecs.VecVec
            := Extract_Samples(grid);
    invcond : Floating_Number;

  begin
    Interpolate1(p,samples,ip,invcond);
    rcond := Round(invcond); Clear(invcond);
  end Multprec_Direct_Interpolate;

  procedure Standard_Experiment
               ( file : in file_type; d : in integer32;
                 cff : in natural32; restab : out Matrix ) is

  -- DESCRIPTION :
  --   Performs the comparison Newton versus direct with standard arithmetic.

    dense,nip,dip : Standard_Complex_Polynomials.Poly;
    nbmon : natural32;
    sys : Standard_Complex_Poly_Systems.Poly_Sys(1..2);
    sols : Standard_Complex_Solutions.Solution_List;
    hyp : Standard_Complex_Vectors.Vector(0..2);
    maxerr,mindst,disnip,disdip,rcond : double_float;
    grid : Array_of_Standard_Sample_Lists(0..d);

  begin
    Add_Symbols;
    Generate_Test_Standard_Polynomial(file,natural32(d),cff,dense,nbmon);
    restab(d,1) := double_float(nbmon);
    Generic_Points(file,dense,sys,sols,hyp);
    Create_Standard_Sample_Grid(file,sys,hyp,sols,grid);
    maxerr := Maximal_Error(grid);
    mindst := Minimal_Distance(grid);
    put(file,"Maximal error on the samples in grid    :");
    put(file,maxerr,3);
    put(file,"  |LOG10| : "); put(file,AbsVal(LOG10(maxerr)),2,3,0);
    new_line(file);
    put(file,"Minimal distance between the samples    :");
    put(file,mindst,3);
    put(file,"  |LOG10| : "); put(file,AbsVal(LOG10(mindst)),2,3,0);
    new_line(file);
    Standard_Newton_Interpolate(file,grid,nip);
    Standard_Direct_Interpolate(file,dense,grid,dip,rcond);
    Deep_Clear(grid);
    put(file,"Estimate of inverse condition number    :");
    put(file,rcond,3);
    restab(d,2) := AbsVal(LOG10(rcond));
    put(file,"  |LOG10| : "); put(file,restab(d,2),2,3,0);
    new_line(file);
    disnip := Standard_Polynomial_Interpolators.Distance(dense,nip);
    disdip := Standard_Polynomial_Interpolators.Distance(dense,dip);
    put(file,"Distance expanded Newton form and exact :");
    put(file,disnip,3);
    restab(d,3) := AbsVal(LOG10(disnip));
    put(file,"  |LOG10| : "); put(file,restab(d,3),2,3,0);
    new_line(file);
    put(file,"Distance direct interpolator and exact  :");
    put(file,disdip,3);
    restab(d,4) := AbsVal(LOG10(disdip));
    put(file,"  |LOG10| : "); put(file,restab(d,4),2,3,0);
    new_line(file);
    new_line(file);
    Standard_Complex_Polynomials.Clear(dense);
    Standard_Complex_Polynomials.Clear(nip);
    Standard_Complex_Polynomials.Clear(dip);
  end Standard_Experiment;

  procedure Embed ( sys : in Standard_Complex_Poly_Systems.Poly_Sys;
                    sols : in Standard_Complex_Solutions.Solution_List;
                    embsys : out Standard_Complex_Poly_Systems.Poly_Sys;
                    embsols : out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Constructs an embedding for the given system and solution list.

    use Standard_Complex_Numbers;
    sysp1 : Standard_Complex_Poly_Systems.Poly_Sys(1..3);

  begin
    Standard_Complex_Polynomials.Copy(sys(1),sysp1(1));
    Standard_Complex_Polynomials.Copy(sys(2),sysp1(3));
    embsys := Embed(sysp1,Random_Vector(1,embsys'last));
    embsols := Add_Component(sols,Create(0.0));
  end Embed;

  procedure Multprec_Experiment
               ( file : in file_type; d : in integer32;
                 cff,size : in natural32; restab : out Matrix ) is

  -- DESCRIPTION :
  --   Performs the comparison Newton versus direct with multi-precision
  --   arithmetic.

    stdense : Standard_Complex_Polynomials.Poly;
    mpdense,nip,dip : Multprec_Complex_Polynomials.Poly;
    nbmon : natural32;
    sys : Standard_Complex_Poly_Systems.Poly_Sys(1..2);
    embsys : Standard_Complex_Poly_Systems.Poly_Sys(1..3);
    sols,embsols : Standard_Complex_Solutions.Solution_List;
    hyp : Standard_Complex_Vectors.Vector(0..2);
    tmp : Floating_Number;
    maxerr,mindst,disnip,disdip,rcond : double_float;
    grid : Array_of_Multprec_Sample_Lists(0..d);

  begin
    Add_Symbols;
    Generate_Test_Standard_Polynomial(file,natural32(d),cff,stdense,nbmon);
    restab(d,1) := double_float(nbmon);
    Generic_Points(file,stdense,sys,sols,hyp);
    Embed(sys,sols,embsys,embsols);
   -- put_line(file,"The embedded polynomial system :");
   -- put_line(file,embsys);
    mpdense := Convert(stdense);
    Set_Size(mpdense,size);
   -- put_line(file,"The multi-precision test polynomial : ");
   -- put_line(file,mpdense);
    Create_Multprec_Sample_Grid(file,size,sys,embsys,hyp,embsols,grid);
    tmp := Maximal_Error(grid); maxerr := Round(tmp); Clear(tmp);
    tmp := Minimal_Distance(grid); mindst := Round(tmp); Clear(tmp);
    put(file,"Maximal error on the samples in grid    :");    
    put(file,maxerr,3);
    put(file,"  |LOG10| : "); put(file,AbsVal(LOG10(maxerr)),2,3,0);
    new_line(file);
    put(file,"Minimal distance between the samples    :");    
    put(file,mindst,3);
    put(file,"  |LOG10| : "); put(file,AbsVal(LOG10(mindst)),2,3,0);
    new_line(file);
    Multprec_Newton_Interpolate(file,grid,nip);
    Multprec_Direct_Interpolate(file,mpdense,grid,dip,rcond);
    Deep_Clear(grid);
    put(file,"Estimate of inverse condition number    :");
    put(file,rcond,3);
    restab(d,2) := AbsVal(LOG10(rcond));
    put(file,"  |LOG10| : "); put(file,restab(d,2),2,3,0);
    new_line(file);
    tmp := Multprec_Polynomial_Interpolators.Distance(mpdense,nip);
    disnip := Round(tmp); Clear(tmp);
    put(file,"Distance expanded Newton form and exact :");
    put(file,disnip,3);
    restab(d,3) := AbsVal(LOG10(disnip));
    put(file,"  |LOG10| : "); put(file,restab(d,3),2,3,0);
    new_line(file);
    tmp := Multprec_Polynomial_Interpolators.Distance(mpdense,dip);
    disdip := Round(tmp); Clear(tmp);
    put(file,"Distance direct interpolator and exact  :");
    put(file,disdip,3);
    restab(d,4) := AbsVal(LOG10(disdip));
    put(file,"  |LOG10| : "); put(file,restab(d,4),2,3,0);
    new_line(file);
    new_line(file);
    Standard_Complex_Polynomials.Clear(stdense);
    Multprec_Complex_Polynomials.Clear(mpdense);
    Multprec_Complex_Polynomials.Clear(nip);
    Multprec_Complex_Polynomials.Clear(dip);
  end Multprec_Experiment;

  procedure Write_Results ( file : in file_type; restab : in Matrix ) is
  begin
    put_line(file,
      "--------------------------------------------------------------");
    put(file,"  d     m ");
    put(file," |LOG10(rcond)| ");
    put(file," |LOG10(NewDf)| ");
    put(file," |LOG10(DirDf)| "); new_line(file);
    put_line(file,
      "--------------------------------------------------------------");
    for i in restab'range(1) loop
      put(file,i,3); put(file,"  ");
      put(file,integer32(restab(i,1)),4); put(file,"      ");
      put(file,restab(i,2),2,3,0);      put(file,"          ");
      put(file,restab(i,3),2,3,0);      put(file,"          ");
      put(file,restab(i,4),2,3,0); new_line(file);
    end loop;
    put_line(file,
      "--------------------------------------------------------------");
  end Write_Results;

  procedure Main is

    file : file_type;
    cff,deci,size : natural32;
    d1,d2 : integer32;

  begin
    new_line;
    put_line("Interpolation experiment on polynomial in two variables.");
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_file(file);
    new_line;
    put_line("A series of polynomials of various degrees will be generated.");
    put("Give lower degree in range : "); get(d1);
    put("Give upper degree in range : "); get(d2);
    new_line;
    put_line("MENU for choice of coefficients :");
    put_line("  0. random complex coefficients of modulus one;");
    put_line("  1. all coefficients are equal to one;");
    put_line("  2. random real coefficients will be generated.");
    put("Give a natural number (0,1, or 2) : "); get(cff);
    declare
      restab : Matrix(d1..d2,1..4);
    begin
      new_line;
      put("Give the number of decimal places (<= 16 is standard) : ");
      get(deci);
      if deci > 16 then
        size := Decimal_to_Size(deci);
        Sampling_Machine.Interactive_Tune_Refiner(size);
        new_line;
        put_line("See the output file for results...");
        new_line;
        for d in d1..d2 loop
          Multprec_Experiment(file,d,cff,size,restab);
        end loop;
      else
        new_line;
        put_line("See the output file for results...");
        new_line;
        for d in d1..d2 loop
          Standard_Experiment(file,d,cff,restab);
        end loop;
      end if;
      Write_Results(file,restab);
    end;
  end Main;

begin
  Main;
end ts_dvdexp;
