with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_VecVecs;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Symbol_Table;                       use Symbol_Table;
with Standard_Natural_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Homotopy;
with Standard_Scaling;                   use Standard_Scaling;
with Total_Degree_Start_Systems;         use Total_Degree_Start_Systems;
with Continuation_Parameters;
with Standard_IncFix_Continuation;       use Standard_IncFix_Continuation;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Standard_Polynomial_Interpolators;
with Sampling_Machine;
with Sample_Points;                      use Sample_Points;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with Irreducible_Components;             use Irreducible_Components;
with Irreducible_Component_Lists;        use Irreducible_Component_Lists;
with Irreducible_Component_Lists_io;     use Irreducible_Component_Lists_io;
with Rectangular_Sample_Grids;           use Rectangular_Sample_Grids;
with Standard_Divided_Differences;       use Standard_Divided_Differences;
with Drivers_to_Component_Creators;      use Drivers_to_Component_Creators;

procedure Sensitivity_of_Factorization is

  procedure Add_Symbols is

  -- DESCRIPTION :
  --   Initializes the symbol table and adds "x" and "y" to it.

    xsb,ysb : Symbol;

  begin
    Symbol_Table.Init(2);
    xsb := (xsb'range => ' ');
    xsb(1) := 'x';
    Symbol_Table.Add(xsb);
    ysb := (ysb'range => ' ');
    ysb(1) := 'y';
    Symbol_Table.Add(ysb);
  end Add_Symbols;

  function Perturb ( p : Poly; eps : double_float ) return Poly is

  -- DESCRIPTION :
  --   Adds random numbers of modulus eps to each coefficient of p.

    res : Poly := Null_Poly;

    procedure Add_Term ( t : in Term; continue : out boolean ) is

      rt : Term;

    begin
      rt.dg := t.dg;
      rt.cf := t.cf + eps*Random1;    
      Add(res,rt);
      continue := true;
    end Add_Term;
    procedure Add_Terms is new Visiting_Iterator(Add_Term);

  begin
    Add_Terms(p);
    return res;
  end Perturb;

  procedure Test_Polynomial
              ( file : in file_type; deg : in Standard_Natural_Vectors.Vector;
                eps : in double_float; p : out Poly ) is

  -- DESCRIPTION :
  --   Returns the product of deg'length polynomials with degrees in deg.
  --   The product is then perturbed randomly with magnitude eps.

    f,prd : Poly;

  begin
    f := Standard_Polynomial_Interpolators.Create(deg(1),2,0);
    put(file,"Factor 1 : "); put_line(file,f);
    Copy(f,prd); Clear(f);
    for i in 2..deg'last loop
      f := Standard_Polynomial_Interpolators.Create(deg(i),2,0);
      put(file,"Factor "); put(file,i,1); put_line(file," :");
      put_line(file,f);
      Mul(prd,f); Clear(f);
    end loop;
    put_line(file,"The product : "); put_line(file,prd);
    p := Perturb(prd,eps);
    Clear(prd);
  end Test_Polynomial;

  procedure Add_Slice ( p : in Poly; sys : out Poly_Sys;
                        hypsli : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Generates a random hyperplane and adds this as slice in the
  --   system on return.  The coefficients of the slice are in hypsli.

    hyp : Poly;
    slice : constant Standard_Complex_Vectors.Vector := Random_Vector(0,2);
    t : Term;

  begin
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

  procedure Generic_Points
               ( file : in file_type; p : in Poly; sys : out Poly_Sys;
                 sols : out Solution_List;
                 hyp : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes the generic points on this polynomial.

    p1,q : Poly_Sys(1..2);
    qsols : Solution_List;
    sli : Standard_Complex_Vectors.Vector(0..2);
    epsxa : constant double_float := 1.0E-13;
    epsfa : constant double_float := 1.0E-13;
    tolsing : constant double_float := 1.0E-8;
    numit : natural32 := 0;
    max : constant natural32 := 5;
   -- cond : double_float;
   -- sccff : Standard_Complex_Vectors.Vector(1..2);
    deflate : boolean := false;

    procedure Continue is
      new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                          Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    Add_Slice(p,p1,sli);
    Scale(p1); -- Scale(p1,2,true,cond,sccff);
    Start_System(p1,q,qsols);
    Standard_Homotopy.Create(p1,q,2,Random1);
    Continuation_Parameters.Tune(2);
    Continue(qsols,false,target=>Create(1.0));
    Reporting_Root_Refiner
      (file,p1,qsols,epsxa,epsfa,tolsing,numit,max,deflate,false);
    Standard_Homotopy.Clear;
    sys := p1; sols := qsols; hyp := sli;
  end Generic_Points;

  procedure Validate_Component
              ( file : in file_type; kc : in natural32;
                c : in Standard_Irreducible_Component ) is

  -- DESCRIPTION :
  --   Creates the Newton interpolating polynomial for the component
  --   and checks the residuals of newly generated test samples.

    sps : constant Standard_Sample_List := Points(c);
    deg : constant integer32 := integer32(Length_Of(sps));
    grid : constant Array_of_Standard_Sample_Lists(0..deg)
         := Create1(sps,natural32(deg));
    cnb : constant Complex_Number := Random1;
    q : constant Newton_Interpolator1 := Create(grid,cnb);
    testsps,testsps_last,tmp : Standard_Sample_List;
    hyps : constant Standard_Complex_VecVecs.VecVec
         := Hyperplane_Sections(Head_Of(sps));
    newsli : Standard_Complex_VecVecs.VecVec(hyps'range);
    eva : Complex_Number;
    abseva : double_float;
    logeva : Standard_Floating_Vectors.Vector(1..deg);

  begin
    put(file,"The maximal error on samples : ");
    put(file,Maximal_Error(grid),3); new_line(file);
    put(file,"Maximal residual of evaluation at grid : ");
    put(file,Maximal_Error(q,grid),3); new_line(file);
    for i in newsli'range loop           -- generating test samples
      newsli(i) := new Standard_Complex_Vectors.Vector'(hyps(i).all);
      newsli(i)(0) := Random1;
    end loop;
    Sample(sps,newsli,testsps,testsps_last);
    put_line(file,"Evaluating the test samples : ");
    tmp := testsps;
    for i in 1..integer32(Length_Of(testsps)) loop
      eva := Eval(q,Sample_Point(Head_Of(tmp)).v);
      abseva := AbsVal(eva);
      if abseva + 1.0 = 1.0
       then logeva(i) := -16.0;
       else logeva(i) := LOG10(AbsVal(eva));
      end if;
      put(file,eva); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
    put(file,"LOG10(res) "); put(file,kc,1);
    put(file," : ");
    for i in logeva'range loop
      put(file,logeva(i),2,3,0);
    end loop;
    new_line(file);
  end Validate_Component;

  procedure Validate_Breakup
              ( file : in file_type;
                p : in Poly_Sys; sols : in solution_List;
                sps : in Standard_Sample_List;
                deco : in Standard_Irreducible_Component_List ) is

    tmp : Standard_Irreducible_Component_List := deco;
    c : Standard_Irreducible_Component;

  begin
    put_line(file,"Validating the breakup : ");
    for i in 1..Length_Of(deco) loop
      put(file,"Validating component number "); put(file,i,1);
      put_line(file," :");
      c := Head_Of(tmp);
      put(file,"The labels : "); put(file,Labels(c)); new_line(file);
     -- put_line(file,"The sample points : ");
     -- put(file,Points(c)); new_line(file);
      Validate_Component(file,i,c);
      tmp := Tail_Of(tmp);
    end loop;
  end Validate_Breakup;

  procedure Monodromy_Breakup
              ( file : in file_type;
                p : in Poly_Sys; sols : in Solution_List; 
                slice : in Standard_Complex_Vectors.Vector;
                nbc,nbit : out natural32 ) is

    timer : Timing_Widget;
    tol : constant double_float := 1.0E-12;
    hyp : Standard_Complex_VecVecs.VecVec(1..1);
    sps : Standard_Sample_List;
    threshold,nit : natural32;
    deco,deco_last : Standard_Irreducible_Component_List;

  begin
    hyp(1) := new Standard_Complex_Vectors.Vector'(slice);
    sps := Create(sols,hyp);
   -- new_line;
    threshold := 10;
   -- put("Give the stabilizing threshold : "); get(threshold);
   -- new_line;
   -- put_line("See the output file for results...");
   -- new_line;
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner;
    tstart(timer);
    Monodromy_Breakup(file,sps,threshold,tol,deco,deco_last,nit);
    tstop(timer);
    new_line(file);
    Write_Summary(file,1,deco);
    new_line(file);
    put_line(file,"The labels of the list of components : ");
    put_labels(file,deco);
    new_line(file);
    put(file,"Number of iterations : "); put(file,nit,1); new_line(file);
    new_line(file);
    print_times(file,timer,"Monodromy Group Action Breakup");
    nbc := Length_Of(deco);
    nbit := nit;
    Validate_Breakup(file,p,sols,sps,deco);
    Sampling_Machine.Clear;
    Clear(sps);
    Clear(deco);
  end Monodromy_Breakup;

  procedure Test_Sensitivity ( file : in file_type;
                               deg : in Standard_Natural_Vectors.Vector ) is

    eps : double_float := 0.0;
    hyp : Standard_Complex_Vectors.Vector(0..2);
    nbc,nit : natural32;
    iters,comps : Standard_Natural_Vectors.Vector(0..14);

  begin
    for i in 0..integer32(14) loop
      eps := 10.0**integer(-i);
      declare
        p : Poly;
        sys : Poly_Sys(1..2);
        sols : Solution_List;
      begin
        Test_Polynomial(file,deg,eps,p);
        Generic_Points(file,p,sys,sols,hyp);
        Monodromy_Breakup(file,sys,sols,hyp,nbc,nit);
      end;
      put(file," eps = ");
      put(file,eps,2,3,3);
      put(file," : "); put(file,nbc,1); 
      put(file," : "); put(file,nit,1); new_line(file);
      comps(i) := nbc;
      iters(i) := nit;
    end loop;
    new_line(file);
    put_line(file,"SUMMARY OF SENSITIVITY EXPERIMENT :");
    put_line(file,"  eps  :  #comps : iterations");
    for i in 0..integer32(14) loop
      eps := 10.0**integer(-i);
      put(file,eps,2,3,3);
      put(file," : "); put(file,comps(i),1);
      put(file," : "); put(file,iters(i),1);
      new_line(file);
    end loop;
  end Test_Sensitivity;

  procedure Main is

    file : file_type;
    n : integer32 := 0;

  begin
    new_line;
    put_line("Multivariate Factorization with Monodromy Group Actions.");
    new_line;
    Add_Symbols;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put("Give the number of factors : "); get(n);
    declare
      deg : Standard_Natural_Vectors.Vector(1..n);
    begin
      put("Give "); put(n,1); put(" degrees for each factor : "); get(deg);
      new_line;
      put_line("See the output file for results...");
      new_line;
      Test_Sensitivity(file,deg);
    end;
  end Main;

begin
  Main;
end Sensitivity_of_Factorization;
