with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;         use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;             use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;       use Standard_Complex_Norms_Equals;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;    use Standard_Complex_Poly_Systems_io;
with Standard_Homotopy;
with Continuation_Parameters;
with Standard_IncFix_Continuation;
with Standard_Root_Refiners;              use Standard_Root_Refiners;
with Standard_to_Multprec_Convertors;     use Standard_to_Multprec_Convertors;
with Planes_and_Polynomials;              use Planes_and_Polynomials;
with Multprec_Floating_Numbers;           use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;
with Multprec_Complex_Number_Tools;       use Multprec_Complex_Number_Tools;
with Multprec_Complex_Vectors;
with Multprec_Complex_Norms_Equals;       use Multprec_Complex_Norms_Equals;
with Multprec_Complex_Matrices;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Jaco_Matrices;
with Multprec_IncFix_Continuation;
with Multprec_Root_Refiners;              use Multprec_Root_Refiners;
with Extrapolate_Solution_Clusters;       use Extrapolate_Solution_Clusters;

package body Multiplicity_Homotopies is

  procedure Write_Summary
                   ( file : in file_type;
                     sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Writes the (err,rco,res) results for every solution in the list.

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      put(file,"solution "); put(file,i,1); put(file," : ");
      put(file,"  err : "); put(file,ls.err,3);
      put(file,"  rco : "); put(file,ls.rco,3);
      put(file,"  res : "); put(file,ls.res,3);
      new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Summary;

  procedure Reconditioning_Homotopy_Continuation
               ( file : in file_type;
                 p,q : in Standard_Complex_Poly_Systems.Poly_Sys;
                 k : in integer32; 
                 a : in Standard_Complex_Numbers.Complex_Number;
                 sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Performs a standard homotopy between p and q starting at the
  --   solutions in the solution list sols.

  -- ON ENTRY :
  --   p         target system;
  --   q         start system;
  --   k         relaxation parameter in the homotopy;
  --   a         accessibility constant;
  --   sols      start solutions leading to the cluster.

  -- ON RETURN :
  --   sols      solutions in cluster that are better conditioned.

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;
    use Standard_IncFix_Continuation;

    target : constant Complex_Number := Create(1.0);

    procedure Cont is
      new Reporting_Continue(Max_Norm,Standard_Homotopy.Eval,
                                      Standard_Homotopy.Diff,
                                      Standard_Homotopy.Diff);

  begin
    Standard_Homotopy.Create(p,q,natural32(k),a);
    Standard_Complex_Solutions.Set_Continuation_Parameter(sols,Create(0.0));
    Continuation_Parameters.Tune(2);
    Continuation_Parameters.tol_path_distance := 0.0;
    Continuation_Parameters.tol_endg_distance := 0.0;
    put_line(file,"The reconditioned cluster : ");
    put(file,Length_Of(sols),1); put(file," ");
    put(file,Head_Of(sols).n,1);
    new_line(file);
    Cont(file,sols,false,target);
    declare
      epsxa : constant double_float := 1.0E-13;
      epsfa : constant double_float := 1.0E-13;
      tolsing : constant double_float := 1.0E-13;
      nit : natural32 := 0;
      max : constant natural32 := 3;
      deflate : boolean := false;
    begin
      Reporting_Root_Refiner
        (file,p,sols,epsxa,epsfa,tolsing,nit,max,deflate,false);
    end;
  end Reconditioning_Homotopy_Continuation;

  procedure Incoming_Homotopy_Continuation
               ( file : in file_type;
                 p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 deci : in natural32; ind : in integer32;
                 ze : in Multprec_Complex_Numbers.Complex_Number;
                 sols : in out Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Performs a homotopy from z = ze to z = 0, by p(ind): z + (1-t)*ze = 0,
  --   for t ranging from zero to one.

    use Multprec_Complex_Numbers,Multprec_Complex_Vectors;
    use Multprec_Complex_Solutions;
    use Multprec_Complex_Poly_SysFun;
    use Multprec_Complex_Matrices,Multprec_Complex_Jaco_Matrices;
    use Multprec_IncFix_Continuation;

    target : constant Complex_Number := Create(integer(1));
    p_eval : constant Eval_Poly_Sys(p'range) := Create(p);
    jac : constant Jaco_Mat(p'range,p'range) := Create(p);
    jac_eval : constant Eval_Jaco_Mat(p'range,p'range) := Create(jac);

    function Hom_Eval ( x : Vector; t : Complex_Number ) return Vector is

      res : Vector(x'range) := Eval(p_eval,x);
      acc : Complex_Number := target - t;

    begin
      Mul(acc,ze);
      Sub(res(ind),acc);
      Clear(acc);
      return res;
    end Hom_Eval;

    function Dummy ( x : Vector; t : Complex_Number ) return Vector is

    -- NOTE : this procedure will not be used with Tune(0) and Tune(2).

      res : constant Vector(x'range) := x;

    begin
      return res;
    end Dummy;

    function Diff_Eval ( x : Vector; t : Complex_Number ) return Matrix is

      res : constant Matrix(x'range,x'range) := Eval(jac_eval,x);

    begin
      return res;
    end Diff_Eval;

    procedure Cont is
      new Reporting_Continue(Max_Norm,Hom_Eval,Dummy,Diff_Eval);

  begin
    Multprec_Complex_Solutions.Set_Continuation_Parameter
      (sols,Create(integer(0)));
    Continuation_Parameters.Tune(2,deci);
    Cont(file,sols,false,target);
  end Incoming_Homotopy_Continuation;

  procedure Reconditioning_Homotopy
               ( file : in file_type;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 q : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 k : in integer32;
                 a : in Standard_Complex_Numbers.Complex_Number;
                 ind : in integer32;
                 ze : out Standard_Complex_Numbers.Complex_Number ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;
    use Standard_Complex_Solutions;

    len : constant natural32 := Length_Of(sols);
    eps : constant double_float := 10.0**integer((-16)/len);
    t : Term;

  begin
    new_line(file);
    put_line(file,"Reconditioning the cluster of solutions ...");
    ze := Random1*Create(eps);
    put(file,"  eps = "); put(file,eps); new_line(file);
    put(file,"  z(eps) = "); put(file,ze); new_line(file);
    t.dg := new Standard_Natural_Vectors.Vector'(p'range => 0);
    t.cf := ze;
    Sub(p(ind),t);
    put_line(file,"The target system :");
    put_line(file,p);
    Reconditioning_Homotopy_Continuation(file,p,q,k,a,sols);
  end Reconditioning_Homotopy;

  procedure Incoming_Homotopy
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 ind : in integer32;
                 ze : in Standard_Complex_Numbers.Complex_Number;
                 mpsols : out Multprec_Complex_Solutions.Solution_List ) is

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Polynomials;

    len : constant natural32
        := Standard_Complex_Solutions.Length_Of(sols);
    deci : constant natural32 := 32*len;
    size : constant natural32 := Decimal_to_Size(deci);
    mp : Multprec_Complex_Poly_Systems.Poly_Sys(p'range) := Convert(p);
    mpze : Multprec_Complex_Numbers.Complex_Number := Create(ze);
    epsxa : constant double_float := 10.0**integer(-deci+6);
    epsfa : constant double_float := 10.0**integer(-deci+3);
    tolsing : constant double_float := 10.0**integer(-8);
    mpepsxa : constant Floating_Number := Create(epsxa);
    mpepsfa : constant Floating_Number := Create(epsfa);
    mptolsing : constant Floating_Number := Create(tolsing);
    max : constant natural32 := 5;
    numit : natural32 := 0;
    zet : Term;
    ten : constant Multprec_Complex_Numbers.Complex_Number
        := Create(integer(10));
    order : constant integer32 := 5;
    frame1 : Multprec_Complex_Solutions.Solution_Array(0..order);
    frame2 : Multprec_Complex_Solutions.Solution_Array(0..order);
    extra1,extra2 : Multprec_Complex_Vectors.Vector(p'range);
    inc : constant natural32 := deci/natural32(order);
    frame_ind : integer32 := 0;

    use Multprec_Complex_Solutions;

  begin
    mpsols := Multprec_Complex_Solutions.Create(sols);
    new_line(file);
    put(file,"Number of decimal places : ");
    put(file,deci,1); new_line(file);
    Set_Size(mp,size);
    Set_Size(mpze,size);
    Multprec_Complex_Solutions.Set_Size(mpsols,size);
    zet.dg := new Standard_Natural_Vectors.Vector'(p'range => 0);
    Copy(mpze,zet.cf); 
    for i in 1..deci loop
     -- Reporting_Root_Refiner
     --   (file,mp,mpsols,mpepsxa,mpepsfa,mptolsing,numit,max,false,false);
      Silent_Root_Refiner
        (mp,mpsols,mpepsxa,mpepsfa,mptolsing,numit,max,false);
      if i mod inc = 0 and frame_ind <= order then
        frame1(frame_ind) := new Solution(p'last);
        frame2(frame_ind) := new Solution(p'last);
        Copy(Head_Of(mpsols).all,frame1(frame_ind).all);
        Copy(Head_Of(Tail_Of(mpsols)).all,frame2(frame_ind).all);
        frame_ind := frame_ind+1;
      end if;
      Add(mp(ind),zet);
      Div(zet.cf,ten);
      Sub(mp(ind),zet);
     -- put_line(file,"The new moving equation : ");
     -- put(file,mp(ind)); new_line(file);
    end loop;
    frame1(order) := new Solution(p'last);
    frame2(order) := new Solution(p'last);
    Copy(Head_Of(mpsols).all,frame1(order).all);
    Copy(Head_Of(Tail_Of(mpsols)).all,frame2(order).all);
    put_line(file,"First Extrapolation : ");
    Extrapolate(file,order,integer32(len),p'last,size,frame1,extra1);
    Multprec_Complex_Vectors.Copy(extra1,Head_Of(mpsols).v);
    put_line(file,"Second Extrapolation : ");
    Extrapolate(file,order,integer32(len),p'last,size,frame2,extra2);
    Multprec_Complex_Vectors.Copy(extra2,Head_Of(Tail_Of(mpsols)).v);
  end Incoming_Homotopy;

  procedure Sampling_Homotopy
               ( file : in file_type; eps : in double_float;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 hyp : out Standard_Complex_VecVecs.VecVec;
                 sols : in out Standard_Complex_Solutions.Solution_List ) is

    slice : Standard_Complex_Vectors.Vector(0..p'last)
          := Polynomial(p(p'last));
    inc : constant Standard_Complex_Numbers.Complex_Number
        := Standard_Complex_Numbers.Create(eps,eps);
    epsxa : constant double_float := 10.0**(-13);
    epsfa : constant double_float := 10.0**(-13);
    tolsing : constant double_float := 10.0**(-8);
    max : constant natural32 := 5;
    nit : natural32 := 0;
    deflate : boolean := false;

  begin
    for i in slice'range loop
      Standard_Complex_Numbers.Add(slice(i),inc);
    end loop;
    p(p'last) := Hyperplane(slice);
   -- put_line(file,p);
   -- put_line(file,"The moving hyperplane : ");
   -- put_line(file,p(p'last));
   -- Reporting_Root_Refiner
   --   (file,p,sols,epsxa,epsfa,tolsing,nit,max,false,false);
    Silent_Root_Refiner(p,sols,epsxa,epsfa,tolsing,nit,max,deflate);
    Write_Summary(file,sols);
    hyp(1) := new Standard_Complex_Vectors.Vector'(slice);
  end Sampling_Homotopy;

end Multiplicity_Homotopies;
