with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;
with Double_Double_Vectors;
with Quad_Double_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Continuation_Parameters;
with Standard_Continuation_Data;
with Standard_Path_Trackers;
with DoblDobl_Continuation_Data;
with DoblDobl_Path_Trackers;
with QuadDobl_Continuation_Data;
with QuadDobl_Path_Trackers;
with Polyhedral_Coefficient_Homotopies;

package body Single_Polyhedral_Trackers is

  procedure Track_Path
               ( mxt : in Standard_Integer_Vectors.Vector;
                 lft : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 nor : in Standard_Floating_Vectors.Link_to_Vector;
                 cff : in Standard_Complex_VecVecs.VecVec;
                 dpw : in Standard_Floating_VecVecs.Link_to_VecVec;
                 cft : in Standard_Complex_VecVecs.Link_to_VecVec;
                 epv : in Exponent_Vectors.Exponent_Vectors_Array;
                 hom : in Standard_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                 ejf : in Standard_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                 jmf : in Standard_Complex_Laur_JacoMats.Mult_Factors;
                 sol : in Standard_Complex_Solutions.Link_to_Solution ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Norms_Equals;
    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;
    use Standard_Complex_Solutions;
    use Standard_Continuation_Data;
    use Standard_Path_Trackers;
    use Polyhedral_Coefficient_Homotopies;

    function Eval ( x : Standard_Complex_Vectors.Vector;
                    t : Complex_Number )
                  return Standard_Complex_Vectors.Vector is
    begin
      Eval(cff,REAL_PART(t),dpw.all,cft.all);
      return Eval(hom,cft.all,x);
    end Eval;

    function dHt ( x : Standard_Complex_Vectors.Vector;
                   t : Standard_Complex_Numbers.Complex_Number )
                 return Standard_Complex_Vectors.Vector is

      res : Standard_Complex_Vectors.Vector(hom'range);
      xtl : constant integer32 := x'last+1;

    begin
      Eval(cff,REAL_PART(t),dpw.all,cft.all);
      for k in res'range loop
        res(k) := Eval(ejf(k,xtl),jmf(k,xtl).all,cft(k).all,x);
      end loop;
      return res;
    end dHt;

    function dHx ( x : Standard_Complex_Vectors.Vector;
                   t : Complex_Number )
                 return Standard_Complex_Matrices.Matrix is

      mt : Standard_Complex_Matrices.Matrix(x'range,x'range);

    begin
      Eval(cff,REAL_PART(t),dpw.all,cft.all);
      for k in mt'range(1) loop
        for L in mt'range(2) loop
          mt(k,L) := Eval(ejf(k,L),jmf(k,L).all,cft(k).all,x);
        end loop;
      end loop;
      return mt;
    end dHx;

    procedure Track_Path_along_Path is
      new Linear_Single_Normal_Silent_Continue(Max_Norm,Eval,dHt,dHx);
    procedure Track_Path_at_End is
      new Linear_Single_Conditioned_Silent_Continue(Max_Norm,Eval,dHt,dHx);

    procedure Main is

      t1 : constant Complex_Number := Create(1.0);
      tol : constant double_float := 1.0E-12;
      pp1 : Continuation_Parameters.Pred_Pars
          := Continuation_Parameters.Create_for_Path;
      pp2 : constant Continuation_Parameters.Pred_Pars
          := Continuation_Parameters.Create_End_Game;
      cp1 : constant Continuation_Parameters.Corr_Pars
          := Continuation_Parameters.Create_for_Path;
      cp2 : constant Continuation_Parameters.Corr_Pars
          := Continuation_Parameters.Create_End_Game;

    begin
      Power_Transform(epv,lft,mxt,nor.all,dpw.all);
      Polyhedral_Coefficient_Homotopies.Scale(dpw.all);
      sol.t := Create(0.0);
      declare
        s : Solu_Info := Shallow_Create(sol);
        w : integer32 := 1;
        v : Standard_Floating_Vectors.Link_to_Vector;
        e : double_float := 0.0;
      begin
        pp1.dist_target := 0.0;
        Track_Path_along_Path(s,t1,tol,false,pp1,cp1);
        Track_Path_at_End(s,t1,tol,false,0,w,v,e,pp2,cp2);
        sol.err := s.cora; sol.rco := s.rcond; sol.res := s.resa;
      end;
    end Main;

  begin
    Main;
  end Track_Path;

  procedure Track_Path
               ( mxt : in Standard_Integer_Vectors.Vector;
                 lft : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 nor : in Standard_Floating_Vectors.Link_to_Vector;
                 cff : in DoblDobl_Complex_VecVecs.VecVec;
                 dpw : in Standard_Floating_VecVecs.Link_to_VecVec;
                 cft : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                 epv : in Exponent_Vectors.Exponent_Vectors_Array;
                 hom : in DoblDobl_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                 ejf : in DoblDobl_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                 jmf : in DoblDobl_Complex_Laur_JacoMats.Mult_Factors;
                 sol : in DoblDobl_Complex_Solutions.Link_to_Solution ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vector_Norms;
    use DoblDobl_Complex_Laur_SysFun;
    use DoblDobl_Complex_Laur_JacoMats;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Continuation_Data;
    use DoblDobl_Path_Trackers;
    use Polyhedral_Coefficient_Homotopies;

    function Eval ( x : DoblDobl_Complex_Vectors.Vector;
                    t : Complex_Number )
                  return DoblDobl_Complex_Vectors.Vector is
    begin
      Eval(cff,REAL_PART(t),dpw.all,cft.all);
      return Eval(hom,cft.all,x);
    end Eval;

    function dHt ( x : DoblDobl_Complex_Vectors.Vector;
                   t : DoblDobl_Complex_Numbers.Complex_Number )
                 return DoblDobl_Complex_Vectors.Vector is

      res : DoblDobl_Complex_Vectors.Vector(hom'range);
      xtl : constant integer32 := x'last+1;

    begin
      Eval(cff,REAL_PART(t),dpw.all,cft.all);
      for k in res'range loop
        res(k) := Eval(ejf(k,xtl),jmf(k,xtl).all,cft(k).all,x);
      end loop;
      return res;
    end dHt;

    function dHx ( x : DoblDobl_Complex_Vectors.Vector;
                   t : Complex_Number )
                 return DoblDobl_Complex_Matrices.Matrix is

      mt : DoblDobl_Complex_Matrices.Matrix(x'range,x'range);

    begin
      Eval(cff,REAL_PART(t),dpw.all,cft.all);
      for k in mt'range(1) loop
        for L in mt'range(2) loop
          mt(k,L) := Eval(ejf(k,L),jmf(k,L).all,cft(k).all,x);
        end loop;
      end loop;
      return mt;
    end dHx;

    procedure Track_Path_along_Path is
      new Linear_Single_Normal_Silent_Continue(Max_Norm,Eval,dHt,dHx);
    procedure Track_Path_at_End is
      new Linear_Single_Conditioned_Silent_Continue(Max_Norm,Eval,dHt,dHx);

    procedure Main is

      zero : constant double_double := create(0.0);
      one : constant double_double := create(1.0);
      t1 : constant Complex_Number := Create(one);
      tol : constant double_float := 1.0E-12;
      pp1 : Continuation_Parameters.Pred_Pars
          := Continuation_Parameters.Create_for_Path;
      pp2 : constant Continuation_Parameters.Pred_Pars
          := Continuation_Parameters.Create_End_Game;
      cp1 : constant Continuation_Parameters.Corr_Pars
          := Continuation_Parameters.Create_for_Path;
      cp2 : constant Continuation_Parameters.Corr_Pars
          := Continuation_Parameters.Create_End_Game;

    begin
      Power_Transform(epv,lft,mxt,nor.all,dpw.all);
      Polyhedral_Coefficient_Homotopies.Scale(dpw.all);
      sol.t := Create(zero);
      declare
        s : Solu_Info := Shallow_Create(sol);
        w : integer32 := 1;
        v : Double_Double_Vectors.Link_to_Vector;
        e : double_double := create(0.0);
      begin
        pp1.dist_target := 0.0;
        Track_Path_along_Path(s,t1,tol,false,pp1,cp1);
        Track_Path_at_End(s,t1,tol,false,0,w,v,e,pp2,cp2);
        sol.err := create(s.cora);
        sol.rco := create(s.rcond);
        sol.res := create(s.resa);
      end;
    end Main;

  begin
    Main;
  end Track_Path;

  procedure Track_Path
               ( mxt : in Standard_Integer_Vectors.Vector;
                 lft : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 nor : in Standard_Floating_Vectors.Link_to_Vector;
                 cff : in QuadDobl_Complex_VecVecs.VecVec;
                 dpw : in Standard_Floating_VecVecs.Link_to_VecVec;
                 cft : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                 epv : in Exponent_Vectors.Exponent_Vectors_Array;
                 hom : in QuadDobl_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                 ejf : in QuadDobl_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                 jmf : in QuadDobl_Complex_Laur_JacoMats.Mult_Factors;
                 sol : in QuadDobl_Complex_Solutions.Link_to_Solution ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Vector_Norms;
    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Continuation_Data;
    use QuadDobl_Path_Trackers;
    use Polyhedral_Coefficient_Homotopies;

    function Eval ( x : QuadDobl_Complex_Vectors.Vector;
                    t : Complex_Number )
                  return QuadDobl_Complex_Vectors.Vector is
    begin
      Eval(cff,REAL_PART(t),dpw.all,cft.all);
      return Eval(hom,cft.all,x);
    end Eval;

    function dHt ( x : QuadDobl_Complex_Vectors.Vector;
                   t : QuadDobl_Complex_Numbers.Complex_Number )
                 return QuadDobl_Complex_Vectors.Vector is

      res : QuadDobl_Complex_Vectors.Vector(hom'range);
      xtl : constant integer32 := x'last+1;

    begin
      Eval(cff,REAL_PART(t),dpw.all,cft.all);
      for k in res'range loop
        res(k) := Eval(ejf(k,xtl),jmf(k,xtl).all,cft(k).all,x);
      end loop;
      return res;
    end dHt;

    function dHx ( x : QuadDobl_Complex_Vectors.Vector;
                   t : Complex_Number )
                 return QuadDobl_Complex_Matrices.Matrix is

      mt : QuadDobl_Complex_Matrices.Matrix(x'range,x'range);

    begin
      Eval(cff,REAL_PART(t),dpw.all,cft.all);
      for k in mt'range(1) loop
        for L in mt'range(2) loop
          mt(k,L) := Eval(ejf(k,L),jmf(k,L).all,cft(k).all,x);
        end loop;
      end loop;
      return mt;
    end dHx;

    procedure Track_Path_along_Path is
      new Linear_Single_Normal_Silent_Continue(Max_Norm,Eval,dHt,dHx);
    procedure Track_Path_at_End is
      new Linear_Single_Conditioned_Silent_Continue(Max_Norm,Eval,dHt,dHx);

    procedure Main is

      zero : constant quad_double := create(0.0);
      one : constant quad_double := create(1.0);
      t1 : constant Complex_Number := Create(one);
      tol : constant double_float := 1.0E-12;
      pp1 : Continuation_Parameters.Pred_Pars
          := Continuation_Parameters.Create_for_Path;
      pp2 : constant Continuation_Parameters.Pred_Pars
          := Continuation_Parameters.Create_End_Game;
      cp1 : constant Continuation_Parameters.Corr_Pars
          := Continuation_Parameters.Create_for_Path;
      cp2 : constant Continuation_Parameters.Corr_Pars
          := Continuation_Parameters.Create_End_Game;

    begin
      Power_Transform(epv,lft,mxt,nor.all,dpw.all);
      Polyhedral_Coefficient_Homotopies.Scale(dpw.all);
      sol.t := Create(zero);
      declare
        s : Solu_Info := Shallow_Create(sol);
        w : integer32 := 1;
        v : Quad_Double_Vectors.Link_to_Vector;
        e : quad_double := create(0.0);
      begin
        pp1.dist_target := 0.0;
        Track_Path_along_Path(s,t1,tol,false,pp1,cp1);
        Track_Path_at_End(s,t1,tol,false,0,w,v,e,pp2,cp2);
        sol.err := create(s.cora);
        sol.rco := create(s.rcond);
        sol.res := create(s.resa);
      end;
    end Main;

  begin
    Main;
  end Track_Path;

end Single_Polyhedral_Trackers;
