with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vector_Norms;      use DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Multprec_Random_Numbers;            use Multprec_Random_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Multprec_Complex_Vectors;
with Multprec_Complex_Matrices;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;
with Multprec_System_and_Solutions_io;
with Process_io;                         use Process_io;
with Continuation_Parameters;
with Standard_Continuation_Data;
with Standard_Predictors;
with DoblDobl_Predictors;
with Standard_Correctors;
with Multprec_Continuation_Data;
with Multprec_Predictors;
with Multprec_Correctors;

procedure ts_preco is

-- DESCRIPTION : test on predictor and corrector methods.

  procedure Standard_Perturb_Solution
              ( sol : in out Standard_Complex_Solutions.Solution;
                permag : in double_float ) is

  -- DESCRIPTION :
  --   Perturb each component of the solution vector with a random
  --   number of the magnitude permag.

    use Standard_Complex_Numbers;

    perfac : Complex_Number;

  begin
    for i in sol.v'range loop
      perfac := Random1;
      perfac := permag*perfac;
      sol.v(i) := sol.v(i) + perfac; 
    end loop;
  end Standard_Perturb_Solution;

  procedure Standard_Perturb_Solutions
              ( sols : in out Standard_Complex_Solutions.Solution_List;
                permag : in double_float ) is

  -- DESCRIPTION :
  --   Perturbs every solution in the list randomly with given
  --   magnitude permag.

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Standard_Perturb_Solution(ls.all,permag);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Standard_Perturb_Solutions;

  procedure Multprec_Perturb_Solution
              ( sol : in out Multprec_Complex_Solutions.Solution;
                size : in natural32; permag : in double_float ) is

  -- DESCRIPTION :
  --   Perturb each component of the solution vector with a random
  --   number of the magnitude permag.

    use Multprec_Complex_Numbers;

    magper : Floating_Number := Create(permag);
    perfac : Complex_Number;

  begin
    for i in sol.v'range loop
      perfac := Random(size);
      Mul(perfac,magper);
      Add(sol.v(i),perfac); 
      Clear(perfac);
    end loop;
    Clear(magper);
  end Multprec_Perturb_Solution;

  procedure Multprec_Perturb_Solutions
              ( sols : in out Multprec_Complex_Solutions.Solution_List;
                size : in natural32; permag : in double_float ) is

  -- DESCRIPTION :
  --   Perturbs every solution in the list randomly with given
  --   magnitude permag.

    use Multprec_Complex_Solutions;

    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Multprec_Perturb_Solution(ls.all,size,permag);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Multprec_Perturb_Solutions;

  procedure Standard_Tangent_Predictor
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                target : in Standard_Complex_Numbers.Complex_Number;
                step : in double_float ) is

  -- DESCRIPTION :
  --   Applies the tangent predictor to the solution of the system.

    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Predictors;

    tol : constant double_float := 10.0**(-13);
    jac : Jaco_Mat(p'range,p'range) := Create(p);
    jac_eval : Eval_Jaco_Mat(p'range,p'range) := Create(jac);

    function Diff_t ( x : Vector; t : Complex_Number ) return Vector is

    -- DESCRIPTION :
    --   By taking ones as derivatives, this assumes the perturbation
    --   to the coefficients is linear.

      res : constant Vector(x'range) := (x'range => Create(1.0));

    begin
      return res;
    end Diff_t;

    function Diff_x ( x : Vector; t : Complex_Number ) return Matrix is

      res : Matrix(x'range,x'range);

    begin
      res := Eval(jac_eval,x);
      return res;
    end Diff_x;

    procedure Predict is
      new Tangent_Single_Real_Predictor(Max_Norm,Diff_t,Diff_x);

  begin
    Predict(sol.v,sol.t,target,step,tol);
    Clear(jac); Clear(jac_eval);
  end Standard_Tangent_Predictor;

  procedure DoblDobl_Tangent_Predictor
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                target : in DoblDobl_Complex_Numbers.Complex_Number;
                step : in double_float ) is

  -- DESCRIPTION :
  --   Applies the tangent predictor to the solution of the system.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Predictors;

    tol : constant double_float := 10.0**(-30);
    jac : Jaco_Mat(p'range,p'range) := Create(p);
    jac_eval : Eval_Jaco_Mat(p'range,p'range) := Create(jac);
    one : constant double_double := create(1.0);

    function Diff_t ( x : Vector; t : Complex_Number ) return Vector is

    -- DESCRIPTION :
    --   By taking ones as derivatives, this assumes the perturbation
    --   to the coefficients is linear.

      res : constant Vector(x'range) := (x'range => Create(one));

    begin
      return res;
    end Diff_t;

    function Diff_x ( x : Vector; t : Complex_Number ) return Matrix is

      res : Matrix(x'range,x'range);

    begin
      res := Eval(jac_eval,x);
      return res;
    end Diff_x;

    procedure Predict is
      new Tangent_Single_Real_Predictor(Max_Norm,Diff_t,Diff_x);

  begin
    Predict(sol.v,sol.t,target,step,tol);
    Clear(jac); Clear(jac_eval);
  end DoblDobl_Tangent_Predictor;

  procedure Multprec_Tangent_Predictor
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sol : in out Multprec_Complex_Solutions.Solution;
                target : in Multprec_Complex_Numbers.Complex_Number;
                step : in Floating_Number; size : in natural32 ) is

  -- DESCRIPTION :
  --   Applies the tangent predictor to the solution of the system.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Vectors;
    use Multprec_Complex_Matrices;
    use Multprec_Complex_Jaco_Matrices;
    use Multprec_Predictors;

    tol : constant Floating_Number := Create(10.0**((-4)*integer(size)));
    jac : Jaco_Mat(p'range,p'range) := Create(p);
    jac_eval : Eval_Jaco_Mat(p'range,p'range) := Create(jac);

    function Diff_t ( x : Vector; t : Complex_Number ) return Vector is

    -- DESCRIPTION :
    --   By taking ones as derivatives, this assumes the perturbation
    --   to the coefficients is linear.

      res : constant Vector(x'range) := (x'range => Create(integer(1)));

    begin
      return res;
    end Diff_t;

    function Diff_x ( x : Vector; t : Complex_Number ) return Matrix is

      res : Matrix(x'range,x'range);

    begin
      res := Eval(jac_eval,x);
      return res;
    end Diff_x;

    procedure Predict is
      new Tangent_Single_Real_Predictor(Max_Norm,Diff_t,Diff_x);

  begin
    Predict(sol.v,sol.t,target,step,tol);
    Clear(jac); Clear(jac_eval);
  end Multprec_Tangent_Predictor;

  procedure Test_Standard_Tangent_Predictor is

  -- DESCRIPTION :
  --   Reads in a polynomial system and a list of solutions,
  --   after setting the parameters, does one predictor step.

    use Standard_Complex_Numbers;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    sols,tmp : Solution_List;
    target : Complex_Number;
    step : double_float := 0.1;

  begin
    Standard_System_and_Solutions_io.get(lp,sols);
    new_line;
    put("Give complex target value : "); get(target);
    put("Give real step size : "); get(step);
    tmp := sols;
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Standard_Tangent_Predictor(lp.all,ls.all,target,step);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    put_line("The predicted solutions : "); put(sols);
  end Test_Standard_Tangent_Predictor;

  procedure Test_DoblDobl_Tangent_Predictor is

  -- DESCRIPTION :
  --   Reads in a polynomial system and a list of solutions,
  --   after setting the parameters, does one predictor step.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    sols,tmp : Solution_List;
    target : Complex_Number;
    step : double_float := 0.1;

  begin
    DoblDobl_System_and_Solutions_io.get(lp,sols);
    new_line;
    put("Give complex target value : "); get(target);
    put("Give real step size : "); get(step);
    tmp := sols;
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        DoblDobl_Tangent_Predictor(lp.all,ls.all,target,step);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    put_line("The predicted solutions : "); put(sols);
  end Test_DoblDobl_Tangent_Predictor;

  procedure Test_QuadDobl_Tangent_Predictor is

  -- DESCRIPTION :
  --   Reads in a polynomial system and a list of solutions,
  --   after setting the parameters, does one predictor step.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    sols,tmp : Solution_List;
    target : Complex_Number;
    step : double_float := 0.1;

  begin
    QuadDobl_System_and_Solutions_io.get(lp,sols);
    new_line;
    put("Give complex target value : "); get(target);
    put("Give real step size : "); get(step);
    tmp := sols;
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
       -- DoblDobl_Tangent_Predictor(lp.all,ls.all,target,step);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    put_line("The predicted solutions : "); put(sols);
  end Test_QuadDobl_Tangent_Predictor;

  procedure Test_Multprec_Tangent_Predictor is

  -- DESCRIPTION :
  --   Reads in a polynomial system and a list of solutions,
  --   after setting the parameters, does one predictor step.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Solutions;

    mplp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    mpsols,tmp : Multprec_Complex_Solutions.Solution_List;
    target : Complex_Number;
    step : Floating_Number;
    size : natural32 := 0;

  begin
    Multprec_System_and_Solutions_io.get(mplp,mpsols);
    new_line;
    put("Give complex target value : "); get(target);
    put("Give real step size : "); get(step);
    put("Give the size of the numbers : "); get(size);
    Set_Size(mplp.all,size);
    Multprec_Complex_Solutions.Set_Size(mpsols,size);
    tmp := mpsols;
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Multprec_Tangent_Predictor(mplp.all,ls.all,target,step,size);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    put_line("The predicted solutions : "); put(mpsols);
  end Test_Multprec_Tangent_Predictor;

  procedure Call_Standard_Corrector 
              ( p_eval : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jac_eval : in Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                cp : in Continuation_Parameters.Corr_Pars; wout : in boolean;
                ls : in out Standard_Complex_Solutions.Link_to_Solution ) is

  -- DESCRIPTION :
  --   Applies the standard corrector to refine the solution of the system p.

    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;
    use Standard_Continuation_Data;
    use Standard_Correctors;

    solinf : Solu_Info := Shallow_Create(ls);

    function Eval ( x : Vector; t : Complex_Number ) return Vector is

      res : constant Vector(x'range)
          := Standard_Complex_Poly_SysFun.Eval(p_eval,x);

    begin
      return res;
    end Eval;

    function Diff ( x : Vector; t : Complex_Number ) return Matrix is

      res : constant Matrix(x'range,x'range)
          := Standard_Complex_Jaco_Matrices.Eval(jac_eval,x);

    begin
      return res;
    end Diff; 

    procedure Silent_Corrector is
      new Affine_Single_Loose_Conditioned_Silent_Corrector
            (Max_Norm,Eval,Diff);
    procedure Reporting_Corrector is
      new Affine_Single_Loose_Conditioned_Reporting_Corrector
            (Max_Norm,Eval,Diff);

  begin
    if wout
     then Process_io.Set_Output_Code(c);
          Reporting_Corrector(Standard_Output,solinf,cp);
     else Silent_Corrector(solinf,cp);
    end if;
  end Call_Standard_Corrector;

  procedure Call_Multprec_Corrector 
              ( p_eval : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                jac_eval : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                cp : in Multprec_Continuation_Data.Corr_Pars;
                wout : in boolean;
                ls : in out Multprec_Complex_Solutions.Link_to_Solution ) is

  -- DESCRIPTION :
  --   Applies the multprec corrector to refine the solution of the system p.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Vectors;
    use Multprec_Complex_Matrices;
    use Multprec_Continuation_Data;
    use Multprec_Correctors;

    solinf : Solu_Info := Shallow_Create(ls);

    function Eval ( x : Vector; t : Complex_Number ) return Vector is

      res : Vector(x'range) := Multprec_Complex_Poly_SysFun.Eval(p_eval,x);

    begin
      return res;
    end Eval;

    function Diff ( x : Vector; t : Complex_Number ) return Matrix is

      res : constant Matrix(x'range,x'range)
          := Multprec_Complex_Jaco_Matrices.Eval(jac_eval,x);

    begin
      return res;
    end Diff; 

    procedure Silent_Corrector is
      new Affine_Single_Loose_Conditioned_Silent_Corrector
            (Max_Norm,Eval,Diff);
    procedure Reporting_Corrector is
      new Affine_Single_Loose_Conditioned_Reporting_Corrector
            (Max_Norm,Eval,Diff);

  begin
    if wout
     then Process_io.Set_Output_Code(c);
          Reporting_Corrector(Standard_Output,solinf,cp);
     else Silent_Corrector(solinf,cp);
    end if;
  end Call_Multprec_Corrector;

  procedure Standard_Corrector
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                cp : in Continuation_Parameters.Corr_Pars; wout : in boolean;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Creates the Jacobian matrix, the appropriate evaluation formats
  --   and calls the standard corrector on every solution in the list.

    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;
    p_eval : Eval_Poly_Sys(p'range) := Create(p);
    jac : Jaco_Mat(p'range,p'range) := Create(p);
    jac_eval : Eval_Jaco_Mat(p'range,p'range) := Create(jac);

  begin
    while not Is_Null(tmp) loop
      declare
        ls : Link_to_Solution := Head_Of(tmp);
      begin
        Call_Standard_Corrector(p_eval,jac_eval,cp,wout,ls);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Clear(p_eval); Clear(jac); Clear(jac_eval);
  end Standard_Corrector;

  procedure Multprec_Corrector
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                cp : in Multprec_Continuation_Data.Corr_Pars;
                wout : in boolean;
                sols : in out Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Creates the Jacobian matrix, the appropriate evaluation formats
  --   and calls the multprec corrector on every solution in the list.

    use Multprec_Complex_Poly_SysFun;
    use Multprec_Complex_Jaco_Matrices;
    use Multprec_Complex_Solutions;

    tmp : Solution_List := sols;
    p_eval : Eval_Poly_Sys(p'range) := Create(p);
    jac : Jaco_Mat(p'range,p'range) := Create(p);
    jac_eval : Eval_Jaco_Mat(p'range,p'range) := Create(jac);

  begin
    while not Is_Null(tmp) loop
      declare
        ls : Link_to_Solution := Head_Of(tmp);
      begin
        Call_Multprec_Corrector(p_eval,jac_eval,cp,wout,ls);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Clear(p_eval); Clear(jac); Clear(jac_eval);
  end Multprec_Corrector;

  procedure Test_Standard_Corrector is

  -- DESCRIPTION : 
  --   Reads in a polynomial system and a list of solutions;
  --   after setting parameters, applies the correct to them.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Solutions;
    use Standard_Continuation_Data;
    use Continuation_Parameters;

    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    permag : double_float := 0.0;
    ans : character;
    wout : boolean;
    cp : Corr_Pars;

  begin
    Standard_System_and_Solutions_io.get(lp,sols);
    new_line;
    put("Give the perturbation magnitude : "); get(permag);
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    wout := (ans = 'y');
    cp.epsrx := 10.0**(-14);
    cp.epsax := 10.0**(-14);
    cp.epsrf := 10.0**(-14);
    cp.epsaf := 10.0**(-14);
    cp.maxit := 4;
    cp.maxtot := cp.maxit*Length_Of(sols);
    Standard_Perturb_Solutions(sols,permag);
    put_line("The perturbed solutions : "); put(sols);
    Standard_Corrector(lp.all,cp,wout,sols);
    put_line("The corrected solutions : "); put(sols);
  end Test_Standard_Corrector;

  procedure Test_Multprec_Corrector is

  -- DESCRIPTION : 
  --   Reads in a polynomial system and a list of solutions;
  --   after setting parameters, applies the correct to them.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Solutions;
    use Multprec_Continuation_Data;

    mplp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    mpsols : Multprec_Complex_Solutions.Solution_List;
    size : natural32 := 0;
    permag : double_float := 0.0;
    ans : character;
    wout : boolean;
    cp : Corr_Pars;

  begin
    Multprec_System_and_Solutions_io.get(mplp,mpsols);
    new_line;
    put("Give the size of the numbers : "); get(size);
    put("Give the perturbation magnitude : "); get(permag);
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    wout := (ans = 'y');
    Set_Size(mplp.all,size);
    Multprec_Complex_Solutions.Set_Size(mpsols,size);
    cp.epsrx := Create(10.0**((-8)*integer(size)));
    cp.epsax := Create(10.0**((-8)*integer(size)));
    cp.epsrf := Create(10.0**((-8)*integer(size)));
    cp.epsaf := Create(10.0**((-8)*integer(size)));
    cp.maxit := 4;
    cp.maxtot := cp.maxit*Multprec_Complex_Solutions.Length_Of(mpsols);
    Multprec_Perturb_Solutions(mpsols,size,permag);
    put_line("The perturbed solutions : "); put(mpsols);
    Multprec_Corrector(mplp.all,cp,wout,mpsols);
    put_line("The corrected solutions : "); put(mpsols);
  end Test_Multprec_Corrector;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing the predictor and corrector methods.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program.");
      put_line("  1. test tangent predictor with standard arithmetic");
      put_line("  2.                        with double double arithmetic");
      put_line("  3.                        with quad double arithmetic");
      put_line("  4.                        with multi-precision arithmetic");
      put_line("  5. test corrector with standard arithmetic");
      put_line("  6.                with multi-precision arithmetic");
      put("Type 0, 1, 2, 3, 4, 5, or 6 to select : ");
      Ask_Alternative(ans,"0123456");
      exit when ans = '0';
      case ans is
        when '1' => Test_Standard_Tangent_Predictor;
        when '2' => Test_DoblDobl_Tangent_Predictor;
        when '3' => Test_QuadDobl_Tangent_Predictor;
        when '4' => Test_Multprec_Tangent_Predictor;
        when '5' => Test_Standard_Corrector;
        when '6' => Test_Multprec_Corrector;
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_preco;
