with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Standard_Complex_Laur_Functions;    use Standard_Complex_Laur_Functions;

package body Polyhedral_Coefficient_Correctors is

  procedure Silent_Newton_Step
              ( hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                x,y : in out Standard_Complex_Vectors.Vector;
                norm_dx,norm_y : out double_float ) is

    jm : Matrix(y'range,x'range);
    nrm : double_float;
    ipvt : Standard_Integer_Vectors.Vector(x'range);
    info : integer32;

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        jm(i,j) := Eval(jacmat(i,j),mulfac(i,j).all,ctm(i).all,x);
      end loop;
    end loop;
    lufac(jm,x'last,ipvt,info);
    Standard_Complex_Vectors.Min(y);
    lusolve(jm,x'last,ipvt,y);
    nrm := Max_Norm(y); norm_dx := nrm;
    Standard_Complex_Vectors.Add(x,y);
    y := Eval(hq,ctm,x);
    nrm := Max_Norm(y); norm_y := nrm;
  end Silent_Newton_Step;

  procedure Silent_Newton_Step
              ( hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                x,y : in out Standard_Complex_Vectors.Vector;
                norm_dx,norm_y,rcond : out double_float ) is

    jm : Matrix(y'range,x'range);
    nrm : double_float;
    ipvt : Standard_Integer_Vectors.Vector(x'range);

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        jm(i,j) := Eval(jacmat(i,j),mulfac(i,j).all,ctm(i).all,x);
      end loop;
    end loop;
    lufco(jm,x'last,ipvt,rcond);
    Standard_Complex_Vectors.Min(y);
    lusolve(jm,x'last,ipvt,y);
    nrm := Max_Norm(y); norm_dx := nrm;
    Standard_Complex_Vectors.Add(x,y);
    y := Eval(hq,ctm,x);
    nrm := Max_Norm(y); norm_y := nrm;
  end Silent_Newton_Step;

  procedure Reporting_Newton_Step
              ( file : in file_type; hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                x,y : in out Standard_Complex_Vectors.Vector;
                norm_dx,norm_y : out double_float ) is

    jm : Matrix(y'range,x'range);
    nrm : double_float;
    ipvt : Standard_Integer_Vectors.Vector(x'range);
    info : integer32;

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        jm(i,j) := Eval(jacmat(i,j),mulfac(i,j).all,ctm(i).all,x);
      end loop;
    end loop;
    lufac(jm,x'last,ipvt,info);
    Standard_Complex_Vectors.Min(y);
    lusolve(jm,x'last,ipvt,y);
    nrm := Max_Norm(y); norm_dx := nrm;
    put(file,"  |dx| ="); put(file,nrm,2); 
    Standard_Complex_Vectors.Add(x,y);
    y := Eval(hq,ctm,x);
    nrm := Max_Norm(y); norm_y := nrm;
    put(file,"  |f(x)| ="); put(file,nrm,2);
  end Reporting_Newton_Step;

  procedure Reporting_Newton_Step
              ( file : in file_type; hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                x,y : in out Standard_Complex_Vectors.Vector;
                norm_dx,norm_y,rcond : out double_float ) is

    jm : Matrix(y'range,x'range);
    nrm : double_float;
    ipvt : Standard_Integer_Vectors.Vector(x'range);

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        jm(i,j) := Eval(jacmat(i,j),mulfac(i,j).all,ctm(i).all,x);
      end loop;
    end loop;
    lufco(jm,x'last,ipvt,rcond);
    Standard_Complex_Vectors.Min(y);
    lusolve(jm,x'last,ipvt,y);
    nrm := Max_Norm(y); norm_dx := nrm;
    put(file,"  |dx| ="); put(file,nrm,2); 
    Standard_Complex_Vectors.Add(x,y);
    y := Eval(hq,ctm,x);
    nrm := Max_Norm(y); norm_y := nrm;
    put(file,"  |f(x)| ="); put(file,nrm,2);
  end Reporting_Newton_Step;

  procedure Silent_Apply_Newton
              ( hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                tol : in double_float; max : in natural32;
                x,y : in out Standard_Complex_Vectors.Vector;
                nit : out natural32; fail : out boolean ) is

    prev_ndx,prev_ny,new_ndx,new_ny : double_float;

  begin
    nit := 1;
    Silent_Newton_Step(hq,ctm,jacmat,mulfac,x,y,prev_ndx,prev_ny);
    fail := (prev_ndx > tol) and (prev_ny > tol);
    while nit < max loop
      nit := nit + 1;
      Silent_Newton_Step(hq,ctm,jacmat,mulfac,x,y,new_ndx,new_ny);
      exit when (new_ndx > prev_ndx) or (new_ny > prev_ny);  -- divergence
      fail := (new_ndx > tol) and (new_ny > tol);
      exit when not fail;
      prev_ndx := new_ndx; prev_ny := new_ny;
    end loop;
  end Silent_Apply_Newton;

  procedure Reporting_Apply_Newton
              ( file : in file_type; hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                tol : in double_float; max : in natural32;
                x,y : in out Standard_Complex_Vectors.Vector;
                nit : out natural32; fail : out boolean ) is

    prev_ndx,prev_ny,new_ndx,new_ny : double_float;

  begin
    nit := 1;
    put(file,nit,3); put(file," :");
    Reporting_Newton_Step(file,hq,ctm,jacmat,mulfac,x,y,prev_ndx,prev_ny);
    new_line(file);
    fail := (prev_ndx > tol) and (prev_ny > tol);
    while nit < max loop
      nit := nit + 1;
      put(file,nit,3); put(file," :");
      Reporting_Newton_Step(file,hq,ctm,jacmat,mulfac,x,y,new_ndx,new_ny);
      exit when (new_ndx > prev_ndx) or (new_ny > prev_ny);  -- divergence
      fail := (new_ndx > tol) and (new_ny > tol);
      exit when not fail;
      prev_ndx := new_ndx; prev_ny := new_ny;
      exit when (nit = max);
      new_line(file);
    end loop;
    if fail
     then put_line(file,"  failure");
     else put_line(file,"  success");
    end if;
  end Reporting_Apply_Newton;

end Polyhedral_Coefficient_Correctors;
