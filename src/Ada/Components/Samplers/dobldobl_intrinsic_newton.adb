with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_Polar;     use DoblDobl_Complex_Numbers_Polar;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Vector_Norms;      use DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with DoblDobl_Complex_QR_Least_Squares;  use DoblDobl_Complex_QR_Least_Squares;
with DoblDobl_Complex_Singular_Values;   use DoblDobl_Complex_Singular_Values;
with DoblDobl_Point_Coordinates;         use DoblDobl_Point_Coordinates;
with Process_io;                         use Process_io;

package body DoblDobl_Intrinsic_Newton is

-- EVALUATION OF JACOBI MATRIX :

  function Affine_Eval ( ejf,p : Matrix ) return Matrix is

  -- NOTE : ejf is the evaluated Jacobi matrix at x,
  --   this routine is auxiliary to the next two Affine_Evals.

    res : Matrix(ejf'range(1),1..p'last(2));
    zero : constant double_double := create(0.0);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Create(zero);
        for k in ejf'range(2) loop
          res(i,j) := res(i,j) + ejf(i,k)*p(k,j);
        end loop;
      end loop;
    end loop;
    return res;
  end Affine_Eval;

  function Affine_Eval
             ( jf : Eval_Jaco_Mat; p : Matrix; x : Vector )
             return Matrix is

    eva : constant Matrix(jf'range(1),jf'range(2)) := Eval(jf,x);

  begin
    return Affine_Eval(eva,p);
  end Affine_Eval;

  function Generic_Affine_Eval ( p : Matrix; x : Vector ) return Matrix is

    eva : constant Matrix := jf(x);

  begin
    return Affine_Eval(eva,p);
  end Generic_Affine_Eval;

  function Projective_Eval
             ( jf : Eval_Jaco_Mat; p : Matrix; x : Vector; k : natural32 )
             return Matrix is

    eva : constant Matrix(jf'range(1),jf'range(2)) := Eval(jf,x);
    res : Matrix(jf'range(1),1..p'last(2));
    zero : constant double_double := create(0.0);

  begin
    for i in res'range(1) loop
      for j in 0..integer32(k)-1 loop
        res(i,j+1) := Create(zero);
        for kk in jf'range(2) loop
          res(i,j+1) := res(i,j+1) + eva(i,kk)*p(kk,j);
        end loop;
      end loop;
      for j in integer32(k)+1..p'last(2) loop
        res(i,j) := Create(zero);
        for kk in jf'range(2) loop
          res(i,j) := res(i,j) + eva(i,kk)*p(kk,j);
        end loop;
      end loop;
    end loop;
    return res;
  end Projective_Eval;

  procedure Add ( x : in out Vector; dx : in Vector; k : in natural32 ) is

  -- DESCRIPTION :
  --   Adds the components of the vector dx to x, skipping the k-th 
  --   entry of x.  If k = 0, then this is just the plain Add.

  -- REQUIRED : x'range = 0..x'last, dx'range = 1..x'last.

  begin
    for i in x'first..integer32(k)-1 loop
      x(i) := x(i) + dx(i+1);
    end loop;
    for i in integer32(k)+1..x'last loop
      x(i) := x(i) + dx(i);
    end loop;
  end Add;

 -- procedure Write_Expanded
 --             ( file : in file_type; p : in Matrix; x : Vector ) is
 --
  -- DESCRIPTION :
  --   Writes the expansion of the vector x wrt the generators p to file.
 --
 --   y : constant Vector := Affine_Expand(x,p);
 --
 -- begin
 --   put_line(file,y);
 -- end Write_Expanded;

 -- function Expanded_Vector_Norm
 --             ( p : Matrix; x : Vector ) return double_double is
 --
  -- DESCRIPTION :
  --   Returns the max norm of the vector x expanded wrt the generators in p.
 --
 --   res : double_double;
 --   y : Vector(p'range(1)) := Affine_Expand(x,p);
 --
 -- begin
 --   for i in y'range loop
 --     y(i) := y(i) - p(i,0);
 --   end loop;
 --   res := Max_Norm(y);
 --   return res;
 -- end Expanded_Vector_Norm;

-- TOP NEWTON ROUTINES :

  procedure Affine_LU_Newton
              ( f : in Poly_Sys;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Affine_LU_Newton(ef,jf,p,x,epsax,epsrx,epsaf,epsrf,
                     incax,incrx,resaf,resrf,nit,maxit,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Affine_LU_Newton;

  procedure Affine_LU_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Affine_LU_Newton(file,ef,jf,p,x,epsax,epsrx,epsaf,epsrf,
                     incax,incrx,resaf,resrf,nit,maxit,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Affine_LU_Newton;

  procedure Projective_LU_Newton
              ( f : in Poly_Sys;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Projective_LU_Newton(ef,jf,p,x,k,epsax,epsrx,epsaf,epsrf,
                         incax,incrx,resaf,resrf,nit,maxit,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Projective_LU_Newton;

  procedure Projective_LU_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Projective_LU_Newton(file,ef,jf,p,x,k,epsax,epsrx,epsaf,epsrf,
                         incax,incrx,resaf,resrf,nit,maxit,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Projective_LU_Newton;

  procedure Affine_LU_Newton
              ( f : in Poly_Sys;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                rco : out double_double; fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Affine_LU_Newton(ef,jf,p,x,epsax,epsrx,epsaf,epsrf,
                     incax,incrx,resaf,resrf,nit,maxit,rco,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Affine_LU_Newton;

  procedure Affine_LU_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                rco : out double_double; fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Affine_LU_Newton(file,ef,jf,p,x,epsax,epsrx,epsaf,epsrf,
                     incax,incrx,resaf,resrf,nit,maxit,rco,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Affine_LU_Newton;

  procedure Projective_LU_Newton
              ( f : in Poly_Sys;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                rco : out double_double; fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Projective_LU_Newton(ef,jf,p,x,k,epsax,epsrx,epsaf,epsrf,
                         incax,incrx,resaf,resrf,nit,maxit,rco,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Projective_LU_Newton;

  procedure Projective_LU_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                rco : out double_double; fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Projective_LU_Newton(file,ef,jf,p,x,k,epsax,epsrx,epsaf,epsrf,
                         incax,incrx,resaf,resrf,nit,maxit,rco,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Projective_LU_Newton;

  procedure Affine_QR_Newton
              ( f : in Poly_Sys;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Affine_QR_Newton(ef,jf,p,x,epsax,epsrx,epsaf,epsrf,
                     incax,incrx,resaf,resrf,nit,maxit,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Affine_QR_Newton;

  procedure Affine_QR_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Affine_QR_Newton(file,ef,jf,p,x,epsax,epsrx,epsaf,epsrf,
                     incax,incrx,resaf,resrf,nit,maxit,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Affine_QR_Newton;

  procedure Projective_QR_Newton
              ( f : in Poly_Sys;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Projective_QR_Newton(ef,jf,p,x,k,epsax,epsrx,epsaf,epsrf,
                         incax,incrx,resaf,resrf,nit,maxit,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Projective_QR_Newton;

  procedure Projective_QR_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Projective_QR_Newton(file,ef,jf,p,x,k,epsax,epsrx,epsaf,epsrf,
                         incax,incrx,resaf,resrf,nit,maxit,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Projective_QR_Newton;

  procedure Affine_SV_Newton
              ( f : in Poly_Sys;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Affine_SV_Newton(ef,jf,p,x,epsax,epsrx,epsaf,epsrf,
                     incax,incrx,resaf,resrf,nit,maxit,sv,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Affine_SV_Newton;

  procedure Affine_SV_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Affine_SV_Newton(file,ef,jf,p,x,epsax,epsrx,epsaf,epsrf,
                     incax,incrx,resaf,resrf,nit,maxit,sv,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Affine_SV_Newton;

  procedure Projective_SV_Newton
              ( f : in Poly_Sys;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Projective_SV_Newton(ef,jf,p,x,k,epsax,epsrx,epsaf,epsrf,
                         incax,incrx,resaf,resrf,nit,maxit,sv,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Projective_SV_Newton;

  procedure Projective_SV_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean ) is

    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

  begin
    Projective_SV_Newton(file,ef,jf,p,x,k,epsax,epsrx,epsaf,epsrf,
                         incax,incrx,resaf,resrf,nit,maxit,sv,fail);
    Clear(ef); Clear(jm); Clear(jf);
  end Projective_SV_Newton;

-- BASIC NEWTON ROUTINES :

  procedure Affine_LU_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    m : Matrix(jf'range(1),x'range);
    z : Vector(p'range(1)) := Affine_Expand(x,p);
    y : Vector(f'range) := Eval(f,z);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    info : integer32;
    nx : double_double;
    prev_resaf : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Affine_Eval(jf,p,z);
      lufac(m,y'last,ipvt,info);
      exit when (info /= 0);
      Min(y);
      lusolve(m,y'last,ipvt,y);
      Add(x,y);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      z := Affine_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if nit > 1 then
        if resaf > prev_resaf
         then return;
        end if;
      end if;
      prev_resaf := resaf;
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Affine_LU_Newton;

  procedure Affine_LU_Newton
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    m : Matrix(jf'range(1),x'range);
    z : Vector(p'range(1)) := Affine_Expand(x,p);
    y : Vector(f'range) := Eval(f,z);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    info : integer32;
    nx : double_double;
    prev_resaf : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Affine_Eval(jf,p,z);
      lufac(m,y'last,ipvt,info);
      exit when (info /= 0);
      Min(y);
      lusolve(m,y'last,ipvt,y);
      Add(x,y);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,nit,3); put(file,"  |dx| ="); put(file,incax,3);
      end if;
      z := Affine_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,"  |f(x)| ="); put(file,resaf,3); new_line(file);
      end if;
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if nit > 1 then
        if resaf > prev_resaf
         then return;
        end if;
      end if;
      prev_resaf := resaf;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Affine_LU_Newton;

  procedure Projective_LU_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    m : Matrix(jf'range(1),1..x'last);
    z : Vector(p'range(1)) := Projective_Expand(x,p);
    y : Vector(f'range) := Eval(f,z);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    info : integer32;
    nx : double_double;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Projective_Eval(jf,p,z,k);
      lufac(m,y'last,ipvt,info);
      exit when (info /= 0);
      Min(y);
      lusolve(m,y'last,ipvt,y);
      Add(x,y,k);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      z := Projective_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Projective_LU_Newton;

  procedure Projective_LU_Newton
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
		p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    m : Matrix(jf'range(1),1..x'last);
    z : Vector(p'range(1)) := Projective_Expand(x,p);
    y : Vector(f'range) := Eval(f,z);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    info : integer32;
    nx : double_double;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Projective_Eval(jf,p,z,k);
      lufac(m,y'last,ipvt,info);
      exit when (info /= 0);
      Min(y);
      lusolve(m,y'last,ipvt,y);
      Add(x,y,k);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,nit,3); put(file,"  |dx| ="); put(file,incax,3);
      end if;
      z := Projective_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c) 
       then put(file,"  |f(x)| ="); put(file,resaf,3); new_line(file);
      end if;
      if nx + one /= one 
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Projective_LU_Newton;

  procedure Affine_LU_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                rco : out double_double; fail : out boolean ) is

    m : Matrix(jf'range(1),x'range);
    z : Vector(p'range(1)) := Affine_Expand(x,p);
    y : Vector(f'range) := Eval(f,z);
    ipvt : Standard_Integer_Vectors.Vector(z'range);
    nx : double_double;
    prev_resaf : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Affine_Eval(jf,p,z);
      lufco(m,y'last,ipvt,rco);
      exit when (one + rco = one);
      Min(y);
      lusolve(m,y'last,ipvt,y);
      Add(x,y);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      z := Affine_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if nit > 1 then
        if resaf > prev_resaf
         then return;
        end if;
      end if;
      if nx + one /= one 
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Affine_LU_Newton;

  procedure Affine_LU_Newton
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                rco : out double_double; fail : out boolean ) is

    m : Matrix(jf'range(1),x'range);
    z : Vector(p'range(1)); -- := Affine_Expand(x,p);
    y : Vector(f'range); --:= Eval(f,z);
    ipvt : Standard_Integer_Vectors.Vector(z'range);
    nx : double_double;
    prev_resaf : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    z := Affine_Expand(x,p);  -- otherwise crash ?!
    y := Eval(f,z);
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Affine_Eval(jf,p,z);
      lufco(m,y'last,ipvt,rco);
      exit when (one + rco = one);
      Min(y);
      lusolve(m,y'last,ipvt,y);
      Add(x,y);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c) then
        put(file,nit,3);
        put(file,"  |dx| ="); put(file,incax,3);
        put(file,"  rco ="); put(file,rco,3);
      end if;
      z := Affine_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,"  |f(x)| ="); put(file,resaf,3); new_line(file);
      end if;
      if nx + one /= one 
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if nit > 1 then
        if resaf > prev_resaf
         then return;
        end if;
      end if;
      prev_resaf := resaf;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
 -- exception 
 --   when others => put_line("Exception raised in Affine_LU_Newton ...");
 --                  raise;
  end Affine_LU_Newton;

  procedure Projective_LU_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                rco : out double_double; fail : out boolean ) is

    m : Matrix(jf'range(1),1..x'last);
    z : Vector(p'range(1)) := Projective_Expand(x,p);
    y : Vector(f'range) := Eval(f,z);
    ipvt : Standard_Integer_Vectors.Vector(z'range);
    nx : double_double;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Projective_Eval(jf,p,z,k);
      lufco(m,y'last,ipvt,rco);
      exit when (one + rco = one);
      Min(y);
      lusolve(m,y'last,ipvt,y);
      Add(x,y,k);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      z := Projective_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Projective_LU_Newton;

  procedure Projective_LU_Newton
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
		p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                rco : out double_double; fail : out boolean ) is

    m : Matrix(jf'range(1),1..x'last);
    z : Vector(p'range(1)) := Projective_Expand(x,p);
    y : Vector(f'range) := Eval(f,z);
    ipvt : Standard_Integer_Vectors.Vector(z'range);
    nx : double_double;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Projective_Eval(jf,p,z,k);
      lufco(m,y'last,ipvt,rco);
      exit when (one + rco = one);
      Min(y);
      lusolve(m,y'last,ipvt,y);
      Add(x,y,k);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c) then
        put(file,nit,3);
        put(file,"  |dx| ="); put(file,incax,3);
        put(file,"  rco ="); put(file,rco,3);
      end if;
      z := Projective_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,"  |f(x)| ="); put(file,resaf,3); new_line(file);
      end if;
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Projective_LU_Newton;

  procedure Affine_QR_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    m : Matrix(jf'range(1),x'range);
    z : Vector(p'range(1)) := Affine_Expand(x,p);
    y : Vector(f'range) := Eval(f,z);
    s : Vector(x'range);
    nx : double_double;
    n1 : constant integer32 := m'last(1);
    n2 : constant integer32 := m'last(2);
    rsd,dum : Vector(1..n1);
    zero : constant double_double := create(0.0);
    qraux : Vector(m'range(2)) := (m'range(2) => Create(zero));
    jpvt : Standard_Integer_Vectors.Vector(m'range(2)) := (m'range(2) => 0);
    info : integer32;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Affine_Eval(jf,p,z); Min(y);
      QRD(m,qraux,jpvt,false);
      QRLS(m,n1,n2,qraux,y,dum,dum,s,rsd,dum,110,info);  
      Add(x,s);
      nx := Max_Norm(x);
      incax := Max_Norm(s);
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      z := Affine_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Affine_QR_Newton;

  procedure Affine_QR_Newton
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    m : Matrix(jf'range(1),x'range);
    z : Vector(p'range(1)) := Affine_Expand(x,p);
    y : Vector(f'range) := Eval(f,z);
    s : Vector(x'range);
    nx : double_double;
    n1 : constant integer32 := m'last(1);
    n2 : constant integer32 := m'last(2);
    rsd,dum : Vector(1..n1);
    zero : constant double_double := create(0.0);
    qraux : Vector(m'range(2)) := (m'range(2) => Create(zero));
    jpvt : Standard_Integer_Vectors.Vector(m'range(2)) := (m'range(2) => 0);
    info : integer32;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Affine_Eval(jf,p,z); Min(y);
      QRD(m,qraux,jpvt,false);
      QRLS(m,n1,n2,qraux,y,dum,dum,s,rsd,dum,110,info);  
      Add(x,s);
      nx := Max_Norm(x);
      incax := Max_Norm(s);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,nit,3); put(file,"  |dx| ="); put(file,incax,3);
      end if;
      z := Affine_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,"  |f(x)| ="); put(file,resaf,3); new_line(file);
      end if;
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Affine_QR_Newton;

  procedure Projective_QR_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    m : Matrix(jf'range(1),1..x'last);
    z : Vector(p'range(1)) := Projective_Expand(x,p);
    y : Vector(f'range) := Eval(f,z);
    s : Vector(1..x'last);
    nx : double_double;
    n1 : constant integer32 := m'last(1);
    n2 : constant integer32 := m'last(2);
    rsd,dum : Vector(1..n1);
    zero : constant double_double := create(0.0);
    qraux : Vector(m'range(2)) := (m'range(2) => Create(zero));
    jpvt : Standard_Integer_Vectors.Vector(m'range(2)) := (m'range(2) => 0);
    info : integer32;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Projective_Eval(jf,p,z,k); Min(y);
      QRD(m,qraux,jpvt,false);
      QRLS(m,n1,n2,qraux,y,dum,dum,s,rsd,dum,110,info);  
      Add(x,s,k);
      nx := Max_Norm(x);
      incax := Max_Norm(s);
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      z := Projective_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Projective_QR_Newton;

  procedure Projective_QR_Newton
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
		p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    m : Matrix(jf'range(1),1..x'last);
    z : Vector(p'range(1)) := Projective_Expand(x,p);
    y : Vector(f'range) := Eval(f,z);
    s : Vector(1..x'last);
    nx : double_double;
    n1 : constant integer32 := m'last(1);
    n2 : constant integer32 := m'last(2);
    rsd,dum : Vector(1..n1);
    zero : constant double_double := create(0.0);
    qraux : Vector(m'range(2)) := (m'range(2) => Create(zero));
    jpvt : Standard_Integer_Vectors.Vector(m'range(2)) := (m'range(2) => 0);
    info : integer32;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Projective_Eval(jf,p,z,k); Min(y);
      QRD(m,qraux,jpvt,false);
      QRLS(m,n1,n2,qraux,y,dum,dum,s,rsd,dum,110,info);  
      Add(x,s,k);
      nx := Max_Norm(x);
      incax := Max_Norm(s);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,nit,3); put(file,"  |dx| ="); put(file,incax,3);
      end if;
      z := Projective_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,"  |f(x)| ="); put(file,resaf,3); new_line(file);
      end if;
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Projective_QR_Newton;

  procedure Affine_SV_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean ) is

    m : Matrix(jf'range(1),x'range);
    z : Vector(p'range(1)) := Affine_Expand(x,p);
    y : Vector(f'range) := Eval(f,z);
    w : Vector(x'range);
    nx : double_double;
    n : constant integer32 := m'last(1);
    q : constant integer32 := m'last(2);
    e : Vector(1..q);
    u : Matrix(1..n,1..n);
    v : Matrix(1..q,1..q);
    info : integer32;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Affine_Eval(jf,p,z); Min(y);
      SVD(m,n,q,sv,e,u,v,11,info);
      w := Solve(u,v,sv,y);
      Add(x,w);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      z := Affine_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Affine_SV_Newton;

  procedure Affine_SV_Newton
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean ) is

    m : Matrix(jf'range(1),x'range);
    z : Vector(p'range(1)) := Affine_Expand(x,p);
    y : Vector(f'range) := Eval(f,z);
    w : Vector(x'range);
    nx : double_double;
    n : constant integer32 := m'last(1);
    q : constant integer32 := m'last(2);
    e : Vector(1..q);
    u : Matrix(1..n,1..n);
    v : Matrix(1..q,1..q);
    info : integer32;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Affine_Eval(jf,p,z); Min(y);
      SVD(m,n,q,sv,e,u,v,11,info);
      w := Solve(u,v,sv,y);
      Add(x,w);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,nit,3); put(file,"  |dx| ="); put(file,incax,3);
      end if;
      z := Affine_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c) then
        put(file,"  rco ="); put(file,Radius(sv(x'last)/sv(sv'first)),3);
        put(file,"  |f(x)| ="); put(file,resaf,3); new_line(file);
      end if;
     -- put_line(file,"singular values : ");
     -- put_line(file,sv);
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Affine_SV_Newton;

  procedure Projective_SV_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean ) is

    m : Matrix(jf'range(1),1..x'last);
    z : Vector(p'range(1)) := Projective_Expand(x,p);
    y : Vector(f'range) := Eval(f,z);
    w : Vector(1..x'last);
    nx : double_double;
    n : constant integer32 := m'last(1);
    q : constant integer32 := m'last(2);
    e : Vector(1..q);
    u : Matrix(1..n,1..n);
    v : Matrix(1..q,1..q);
    info : integer32;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Projective_Eval(jf,p,z,k); Min(y);
      SVD(m,n,q,sv,e,u,v,11,info);
      w := Solve(u,v,sv,y);
      Add(x,w,k);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      z := Projective_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Projective_SV_Newton;

  procedure Projective_SV_Newton
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
		p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean ) is

    m : Matrix(jf'range(1),1..x'last);
    z : Vector(p'range(1)) := Projective_Expand(x,p);
    y : Vector(f'range) := Eval(f,z);
    w : Vector(1..x'last);
    nx : double_double;
    n : constant integer32 := m'last(1);
    q : constant integer32 := m'last(2);
    e : Vector(1..q);
    u : Matrix(1..n,1..n);
    v : Matrix(1..q,1..q);
    info : integer32;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Projective_Eval(jf,p,z,k); Min(y);
      SVD(m,n,q,sv,e,u,v,11,info);
      w := Solve(u,v,sv,y);
      Add(x,w,k);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,nit,3); put(file,"  |dx| ="); put(file,incax,3);
      end if;
      z := Projective_Expand(x,p);
      y := Eval(f,z);
      resaf := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c) then
        put(file,"  rco ="); put(file,Radius(sv(x'last)/sv(sv'first)),3);
        put(file,"  |f(x)| ="); put(file,resaf,3); new_line(file);
      end if;
     -- put_line(file,"singular values : ");
     -- put_line(file,sv);
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Projective_SV_Newton;

-- GENERIC VERSIONS :

  procedure Silent_Affine_LU_Newton
              ( n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    m : Matrix(1..integer32(n),x'range);
    z : Vector(p'range(1)) := Affine_Expand(x,p);
    y : Vector(1..integer32(n)) := f(z);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    info : integer32;
    nx : double_double;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

    function Jacobi_Eval is new Generic_Affine_Eval(jf);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Jacobi_Eval(p,z);
      lufac(m,y'last,ipvt,info);
      exit when (info /= 0);
      Min(y);
      lusolve(m,y'last,ipvt,y);
      Add(x,y);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      z := Affine_Expand(x,p);
      y := f(z);
      resaf := Max_Norm(y);
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Silent_Affine_LU_Newton;

  procedure Reporting_Affine_LU_Newton
              ( file : in file_type;
                n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    m : Matrix(1..integer32(n),x'range);
    z : Vector(p'range(1)) := Affine_Expand(x,p);
    y : Vector(1..integer32(n)) := f(z);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    info : integer32;
    nx : double_double;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

    function Jacobi_Eval is new Generic_Affine_Eval(jf);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Jacobi_Eval(p,z);
      lufac(m,y'last,ipvt,info);
      exit when (info /= 0);
      Min(y);
      lusolve(m,y'last,ipvt,y);
      Add(x,y);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,nit,3); put(file,"  |dx| ="); put(file,incax,3);
      end if;
      z := Affine_Expand(x,p);
      y := f(z);
      resaf := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,"  |f(x)| ="); put(file,resaf,3); new_line(file);
      end if;
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Reporting_Affine_LU_Newton;

  procedure Silent_Affine_LU_RCO_Newton
              ( n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                rco : out double_double; fail : out boolean ) is

    m : Matrix(1..integer32(n),x'range);
    z : Vector(p'range(1)) := Affine_Expand(x,p);
    y : Vector(1..integer32(n)) := f(z);
    ipvt : Standard_Integer_Vectors.Vector(z'range);
    nx : double_double;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

    function Jacobi_Eval is new Generic_Affine_Eval(jf);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Jacobi_Eval(p,z);
      lufco(m,y'last,ipvt,rco);
      exit when (one + rco = one);
      Min(y);
      lusolve(m,y'last,ipvt,y);
      Add(x,y);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      z := Affine_Expand(x,p);
      y := f(z);
      resaf := Max_Norm(y);
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Silent_Affine_LU_RCO_Newton;

  procedure Reporting_Affine_LU_RCO_Newton
              ( file : in file_type;
                n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                rco : out double_double; fail : out boolean ) is

    m : Matrix(1..integer32(n),x'range);
    z : Vector(p'range(1)) := Affine_Expand(x,p);
    y : Vector(1..integer32(n)) := f(z);
    ipvt : Standard_Integer_Vectors.Vector(z'range);
    nx : double_double;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

    function Jacobi_Eval is new Generic_Affine_Eval(jf);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Jacobi_Eval(p,z);
      lufco(m,y'last,ipvt,rco);
      exit when (one + rco = one);
      Min(y);
      lusolve(m,y'last,ipvt,y);
      Add(x,y);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c) then
        put(file,nit,3);
        put(file,"  |dx| ="); put(file,incax,3);
        put(file,"  rco ="); put(file,rco,3);
      end if;
      z := Affine_Expand(x,p);
      y := f(z);
      resaf := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c) 
       then put(file,"  |f(x)| ="); put(file,resaf,3); new_line(file);
      end if;
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Reporting_Affine_LU_RCO_Newton;

  procedure Silent_Affine_QR_Newton
              ( n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    m : Matrix(1..integer32(n),x'range);
    z : Vector(p'range(1)) := Affine_Expand(x,p);
    y : Vector(1..integer32(n)) := f(z);
    s : Vector(x'range);
    nx : double_double;
    n1 : constant integer32 := m'last(1);
    n2 : constant integer32 := m'last(2);
    rsd,dum : Vector(1..n1);
    zero : constant double_double := create(0.0);
    qraux : Vector(m'range(2)) := (m'range(2) => Create(zero));
    jpvt : Standard_Integer_Vectors.Vector(m'range(2)) := (m'range(2) => 0);
    info : integer32;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

    function Jacobi_Eval is new Generic_Affine_Eval(jf);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Jacobi_Eval(p,z); Min(y);
      QRD(m,qraux,jpvt,false);
      QRLS(m,n1,n2,qraux,y,dum,dum,s,rsd,dum,110,info);  
      Add(x,s);
      nx := Max_Norm(x);
      incax := Max_Norm(s);
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      z := Affine_Expand(x,p);
      y := f(z);
      resaf := Max_Norm(y);
      if nx + one /= one 
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Silent_Affine_QR_Newton;

  procedure Reporting_Affine_QR_Newton
              ( file : in file_type;
                n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                fail : out boolean ) is

    m : Matrix(1..integer32(n),x'range);
    z : Vector(p'range(1)) := Affine_Expand(x,p);
    y : Vector(1..integer32(n)) := f(z);
    s : Vector(x'range);
    nx : double_double;
    n1 : constant integer32 := m'last(1);
    n2 : constant integer32 := m'last(2);
    rsd,dum : Vector(1..n1);
    zero : constant double_double := create(0.0);
    qraux : Vector(m'range(2)) := (m'range(2) => Create(zero));
    jpvt : Standard_Integer_Vectors.Vector(m'range(2)) := (m'range(2) => 0);
    info : integer32;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

    function Jacobi_Eval is new Generic_Affine_Eval(jf);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Jacobi_Eval(p,z); Min(y);
      QRD(m,qraux,jpvt,false);
      QRLS(m,n1,n2,qraux,y,dum,dum,s,rsd,dum,110,info);  
      Add(x,s);
      nx := Max_Norm(x);
      incax := Max_Norm(s);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,nit,3); put(file,"  |dx| ="); put(file,incax,3);
      end if;
      z := Affine_Expand(x,p);
      y := f(z);
      resaf := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,"  |f(x)| ="); put(file,resaf,3); new_line(file);
      end if;
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Reporting_Affine_QR_Newton;

  procedure Silent_Affine_SV_Newton
              ( n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean ) is

    m : Matrix(1..integer32(n),x'range);
    z : Vector(p'range(1)) := Affine_Expand(x,p);
    y : Vector(1..integer32(n)) := f(z);
    w : Vector(x'range);
    nx : double_double;
    q : constant integer32 := m'last(2);
    e : Vector(1..q);
    u : Matrix(1..integer32(n),1..integer32(n));
    v : Matrix(1..q,1..q);
    info : integer32;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

    function Jacobi_Eval is new Generic_Affine_Eval(jf);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Jacobi_Eval(p,z); Min(y);
      SVD(m,integer32(n),q,sv,e,u,v,11,info);
      w := Solve(u,v,sv,y);
      Add(x,w);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      z := Affine_Expand(x,p);
      y := f(z);
      resaf := Max_Norm(y);
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Silent_Affine_SV_Newton;

  procedure Reporting_Affine_SV_Newton
              ( file : in file_type;
                n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean ) is

    m : Matrix(1..integer32(n),x'range);
    z : Vector(p'range(1)) := Affine_Expand(x,p);
    y : Vector(1..integer32(n)) := f(z);
    w : Vector(x'range);
    nx : double_double;
    q : constant integer32 := m'last(2);
    e : Vector(1..q);
    u : Matrix(1..integer32(n),1..integer32(n));
    v : Matrix(1..q,1..q);
    info : integer32;
    prev_incax : double_double := create(1.0E+8);
    one : constant double_double := create(1.0);

    function Jacobi_Eval is new Generic_Affine_Eval(jf);

  begin
    fail := true;
    nit := 0;
    while nit < maxit loop
      nit := nit + 1;
      m := Jacobi_Eval(p,z); Min(y);
      SVD(m,integer32(n),q,sv,e,u,v,11,info);
      w := Solve(u,v,sv,y);
      Add(x,w);
      nx := Max_Norm(x);
      incax := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c)
       then put(file,nit,3); put(file,"  |dx| ="); put(file,incax,3);
      end if;
      z := Affine_Expand(x,p);
      y := f(z);
      resaf := Max_Norm(y);
      if Process_io.Contains_Output_Code(Process_io.c) then
        put(file,"  rco ="); put(file,Radius(sv(x'last)/sv(sv'first)),3);
        put(file,"  |f(x)| ="); put(file,resaf,3); new_line(file);
      end if;
     -- put_line(file,"singular values : ");
     -- put_line(file,sv);
      if nx + one /= one
       then incrx := incax/nx; resrf := resaf/nx;
       else incrx := one;      resrf := one;
      end if;
      if nit > 1 then
        if incax > prev_incax
         then return;
        end if;
      end if;
      prev_incax := incax;
      if (((incax <= epsax) or else (resaf <= epsaf)) and then
          ((incrx <= epsrx) or else (resrf <= epsrf)))
       then fail := false; exit;
      end if;
    end loop;
  end Reporting_Affine_SV_Newton;

end DoblDobl_Intrinsic_Newton;
