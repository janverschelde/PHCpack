with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;
with Standard_Random_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Series;
with Standard_Complex_Series_Vectors;
with Standard_CSeries_Polynomials;
with Standard_CSeries_Poly_Functions;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_CSeries_Polynomials;
with DoblDobl_CSeries_Poly_Functions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_CSeries_Polynomials;
with QuadDobl_CSeries_Poly_Functions;
with QuadDobl_Speelpenning_Convolutions;
with Series_Polynomial_Gradients;        use Series_Polynomial_Gradients;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Hessian_Convolution_Circuits;
with Standard_Complex_Circuits;
with DoblDobl_Complex_Circuits;
with QuadDobl_Complex_Circuits;

procedure ts_hesscnv is

-- DESCRIPTION :
--   Tests the computation of the second derivatives of a polynomial
--   stored as a convolution circuit.

  procedure Standard_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power
  --   to generate a random circuit for testing in double precision.

    use Standard_Complex_Series;
    use Standard_Complex_Series_Vectors;
    use Standard_CSeries_Polynomials;
    use Standard_Speelpenning_Convolutions;

    c : constant Circuit(nbr,dim,dim-1,dim-2)
      := Standard_Random_Convolution_Circuit(dim,deg,nbr,pwr);
    p : constant Poly := Standard_Polynomial(c);
    q,h : Poly;
    d : Standard_Complex_Numbers.Complex_Number;
    x : constant Standard_Complex_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    sx : Standard_Complex_Series_Vectors.Vector(1..dim);
    idx1,idx2 : integer32 := 0;
    eva : Standard_Complex_Series.Link_to_Series;

  begin
    put_line("The exponents : ");
    for i in c.xps'range loop
      put(c.xps(i).all); new_line;
    end loop;
    put("Give the first index : "); get(idx1);
    put("Give the second index : "); get(idx2);
    d := Diff(c,x,idx1,idx2);
    put("The value of the second derivative w.r.t. ");
    put(idx1,1); put(" and "); put(idx2,1); put_line(" :");
    put(d); new_line;
    q := Diff(p,idx1);
    h := Diff(q,idx2);
    for k in sx'range loop
      sx(k) := new Series'(Create(x(k),deg));
    end loop;
    eva := Standard_CSeries_Poly_Functions.Eval(h,sx);
    put_line("The value computed via symbolic differentiation :");
    put(eva.cff(0)); new_line;
  end Standard_Test;

  function Make ( c : Standard_Speelpenning_Convolutions.Circuit )
                return Standard_Complex_Circuits.Circuit is

  -- DESCRIPTION :
  --   Returns the equivalent circuit, using the leading coefficients
  --   in c as the coefficients in the returned circuit.

    res : Standard_Complex_Circuits.Circuit(c.nbr)
        := Standard_Complex_Circuits.Allocate(c.nbr,c.dim);

    use Standard_Complex_Vectors;

  begin
    res.dim := c.dim;
    res.xps := c.xps;
    res.idx := c.idx;
    res.fac := c.fac;
    for k in c.cff'range loop
      res.cff(k) := c.cff(k)(0);
    end loop;
    if c.cst = null
     then res.cst := Standard_Complex_Numbers.Create(0.0);
     else res.cst := c.cst(0);
    end if;
    return res;
  end Make;

  function Make ( c : DoblDobl_Speelpenning_Convolutions.Circuit )
                return DoblDobl_Complex_Circuits.Circuit is

  -- DESCRIPTION :
  --   Returns the equivalent circuit, using the leading coefficients
  --   in c as the coefficients in the returned circuit.

    res : DoblDobl_Complex_Circuits.Circuit(c.nbr)
        := DoblDobl_Complex_Circuits.Allocate(c.nbr,c.dim);

    use DoblDobl_Complex_Vectors;

  begin
    res.dim := c.dim;
    res.xps := c.xps;
    res.idx := c.idx;
    res.fac := c.fac;
    for k in c.cff'range loop
      res.cff(k) := c.cff(k)(0);
    end loop;
    if c.cst = null
     then res.cst := DoblDobl_Complex_Numbers.Create(integer(0));
     else res.cst := c.cst(0);
    end if;
    return res;
  end Make;

  function Make ( c : QuadDobl_Speelpenning_Convolutions.Circuit )
                return QuadDobl_Complex_Circuits.Circuit is

  -- DESCRIPTION :
  --   Returns the equivalent circuit, using the leading coefficients
  --   in c as the coefficients in the returned circuit.

    res : QuadDobl_Complex_Circuits.Circuit(c.nbr)
        := QuadDobl_Complex_Circuits.Allocate(c.nbr,c.dim);

    use QuadDobl_Complex_Vectors;

  begin
    res.dim := c.dim;
    res.xps := c.xps;
    res.idx := c.idx;
    res.fac := c.fac;
    for k in c.cff'range loop
      res.cff(k) := c.cff(k)(0);
    end loop;
    if c.cst = null
     then res.cst := QuadDobl_Complex_Numbers.Create(integer(0));
     else res.cst := c.cst(0);
    end if;
    return res;
  end Make;

  procedure Standard_SVD_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power
  --   to generate a random circuit for testing in double precision.
  --   This test checks whether the singular values computed by the
  --   new Complex_Circuits match the singular values computed by
  --   the procedures in the Hession_Convolution_Circuits.

    use Standard_Speelpenning_Convolutions;

    c : constant Standard_Speelpenning_Convolutions.Circuit
                   (nbr,dim,dim-1,dim-2)
      := Standard_Random_Convolution_Circuit(dim,deg,nbr,pwr);
    c2 : constant Standard_Complex_Circuits.Circuit := Make(c);
    x : constant Standard_Complex_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    xv : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(x);
    y : constant Standard_Complex_Vectors.Vector(0..dim)
      := (0..dim => Standard_Complex_Numbers.Create(0.0));
    yd : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(y);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim) := (1..dim => pwr);
    pwt : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Standard_Complex_Circuits.Allocate(mxe);
    A,U,V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    e : Standard_Complex_Vectors.Vector(1..dim);
    s,s2 : Standard_Complex_Vectors.Vector(1..dim+1);

  begin
    Hessian_Convolution_Circuits.Singular_Values(c,x,A,U,V,e,s);
    put_line("The singular values :"); put_line(s(1..dim));
    Standard_Complex_Circuits.Power_Table(mxe,xv,pwt);
    Standard_Complex_Circuits.Singular_Values(c2,xv,yd,pwt,A,U,V,e,s2);
    put_line("Recomputed singular values :"); put_line(s2(1..dim));
  end Standard_SVD_Test;

  procedure DoblDobl_SVD_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power to
  --   generate a random circuit for testing in double double precision.
  --   This test checks whether the singular values computed by the
  --   new Complex_Circuits match the singular values computed by
  --   the procedures in the Hession_Convolution_Circuits.

    use DoblDobl_Speelpenning_Convolutions;

    c : constant DoblDobl_Speelpenning_Convolutions.Circuit
                   (nbr,dim,dim-1,dim-2)
      := DoblDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
    c2 : constant DoblDobl_Complex_Circuits.Circuit := Make(c);
    x : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    xv : constant DoblDobl_Complex_Vectors.Link_to_Vector
       := new DoblDobl_Complex_Vectors.Vector'(x);
    y : constant DoblDobl_Complex_Vectors.Vector(0..dim)
      := (0..dim => DoblDobl_Complex_Numbers.Create(integer(0)));
    yd : constant DoblDobl_Complex_Vectors.Link_to_Vector
       := new DoblDobl_Complex_Vectors.Vector'(y);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim) := (1..dim => pwr);
    pwt : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
        := DoblDobl_Complex_Circuits.Allocate(mxe);
    A,U,V : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : DoblDobl_Complex_Vectors.Vector(1..dim);
    s,s2 : DoblDobl_Complex_Vectors.Vector(1..dim+1);

  begin
    Hessian_Convolution_Circuits.Singular_Values(c,x,A,U,V,e,s);
    put_line("The singular values :"); put_line(s(1..dim));
    DoblDobl_Complex_Circuits.Power_Table(mxe,xv,pwt);
    DoblDobl_Complex_Circuits.Singular_Values(c2,xv,yd,pwt,A,U,V,e,s2);
    put_line("Recomputed singular values :"); put_line(s2(1..dim));
  end DoblDobl_SVD_Test;

  procedure QuadDobl_SVD_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power to
  --   generate a random circuit for testing in double double precision.
  --   This test checks whether the singular values computed by the
  --   new Complex_Circuits match the singular values computed by
  --   the procedures in the Hession_Convolution_Circuits.

    use QuadDobl_Speelpenning_Convolutions;

    c : constant QuadDobl_Speelpenning_Convolutions.Circuit
                   (nbr,dim,dim-1,dim-2)
      := QuadDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
    c2 : constant QuadDobl_Complex_Circuits.Circuit := Make(c);
    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    xv : constant QuadDobl_Complex_Vectors.Link_to_Vector
       := new QuadDobl_Complex_Vectors.Vector'(x);
    y : constant QuadDobl_Complex_Vectors.Vector(0..dim)
      := (0..dim => QuadDobl_Complex_Numbers.Create(integer(0)));
    yd : constant QuadDobl_Complex_Vectors.Link_to_Vector
       := new QuadDobl_Complex_Vectors.Vector'(y);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim) := (1..dim => pwr);
    pwt : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
        := QuadDobl_Complex_Circuits.Allocate(mxe);
    A,U,V : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : QuadDobl_Complex_Vectors.Vector(1..dim);
    s,s2 : QuadDobl_Complex_Vectors.Vector(1..dim+1);

  begin
    Hessian_Convolution_Circuits.Singular_Values(c,x,A,U,V,e,s);
    put_line("The singular values :"); put_line(s(1..dim));
    QuadDobl_Complex_Circuits.Power_Table(mxe,xv,pwt);
    QuadDobl_Complex_Circuits.Singular_Values(c2,xv,yd,pwt,A,U,V,e,s2);
    put_line("Recomputed singular values :"); put_line(s2(1..dim));
  end QuadDobl_SVD_Test;

  procedure DoblDobl_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power
  --   to generate a random circuit for testing in double double precision.

    use DoblDobl_Complex_Series;
    use DoblDobl_Complex_Series_Vectors;
    use DoblDobl_CSeries_Polynomials;
    use DoblDobl_Speelpenning_Convolutions;

    c : constant Circuit(nbr,dim,dim-1,dim-2)
      := DoblDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
    p : constant Poly := DoblDobl_Polynomial(c);
    q,h : Poly;
    d : DoblDobl_Complex_Numbers.Complex_Number;
    x : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    sx : DoblDobl_Complex_Series_Vectors.Vector(1..dim);
    idx1,idx2 : integer32 := 0;
    eva : DoblDobl_Complex_Series.Link_to_Series;

  begin
    put_line("The exponents : ");
    for i in c.xps'range loop
      put(c.xps(i).all); new_line;
    end loop;
    put("Give the first index : "); get(idx1);
    put("Give the second index : "); get(idx2);
    d := Diff(c,x,idx1,idx2);
    put("The value of the second derivative w.r.t. ");
    put(idx1,1); put(" and "); put(idx2,1); put_line(" :");
    put(d); new_line;
    q := Diff(p,idx1);
    h := Diff(q,idx2);
    for k in sx'range loop
      sx(k) := new Series'(Create(x(k),deg));
    end loop;
    eva := DoblDobl_CSeries_Poly_Functions.Eval(h,sx);
    put_line("The value computed via symbolic differentiation :");
    put(eva.cff(0)); new_line;
  end DoblDobl_Test;

  procedure QuadDobl_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power
  --   to generate a random circuit for testing in quad double precision.

    use QuadDobl_Complex_Series;
    use QuadDobl_Complex_Series_Vectors;
    use QuadDobl_CSeries_Polynomials;
    use QuadDobl_Speelpenning_Convolutions;

    c : constant Circuit(nbr,dim,dim-1,dim-2)
      := QuadDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
    p : constant Poly := QuadDobl_Polynomial(c);
    q,h : Poly;
    d : QuadDobl_Complex_Numbers.Complex_Number;
    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    sx : QuadDobl_Complex_Series_Vectors.Vector(1..dim);
    idx1,idx2 : integer32 := 0;
    eva : QuadDobl_Complex_Series.Link_to_Series;

  begin
    put_line("The exponents : ");
    for i in c.xps'range loop
      put(c.xps(i).all); new_line;
    end loop;
    put("Give the first index : "); get(idx1);
    put("Give the second index : "); get(idx2);
    d := Diff(c,x,idx1,idx2);
    put("The value of the second derivative w.r.t. ");
    put(idx1,1); put(" and "); put(idx2,1); put_line(" :");
    put(d); new_line;
    q := Diff(p,idx1);
    h := Diff(q,idx2);
    for k in sx'range loop
      sx(k) := new Series'(Create(x(k),deg));
    end loop;
    eva := QuadDobl_CSeries_Poly_Functions.Eval(h,sx);
    put_line("The value computed via symbolic differentiation :");
    put(eva.cff(0)); new_line;
  end QuadDobl_Test;

  procedure Standard_Test_System ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power
  --   to generate random circuits for testing in double precision.
  --   Computes the Hessian for every polynomial in the circuit
  --   and compares with the evaluation at the symbolic derivative.

    use Standard_Complex_Numbers;
    use Standard_Complex_Matrices;
    use Standard_Complex_Series;
    use Standard_Complex_Series_Vectors;
    use Standard_CSeries_Polynomials;
    use Standard_Speelpenning_Convolutions;

    c : constant Circuits(1..dim)
      := Standard_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    h,ph : Matrix(1..dim,1..dim);
    x : constant Standard_Complex_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    p,dp,d2p : Poly;
    sx : Standard_Complex_Series_Vectors.Vector(1..dim);
    eva : Standard_Complex_Series.Link_to_Series;
    err : Complex_Number;
    val,errsum : double_float := 0.0;

  begin
    for k in sx'range loop
      sx(k) := new Series'(Create(x(k),deg));
    end loop;
    for k in 1..dim loop
      h := Hessian_Convolution_Circuits.Hessian(c(k),x);
      p := Standard_Polynomial(c(k).all);
      for i in 1..dim loop
        dp := Diff(p,i);
        for j in i..dim loop
          d2p := Diff(dp,j);
          eva := Standard_CSeries_Poly_Functions.Eval(d2p,sx);
          ph(i,j) := eva.cff(0); ph(j,i) := ph(i,j);
          Standard_CSeries_Polynomials.Clear(d2p);
	  err := ph(i,j) - h(i,j); val := AbsVal(err);
          errsum := errsum + val;
        end loop;
        Standard_CSeries_Polynomials.Clear(dp);
      end loop;
      Standard_CSeries_Polynomials.Clear(p);
    end loop;
    put("Sum of errors :"); put(errsum,3); new_line;
  end Standard_Test_System;

  procedure DoblDobl_Test_System ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power
  --   to generate random circuits for testing in double double precision.
  --   Computes the Hessian for every polynomial in the circuit
  --   and compares with the evaluation at the symbolic derivative.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Series;
    use DoblDobl_Complex_Series_Vectors;
    use DoblDobl_CSeries_Polynomials;
    use DoblDobl_Speelpenning_Convolutions;

    c : constant Circuits(1..dim)
      := DoblDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    h,ph : Matrix(1..dim,1..dim);
    x : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    p,dp,d2p : Poly;
    sx : DoblDobl_Complex_Series_Vectors.Vector(1..dim);
    eva : DoblDobl_Complex_Series.Link_to_Series;
    err : Complex_Number;
    val,errsum : double_double := create(0.0);

  begin
    for k in sx'range loop
      sx(k) := new Series'(Create(x(k),deg));
    end loop;
    for k in 1..dim loop
      h := Hessian_Convolution_Circuits.Hessian(c(k),x);
      p := DoblDobl_Polynomial(c(k).all);
      for i in 1..dim loop
        dp := Diff(p,i);
        for j in i..dim loop
          d2p := Diff(dp,j);
          eva := DoblDobl_CSeries_Poly_Functions.Eval(d2p,sx);
          ph(i,j) := eva.cff(0); ph(j,i) := ph(i,j);
          DoblDobl_CSeries_Polynomials.Clear(d2p);
	  err := ph(i,j) - h(i,j); val := AbsVal(err);
          errsum := errsum + val;
        end loop;
        DoblDobl_CSeries_Polynomials.Clear(dp);
      end loop;
      DoblDobl_CSeries_Polynomials.Clear(p);
    end loop;
    put("Sum of errors : "); put(errsum,3); new_line;
  end DoblDobl_Test_System;

  procedure QuadDobl_Test_System ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power
  --   to generate random circuits for testing in quad double precision.
  --   Computes the Hessian for every polynomial in the circuit
  --   and compares with the evaluation at the symbolic derivative.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Series;
    use QuadDobl_Complex_Series_Vectors;
    use QuadDobl_CSeries_Polynomials;
    use QuadDobl_Speelpenning_Convolutions;

    c : constant Circuits(1..dim)
      := QuadDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    h,ph : Matrix(1..dim,1..dim);
    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    p,dp,d2p : Poly;
    sx : QuadDobl_Complex_Series_Vectors.Vector(1..dim);
    eva : QuadDobl_Complex_Series.Link_to_Series;
    err : Complex_Number;
    val,errsum : quad_double := create(0.0);

  begin
    for k in sx'range loop
      sx(k) := new Series'(Create(x(k),deg));
    end loop;
    for k in 1..dim loop
      h := Hessian_Convolution_Circuits.Hessian(c(k),x);
      p := QuadDobl_Polynomial(c(k).all);
      for i in 1..dim loop
        dp := Diff(p,i);
        for j in i..dim loop
          d2p := Diff(dp,j);
          eva := QuadDobl_CSeries_Poly_Functions.Eval(d2p,sx);
          ph(i,j) := eva.cff(0); ph(j,i) := ph(i,j);
          QuadDobl_CSeries_Polynomials.Clear(d2p);
	  err := ph(i,j) - h(i,j); val := AbsVal(err);
          errsum := errsum + val;
        end loop;
        QuadDobl_CSeries_Polynomials.Clear(dp);
      end loop;
      QuadDobl_CSeries_Polynomials.Clear(p);
    end loop;
    put("Sum of errors : "); put(errsum,3); new_line;
  end QuadDobl_Test_System;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension, degree, number of terms,
  --   and largest power to generate a random circuit.

    dim,deg,nbr,pwr : integer32 := 0;
    ans,precision : character;

  begin
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    put("Give the number of terms : "); get(nbr);
    put("Give the largest power : "); get(pwr);
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(precision,"012");
    new_line;
    put("Test system ? (y/n) ");  Ask_Yes_or_No(ans);
    if ans = 'y' then
      case precision is
        when '0' => Standard_Test_System(dim,deg,nbr,pwr);
        when '1' => DoblDobl_Test_System(dim,deg,nbr,pwr);
        when '2' => QuadDobl_Test_System(dim,deg,nbr,pwr);
        when others => null;
      end case;
    else
      new_line;
      put("Test Singular Value Decomposition ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        case precision is
          when '0' => Standard_SVD_Test(dim,deg,nbr,pwr);
          when '1' => DoblDobl_SVD_Test(dim,deg,nbr,pwr);
          when '2' => QuadDobl_SVD_Test(dim,deg,nbr,pwr);
          when others => null;
        end case;
      else
        case precision is
          when '0' => Standard_Test(dim,deg,nbr,pwr);
          when '1' => DoblDobl_Test(dim,deg,nbr,pwr);
          when '2' => QuadDobl_Test(dim,deg,nbr,pwr);
          when others => null;
        end case;
      end if;
    end if;
  end Main;

begin
  Main;
end ts_hesscnv;
