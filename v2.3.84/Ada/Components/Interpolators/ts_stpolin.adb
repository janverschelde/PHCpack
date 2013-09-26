with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Standard_Polynomial_Interpolators;  use Standard_Polynomial_Interpolators;

procedure ts_stpolin is

-- DESCRIPTION :
--   This procedure tests interpolating hypersurfaces.

  procedure Test_Sampler
              ( p : in Poly; m : in natural32; luf : in boolean ) is

  -- DESCRIPTION :
  --   Samples m points of the polynomial p and interpolates with
  --   LU factorization if luf is true.

    vv : constant Standard_Complex_VecVecs.VecVec := Sample(p,m);
    rcond : double_float;
    ip : Poly;

  begin
   -- for i in vv'range loop
   --   put("Point no. "); put(i,1); put_line(" has coordinates :");
   --   put_line(vv(i));
   --   put("-> evaluation : "); put(Eval(p,vv(i).all)); new_line;
   -- end loop;
    if luf then
      Interpolate(p,vv,ip,rcond);
      put("rcond : "); put(rcond,3); new_line;
    else
      ip := Interpolate(p,vv);
    end if;
    put_line("the interpolated polynomial : "); put_line(ip);
    for i in vv'range loop
      put("Eval at "); put(natural32(i),1); put(" : ");
      put(Eval(ip,vv(i).all)); new_line;
    end loop;
    put("The distance between original and interpolated polynomial :");
    put(Distance(p,ip),3); new_line;
  end Test_Sampler;

  procedure Main is

    d,n,m,cff,method : natural32 := 0;
    luf : boolean;
    p : Poly;

  begin
    new_line;
    put_line("Multivariate Polynomial Interpolation in standard arithmetic");
    new_line;
    put("Give the degree : "); get(d);
    put("Give the number of variables : "); get(n);
    put_line("Choose one of the following : ");
    put_line("  1. all coefficients equal to one;");
    put_line("  2. random real coefficients;");
    put_line("  3. random complex coefficients.");
    put("Type 1,2, or 3 to make your choice : "); get(cff);
    put_line("Which method you want to solve the interpolation conditions : ");
    put_line("  1. LU factorization with condition number estimate;");
    put_line("  2. QR factorization followed by least squares.");
    put("Type 1 or 2 to select : "); get(method);
    luf := (method = 1);
    p := Create(d,n,cff);
    put("A random equation of a hypersurface of degree "); put(d,1);
    put(" in "); put(n,1); put_line(" variables :");
    put_line(p);
    new_line;
    put("The number of terms : "); put(Number_of_Terms(d,n),1); new_line;
    put("Give the number of points to sample : "); get(m);
    Test_Sampler(p,m,luf);
  end Main;

begin
  Main;
end ts_stpolin;
