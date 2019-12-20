with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Integer_VecVecs_io;
with Standard_Random_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;         use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;         use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;         use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with Exponent_Indices;
with Standard_Dense_Series;
with Standard_Dense_Series_Vectors;
with DoblDobl_Dense_Series;
with DoblDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series;
with QuadDobl_Dense_Series_Vectors;
with Standard_Series_Polynomials;
with Standard_Series_Poly_Functions;
with DoblDobl_Series_Polynomials;
with DoblDobl_Series_Poly_Functions;
with QuadDobl_Series_Polynomials;
with QuadDobl_Series_Poly_Functions;
with Standard_Random_Series;
with DoblDobl_Random_Series;
with QuadDobl_Random_Series;
with Series_and_Polynomials_io;           use Series_and_Polynomials_io;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Series_Polynomial_Gradients;         use Series_Polynomial_Gradients;

procedure ts_speelcnv is

-- DESCRIPTION :
--   Tests the evaluation of the gradient of a polynomial in many variables,
--   in a power series of some fixed degree.

  function Random_Exponents
             ( dim,nbr : integer32 ) 
             return Standard_Integer_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Generates as many random 0/1 vectors of dimension dim
  --   as the value of nbr.  All vectors have a nonzero sum.

    res : Standard_Integer_VecVecs.VecVec(1..nbr);
    nz : integer32;

  begin
    for i in 1..nbr loop
      loop
        declare
          xp : constant Standard_Integer_Vectors.Vector(1..dim)
             := Standard_Random_Vectors.Random_Vector(1,dim,0,1);
        begin
          nz := Standard_Integer_Vectors.Sum(xp);
          if nz > 0
           then res(i) := new Standard_Integer_Vectors.Vector'(xp);
          end if;
        end;
        exit when (nz > 0);
      end loop;
    end loop;
    return res;
  end Random_Exponents;

  procedure Standard_Test ( dim,deg,nbr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and tests the
  --   evaluation and differentiation in standard double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products.

    use Standard_Speelpenning_Convolutions;

    xps : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Random_Exponents(dim,nbr);
    idx : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Exponent_Indices.Exponent_Index(xps);
    polcff : constant Standard_Dense_Series_Vectors.Vector(1..nbr)
           := Standard_Random_Series.Random_Series_Vector(1,nbr,deg);
    pol : constant Standard_Series_Polynomials.Poly
       -- := Standard_Polynomial(dim,deg,idx); -- all coefficients are one
        := Standard_Polynomial(dim,idx,polcff);
    x : constant Standard_Dense_Series_Vectors.Vector(1..dim)
      := Standard_Random_Series.Random_Series_Vector(1,dim,deg);
    y : Standard_Dense_Series.Series;
    grad : Standard_Dense_Series_Vectors.Vector(1..dim);
    xcff : constant Standard_Complex_VecVecs.VecVec(1..dim)
         := Standard_Series_Coefficients(x);
    pcff : constant Standard_Complex_VecVecs.VecVec(1..nbr)
         := Standard_Series_Coefficients(polcff);
    forward : Standard_Complex_VecVecs.VecVec(1..dim-1)
            := Allocate_Coefficients(dim-1,deg);
    backward : Standard_Complex_VecVecs.VecVec(1..dim-2)
             := Allocate_Coefficients(dim-2,deg);
    cross : Standard_Complex_VecVecs.VecVec(1..dim-2)
          := Allocate_Coefficients(dim-2,deg);
    ygrad : Standard_Complex_VecVecs.VecVec(1..dim+1)
          := Allocate_Coefficients(dim+1,deg);
    work : constant Standard_Complex_Vectors.Link_to_Vector
         := Allocate_Coefficients(deg);

  begin
    put_line("Some random exponents :");
    Standard_Integer_VecVecs_io.put(xps);
    put_line("its exponent indices :");
    Standard_Integer_VecVecs_io.put(idx);
    put_line("the polynomial :"); put(pol); new_line;
    y := Standard_Series_Poly_Functions.Eval(pol,x);
   -- Speel(idx,xcff,forward,backward,cross,ygrad); -- if all coefficients one
    Speel(idx,pcff,xcff,forward,backward,cross,ygrad,work);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    put_line("The coefficient vector of the value of the product :");
    put_line(ygrad(ygrad'last));
    grad := Standard_Gradient(pol,x);
    for k in grad'range loop
      put("derivative "); put(k,1); put_line(" :");
      put(grad(k)); new_line;
      put("The coefficient vector of derivative ");
      put(k,1); put_line(" :"); put_line(ygrad(k));
    end loop;
  end Standard_Test;

  procedure DoblDobl_Test ( dim,deg,nbr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and tests the
  --   evaluation and differentiation in double double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products.

    use DoblDobl_Speelpenning_Convolutions;

    xps : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Random_Exponents(dim,nbr);
    idx : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Exponent_Indices.Exponent_Index(xps);
    polcff : constant DoblDobl_Dense_Series_Vectors.Vector(1..nbr)
           := DoblDobl_Random_Series.Random_Series_Vector(1,nbr,deg);
    pol : constant DoblDobl_Series_Polynomials.Poly
       -- := DoblDobl_Polynomial(dim,deg,idx); -- all coefficients are one
        := DoblDobl_Polynomial(dim,idx,polcff);
    x : constant DoblDobl_Dense_Series_Vectors.Vector(1..dim)
      := DoblDobl_Random_Series.Random_Series_Vector(1,dim,deg);
    y : DoblDobl_Dense_Series.Series;
    grad : DoblDobl_Dense_Series_Vectors.Vector(1..dim);
    xcff : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
         := DoblDobl_Series_Coefficients(x);
    pcff : constant DoblDobl_Complex_VecVecs.VecVec(1..nbr)
         := DOblDobl_Series_Coefficients(polcff);
    forward : DoblDobl_Complex_VecVecs.VecVec(1..dim-1)
            := Allocate_Coefficients(dim-1,deg);
    backward : DoblDobl_Complex_VecVecs.VecVec(1..dim-2)
             := Allocate_Coefficients(dim-2,deg);
    cross : DoblDobl_Complex_VecVecs.VecVec(1..dim-2)
          := Allocate_Coefficients(dim-2,deg);
    ygrad : DoblDobl_Complex_VecVecs.VecVec(1..dim+1)
          := Allocate_Coefficients(dim+1,deg);
    work : constant DoblDobl_Complex_Vectors.Link_to_Vector
         := Allocate_Coefficients(deg);

  begin
    put_line("Some random exponents :");
    Standard_Integer_VecVecs_io.put(xps);
    put_line("its exponent indices :");
    Standard_Integer_VecVecs_io.put(idx);
    put_line("the polynomial :"); put(pol); new_line;
    y := DoblDobl_Series_Poly_Functions.Eval(pol,x);
   -- Speel(idx,xcff,forward,backward,cross,ygrad); -- if all coefficients one
    Speel(idx,pcff,xcff,forward,backward,cross,ygrad,work);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    put_line("The coefficient vector of the value of the product :");
    put_line(ygrad(ygrad'last));
    grad := DoblDobl_Gradient(pol,x);
    for k in grad'range loop
      put("derivative "); put(k,1); put_line(" :");
      put(grad(k)); new_line;
      put("The coefficient vector of derivative ");
      put(k,1); put_line(" :"); put_line(ygrad(k));
    end loop;
  end DoblDobl_Test;

  procedure QuadDobl_Test ( dim,deg,nbr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and tests the
  --   evaluation and differentiation in quad double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products.

    use QuadDobl_Speelpenning_Convolutions;

    xps : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Random_Exponents(dim,nbr);
    idx : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Exponent_Indices.Exponent_Index(xps);
    polcff : constant QuadDobl_Dense_Series_Vectors.Vector(1..nbr)
           := QuadDobl_Random_Series.Random_Series_Vector(1,nbr,deg);
    pol : constant QuadDobl_Series_Polynomials.Poly
       -- := QuadDobl_Polynomial(dim,deg,idx); -- all coefficients are one
        := QuadDobl_Polynomial(dim,idx,polcff);
    x : constant QuadDobl_Dense_Series_Vectors.Vector(1..dim)
      := QuadDobl_Random_Series.Random_Series_Vector(1,dim,deg);
    y : QuadDobl_Dense_Series.Series;
    grad : QuadDobl_Dense_Series_Vectors.Vector(1..dim);
    xcff : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
         := QuadDobl_Series_Coefficients(x);
    pcff : constant QuadDobl_Complex_VecVecs.VecVec(1..nbr)
         := QuadDobl_Series_Coefficients(polcff);
    forward : QuadDobl_Complex_VecVecs.VecVec(1..dim-1)
            := Allocate_Coefficients(dim-1,deg);
    backward : QuadDobl_Complex_VecVecs.VecVec(1..dim-2)
             := Allocate_Coefficients(dim-2,deg);
    cross : QuadDobl_Complex_VecVecs.VecVec(1..dim-2)
          := Allocate_Coefficients(dim-2,deg);
    ygrad : QuadDobl_Complex_VecVecs.VecVec(1..dim+1)
          := Allocate_Coefficients(dim+1,deg);
    work : constant QuadDobl_Complex_Vectors.Link_to_Vector
         := Allocate_Coefficients(deg);

  begin
    put_line("Some random exponents :");
    Standard_Integer_VecVecs_io.put(xps);
    put_line("its exponent indices :");
    Standard_Integer_VecVecs_io.put(idx);
    put_line("the polynomial :"); put(pol); new_line;
    y := QuadDobl_Series_Poly_Functions.Eval(pol,x);
   -- Speel(idx,xcff,forward,backward,cross,ygrad); -- if all coefficients one
    Speel(idx,pcff,xcff,forward,backward,cross,ygrad,work);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    put_line("The coefficient vector of the value of the product :");
    put_line(ygrad(ygrad'last));
    grad := QuadDobl_Gradient(pol,x);
    for k in grad'range loop
      put("derivative "); put(k,1); put_line(" :");
      put(grad(k)); new_line;
      put("The coefficient vector of derivative ");
      put(k,1); put_line(" :"); put_line(ygrad(k));
    end loop;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree, the dimension,
  --   the number of monomials, and the precision.  Then runs the tests.

    dim,deg,nbr : integer32 := 0;
    precision : character;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    put("Give the number of monomials : "); get(nbr);
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(precision,"012");
    new_line;
    case precision is
      when '0' => Standard_Test(dim,deg,nbr);
      when '1' => DoblDobl_Test(dim,deg,nbr);
      when '2' => QuadDobl_Test(dim,deg,nbr);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_speelcnv;
