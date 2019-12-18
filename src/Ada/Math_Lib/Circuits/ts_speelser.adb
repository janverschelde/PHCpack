with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;         use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Dense_Series;
with Standard_Dense_Series_Vectors;
with Standard_Random_Series;
with Standard_Series_Polynomials;
with Series_and_Polynomials_io;           use Series_and_Polynomials_io;
with Standard_Series_Poly_Functions;

procedure ts_speelser is

-- DESCRIPTION :
--   A vectorized version of the computation of Speelpenning products
--   for truncated power series.

  function Standard_Product
             ( dim,deg : in integer32 )
             return Standard_Series_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the product of the first dim variables,
  --   as a polynomial where the coefficients are truncated power series,
  --   truncated to degree deg, with double precision coefficients.

    res : Standard_Series_Polynomials.Poly;
    trm : Standard_Series_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
    trm.dg(1) := 1;
    trm.cf := Standard_Dense_Series.Create(1.0,deg);
    res := Standard_Series_Polynomials.Create(trm);
    Standard_Series_Polynomials.Clear(trm);
    for i in 2..dim loop
      declare
        trmtwo : Standard_Series_Polynomials.Term;
      begin
        trmtwo.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
        trmtwo.dg(i) := 1;
        trmtwo.cf := Standard_Dense_Series.Create(1.0,deg);
        Standard_Series_Polynomials.Mul(res,trmtwo);
        Standard_Series_Polynomials.Clear(trmtwo);
      end;
    end loop;
    return res;
  end Standard_Product;

  function Standard_Allocate_Coefficients
             ( dim,deg : integer32 )
             return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns allocated space for the coefficients of the series
  --   truncated to degree deg.  The vector on return has range 1..dim.

    res : Standard_Complex_VecVecs.VecVec(1..dim);

  begin
    for k in 1..dim loop
      declare
        cff : constant Standard_Complex_Vectors.Vector(0..deg)
            := (0..deg => Standard_Complex_Numbers.Create(0.0));
      begin
        res(k) := new Standard_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Standard_Allocate_Coefficients;

  procedure Standard_Multiply
              ( first,second : in Standard_Complex_Vectors.Link_to_Vector;
                product : in Standard_Complex_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Convolutes the vectors first and second into the product,
  --   corresponding to the multiplication of two series of the same degree.

  -- REQUIRED : first'last = second'last = product'last.

    deg : constant integer32 := first'last;

  begin
    product(0) := first(0)*second(0);
    for k in 1..deg loop
      product(k) := first(0)*second(k);
      for i in 1..k loop
        product(k) := product(k) + first(i)*second(k-i);
      end loop;
    end loop;
  end Standard_Multiply;

  procedure Standard_Evaluate_Product
              ( x : in Standard_Complex_VecVecs.VecVec;
                forward : in out Standard_Complex_VecVecs.VecVec;
                backward : in out Standard_Complex_VecVecs.VecVec;
                cross : in out Standard_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Inline computation of the coefficients of the product
  --   of the series in x, which are all of degree d.
  --   The range of x starts at 1 and ends at 2 or a larger value.

  -- ON ENTRY :
  --   x        a vector with range starting at 1 and ending at 2 or higher;
  --   forward  has space allocated for x'last-1 coefficient vectors,
  --            for series of degree d;
  --   backward has space reserved for x'last-2 coefficient vectors,
  --            for series of degree d.
  --   cross    has space reserved for x'last-2 coefficient vectors,
  --            for series of degree d.

  -- ON RETURN :
  --   forward  accumulates the forward products,
  --            forward(x'last-1) holds the coefficients for the product
  --            of all coefficients in x;
  --   backward accumulates the backward products,
  --            backward(x'last-2) holds the coefficients of the first
  --            partial derivative of the product;
  --   cross    stores the cross products, cross(k) contains the
  --            coefficients of the partial derivative w.r.t. k+1.

  begin
    Standard_Multiply(x(1),x(2),forward(1));
    for k in 3..x'last loop
      Standard_Multiply(forward(k-2),x(k),forward(k-1));
    end loop;
    if x'last > 2 then
      Standard_Multiply(x(x'last),x(x'last-1),backward(1));
      for k in 2..x'last-2 loop
        Standard_Multiply(backward(k-1),x(x'last-k),backward(k));
      end loop;
      Standard_Multiply(x(1),backward(x'last-3),cross(1));
      for k in 2..x'last-3 loop
        Standard_Multiply(forward(k-1),backward(x'last-2-k),cross(k));
      end loop;
      Standard_Multiply(forward(x'last-3),x(x'last),cross(x'last-2));
    end if;
  end Standard_Evaluate_Product;

  function Standard_Gradient
             ( p : Standard_Series_Polynomials.Poly;
               x : Standard_Dense_Series_Vectors.Vector )
             return Standard_Dense_Series_Vectors.Vector is

  -- DESCRIPTION :
  --   Evaluates the gradient of p at x.

    res : Standard_Dense_Series_Vectors.Vector(x'range);

  begin
    for k in x'range loop
      declare
        dpk : Standard_Series_Polynomials.Poly
            := Standard_Series_Polynomials.Diff(p,k);
        val : constant Standard_Dense_Series.Series
            := Standard_Series_Poly_Functions.Eval(dpk,x);
      begin
        Standard_Series_Polynomials.Clear(dpk);
        res(k) := val;
      end;
    end loop;
    return res;
  end Standard_Gradient;

  function Standard_Series_Coefficients
             ( s : Standard_Dense_Series_Vectors.Vector )
             return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns the coefficients of the series in the vector of vectors.
  --   The range of the k-th vector is 0..s(k).deg.

    res : Standard_Complex_VecVecs.VecVec(s'range);

  begin
    for k in s'range loop
      declare
        cff : constant Standard_Complex_Vectors.Vector(0..s(k).deg)
            := s(k).cff(0..s(k).deg);
      begin
        res(k) := new Standard_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Standard_Series_Coefficients;

  procedure Standard_Test ( dim,deg : in integer32 ) is

  -- DESCRIPTION :
  --   Generates random coefficient vectors for as many series as
  --   the value of dim, series truncated to degree deg.

    prd : constant Standard_Series_Polynomials.Poly
        := Standard_Product(dim,deg);
    x : constant Standard_Dense_Series_Vectors.Vector(1..dim)
      := Standard_Random_Series.Random_Series_Vector(1,dim,deg);
    xcff : constant Standard_Complex_VecVecs.VecVec(1..dim)
         := Standard_Series_Coefficients(x);
    y : Standard_Dense_Series.Series;
    grad : Standard_Dense_Series_Vectors.Vector(1..dim);
    forward : Standard_Complex_VecVecs.VecVec(1..dim-1)
            := Standard_Allocate_Coefficients(dim-1,deg);
    backward : Standard_Complex_VecVecs.VecVec(1..dim-2)
             := Standard_Allocate_Coefficients(dim-2,deg);
    cross : Standard_Complex_VecVecs.VecVec(1..dim-2)
          := Standard_Allocate_Coefficients(dim-2,deg);

  begin
    put(dim,1); put(" random series of degree "); put(deg,1);
    put_line(" :"); put(x);
    put("The product of variables : ");
    put(prd); new_line;
    y := Standard_Series_Poly_Functions.Eval(prd,x);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    Standard_Evaluate_Product(xcff,forward,backward,cross);
    put_line("The coefficients of the product : ");
    put_line(forward(dim-1));
    grad := Standard_Gradient(prd,x);
    put_line("the last derivative :"); put(grad(dim)); new_line;
    put_line("coefficients of the last derivative :");
    put_line(forward(dim-2));
    if dim > 2 then
      put_line("the first derivative :"); put(grad(1)); new_line;
      put_line("coefficients of the first derivative :");
      put_line(backward(dim-2));
      for k in 2..dim-1 loop
        put("derivative "); put(k,1); put_line(" :");
        put(grad(k)); new_line;
        put("coefficients of derivative "); put(k,1); put_line(" :");
        put_line(cross(k-1));
      end loop;
    end if;
  end Standard_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree and the dimension.

    dim,deg : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    Standard_Test(dim,deg);
  end Main;

begin
  Main;
end ts_speelser;
