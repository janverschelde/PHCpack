with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;         use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;         use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Random_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;         use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;         use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with Standard_Dense_Series;
with Standard_Dense_Series_Vectors;
with DoblDobl_Dense_Series;
with DoblDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series;
with QuadDobl_Dense_Series_Vectors;
with Standard_Random_Series;
with DoblDobl_Random_Series;
with QuadDobl_Random_Series;
with Standard_Series_Polynomials;
with Standard_Series_Poly_Functions;
with Series_and_Polynomials_io;           use Series_and_Polynomials_io;
with DoblDobl_Series_Polynomials;
with DoblDobl_Series_Poly_Functions;
with QuadDobl_Series_Polynomials;
with QuadDobl_Series_Poly_Functions;
with Exponent_Indices;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

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
  --   truncated to degree deg, with standard double precision coefficients.

    res : Standard_Series_Polynomials.Poly;
    trm : Standard_Series_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 1);
    trm.cf := Standard_Dense_Series.Create(1.0,deg);
    res := Standard_Series_Polynomials.Create(trm);
    Standard_Series_Polynomials.Clear(trm);
    return res;
  end Standard_Product;

  function DoblDobl_Product
             ( dim,deg : in integer32 )
             return DoblDobl_Series_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the product of the first dim variables,
  --   as a polynomial where the coefficients are truncated power series,
  --   truncated to degree deg, with double double precision coefficients.

    res : DoblDobl_Series_Polynomials.Poly;
    trm : DoblDobl_Series_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 1);
    trm.cf := DoblDobl_Dense_Series.Create(1.0,deg);
    res := DoblDobl_Series_Polynomials.Create(trm);
    DoblDobl_Series_Polynomials.Clear(trm);
    return res;
  end DoblDobl_Product;

  function QuadDobl_Product
             ( dim,deg : in integer32 )
             return QuadDobl_Series_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the product of the first dim variables,
  --   as a polynomial where the coefficients are truncated power series,
  --   truncated to degree deg, with quad double precision coefficients.

    res : QuadDobl_Series_Polynomials.Poly;
    trm : QuadDobl_Series_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 1);
    trm.cf := QuadDobl_Dense_Series.Create(1.0,deg);
    res := QuadDobl_Series_Polynomials.Create(trm);
    QuadDobl_Series_Polynomials.Clear(trm);
    return res;
  end QuadDobl_Product;

  function Standard_Product
             ( deg : in integer32;
               xp : Standard_Integer_Vectors.Vector )
             return Standard_Series_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the product of the first dim variables, dim = xp'last,
  --   using the exponent vector in xp,
  --   as a polynomial where the coefficients are truncated power series,
  --   truncated to degree deg, with standard double precision coefficients.

    res : Standard_Series_Polynomials.Poly;
    trm : Standard_Series_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector(1..xp'last);
    for i in xp'range loop
      trm.dg(i) := natural32(xp(i));
    end loop;
    trm.cf := Standard_Dense_Series.Create(1.0,deg);
    res := Standard_Series_Polynomials.Create(trm);
    Standard_Series_Polynomials.Clear(trm);
    return res;
  end Standard_Product;

  function DoblDobl_Product
             ( deg : in integer32;
               xp : Standard_Integer_Vectors.Vector )
             return DoblDobl_Series_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the product of the first dim variables, dim = xp'last,
  --   using the exponent vector in xp,
  --   as a polynomial where the coefficients are truncated power series,
  --   truncated to degree deg, with double double precision coefficients.

    res : DoblDobl_Series_Polynomials.Poly;
    trm : DoblDobl_Series_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector(1..xp'last);
    for i in xp'range loop
      trm.dg(i) := natural32(xp(i));
    end loop;
    trm.cf := DoblDobl_Dense_Series.Create(1.0,deg);
    res := DoblDobl_Series_Polynomials.Create(trm);
    DoblDobl_Series_Polynomials.Clear(trm);
    return res;
  end DoblDobl_Product;

  function QuadDobl_Product
             ( deg : in integer32;
               xp : Standard_Integer_Vectors.Vector )
             return QuadDobl_Series_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the product of the first dim variables, dim = xp'last,
  --   using the exponent vector in xp,
  --   as a polynomial where the coefficients are truncated power series,
  --   truncated to degree deg, with quad double precision coefficients.

    res : QuadDobl_Series_Polynomials.Poly;
    trm : QuadDobl_Series_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector(1..xp'last);
    for i in xp'range loop
      trm.dg(i) := natural32(xp(i));
    end loop;
    trm.cf := QuadDobl_Dense_Series.Create(1.0,deg);
    res := QuadDobl_Series_Polynomials.Create(trm);
    QuadDobl_Series_Polynomials.Clear(trm);
    return res;
  end QuadDobl_Product;

  function Standard_Gradient
             ( p : Standard_Series_Polynomials.Poly;
               x : Standard_Dense_Series_Vectors.Vector )
             return Standard_Dense_Series_Vectors.Vector is

  -- DESCRIPTION :
  --   Evaluates the gradient of p at x, for testing purposes,
  --   in standard double precision.

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

  function DoblDobl_Gradient
             ( p : DoblDobl_Series_Polynomials.Poly;
               x : DoblDobl_Dense_Series_Vectors.Vector )
             return DoblDobl_Dense_Series_Vectors.Vector is

  -- DESCRIPTION :
  --   Evaluates the gradient of p at x, for testing purposes,
  --   in double double precision.

    res : DoblDobl_Dense_Series_Vectors.Vector(x'range);

  begin
    for k in x'range loop
      declare
        dpk : DoblDobl_Series_Polynomials.Poly
            := DoblDobl_Series_Polynomials.Diff(p,k);
        val : constant DoblDobl_Dense_Series.Series
            := DoblDobl_Series_Poly_Functions.Eval(dpk,x);
      begin
        DoblDobl_Series_Polynomials.Clear(dpk);
        res(k) := val;
      end;
    end loop;
    return res;
  end DoblDobl_Gradient;

  function QuadDobl_Gradient
             ( p : QuadDobl_Series_Polynomials.Poly;
               x : QuadDobl_Dense_Series_Vectors.Vector )
             return QuadDobl_Dense_Series_Vectors.Vector is

  -- DESCRIPTION :
  --   Evaluates the gradient of p at x, for testing purposes,
  --   in quad double precision.

    res : QuadDobl_Dense_Series_Vectors.Vector(x'range);

  begin
    for k in x'range loop
      declare
        dpk : QuadDobl_Series_Polynomials.Poly
            := QuadDobl_Series_Polynomials.Diff(p,k);
        val : constant QuadDobl_Dense_Series.Series
            := QuadDobl_Series_Poly_Functions.Eval(dpk,x);
      begin
        QuadDobl_Series_Polynomials.Clear(dpk);
        res(k) := val;
      end;
    end loop;
    return res;
  end QuadDobl_Gradient;

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

  function DoblDobl_Series_Coefficients
             ( s : DoblDobl_Dense_Series_Vectors.Vector )
             return DoblDobl_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns the coefficients of the series in the vector of vectors.
  --   The range of the k-th vector is 0..s(k).deg.

    res : DoblDobl_Complex_VecVecs.VecVec(s'range);

  begin
    for k in s'range loop
      declare
        cff : constant DoblDobl_Complex_Vectors.Vector(0..s(k).deg)
            := s(k).cff(0..s(k).deg);
      begin
        res(k) := new DoblDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end DoblDobl_Series_Coefficients;

  function QuadDobl_Series_Coefficients
             ( s : QuadDobl_Dense_Series_Vectors.Vector )
             return QuadDobl_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns the coefficients of the series in the vector of vectors.
  --   The range of the k-th vector is 0..s(k).deg.

    res : QuadDobl_Complex_VecVecs.VecVec(s'range);

  begin
    for k in s'range loop
      declare
        cff : constant QuadDobl_Complex_Vectors.Vector(0..s(k).deg)
            := s(k).cff(0..s(k).deg);
      begin
        res(k) := new QuadDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end QuadDobl_Series_Coefficients;

  procedure Standard_Test ( dim,deg : in integer32 ) is

  -- DESCRIPTION :
  --   Generates random coefficient vectors for as many series as
  --   the value of dim, series truncated to degree deg,
  --   and runs a test in standard double precision.

    use Standard_Speelpenning_Convolutions;

    prd : constant Standard_Series_Polynomials.Poly
        := Standard_Product(dim,deg);
    x : constant Standard_Dense_Series_Vectors.Vector(1..dim)
      := Standard_Random_Series.Random_Series_Vector(1,dim,deg);
    xcff : constant Standard_Complex_VecVecs.VecVec(1..dim)
         := Standard_Series_Coefficients(x);
    y : Standard_Dense_Series.Series;
    grad : Standard_Dense_Series_Vectors.Vector(1..dim);
    forward : Standard_Complex_VecVecs.VecVec(1..dim-1)
            := Allocate_Coefficients(dim-1,deg);
    backward : Standard_Complex_VecVecs.VecVec(1..dim-2)
             := Allocate_Coefficients(dim-2,deg);
    cross : Standard_Complex_VecVecs.VecVec(1..dim-2)
          := Allocate_Coefficients(dim-2,deg);

  begin
    put(dim,1); put(" random series of degree "); put(deg,1);
    put_line(" :"); put(x);
    put("The product of variables : ");
    put(prd); new_line;
    y := Standard_Series_Poly_Functions.Eval(prd,x);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    Speel(xcff,forward,backward,cross);
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

  procedure DoblDobl_Test ( dim,deg : in integer32 ) is

  -- DESCRIPTION :
  --   Generates random coefficient vectors for as many series as
  --   the value of dim, series truncated to degree deg,
  --   and runs a test in double double precision.

    use DoblDobl_Speelpenning_Convolutions;

    prd : constant DoblDobl_Series_Polynomials.Poly
        := DoblDobl_Product(dim,deg);
    x : constant DoblDobl_Dense_Series_Vectors.Vector(1..dim)
      := DoblDobl_Random_Series.Random_Series_Vector(1,dim,deg);
    xcff : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
         := DoblDobl_Series_Coefficients(x);
    y : DoblDobl_Dense_Series.Series;
    grad : DoblDobl_Dense_Series_Vectors.Vector(1..dim);
    forward : DoblDobl_Complex_VecVecs.VecVec(1..dim-1)
            := Allocate_Coefficients(dim-1,deg);
    backward : DoblDobl_Complex_VecVecs.VecVec(1..dim-2)
             := Allocate_Coefficients(dim-2,deg);
    cross : DoblDobl_Complex_VecVecs.VecVec(1..dim-2)
          := Allocate_Coefficients(dim-2,deg);

  begin
    put(dim,1); put(" random series of degree "); put(deg,1);
    put_line(" :"); put(x);
    put("The product of variables : ");
    put(prd); new_line;
    y := DoblDobl_Series_Poly_Functions.Eval(prd,x);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    Speel(xcff,forward,backward,cross);
    put_line("The coefficients of the product : ");
    put_line(forward(dim-1));
    grad := DoblDobl_Gradient(prd,x);
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
  end DoblDobl_Test;

  procedure QuadDobl_Test ( dim,deg : in integer32 ) is

  -- DESCRIPTION :
  --   Generates random coefficient vectors for as many series as
  --   the value of dim, series truncated to degree deg,
  --   and runs a test in quad double precision.

    use QuadDobl_Speelpenning_Convolutions;

    prd : constant QuadDobl_Series_Polynomials.Poly
        := QuadDobl_Product(dim,deg);
    x : constant QuadDobl_Dense_Series_Vectors.Vector(1..dim)
      := QuadDobl_Random_Series.Random_Series_Vector(1,dim,deg);
    xcff : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
         := QuadDobl_Series_Coefficients(x);
    y : QuadDobl_Dense_Series.Series;
    grad : QuadDobl_Dense_Series_Vectors.Vector(1..dim);
    forward : QuadDobl_Complex_VecVecs.VecVec(1..dim-1)
            := Allocate_Coefficients(dim-1,deg);
    backward : QuadDobl_Complex_VecVecs.VecVec(1..dim-2)
             := Allocate_Coefficients(dim-2,deg);
    cross : QuadDobl_Complex_VecVecs.VecVec(1..dim-2)
          := Allocate_Coefficients(dim-2,deg);

  begin
    put(dim,1); put(" random series of degree "); put(deg,1);
    put_line(" :"); put(x);
    put("The product of variables : ");
    put(prd); new_line;
    y := QuadDobl_Series_Poly_Functions.Eval(prd,x);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    Speel(xcff,forward,backward,cross);
    put_line("The coefficients of the product : ");
    put_line(forward(dim-1));
    grad := QuadDobl_Gradient(prd,x);
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
  end QuadDobl_Test;

  procedure Standard_Indexed_Test
              ( dim,deg,nz : in integer32;
                xp : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Generates random coefficient vectors for as many series as
  --   the value of dim, series truncated to degree deg,
  --   and runs a test in standard double precision,
  --   for a random exponent vector xp of zeros and ones,
  --   where the number of nonzeros nz > 2.

    use Standard_Speelpenning_Convolutions;

    idx : constant Standard_Integer_Vectors.Vector(1..nz)
        := Exponent_Indices.Exponent_Index(xp);
    prd : constant Standard_Series_Polynomials.Poly
        := Standard_Product(deg,xp);
    x : constant Standard_Dense_Series_Vectors.Vector(1..dim)
      := Standard_Random_Series.Random_Series_Vector(1,dim,deg);
    xcff : constant Standard_Complex_VecVecs.VecVec(1..dim)
         := Standard_Series_Coefficients(x);
    y : Standard_Dense_Series.Series;
    grad : Standard_Dense_Series_Vectors.Vector(1..dim);
    forward : Standard_Complex_VecVecs.VecVec(1..nz-1)
            := Allocate_Coefficients(nz-1,deg);
    backward : Standard_Complex_VecVecs.VecVec(1..nz-2)
             := Allocate_Coefficients(nz-2,deg);
    cross : Standard_Complex_VecVecs.VecVec(1..nz-2)
          := Allocate_Coefficients(nz-2,deg);

  begin
    put("its exponent index : "); put(idx); new_line;
    put("The product : "); put(prd); new_line;
    put(dim,1); put(" random series of degree "); put(deg,1);
    put_line(" :"); put(x);
    put("The product of variables : ");
    put(prd); new_line;
    Speel(xcff,idx,forward,backward,cross);
    y := Standard_Series_Poly_Functions.Eval(prd,x);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    put_line("The coefficients of the product :");
    put_line(forward(idx'last-1));
    grad := Standard_Gradient(prd,x);
    put("derivative "); put(idx(idx'first),1); put_line(" :");
    put(grad(idx(idx'first))); new_line;
    put("coefficients of derivative "); put(idx(idx'first),1);
    put_line(" :"); put_line(backward(idx'last-2));
    for k in idx'first+1..idx'last-1 loop
      put("derivative "); put(idx(k),1); put_line(" :");
      put(grad(idx(k))); new_line;
      put("coefficients of derivative "); put(idx(k),1);
      put_line(" :"); put_line(cross(k-1));
    end loop;
    put("derivative "); put(idx(idx'last),1); put_line(" :");
    put(grad(idx(idx'last))); new_line;
    put("coefficients of derivative "); put(idx(idx'last),1);
    put_line(" :"); put_line(forward(idx'last-2));
  end Standard_Indexed_Test;

  procedure DoblDobl_Indexed_Test
              ( dim,deg,nz : in integer32;
                xp : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Generates random coefficient vectors for as many series as
  --   the value of dim, series truncated to degree deg,
  --   and runs a test in double double precision,
  --   for a random exponent vector xp of zeros and ones,
  --   where the number of nonzeros nz > 2.

    use DoblDobl_Speelpenning_Convolutions;

    idx : constant Standard_Integer_Vectors.Vector(1..nz)
        := Exponent_Indices.Exponent_Index(xp);
    prd : constant DoblDobl_Series_Polynomials.Poly
        := DoblDobl_Product(deg,xp);
    x : constant DoblDobl_Dense_Series_Vectors.Vector(1..dim)
      := DoblDobl_Random_Series.Random_Series_Vector(1,dim,deg);
    xcff : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
         := DoblDobl_Series_Coefficients(x);
    y : DoblDobl_Dense_Series.Series;
    grad : DoblDobl_Dense_Series_Vectors.Vector(1..dim);
    forward : DoblDobl_Complex_VecVecs.VecVec(1..nz-1)
            := Allocate_Coefficients(nz-1,deg);
    backward : DoblDobl_Complex_VecVecs.VecVec(1..nz-2)
             := Allocate_Coefficients(nz-2,deg);
    cross : DoblDobl_Complex_VecVecs.VecVec(1..nz-2)
          := Allocate_Coefficients(nz-2,deg);

  begin
    put("its exponent index : "); put(idx); new_line;
    put("The product : "); put(prd); new_line;
    put(dim,1); put(" random series of degree "); put(deg,1);
    put_line(" :"); put(x);
    put("The product of variables : ");
    put(prd); new_line;
    Speel(xcff,idx,forward,backward,cross);
    y := DoblDobl_Series_Poly_Functions.Eval(prd,x);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    put_line("The coefficients of the product :");
    put_line(forward(idx'last-1));
    grad := DoblDobl_Gradient(prd,x);
    put("derivative "); put(idx(idx'first),1); put_line(" :");
    put(grad(idx(idx'first))); new_line;
    put("coefficients of derivative "); put(idx(idx'first),1);
    put_line(" :"); put_line(backward(idx'last-2));
    for k in idx'first+1..idx'last-1 loop
      put("derivative "); put(idx(k),1); put_line(" :");
      put(grad(idx(k))); new_line;
      put("coefficients of derivative "); put(idx(k),1);
      put_line(" :"); put_line(cross(k-1));
    end loop;
    put("derivative "); put(idx(idx'last),1); put_line(" :");
    put(grad(idx(idx'last))); new_line;
    put("coefficients of derivative "); put(idx(idx'last),1);
    put_line(" :"); put_line(forward(idx'last-2));
  end DoblDobl_Indexed_Test;

  procedure QuadDobl_Indexed_Test
              ( dim,deg,nz : in integer32;
                xp : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Generates random coefficient vectors for as many series as
  --   the value of dim, series truncated to degree deg,
  --   and runs a test in quad double precision,
  --   for a random exponent vector xp of zeros and ones,
  --   where the number of nonzeros nz > 2.

    use QuadDobl_Speelpenning_Convolutions;

    idx : constant Standard_Integer_Vectors.Vector(1..nz)
        := Exponent_Indices.Exponent_Index(xp);
    prd : constant QuadDobl_Series_Polynomials.Poly
        := QuadDobl_Product(deg,xp);
    x : constant QuadDobl_Dense_Series_Vectors.Vector(1..dim)
      := QuadDobl_Random_Series.Random_Series_Vector(1,dim,deg);
    xcff : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
         := QuadDobl_Series_Coefficients(x);
    y : QuadDobl_Dense_Series.Series;
    grad : QuadDobl_Dense_Series_Vectors.Vector(1..dim);
    forward : QuadDobl_Complex_VecVecs.VecVec(1..nz-1)
            := Allocate_Coefficients(nz-1,deg);
    backward : QuadDobl_Complex_VecVecs.VecVec(1..nz-2)
             := Allocate_Coefficients(nz-2,deg);
    cross : QuadDobl_Complex_VecVecs.VecVec(1..nz-2)
          := Allocate_Coefficients(nz-2,deg);

  begin
    put("its exponent index : "); put(idx); new_line;
    put("The product : "); put(prd); new_line;
    put(dim,1); put(" random series of degree "); put(deg,1);
    put_line(" :"); put(x);
    put("The product of variables : ");
    put(prd); new_line;
    Speel(xcff,idx,forward,backward,cross);
    y := QuadDobl_Series_Poly_Functions.Eval(prd,x);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    put_line("The coefficients of the product :");
    put_line(forward(idx'last-1));
    grad := QuadDobl_Gradient(prd,x);
    put("derivative "); put(idx(idx'first),1); put_line(" :");
    put(grad(idx(idx'first))); new_line;
    put("coefficients of derivative "); put(idx(idx'first),1);
    put_line(" :"); put_line(backward(idx'last-2));
    for k in idx'first+1..idx'last-1 loop
      put("derivative "); put(idx(k),1); put_line(" :");
      put(grad(idx(k))); new_line;
      put("coefficients of derivative "); put(idx(k),1);
      put_line(" :"); put_line(cross(k-1));
    end loop;
    put("derivative "); put(idx(idx'last),1); put_line(" :");
    put(grad(idx(idx'last))); new_line;
    put("coefficients of derivative "); put(idx(idx'last),1);
    put_line(" :"); put_line(forward(idx'last-2));
  end QuadDobl_Indexed_Test;

  procedure Indexed_Test
              ( dim,deg : in integer32; precision : in character ) is

  -- DESCRIPTION :
  --   Generates random coefficient vectors for as many series as
  --   the value of dim, series truncated to degree deg,
  --   and runs a test in double, double double, or quad double precision
  --   (if precision is '0', '1', or '2' respectively),
  --   for a random exponent vector of zeros and ones.

    xp : Standard_Integer_Vectors.Vector(1..dim);
    nz : integer32;

  begin
    loop
      xp := Standard_Random_Vectors.Random_Vector(1,dim,0,1);
      nz :=  Standard_Integer_Vectors.Sum(xp);
      put("some random exponents : "); put(xp);
      put("  #nonzeros : "); put(nz,1); new_line;
      exit when (nz > 2);
    end loop;
    case precision is
      when '0' => Standard_Indexed_Test(dim,deg,nz,xp);
      when '1' => DoblDobl_Indexed_Test(dim,deg,nz,xp);
      when '2' => QuadDobl_Indexed_Test(dim,deg,nz,xp);
      when others => null;
    end case;
  end Indexed_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree, the dimension,
  --   and the precision of the coefficients.

    dim,deg : integer32 := 0;
    precision,ans : character;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(precision,"012");
    new_line;
    put("Test indexed version ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    if ans = 'y' then
      Indexed_Test(dim,deg,precision);
    else
      case precision is
        when '0' => Standard_Test(dim,deg);
        when '1' => DoblDobl_Test(dim,deg);
        when '2' => QuadDobl_Test(dim,deg);
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_speelser;
