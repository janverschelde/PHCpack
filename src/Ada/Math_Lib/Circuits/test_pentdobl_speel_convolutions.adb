with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Random_Vectors;
with PentDobl_Complex_Vectors_io;        use PentDobl_Complex_Vectors_io;
with PentDobl_Complex_VecVecs;
with PentDobl_Complex_Series;
with PentDobl_Complex_Series_io;         use PentDobl_Complex_Series_io;
with PentDobl_Complex_Series_Vectors;
with PentDobl_Complex_Series_Vectors_io; use PentDobl_Complex_Series_Vectors_io;
with PentDobl_Random_Series_Vectors;
with PentDobl_CSeries_Polynomials;
with PentDobl_CSeries_Polynomials_io;    use PentDobl_CSeries_Polynomials_io;
with PentDobl_CSeries_Poly_Functions;
with Series_Coefficient_Vectors;         use Series_Coefficient_Vectors;
with Exponent_Indices;
with PentDobl_Speelpenning_Convolutions;
with Series_Polynomial_Gradients;        use Series_Polynomial_Gradients;

package body Test_PentDobl_Speel_Convolutions is

  procedure PentDobl_Test ( dim,deg : in integer32 ) is

    use PentDobl_Speelpenning_Convolutions;

    prd : constant PentDobl_CSeries_Polynomials.Poly
        := PentDobl_Product(dim,deg);
    x : constant PentDobl_Complex_Series_Vectors.Vector(1..dim)
      := PentDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant PentDobl_Complex_VecVecs.VecVec(1..dim)
         := PentDobl_Series_Coefficients(x);
    y : PentDobl_Complex_Series.Link_to_Series;
    grad : PentDobl_Complex_Series_Vectors.Vector(1..dim);
    forward : constant PentDobl_Complex_VecVecs.VecVec(1..dim-1)
            := Allocate_Coefficients(dim-1,deg);
    backward : constant PentDobl_Complex_VecVecs.VecVec(1..dim-2)
             := Allocate_Coefficients(dim-2,deg);
    cross : constant PentDobl_Complex_VecVecs.VecVec(1..dim-2)
          := Allocate_Coefficients(dim-2,deg);

  begin
    put(dim,1); put(" random series of degree "); put(deg,1);
    put_line(" :"); put(x);
    put("The product of variables : ");
    put(prd); new_line;
    y := PentDobl_CSeries_Poly_Functions.Eval(prd,x);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    Speel(xcff,forward,backward,cross);
    put_line("The coefficients of the product : ");
    put_line(forward(dim-1));
    grad := PentDobl_Gradient(prd,x);
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
  end PentDobl_Test;

  procedure PentDobl_Indexed_Test
              ( dim,deg,nz : in integer32;
                xp : in Standard_Integer_Vectors.Vector ) is

    use PentDobl_Speelpenning_Convolutions;

    idx : constant Standard_Integer_Vectors.Vector(1..nz)
        := Exponent_Indices.Exponent_Index(xp);
    prd : constant PentDobl_CSeries_Polynomials.Poly
        := PentDobl_Product(deg,xp);
    x : constant PentDobl_Complex_Series_Vectors.Vector(1..dim)
      := PentDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant PentDobl_Complex_VecVecs.VecVec(1..dim)
         := PentDobl_Series_Coefficients(x);
    y : PentDobl_Complex_Series.Link_to_Series;
    grad : PentDobl_Complex_Series_Vectors.Vector(1..dim);
    forward : constant PentDobl_Complex_VecVecs.VecVec(1..nz-1)
            := Allocate_Coefficients(nz-1,deg);
    backward : constant PentDobl_Complex_VecVecs.VecVec(1..nz-2)
             := Allocate_Coefficients(nz-2,deg);
    cross : constant PentDobl_Complex_VecVecs.VecVec(1..nz-2)
          := Allocate_Coefficients(nz-2,deg);

  begin
    put("its exponent index : "); put(idx); new_line;
    put("The product : "); put(prd); new_line;
    put(dim,1); put(" random series of degree "); put(deg,1);
    put_line(" :"); put(x);
    put("The product of variables : ");
    put(prd); new_line;
    Speel(xcff,idx,forward,backward,cross);
    y := PentDobl_CSeries_Poly_Functions.Eval(prd,x);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    put_line("The coefficients of the product :");
    put_line(forward(idx'last-1));
    grad := PentDobl_Gradient(prd,x);
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
  end PentDobl_Indexed_Test;

  procedure Indexed_Test ( dim,deg : in integer32 ) is

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
    PentDobl_Indexed_Test(dim,deg,nz,xp);
  end Indexed_Test;

  procedure Main is

    dim,deg : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    put("Test indexed version ? (y/n) "); Ask_Yes_or_No(ans);
    new_line;
    if ans = 'y'
     then Indexed_Test(dim,deg);
     else PentDobl_Test(dim,deg);
    end if;
  end Main;

end Test_PentDobl_Speel_Convolutions;
