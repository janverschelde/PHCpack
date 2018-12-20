with text_io;                           use text_io;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Complex_Monomials;
with DoblDobl_Complex_Monomials;
with QuadDobl_Complex_Monomials;
with Random_Monomials;

package body Random_Monomial_Vectors is

  function Standard_Random_Monomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return Standard_Monomial_Vectors.Monomial_Vector is

    res : Standard_Monomial_Vectors.Monomial_Vector(1..size);

  begin
    if expmax = 0 then
      put_line("The upper bound on the exponent, expmax must be positive!");
    else
      for i in 1..size loop
        res(i) := Random_Monomials.Standard_Random_Monomial(dim,expmax,verbose);
      end loop;
    end if;
    return res;
  end Standard_Random_Monomial_Vector;

  function DoblDobl_Random_Monomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return DoblDobl_Monomial_Vectors.Monomial_Vector is

    res : DoblDobl_Monomial_Vectors.Monomial_Vector(1..size);

  begin
    if expmax = 0 then
      put_line("The upper bound on the exponent, expmax must be positive!");
    else
      for i in 1..size loop
        res(i) := Random_Monomials.DoblDobl_Random_Monomial(dim,expmax,verbose);
      end loop;
    end if;
    return res;
  end DoblDobl_Random_Monomial_Vector;

  function QuadDobl_Random_Monomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return QuadDobl_Monomial_Vectors.Monomial_Vector is

    res : QuadDobl_Monomial_Vectors.Monomial_Vector(1..size);

  begin
    if expmax = 0 then
      put_line("The upper bound on the exponent, expmax must be positive!");
    else
      for i in 1..size loop
        res(i) := Random_Monomials.QuadDobl_Random_Monomial(dim,expmax,verbose);
      end loop;
    end if;
    return res;
  end QuadDobl_Random_Monomial_Vector;

  function Standard_Random_Polynomial
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return Standard_Monomial_Vectors.Polynomial is

    res : Standard_Monomial_Vectors.Polynomial(dim,size);

  begin
    res.cff0 := Standard_Random_Numbers.Random1;
    res.mons := Standard_Random_Monomial_Vector(size,dim,expmax,verbose);
    Standard_Monomial_Vectors.Power_Update(res);
    return res;
  end Standard_Random_Polynomial;

  function Standard_Random_Polynomial
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return Standard_Monomial_Vectors.Link_to_Polynomial is

    res : Standard_Monomial_Vectors.Link_to_Polynomial;
    rnd : constant Standard_Monomial_Vectors.Polynomial(dim,size)
        := Standard_Random_Polynomial(size,dim,expmax,verbose);

  begin
    res := new Standard_Monomial_Vectors.Polynomial'(rnd);
    return res;
  end Standard_Random_Polynomial;

  function DoblDobl_Random_Polynomial
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return DoblDobl_Monomial_Vectors.Polynomial is

    res : DoblDobl_Monomial_Vectors.Polynomial(dim,size);

  begin
    res.cff0 := DoblDobl_Random_Numbers.Random1;
    res.mons := DoblDobl_Random_Monomial_Vector(size,dim,expmax,verbose);
    DoblDobl_Monomial_Vectors.Power_Update(res);
    return res;
  end DoblDobl_Random_Polynomial;

  function DoblDobl_Random_Polynomial
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return DoblDobl_Monomial_Vectors.Link_to_Polynomial is

    res : DoblDobl_Monomial_Vectors.Link_to_Polynomial;
    rnd : constant DoblDobl_Monomial_Vectors.Polynomial(dim,size)
        := DoblDobl_Random_Polynomial(size,dim,expmax,verbose);

  begin
    res := new DoblDobl_Monomial_Vectors.Polynomial'(rnd);
    return res;
  end DoblDobl_Random_Polynomial;

  function QuadDobl_Random_Polynomial
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return QuadDobl_Monomial_Vectors.Polynomial is

    res : QuadDobl_Monomial_Vectors.Polynomial(dim,size);

  begin
    res.cff0 := QuadDobl_Random_Numbers.Random1;
    res.mons := QuadDobl_Random_Monomial_Vector(size,dim,expmax,verbose);
    QuadDobl_Monomial_Vectors.Power_Update(res);
    return res;
  end QuadDobl_Random_Polynomial;

  function QuadDobl_Random_Polynomial
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return QuadDobl_Monomial_Vectors.Link_to_Polynomial is

    res : QuadDobl_Monomial_Vectors.Link_to_Polynomial;
    rnd : constant QuadDobl_Monomial_Vectors.Polynomial(dim,size)
        := QuadDobl_Random_Polynomial(size,dim,expmax,verbose);

  begin
    res := new QuadDobl_Monomial_Vectors.Polynomial'(rnd);
    return res;
  end QuadDobl_Random_Polynomial;

end Random_Monomial_Vectors;
