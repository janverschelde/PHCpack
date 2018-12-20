with text_io;                           use text_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with DoblDobl_Random_Numbers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with QuadDobl_Random_Numbers;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;

package body Random_Monomials is

  function Random_Exponents
             ( dim : integer32; expmax : natural32 )  
             return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(1..dim);
    rnd : integer32;

  begin
    for i in 1..dim loop
      rnd := Standard_Random_Numbers.Random(0,integer32(expmax));
      res(i) := natural32(rnd);
    end loop;
    return res;
  end Random_Exponents;

  function Standard_Random_Monomial
             ( dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return Standard_Complex_Monomials.Monomial is

    res : Standard_Complex_Monomials.Monomial;
    exp : Standard_Natural_Vectors.Vector(1..dim);
    cff : constant Standard_Complex_Numbers.Complex_Number
        := Standard_Random_Numbers.Random1;

  begin
    if expmax = 0 then
      put_line("The upper bound on the exponent, expmax must be positive!");
    else
      loop
        exp := Random_Exponents(dim,expmax);
        if verbose then
          put("Coefficient : "); put(cff); new_line;
          put("Exponents : "); put(exp); new_line;
        end if;
        res := Standard_Complex_Monomials.Create(cff,exp);
        exit when Standard_Complex_Monomials.Degree(res) > 0;
      end loop;
    end if;
    return res;
  end Standard_Random_Monomial;

  function Standard_Random_Monomial
             ( dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return Standard_Complex_Monomials.Link_to_Monomial is

    res : Standard_Complex_Monomials.Link_to_Monomial;
    rnd : constant Standard_Complex_Monomials.Monomial
        := Standard_Random_Monomial(dim,expmax,verbose);

  begin
    res := new Standard_Complex_Monomials.Monomial'(rnd);
    return res;
  end Standard_Random_Monomial;

  function DoblDobl_Random_Monomial
             ( dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return DoblDobl_Complex_Monomials.Monomial is

    res : DoblDobl_Complex_Monomials.Monomial;
    exp : Standard_Natural_Vectors.Vector(1..dim);
    cff : constant DoblDobl_Complex_Numbers.Complex_Number
        := DoblDobl_Random_Numbers.Random1;

  begin
    if expmax = 0 then
      put_line("The upper bound on the exponent, expmax must be positive!");
    else
      loop
        exp := Random_Exponents(dim,expmax);
        if verbose then
          put("Coefficient : "); put(cff); new_line;
          put("Exponents : "); put(exp); new_line;
        end if;
        res := DoblDobl_Complex_Monomials.Create(cff,exp);
        exit when DoblDobl_Complex_Monomials.Degree(res) > 0;
      end loop;
    end if;
    return res;
  end DoblDobl_Random_Monomial;

  function DoblDobl_Random_Monomial
             ( dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return DoblDobl_Complex_Monomials.Link_to_Monomial is

    res : DoblDobl_Complex_Monomials.Link_to_Monomial;
    rnd : constant DoblDobl_Complex_Monomials.Monomial
        := DoblDobl_Random_Monomial(dim,expmax,verbose);

  begin
    res := new DoblDobl_Complex_Monomials.Monomial'(rnd);
    return res;
  end DoblDobl_Random_Monomial;

  function QuadDobl_Random_Monomial
             ( dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return QuadDobl_Complex_Monomials.Monomial is

    res : QuadDobl_Complex_Monomials.Monomial;
    exp : Standard_Natural_Vectors.Vector(1..dim);
    cff : constant QuadDobl_Complex_Numbers.Complex_Number
        := QuadDobl_Random_Numbers.Random1;

  begin
    if expmax = 0 then
      put_line("The upper bound on the exponent, expmax must be positive!");
    else
      loop
        exp := Random_Exponents(dim,expmax);
        if verbose then
          put("Coefficient : "); put(cff); new_line;
          put("Exponents : "); put(exp); new_line;
        end if;
        res := QuadDobl_Complex_Monomials.Create(cff,exp);
        exit when QuadDobl_Complex_Monomials.Degree(res) > 0;
      end loop;
    end if;
    return res;
  end QuadDobl_Random_Monomial;

  function QuadDobl_Random_Monomial
             ( dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return QuadDobl_Complex_Monomials.Link_to_Monomial is

    res : QuadDobl_Complex_Monomials.Link_to_Monomial;
    rnd : constant QuadDobl_Complex_Monomials.Monomial
        := QuadDobl_Random_Monomial(dim,expmax,verbose);

  begin
    res := new QuadDobl_Complex_Monomials.Monomial'(rnd);
    return res;
  end QuadDobl_Random_Monomial;

end Random_Monomials;
