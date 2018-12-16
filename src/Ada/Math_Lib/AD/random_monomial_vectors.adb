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
    for i in 1..size loop
      declare
        mon : constant Standard_Complex_Monomials.Monomial
            := Random_Monomials.Standard_Random_Monomial(dim,expmax,verbose);
      begin
        res(i) := new Standard_Complex_Monomials.Monomial'(mon);
      end;
    end loop;
    return res;
  end Standard_Random_Monomial_Vector;

  function DoblDobl_Random_Monomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return DoblDobl_Monomial_Vectors.Monomial_Vector is

    res : DoblDobl_Monomial_Vectors.Monomial_Vector(1..size);

  begin
    for i in 1..size loop
      declare
        mon : constant DoblDobl_Complex_Monomials.Monomial
            := Random_Monomials.DoblDobl_Random_Monomial(dim,expmax,verbose);
      begin
        res(i) := new DoblDobl_Complex_Monomials.Monomial'(mon);
      end;
    end loop;
    return res;
  end DoblDobl_Random_Monomial_Vector;

  function QuadDobl_Random_Monomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return QuadDobl_Monomial_Vectors.Monomial_Vector is

    res : QuadDobl_Monomial_Vectors.Monomial_Vector(1..size);

  begin
    for i in 1..size loop
      declare
        mon : constant QuadDobl_Complex_Monomials.Monomial
            := Random_Monomials.QuadDobl_Random_Monomial(dim,expmax,verbose);
      begin
        res(i) := new QuadDobl_Complex_Monomials.Monomial'(mon);
      end;
    end loop;
    return res;
  end QuadDobl_Random_Monomial_Vector;
 
end Random_Monomial_Vectors;
