with Standard_Monomial_Vectors;
with DoblDobl_Monomial_Vectors;
with QuadDobl_Monomial_Vectors;
with Random_Monomial_Vectors;           use Random_Monomial_Vectors;

package body Random_Polynomial_Vectors is

  function Standard_Random_Polynomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return Standard_Polynomial_Vectors.Polynomial_Vector is

    res : Standard_Polynomial_Vectors.Polynomial_Vector(1..size);

  begin
    for i in 1..size loop
      declare
        pol : constant Standard_Monomial_Vectors.Monomial_Vector(1..size)
            := Standard_Random_Monomial_Vector(size,dim,expmax,verbose);
      begin
        res(i) := new Standard_Monomial_Vectors.Monomial_Vector'(pol);
      end;
    end loop;
    return res;
  end Standard_Random_Polynomial_Vector;

  function DoblDobl_Random_Polynomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return DoblDobl_Polynomial_Vectors.Polynomial_Vector is

    res : DoblDobl_Polynomial_Vectors.Polynomial_Vector(1..size);

  begin
    for i in 1..size loop
      declare
        pol : constant DoblDobl_Monomial_Vectors.Monomial_Vector(1..size)
            := DoblDobl_Random_Monomial_Vector(size,dim,expmax,verbose);
      begin
        res(i) := new DoblDobl_Monomial_Vectors.Monomial_Vector'(pol);
      end;
    end loop;
    return res;
  end DoblDobl_Random_Polynomial_Vector;

  function QuadDobl_Random_Polynomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return QuadDobl_Polynomial_Vectors.Polynomial_Vector is

    res : QuadDobl_Polynomial_Vectors.Polynomial_Vector(1..size);

  begin
    for i in 1..size loop
      declare
        pol : constant QuadDobl_Monomial_Vectors.Monomial_Vector(1..size)
            := QuadDobl_Random_Monomial_Vector(size,dim,expmax,verbose);
      begin
        res(i) := new QuadDobl_Monomial_Vectors.Monomial_Vector'(pol);
      end;
    end loop;
    return res;
  end QuadDobl_Random_Polynomial_Vector;

end Random_Polynomial_Vectors;
