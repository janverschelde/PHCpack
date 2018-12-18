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
      res(i) := Standard_Random_Polynomial(size,dim,expmax,verbose);
    end loop;
    return res;
  end Standard_Random_Polynomial_Vector;

  function Standard_Random_Polynomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return Standard_Polynomial_Vectors.Link_to_Polynomial_Vector is

    res : Standard_Polynomial_Vectors.Link_to_Polynomial_Vector;
    rnd : constant Standard_Polynomial_Vectors.Polynomial_Vector(1..size)
        := Standard_Random_Polynomial_Vector(size,dim,expmax,verbose);

  begin
    res := new Standard_Polynomial_Vectors.Polynomial_Vector'(rnd);
    return res;
  end Standard_Random_Polynomial_Vector;

  function DoblDobl_Random_Polynomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return DoblDobl_Polynomial_Vectors.Polynomial_Vector is

    res : DoblDobl_Polynomial_Vectors.Polynomial_Vector(1..size);

  begin
    for i in 1..size loop
      res(i) := DoblDobl_Random_Polynomial(size,dim,expmax,verbose);
    end loop;
    return res;
  end DoblDobl_Random_Polynomial_Vector;

  function DoblDobl_Random_Polynomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return DoblDobl_Polynomial_Vectors.Link_to_Polynomial_Vector is

    res : DoblDobl_Polynomial_Vectors.Link_to_Polynomial_Vector;
    rnd : constant DoblDobl_Polynomial_Vectors.Polynomial_Vector(1..size)
        := DoblDobl_Random_Polynomial_Vector(size,dim,expmax,verbose);

  begin
    res := new DoblDobl_Polynomial_Vectors.Polynomial_Vector'(rnd);
    return res;
  end DoblDobl_Random_Polynomial_Vector;

  function QuadDobl_Random_Polynomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return QuadDobl_Polynomial_Vectors.Polynomial_Vector is

    res : QuadDobl_Polynomial_Vectors.Polynomial_Vector(1..size);

  begin
    for i in 1..size loop
      res(i) := QuadDobl_Random_Polynomial(size,dim,expmax,verbose);
    end loop;
    return res;
  end QuadDobl_Random_Polynomial_Vector;

  function QuadDobl_Random_Polynomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return QuadDobl_Polynomial_Vectors.Link_to_Polynomial_Vector is

    res : QuadDobl_Polynomial_Vectors.Link_to_Polynomial_Vector;
    rnd : constant QuadDobl_Polynomial_Vectors.Polynomial_Vector(1..size)
        := QuadDobl_Random_Polynomial_Vector(size,dim,expmax,verbose);

  begin
    res := new QuadDobl_Polynomial_Vectors.Polynomial_Vector'(rnd);
    return res;
  end QuadDobl_Random_Polynomial_Vector;

end Random_Polynomial_Vectors;
