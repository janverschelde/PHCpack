with text_io;                           use text_io;
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
    if expmax = 0 then
      put_line("The upper bound on the exponent, expmax must be positive!");
    else
      for i in 1..size loop
        res(i) := Standard_Random_Polynomial(size,dim,expmax,verbose);
      end loop;
    end if;
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
    if expmax = 0 then
      put_line("The upper bound on the exponent, expmax must be positive!");
    else
      for i in 1..size loop
        res(i) := DoblDobl_Random_Polynomial(size,dim,expmax,verbose);
      end loop;
    end if;
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
    if expmax = 0 then
      put_line("The upper bound on the exponent, expmax must be positive!");
    else
      for i in 1..size loop
        res(i) := QuadDobl_Random_Polynomial(size,dim,expmax,verbose);
      end loop;
    end if;
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

  function Standard_Random_System
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return Standard_Polynomial_Vectors.System is

    res : Standard_Polynomial_Vectors.System(dim,size);

  begin
    res.pols := Standard_Random_Polynomial_Vector(size,dim,expmax,verbose);
    Standard_Polynomial_Vectors.Power_Update(res);
    return res;
  end Standard_Random_System;

  function Standard_Random_System
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return Standard_Polynomial_Vectors.Link_to_System is

    res : Standard_Polynomial_Vectors.Link_to_System;
    rnd : constant Standard_Polynomial_Vectors.System(dim,size)
        := Standard_Random_System(size,dim,expmax,verbose);

  begin
    res := new Standard_Polynomial_Vectors.System'(rnd);
    return res;
  end Standard_Random_System;

  function DoblDobl_Random_System
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return DoblDobl_Polynomial_Vectors.System is

    res : DoblDobl_Polynomial_Vectors.System(dim,size);

  begin
    res.pols := DoblDobl_Random_Polynomial_Vector(size,dim,expmax,verbose);
    DoblDobl_Polynomial_Vectors.Power_Update(res);
    return res;
  end DoblDobl_Random_System;

  function DoblDobl_Random_System
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return DoblDobl_Polynomial_Vectors.Link_to_System is

    res : DoblDobl_Polynomial_Vectors.Link_to_System;
    rnd : constant DoblDobl_Polynomial_Vectors.System(dim,size)
        := DoblDobl_Random_System(size,dim,expmax,verbose);

  begin
    res := new DoblDobl_Polynomial_Vectors.System'(rnd);
    return res;
  end DoblDobl_Random_System;

  function QuadDobl_Random_System
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return QuadDobl_Polynomial_Vectors.System is

    res : QuadDobl_Polynomial_Vectors.System(dim,size);

  begin
    res.pols := QuadDobl_Random_Polynomial_Vector(size,dim,expmax,verbose);
    QuadDobl_Polynomial_Vectors.Power_Update(res);
    return res;
  end QuadDobl_Random_System;

  function QuadDobl_Random_System
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return QuadDobl_Polynomial_Vectors.Link_to_System is

    res : QuadDobl_Polynomial_Vectors.Link_to_System;
    rnd : constant QuadDobl_Polynomial_Vectors.System(dim,size)
        := QuadDobl_Random_System(size,dim,expmax,verbose);

  begin
    res := new QuadDobl_Polynomial_Vectors.System'(rnd);
    return res;
  end QuadDobl_Random_System;

end Random_Polynomial_Vectors;
