with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Complex_Numbers;

package body Real_Powered_Homotopy is

  function Is_Linear ( p : Standard_Complex_Laurentials.Poly )
                     return boolean is

    res : boolean := true;

    procedure Is_Term_Linear ( t : in Standard_Complex_Laurentials.Term;
                               continue : out boolean ) is
   
      sumdeg : integer32 := 0;

    begin
      for i in t.dg'range loop
        if t.dg(i) < 0 then
          res := false;
        elsif t.dg(i) /= 0 and t.dg(i) /= 1 then
          res := false;
        else
          sumdeg := sumdeg + t.dg(i);
          if sumdeg > 1
           then res := false;
          end if;
        end if;
        exit when (not res);
      end loop;
      continue := res;
    end Is_Term_Linear;

    procedure Are_Terms_Linear is
      new Standard_Complex_Laurentials.Visiting_Iterator(Is_Term_Linear);
   
  begin
    Are_Terms_Linear(p);
    return res;
  end Is_Linear;

  function Is_Linear ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                     return boolean is
  begin
    for i in p'range loop
      if not Is_Linear(p(i))
       then return false;
      end if;
    end loop;
    return true;
  end Is_Linear;

  function One_Index ( d : Standard_Complex_Laurentials.Degrees ) 
                     return integer32 is

  -- DESCRIPTION :
  --   Returns -1 if d is not a linear term,
  --   otherwise returns the index in d where the unique 1 appears,
  --   or returns 0 if d consists only of zeros.

    res : integer32 := 0; -- constant as default

  begin
    for i in d'range loop
      if d(i) < 0 then
        return -1;
      elsif d(i) > 1 then
        return -1;
      elsif d(i) = 1 then
        if res /= 0
         then return -1; -- found x(res)*x(i) term
         else res := i;  -- first time linear term
        end if;
      end if;
    end loop;
    return res;
  end One_Index;

  procedure Get_Constant_Row
              ( q : in Standard_Complex_Laurentials.Poly;
                c : in Standard_Complex_VecVecs.Link_to_VecVec;
                rowidx : in integer32;
                A : out Standard_Complex_Matrices.Matrix;
                b : out Standard_Complex_Numbers.Complex_Number;
                vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Given a real powered homotopy polynomial q with corresponding
  --   coefficients in c, extracts the constant coefficient for the
  --   row with index rowidx.

    cffidx : integer32 := 0; -- runs over the coefficients in c

    procedure Visit ( trm : in Standard_Complex_Laurentials.Term;
                      continue : out boolean ) is

      colidx : constant integer32 := One_Index(trm.dg);

    begin
      cffidx := cffidx + 1;
      if colidx = 0 then
        b := c(cffidx)(0);
      elsif colidx /= -1 then
        A(rowidx,colidx) := c(cffidx)(0);
      end if;
      continue := true;
    end Visit;

    procedure Visit_Terms is 
      new Standard_Complex_Laurentials.Visiting_Iterator(Visit);

  begin
    if vrblvl > 0 then
      put_line("-> in Real_Powered_Homotopy.get_constant_row ...");
    end if;
    b := Standard_Complex_Numbers.Create(0.0);
    for j in A'range(2) loop
      A(rowidx,j) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    Visit_Terms(q);
  end Get_Constant_Row;

  procedure Get_Constant_Coefficients
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                c : in Standard_Complex_VecVecs.Array_of_VecVecs;
                A : out Standard_Complex_Matrices.Matrix;
                b : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in Real_Powered_Homotopy.get_constant_coefficients ...");
    end if;
    for i in q'range loop
      Get_Constant_Row(q(i),c(i),i,A,b(i),vrblvl-1);
    end loop;
  end Get_Constant_Coefficients;

end Real_Powered_Homotopy;
