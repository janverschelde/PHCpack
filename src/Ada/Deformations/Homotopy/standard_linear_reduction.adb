with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;

package body Standard_Linear_Reduction is

  procedure Pop_First_Term ( p : in out Poly; t : in out Term ) is

    procedure First_Term ( tt : in Term; continue : out boolean ) is
    begin
      Copy(tt,t);
      continue := false;
    end First_Term;
    procedure Get_First_Term is new Visiting_Iterator(First_Term);

  begin
    t.cf := Create(0.0);
    Get_First_Term(p);
    if t.cf /= Create(0.0)
     then Sub(p,t);
    end if;
  end Pop_First_Term;

  procedure Leading_Terms ( p : in out Poly_Sys; ta : in out Terms_Array ) is
  begin
    for i in p'range loop
      Clear(ta(i));
      Pop_First_Term(p(i),ta(i));
    end loop;
  end Leading_Terms;

  procedure Find_Max ( ta : in Terms_Array; index : in out Boolean_Array;
                       stop : in out boolean ) is

    res : Degrees := new Standard_Natural_Vectors.Vector'(ta(1).dg'range => 0);

  begin
    stop := true;
    for i in ta'range loop
      if ta(i).cf /= Create(0.0) then
        if ta(i).dg > res then
          res.all := ta(i).dg.all;
          index(1..(i-1)) := (1..(i-1) => false);
          index(i) := true;
          stop := false;
        elsif Equal(ta(i).dg,res) then
          index(i) := true;
          stop := false;
        end if;
      end if;
    end loop;
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(res));
  end Find_Max;

  procedure Update ( p : in out Poly_Sys; n : in natural32;
                     ta : in out Terms_Array; da : in out Degrees_Array;
                     nda,cnt : in out natural32; mat : in out Matrix;
                     stop : in out boolean ) is

    index : integer32;
    is_max : Boolean_Array(1..integer32(n)) := (1..integer32(n) => false);

  begin
    Find_Max(ta,is_max,stop);
    if not stop then
      -- Get the next Degrees for in the Degrees_Array
      for i in is_max'range loop
        if is_max(i)
         then index := i;  exit;
        end if;
      end loop;
      nda := nda + 1;
      Standard_Natural_Vectors.Copy
        (Standard_Natural_Vectors.Link_to_Vector(ta(index).dg),
         Standard_Natural_Vectors.Link_to_Vector(da(integer32(nda))));
      -- Fill in the matrix and update the window :
      for i in is_max'range loop
        if is_max(i) then
          mat(i,integer32(nda)) := ta(i).cf;
          Pop_First_Term(p(i),ta(i));
          cnt := cnt+1;
        else 
          mat(i,integer32(nda)) := Create(0.0);
        end if;
       end loop;
    end if;
  end Update;

  procedure Coefficient_Matrix
                ( p : in Poly_Sys; mat : in out Matrix;
                  da : in out Degrees_Array; nda : in out natural32;
                  diagonal : in out boolean ) is

    work : Poly_Sys(p'range);
    stop : boolean := false;
    Window : Terms_Array(p'range);
    cnt : natural32 := 0;

  begin
    Copy(p,work);
    Leading_Terms(work,Window);
    diagonal := false;
    while not stop and not diagonal loop
      Update(work,natural32(p'last),Window,da,nda,cnt,mat,stop);
      if (cnt = natural32(p'last)) and then (cnt = nda)
       then diagonal := true;
      end if;
      exit when diagonal;
    end loop;
    Clear(work);
  end Coefficient_Matrix;

  procedure Coefficient_Matrix
                ( p : in Poly_Sys; mat : in out Matrix;
                  da : in out Degrees_Array; nda : in out natural32 ) is

    work : Poly_Sys(p'range);
    stop : boolean := false;
    Window : Terms_Array(p'range);
    cnt : natural32 := 0;

  begin
    Copy(p,work);
    Leading_Terms(work,Window);
    while not stop loop
      Update(work,natural32(p'last),Window,da,nda,cnt,mat,stop);
    end loop;
    Clear(work);
  end Coefficient_Matrix;

  procedure Make_Polynomial_System
                    ( p : in out Poly_Sys; mat : in Matrix;
                      da : in Degrees_Array; nda : in natural32;
                      inconsistent,infinite : out boolean ) is

    t : Term;

  begin
    inconsistent := false;
    infinite := false;
    Clear(p);
    for i in p'range loop
      p(i) := Null_Poly;
      for j in 1..integer32(nda) loop
        if AbsVal(mat(i,j)) > mach_eps then
          t.dg := new Standard_Natural_Vectors.Vector'(da(j).all); 
          t.cf := mat(i,j);
          Add(p(i),t);
          Clear(t);
        end if;
      end loop;
      if p(i) = Null_Poly then
        infinite := true;
      elsif Degree(p(i)) = 0 then
        inconsistent := true;
      end if; 
    end loop;
  end Make_Polynomial_System;

  function Sum_Number_of_Terms ( p : Poly_Sys ) return natural32 is

  -- DESCRIPTION :
  --   Returns the sum of the number of terms of every polynomial in p.

    sum : natural32 := 0;

  begin
    for i in p'range loop
      sum := sum + Number_of_Terms(p(i));
    end loop;
    return sum;
  end Sum_Number_of_Terms;

end Standard_Linear_Reduction;
