with unchecked_deallocation;
with Double_Taylor_Developments;         use Double_Taylor_Developments;

package body Double_Taylor_Homotopies is

  function Make ( deg : integer32; alpha,point : double_float;
                  cst : Standard_Complex_Numbers.Complex_Number;
                  exp : Standard_Integer_Vectors.Vector )
                return Taylor_Monomial is

    dim : constant integer32 := exp'last;
    res : Taylor_Monomial(dim,deg);

  begin
    res.cst := cst;
    res.pwr := alpha;
    if alpha = 0.0 then
      res.cff := (0..deg => 0.0);
      res.cff(0) := 1.0;
    else  
      res.cff := Double_Taylor_Coefficients(deg,alpha,point);
    end if;
    res.exp := exp;
    return res;
  end Make;

  function Make ( deg : integer32; alpha,point : double_float;
                  cst : Standard_Complex_Numbers.Complex_Number;
                  exp : Standard_Integer_Vectors.Vector )
                return Link_to_Taylor_Monomial is

    res : constant Link_to_Taylor_Monomial
        := new Taylor_Monomial'(Make(deg,alpha,point,cst,exp));

  begin
    return res;
  end Make;

  procedure Clear ( tm : in out Link_to_Taylor_Monomial ) is

    procedure free is
      new unchecked_deallocation(Taylor_Monomial,Link_to_Taylor_Monomial);

  begin
    if tm /= null 
     then free(tm);
    end if;
  end Clear;

  procedure Clear ( tmv : in out Taylor_Monomial_Vector ) is
  begin
    for i in tmv'range loop
      Clear(tmv(i));
    end loop;
  end Clear;

  procedure Clear ( tmv : in out Link_to_Taylor_Monomial_Vector ) is

    procedure free is
      new unchecked_deallocation(Taylor_Monomial_Vector,
                                 Link_to_Taylor_Monomial_Vector);

  begin
    if tmv /= null then
      for i in tmv'range loop
        Clear(tmv(i));
      end loop;
      free(tmv);
    end if;
  end Clear;

  procedure Clear ( th : in out Taylor_Homotopy ) is
  begin
    for i in th'range loop
      Clear(th(i));
    end loop;
  end Clear;

  procedure Clear ( th : in out Link_to_Taylor_Homotopy ) is

    procedure free is
      new unchecked_deallocation(Taylor_Homotopy,Link_to_Taylor_Homotopy);

  begin
    if th /= null then
      for i in th'range loop
        Clear(th(i));
      end loop;
      free(th);
    end if;
  end Clear;

end Double_Taylor_Homotopies;
