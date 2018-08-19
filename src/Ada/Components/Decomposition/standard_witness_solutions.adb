with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package body Standard_Witness_Solutions is

-- DATA :

  nbrvariables : integer32; -- number of variables is the ambient dimension
  topdimension : integer32; -- top dimension of the start of the cascade

  embpolysys : Standard_Complex_Poly_Systems.Link_to_Array_of_Poly_Sys;
  emblaursys : Standard_Complex_Laur_Systems.Link_to_Array_of_Laur_Sys;

-- CONSTRUCTORS :

  procedure Initialize ( nbvar,topdim : in natural32 ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Laur_Systems;

  begin
    nbrvariables := integer32(nbvar);
    topdimension := integer32(topdim);
    embpolysys := new Array_of_Poly_Sys(0..topdimension);
    emblaursys := new Array_of_Laur_Sys(0..topdimension);
  end Initialize;

  procedure Save_Embedded_System
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                k : in natural32 ) is

    use Standard_Complex_Poly_Systems;

    q : Poly_Sys(p'range);

  begin
    Copy(p,q);
    embpolysys(integer32(k)) := new Poly_Sys'(q);
  end Save_Embedded_System;

  procedure Save_Embedded_System
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                k : in natural32 ) is

    use Standard_Complex_Laur_Systems;

    q : Laur_Sys(p'range);

  begin
    Copy(p,q);
    emblaursys(integer32(k)) := new Laur_Sys'(q);
  end Save_Embedded_System;

-- SELECTORS :

  function Load_Embedded_System
              ( k : natural32 )
              return Standard_Complex_Poly_Systems.Link_to_Poly_Sys is

    res : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
        := embpolysys(integer32(k));

  begin
    return res;
  end Load_Embedded_System;

  function Load_Embedded_System
              ( k : natural32 )
              return Standard_Complex_Laur_Systems.Link_to_Laur_Sys is

    res : constant Standard_Complex_Laur_Systems.Link_to_Laur_Sys
        := emblaursys(integer32(k));

  begin
    return res;
  end Load_Embedded_System;

-- DESTRUCTOR :

  procedure Clear is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Laur_Systems;

  begin
    nbrvariables := -1;
    topdimension := -1;
    if embpolysys /= null
     then Clear(embpolysys);
    end if;
    if emblaursys /= null
     then Clear(emblaursys);
    end if;
  end Clear;

end Standard_Witness_Solutions;
