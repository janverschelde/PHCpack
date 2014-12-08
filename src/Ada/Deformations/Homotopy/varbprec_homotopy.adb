with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Numbers;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Strings;
with Standard_Homotopy;

package body Varbprec_Homotopy is

-- INTERNAL DATA :

  start,target : Link_to_Array_of_Strings;
  st_start,st_target : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
  gamma : Standard_Complex_Numbers.Complex_Number;
  exp4t : natural32;
  standard_homotopy_initialized : boolean;

-- AUXILIARY CONSTRUCTOR :

  procedure Initialize_Standard_Homotopy is

  -- DESCRIPTION :
  --   Initializes the standard homotopy, using the string representations
  --   for the start and target polynomial systems.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_Strings;

    nvr : natural32;

  begin
    Standard_Homotopy.Clear;
    if start /= null then
      nvr := natural32(start'last);
      if Symbol_Table.Number < nvr + 1
       then Symbol_Table.Init(nvr+1);
      end if;
      st_start := new Poly_Sys'(Parse(nvr,start.all));
      if target /= null
       then st_target := new Poly_Sys'(Parse(nvr,target.all));
      end if;
      Standard_Homotopy.Create(st_target.all,st_start.all,exp4t,gamma);
      standard_homotopy_initialized := true;
    end if;
  end Initialize_Standard_Homotopy;

-- CONSTRUCTORS :

  procedure Create ( p,q : in Link_to_Array_of_Strings; k : in natural32;
                     a : in Standard_Complex_Numbers.Complex_Number ) is
  begin
    gamma := a;
    exp4t := k;
    if start /= null
     then Clear(start);
    end if;
    start := new Array_of_Strings(q'range);
    for i in p'range loop
      start(i) := new string'(q(i).all);
    end loop;
    if target /= null
     then Clear(target);
    end if;
    target := new Array_of_Strings(p'range);
    for i in p'range loop
      target(i) := new string'(p(i).all);
    end loop;
  end Create;

-- SELECTORS :

  function Standard_Homotopy_System
             return Standard_Complex_Poly_Systems.Link_to_Poly_Sys is

    use Standard_Complex_Poly_Systems;
    res : Link_to_Poly_Sys;

  begin
    if not standard_homotopy_initialized
     then Initialize_Standard_Homotopy;
    end if;
    if standard_homotopy_initialized
     then res := new Poly_Sys'(Standard_Homotopy.Homotopy_System);
    end if;
    return res;
  end Standard_Homotopy_System;

  function Eval ( x : Standard_Complex_Vectors.Vector;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(x'range);

  begin
    if not standard_homotopy_initialized
     then Initialize_Standard_Homotopy;
    end if;
    if standard_homotopy_initialized
     then res := Standard_Homotopy.Eval(x,t);
    end if;
    return res;
  end Eval;

  function Diff ( x : Standard_Complex_Vectors.Vector;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(x'range);

  begin
    if not standard_homotopy_initialized
     then Initialize_Standard_Homotopy;
    end if;
    if standard_homotopy_initialized
     then res := Standard_Homotopy.Diff(x,t);
    end if;
    return res;
  end Diff;

  function Diff ( x : Standard_Complex_Vectors.Vector;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(x'range,x'range);

  begin
    if not standard_homotopy_initialized
     then Initialize_Standard_Homotopy;
    end if;
    if standard_homotopy_initialized
     then res := Standard_Homotopy.Diff(x,t);
    end if;
    return res;
  end Diff;

-- DESTRUCTOR :

  procedure Clear is
  begin
    if start /= null
     then Clear(start);
    end if;
    if target /= null
     then Clear(target);
    end if;
    if standard_homotopy_initialized then
      Standard_Homotopy.Clear;
      Standard_Complex_Poly_Systems.Clear(st_start);
      Standard_Complex_Poly_Systems.Clear(st_target);
    end if;
  end Clear;

begin
  start := null;
  target := null;
  standard_homotopy_initialized := false;
end Varbprec_Homotopy;
