with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Strings;
with Standard_Complex_Solutions;
with Standard_Solution_Strings;
with Solution_String_Splitters;
with Continuation_Parameters;
with Varbprec_Homotopy;
with Standard_Path_Tracker;
with Varbprec_Corrector_Steps;

package body Varbprec_Path_Tracker is

-- INTERNAL DATA :

  current : Link_to_String;

-- STANDARD DATA : the first approximation will always be done
--   first with the standard path tracker, so we initialize it:

  procedure Standard_Initialize_Homotopy
              ( p,q : in Link_to_Array_of_Strings; k : in natural32;
                gamma : in Standard_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   Parses the strings p and q into target and start systems
  --   for use in the initialization of the standard path tracker.

    use Standard_Complex_Poly_Systems;

    start,target : Poly_Sys(integer32(p'first)..integer32(p'last));
    nvr : constant natural32 := natural32(p'last);
    lp,lq : Link_to_Poly_Sys;

  begin
    if Symbol_Table.Number < nvr
     then Symbol_Table.Init(nvr);
    end if;
    target := Standard_Complex_Poly_Strings.Parse(nvr,p.all);
    start := Standard_Complex_Poly_Strings.Parse(nvr,q.all);
    lp := new Poly_Sys'(target);
    lq := new Poly_Sys'(start);
    Standard_Path_Tracker.Init(lp,lq,gamma,k);
  end Standard_Initialize_Homotopy;

  procedure Standard_Initialize_Homotopy 
               ( h : in Link_to_Array_of_Strings; txk : in integer32 ) is

  -- DESCRIPTION :
  --   Parses the string for a natural parameter homotopy
  --   with parameter txk and then initializes the standard path tracker.

    use Standard_Complex_Poly_Systems;

    hom : Poly_Sys(integer32(h'first)..integer32(h'last));
    nvr : constant natural32 := natural32(h'last) + 1;
    lh : Link_to_Poly_Sys;

  begin
    if Symbol_Table.Number < nvr
     then Symbol_Table.Init(nvr);
    end if;
    hom := Standard_Complex_Poly_Strings.Parse(nvr,h.all);
    lh := new Poly_Sys'(hom);
    Standard_Path_Tracker.Init(lh,txk);
  end Standard_Initialize_Homotopy;

  procedure Standard_Initialize_Solution
               ( s : in Link_to_String; nv : in natural32 ) is

  -- DESCRIPTION :
  --   Parses the string s into a solution in nv variables.
  --   If no failure, then the standard path tracker is initialized.

    sol : Standard_Complex_Solutions.Solution(integer32(nv));
    ls : Standard_Complex_Solutions.Link_to_Solution;
    p : integer := s'first;
    fail : boolean;

  begin
    Standard_Solution_Strings.Parse(s.all,p,nv,sol,fail);
    if not fail then
      ls := new Standard_Complex_Solutions.Solution'(sol);
      Standard_Path_Tracker.Init(ls);
    end if;
  end Standard_Initialize_Solution;

-- CONSTRUCTORS :

  procedure Init ( s : in Link_to_String; nv : in natural32 ) is
  begin
    current := s;
    Standard_Initialize_Solution(s,nv);
  end Init;

  procedure Init ( p,q : in Link_to_Array_of_Strings;
                   fixed_gamma : in boolean ) is

    re,im : double_float;
    gamma : Standard_Complex_Numbers.Complex_Number;

  begin
    if fixed_gamma then
      re := 0.57670012968461137;
      im := 0.8169559109411918;
      gamma := Standard_Complex_Numbers.Create(re,im);
    else
      gamma := Standard_Random_Numbers.Random1;
    end if;
    Init(p,q,2,gamma);
    Standard_Initialize_Homotopy(p,q,2,gamma);
  end Init;

  procedure Init ( p,q : in Link_to_Array_of_Strings;
                   fixed_gamma : in boolean; s : in Link_to_String ) is
  begin
    Init(p,q,fixed_gamma);
    Init(s,natural32(p'last));
  end Init;

  procedure Init ( p,q : in Link_to_Array_of_Strings; k : in natural32;
                   gamma : in Standard_Complex_Numbers.Complex_Number ) is
  begin
    Init(p,q,k,0,gamma);
  end Init;

  procedure Init ( p,q : in Link_to_Array_of_Strings; k : in natural32;
                   gamma : in Standard_Complex_Numbers.Complex_Number;
                   s : in Link_to_String ) is
  begin
    Init(p,q,k,0,gamma,s);
  end Init;

  procedure Init ( p,q : in Link_to_Array_of_Strings; k,cp : in natural32;
                   gamma : in Standard_Complex_Numbers.Complex_Number ) is
  begin
    Varbprec_Homotopy.Clear;
    Varbprec_Homotopy.Create(p,q,k,gamma);
    Standard_Initialize_Homotopy(p,q,k,gamma);
    Continuation_Parameters.Tune(cp);
  end Init;

  procedure Init ( p,q : in Link_to_Array_of_Strings; k,cp : in natural32;
                   gamma : in Standard_Complex_Numbers.Complex_Number;
                   s : in Link_to_String ) is
  begin
    Init(p,q,k,cp,gamma);
    Init(s,natural32(p'last));
  end Init;

  procedure Init ( h : in Link_to_Array_of_Strings; txk : in integer32 ) is
  begin
    Varbprec_Homotopy.Clear;
    Varbprec_Homotopy.Create(h,txk);
    Standard_Initialize_Homotopy(h,txk);
  end Init;

  procedure Init ( h : in Link_to_Array_of_Strings; txk : in integer32;
                   s : in Link_to_String ) is
  begin
    Varbprec_Homotopy.Clear;
    Varbprec_Homotopy.Create(h,txk);
    Standard_Initialize_Homotopy(h,txk);
    Init(s,natural32(h'last));
  end Init;

-- SELECTORS :

  function get_current return Link_to_String is
  begin
    return current;
  end get_current;

  function get_next ( want,maxprc,maxitr : natural32; output : boolean )
                    return Link_to_String is

    result,coords : Link_to_String;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    previous_t,t : Standard_Complex_Numbers.Complex_Number;
    m,loss : integer32;
    fail : boolean;
    err,rco,res : double_float;

    use Solution_String_Splitters,Varbprec_Corrector_Steps;

  begin
    previous_t := Standard_Path_Tracker.get_current.t;
    ls := Standard_Path_Tracker.get_next;
    if Standard_Complex_Numbers.Equal(previous_t,ls.t) then
      result := current;
    else
      result := new string'(Standard_Solution_Strings.Write(ls.all));
      Split_Coordinates(result.all,m,t,coords,fail);
      put_line("The coordinates : " & coords.all);
      loss := Estimate_Loss_for_Polynomial_Homotopy(coords.all,t,maxprc);
      put("-> estimated loss : "); put(loss,1); new_line;
      if output then
        Newton_Steps_on_Polynomial_Homotopy
          (standard_output,coords,t,want,maxprc,maxitr,loss,err,rco,res);
      else
        Newton_Steps_on_Polynomial_Homotopy
          (coords,t,want,maxprc,maxitr,loss,err,rco,res);
      end if;
      String_Splitters.Clear(result);
      declare
        n : constant integer32 := ls.v'last;
        r : constant string
          := Standard_Solution_Strings.Write(t,n,m,coords.all,err,rco,res);
      begin
        result := new string'(r);
      end;
      String_Splitters.Clear(coords);
      String_Splitters.Clear(current);
      current := result;
    end if;
    return result;
  end get_next;

  function get_next ( target_t : Standard_Complex_Numbers.Complex_Number;
                      want,maxprc,maxitr : natural32; output : boolean )
                    return Link_to_String is

    result,coords : Link_to_String;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    previous_t,t : Standard_Complex_Numbers.Complex_Number;
    m,loss : integer32;
    fail : boolean;
    err,rco,res : double_float;

    use Solution_String_Splitters,Varbprec_Corrector_Steps;

  begin
    previous_t := Standard_Path_Tracker.get_current.t;
    ls := Standard_Path_Tracker.get_next(target_t);
    if Standard_Complex_Numbers.Equal(previous_t,ls.t) then
      result := current;
    else
      result := new string'(Standard_Solution_Strings.Write(ls.all));
      Split_Coordinates(result.all,m,t,coords,fail);
      put_line("The coordinates : " & coords.all);
      loss := Estimate_Loss_for_Polynomial_Homotopy(coords.all,t,maxprc);
      put("-> estimated loss : "); put(loss,1); new_line;
      if output then
        Newton_Steps_on_Polynomial_Homotopy
          (standard_output,coords,t,want,maxprc,maxitr,loss,err,rco,res);
      else
        Newton_Steps_on_Polynomial_Homotopy
          (coords,t,want,maxprc,maxitr,loss,err,rco,res);
      end if;
      String_Splitters.Clear(result);
      declare
        n : constant integer32 := ls.v'last;
        r : constant string
          := Standard_Solution_Strings.Write(t,n,m,coords.all,err,rco,res);
      begin
        result := new string'(r);
      end;
      String_Splitters.Clear(coords);
      String_Splitters.Clear(current);
      current := result;
    end if;
    return result;
  end get_next;

-- DESTRUCTOR :

  procedure Clear is
  begin
    Varbprec_Homotopy.Clear;
  end Clear;

begin
  current := null;
end Varbprec_Path_Tracker;
