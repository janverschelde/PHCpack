with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Continuation_Parameters;
with Varbprec_Homotopy;

package body Varbprec_Path_Tracker is

-- INTERNAL DATA :

  current : Link_to_String;

-- CONSTRUCTORS :

  procedure Init ( s : in Link_to_String ) is
  begin
    current := s;
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
  end Init;

  procedure Init ( p,q : in Link_to_Array_of_Strings;
                   fixed_gamma : in boolean; s : in Link_to_String ) is
  begin
    Init(p,q,fixed_gamma);
    Init(s);
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
    Continuation_Parameters.Tune(cp);
  end Init;

  procedure Init ( p,q : in Link_to_Array_of_Strings; k,cp : in natural32;
                   gamma : in Standard_Complex_Numbers.Complex_Number;
                   s : in Link_to_String ) is
  begin
    Init(p,q,k,cp,gamma);
    Init(s);
  end Init;

  procedure Init ( h : in Link_to_Array_of_Strings; txk : in integer32 ) is
  begin
    Varbprec_Homotopy.Clear;
   -- Varbprec_Homotopy.Create(h,txk);
  end Init;

  procedure Init ( h : in Link_to_Array_of_Strings; txk : in integer32;
                   s : in Link_to_String ) is
  begin
    Varbprec_Homotopy.Clear;
   -- Varbprec_Homotopy.Create(h.all,txk);
    Init(s);
  end Init;

-- DESTRUCTOR :

  procedure Clear is
  begin
    Varbprec_Homotopy.Clear;
  end Clear;

begin
  current := null;
end Varbprec_Path_Tracker;
