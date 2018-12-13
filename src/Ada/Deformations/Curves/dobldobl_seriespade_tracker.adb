with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_cv;        use DoblDobl_Complex_Numbers_cv;
with DoblDobl_Homotopy;

package body DoblDobl_SeriesPade_Tracker is

-- INTERNAL DATA :

  homconpars : Homotopy_Continuation_Parameters.Link_to_Parameters;
  current : Link_to_Solution;

-- CONSTRUCTORS :

  procedure Init ( pars : in Homotopy_Continuation_Parameters.Parameters ) is
  begin
    homconpars := new Homotopy_Continuation_Parameters.Parameters'(pars);
  end Init;

  procedure Init ( p,q : in Link_to_Poly_Sys ) is

    tpow : constant natural32 := 2;
    d_gamma : constant Standard_Complex_Numbers.Complex_Number
            := homconpars.gamma;
    dd_gamma : constant DoblDobl_Complex_Numbers.Complex_Number
             := Standard_to_DoblDobl_Complex(d_gamma);

  begin
    DoblDobl_Homotopy.Create(p.all,q.all,tpow,dd_gamma);
  end Init;

  procedure Init ( s : in Link_to_Solution ) is
  begin
    current := s;
  end Init;

-- SELECTORS :

  function Get_Parameters
    return Homotopy_Continuation_Parameters.Link_to_Parameters is

  begin
    return homconpars;
  end Get_Parameters;

-- DESTRUCTOR :

  procedure Clear is
  begin
    Homotopy_Continuation_Parameters.Clear(homconpars);
  end Clear;

end DoblDobl_SeriesPade_Tracker;
