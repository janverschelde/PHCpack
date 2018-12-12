package body Standard_SeriesPade_Tracker is

-- INTERNAL DATA :

  homconpars : Homotopy_Continuation_Parameters.Link_to_Parameters;

-- CONSTRUCTORS :

  procedure Init ( pars : in Homotopy_Continuation_Parameters.Parameters ) is
  begin
    homconpars := new Homotopy_Continuation_Parameters.Parameters'(pars);
  end Init;

  procedure Init ( p,q : in Link_to_Poly_Sys ) is
  begin
    null;
  end Init;

  procedure Init ( s : in Link_to_Solution ) is
  begin
    null;
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

end Standard_SeriesPade_Tracker;
