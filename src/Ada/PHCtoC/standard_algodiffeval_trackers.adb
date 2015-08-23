package body Standard_AlgoDiffEval_Trackers is

  procedure Standard_ADE_Newton ( verbose : integer32 ) is

    return_of_call : integer32;

    function newton ( v : integer32 ) return integer32;
    pragma import(C, newton, "adenewton_d");

  begin
    return_of_call := newton(verbose);
  end Standard_ADE_Newton;

  procedure Standard_ADE_Track_One ( verbose : integer32 ) is

    return_of_call : integer32;

    function one_track ( v : integer32 ) return integer32;
    pragma import(C, one_track, "adeonepath_d");

  begin
    return_of_call := one_track(verbose);
  end Standard_ADE_Track_One;

  procedure Standard_ADE_Track_Many ( verbose : integer32 ) is

    return_of_call : integer32;

    function many_track ( v : integer32 ) return integer32;
    pragma import(C, many_track, "ademanypaths_d");

  begin
    return_of_call := many_track(verbose);
  end Standard_ADE_Track_Many;

end Standard_AlgoDiffEval_Trackers; 
