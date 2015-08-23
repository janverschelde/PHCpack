package body DoblDobl_AlgoDiffEval_Trackers is

  procedure DoblDobl_ADE_Newton ( verbose : integer32 ) is

    return_of_call : integer32;

    function newton ( v : integer32 ) return integer32;
    pragma import(C, newton, "adenewton_dd");

  begin
    return_of_call := newton(verbose);
  end DoblDobl_ADE_Newton;

  procedure DoblDobl_ADE_Track_One ( verbose : integer32 ) is

    return_of_call : integer32;

    function one_track ( v : integer32 ) return integer32;
    pragma import(C, one_track, "adeonepath_dd");

  begin
    return_of_call := one_track(verbose);
  end DoblDobl_ADE_Track_One;

  procedure DoblDobl_ADE_Track_Many ( verbose : integer32 ) is

    return_of_call : integer32;

    function many_track ( v : integer32 ) return integer32;
    pragma import(C, many_track, "ademanypaths_dd");

  begin
    return_of_call := many_track(verbose);
  end DoblDobl_ADE_Track_Many;

end DoblDobl_AlgoDiffEval_Trackers; 
