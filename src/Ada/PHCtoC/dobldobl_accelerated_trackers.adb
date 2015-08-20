package body DoblDobl_Accelerated_Trackers is

  procedure DoblDobl_GPU_Newton ( execmode,verbose : integer32 ) is

    return_of_call : integer32;

    function newton ( m,v : integer32 ) return integer32;
    pragma import(C, newton, "gpunewton_dd");

  begin
    return_of_call := newton(execmode,verbose);
  end DoblDobl_GPU_Newton;

  procedure DoblDobl_GPU_Track_One ( execmode,verbose : integer32 ) is

    return_of_call : integer32;

    function one_track ( m,v : integer32 ) return integer32;
    pragma import(C, one_track, "gpuonepath_dd");

  begin
    return_of_call := one_track(execmode,verbose);
  end DoblDobl_GPU_Track_One;

  procedure DoblDobl_GPU_Track_Many ( execmode,verbose : integer32 ) is

    return_of_call : integer32;

    function many_track ( m,v : integer32 ) return integer32;
    pragma import(C, many_track, "gpumanypaths_dd");

  begin
    return_of_call := many_track(execmode,verbose);
  end DoblDobl_GPU_Track_Many;

end DoblDobl_Accelerated_Trackers; 
