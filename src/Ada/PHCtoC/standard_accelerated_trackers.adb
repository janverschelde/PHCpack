package body Standard_Accelerated_Trackers is

  procedure Standard_GPU_Newton ( execmode,verbose : integer32 ) is

    return_of_call : integer32;

    function newton ( m,v : integer32 ) return integer32;
    pragma import(C, newton, "gpunewton_d");

  begin
    return_of_call := newton(execmode,verbose);
  end Standard_GPU_Newton;

  procedure Standard_GPU_Track_One ( execmode,verbose : integer32 ) is

    return_of_call : integer32;

    function one_track ( m,v : integer32 ) return integer32;
    pragma import(C, one_track, "gpuonepath_d");

  begin
    return_of_call := one_track(execmode,verbose);
  end Standard_GPU_Track_One;

  procedure Standard_GPU_Track_Many ( execmode,verbose : integer32 ) is

    return_of_call : integer32;

    function many_track ( m,v : integer32 ) return integer32;
    pragma import(C, many_track, "gpumanypaths_d");

  begin
    return_of_call := many_track(execmode,verbose);
  end Standard_GPU_Track_Many;

end Standard_Accelerated_Trackers; 
