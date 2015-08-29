with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package body Standard_Accelerated_Trackers is

  procedure Standard_GPU_Newton ( execmode,verbose : in integer32 ) is

    return_of_call : integer32;

    function newton ( m,v : integer32 ) return integer32;
    pragma import(C, newton, "gpunewton_d");

  begin
    return_of_call := newton(execmode,verbose);
  end Standard_GPU_Newton;

  procedure Standard_GPU_Track_One
              ( execmode,verbose : in integer32;
                gamma : in Complex_Number ) is

    return_of_call : integer32;
    regam : constant double_float := REAL_PART(gamma);
    imgam : constant double_float := IMAG_PART(gamma);

    function one_track ( m,v : integer32;
                         r,i : double_float ) return integer32;
    pragma import(C, one_track, "gpuonepath_d");

  begin
    return_of_call := one_track(execmode,verbose,regam,imgam);
  end Standard_GPU_Track_One;

  procedure Standard_GPU_Track_Many
              ( execmode,verbose : in integer32;
                gamma : in Complex_Number ) is

    return_of_call : integer32;
    regam : constant double_float := REAL_PART(gamma);
    imgam : constant double_float := IMAG_PART(gamma);

    function many_track ( m,v : integer32;
                          r,i : double_float ) return integer32;
    pragma import(C, many_track, "gpumanypaths_d");

  begin
    return_of_call := many_track(execmode,verbose,regam,imgam);
  end Standard_GPU_Track_Many;

end Standard_Accelerated_Trackers; 
