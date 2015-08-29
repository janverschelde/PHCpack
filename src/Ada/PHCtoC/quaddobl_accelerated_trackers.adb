with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package body QuadDobl_Accelerated_Trackers is

  procedure QuadDobl_GPU_Newton ( execmode,verbose : in integer32 ) is

    return_of_call : integer32;

    function newton ( m,v : integer32 ) return integer32;
    pragma import(C, newton, "gpunewton_qd");

  begin
    return_of_call := newton(execmode,verbose);
  end QuadDobl_GPU_Newton;

  procedure QuadDobl_GPU_Track_One
             ( execmode,verbose : in integer32; gamma : in Complex_Number ) is

    return_of_call : integer32;
    regamma : constant double_float := REAL_PART(gamma);
    imgamma : constant double_float := IMAG_PART(gamma);

    function one_track ( m,v : integer32;
                         r,i : double_float ) return integer32;
    pragma import(C, one_track, "gpuonepath_qd");

  begin
    return_of_call := one_track(execmode,verbose,regamma,imgamma);
  end QuadDobl_GPU_Track_One;

  procedure QuadDobl_GPU_Track_Many
             ( execmode,verbose : in integer32; gamma : in Complex_Number ) is

    return_of_call : integer32;
    regamma : constant double_float := REAL_PART(gamma);
    imgamma : constant double_float := IMAG_PART(gamma);

    function many_track ( m,v : integer32;
                         r,i : double_float ) return integer32;
    pragma import(C, many_track, "gpumanypaths_qd");

  begin
    return_of_call := many_track(execmode,verbose,regamma,imgamma);
  end QuadDobl_GPU_Track_Many;

end QuadDobl_Accelerated_Trackers; 
