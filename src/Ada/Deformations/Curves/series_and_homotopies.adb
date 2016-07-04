with Series_and_Polynomials;

package body Series_and_Homotopies is

  function Create ( h : in Standard_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return Standard_Series_Poly_Systems.Poly_Sys is

    res : constant Standard_Series_Poly_Systems.Poly_Sys
        := Series_and_Polynomials.System_to_Series_System(h,idx,verbose);

  begin
    return res;
  end Create;

end Series_and_Homotopies;
