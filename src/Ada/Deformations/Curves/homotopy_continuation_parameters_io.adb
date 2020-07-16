with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Numbers_io;

package body Homotopy_Continuation_Parameters_io is

  procedure put ( pars : in Parameters ) is
  begin
    put(standard_output,pars);
  end put;

  procedure put ( file : in file_type; pars : in Parameters ) is
  begin
    put_line(file,"Values of the HOMOTOPY CONTINUATION PARAMETERS :");
    put(file," 1. gamma : "); put(file,pars.gamma); new_line(file);
    put(file," 2. degree of numerator of Pade approximant    : ");
    put(file,pars.numdeg,1); new_line(file);
    put(file," 3. degree of denominator of Pade approximant  : ");
    put(file,pars.dendeg,1); new_line(file);
    put(file," 4. maximum step size                          :");
    put(file,pars.maxsize,2); new_line(file);
    put(file," 5. minimum step size                          :");
    put(file,pars.minsize,2); new_line(file);
    put(file," 6. multiplication factor of the pole radius   :");
    put(file,pars.pbeta,2); new_line(file);
    put(file," 7. multiplication factor for the curvature    :");
    put(file,pars.cbeta,2); new_line(file);
    put(file," 8. tolerance on the residual of the predictor :");
    put(file,pars.alpha,2); new_line(file);
    put(file," 9. tolerance on the residual of the corrector :");
    put(file,pars.tolres,2); new_line(file);
    put(file,"10. tolerance on zero series coefficients      :");
    put(file,pars.epsilon,2); new_line(file);
    put(file,"11. maximum number of corrector steps          : ");
    put(file,pars.corsteps,1); new_line(file);
    put(file,"12. maximum number of steps on a path          : ");
    put(file,pars.maxsteps,1); new_line(file);
  end put;

  procedure Prompt_for_Selection ( nbr : out natural32 ) is
  begin
    loop
      put("Type a number to change a value, or 0 to exit : ");
      Numbers_io.Read_Natural(nbr);
      exit when (nbr < 13);
      put_line("Your number should be 12 or less.  Please try again.");
    end loop;
  end Prompt_for_Selection;

  procedure Prompt_for_Parameter
              ( pars : in out Parameters; nbr : in natural32 ) is

    nat : natural32;
    pos : positive;
    regam,imgam : double_float;

  begin
    case nbr is
      when 1 =>
        put("-> give the real part of the new gamma : ");
        Numbers_io.Read_Double_Float(regam);
        put("-> give the imaginary part of the new gamma : ");
        Numbers_io.Read_Double_Float(imgam);
        pars.gamma := Standard_Complex_Numbers.Create(regam,imgam);
      when 2 =>
        put("-> give a new numerator degree for the Pade approximant : ");
        Numbers_io.Read_Positive(pos);
        pars.numdeg := natural32(pos);
      when 3 =>
        put("-> give a new denominator degree for the Pade approximant : ");
        Numbers_io.Read_Positive(pos);
        pars.dendeg := natural32(pos);
      when 4 =>
        put("-> give a new value for the maximum step size : ");
        Numbers_io.Read_Positive_Float(pars.maxsize);
      when 5 =>
        put("-> give a new value for the minimum step size  : ");
        Numbers_io.Read_Positive_Float(pars.minsize);
      when 6 =>
        put("-> give a new multiplication factor for the pole radius : ");
        Numbers_io.Read_Positive_Float(pars.pbeta);
      when 7 =>
        put("-> give a new multiplication factor for the curvature : ");
        Numbers_io.Read_Positive_Float(pars.cbeta);
      when 8 =>
        put("-> give a new tolerance on the predictor residual : ");
        Numbers_io.Read_Positive_Float(pars.alpha);
      when 9 =>
        put("-> give a new tolerance on the corrector residual : ");
        Numbers_io.Read_Positive_Float(pars.tolres);
      when 10 =>
        put("-> give a new tolerance on a zero series coefficient : ");
        Numbers_io.Read_Positive_Float(pars.epsilon);
      when 11 =>
        put("-> give a new maximum number of corrector steps : ");
        Numbers_io.Read_Natural(nat);
        pars.corsteps := nat;
      when 12 =>
        put("-> give a new maximum number of steps on a path : ");
        Numbers_io.Read_Positive(pos);
        pars.maxsteps := natural32(pos);
      when others => null;
    end case;
  end Prompt_for_Parameter;

  procedure Tune ( pars : in out Parameters ) is

    nbr : natural32;

  begin
    loop
      put(pars);
      Prompt_for_Selection(nbr);
      exit when (nbr = 0);
      Prompt_for_Parameter(pars,nbr);
    end loop;
  end Tune;

end Homotopy_Continuation_Parameters_io;
