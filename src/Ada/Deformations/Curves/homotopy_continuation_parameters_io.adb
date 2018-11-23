with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
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
    put(file," 2. degree of numerator of Pade approximant   : ");
    put(file,pars.numdeg,1); new_line(file);
    put(file," 3. degree of denominator of Pade approximant : ");
    put(file,pars.dendeg,1); new_line(file);
    put(file," 4. maximum step size                         :");
    put(file,pars.maxsize,2); new_line(file);
    put(file," 5. minimum step size                         :");
    put(file,pars.minsize,2); new_line(file);
    put(file," 6. multiplication factor of the pole radius  :");
    put(file,pars.beta,2); new_line(file);
    put(file," 7. tolerance on the residual                 :");
    put(file,pars.alpha,2); new_line(file);
    put(file," 8. tolerance on zero series coefficient      :");
    put(file,pars.tolcff,2); new_line(file);
    put(file," 9. number of corrector steps                 : ");
    put(file,pars.corsteps,1); new_line(file);
    put(file,"10. maximum steps on a path                   : ");
    put(file,pars.maxsteps,1); new_line(file);
  end put;

  procedure Prompt_for_Selection ( nbr : out natural32 ) is
  begin
    loop
      put("Type a number to change a value, or 0 to exit : ");
      Numbers_io.Read_Natural(nbr);
      exit when (nbr < 11);
      put_line("Your number should be 10 or less.  Please try again.");
    end loop;
  end Prompt_for_Selection;

  procedure Prompt_for_Parameter
              ( pars : in out Parameters; nbr : in natural32 ) is

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
        put("-> give a new degree of numerator of the Pade approximant : ");
        Numbers_io.Read_Positive(pos);
        pars.numdeg := natural32(pos);
      when 3 =>
        put("-> give a new degree of denominator of the Pade approximant : ");
        Numbers_io.Read_Positive(pos);
        pars.dendeg := natural32(pos);
      when 4 =>
        put("-> give a new value for the maximum step size : ");
        Numbers_io.Read_Positive_Float(pars.maxsize);
      when 5 =>
        put("-> give a new value for the minimum step size  : ");
        Numbers_io.Read_Positive_Float(pars.minsize);
      when 6 =>
        put("-> give a new value for the multiplication factor : ");
        Numbers_io.Read_Positive_Float(pars.beta);
      when 7 =>
        put("-> give a new value for the tolerance on the residual : ");
        Numbers_io.Read_Positive_Float(pars.alpha);
      when 8 =>
        put("-> give a new tolerance on a zero series coefficient : ");
        Numbers_io.Read_Positive_Float(pars.tolcff);
      when 9 =>
        put("-> give a new value for the number of corrector steps : ");
        Numbers_io.Read_Positive(pos);
        pars.corsteps := natural32(pos);
      when 10 =>
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
