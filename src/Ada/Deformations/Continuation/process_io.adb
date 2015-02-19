with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;

package body Process_io is

  out_code : output_code;

  procedure Write_path ( file : in file_type; n : in natural32 ) is
  begin
    if out_code /= nil then
      put(file,"*****       path");
      put(file,n);
      put(file,"       *****");
      new_line(file);
    end if;
  end Write_path;

  procedure Write_block ( file : in file_type; n : in natural32 ) is
  begin
    put(file,"#####       block");
    put(file,n);
    put_line(file,"       #####");
  end Write_block;

  procedure Set_Output_Code ( u : in output_code ) is
  begin
    out_code := u;
  end Set_output_code;

  function Get_Output_Code return output_code is
  begin
    return out_code;
  end Get_Output_Code;

  function Contains_Output_Code ( u : output_code ) return boolean is

    res : boolean;

  begin
    if u = p then
      res := (out_code = p) or (out_code = sp)
          or (out_code = pc) or (out_code = spc);
    elsif u = c then
      res := (out_code = c) or (out_code = sc)
          or (out_code = pc) or (out_code = spc);
    elsif u = s then
      res := (out_code = s) or (out_code = sc)
          or (out_code = pc) or (out_code = spc);
    end if;
    return res;
  end Contains_Output_Code;

  procedure sWrite ( file : in file_type;
                     sol : in Standard_Complex_Solutions.Solution ) is
  begin
    if (out_code = s) or (out_code = sp) 
       or (out_code = sc) or (out_code = spc)
     then put(file,sol); new_line(file);
    end if;
  end sWrite;

  procedure sWrite ( file : in file_type;
                     sol : in DoblDobl_Complex_Solutions.Solution ) is
  begin
    if (out_code = s) or (out_code = sp) 
       or (out_code = sc) or (out_code = spc)
     then put(file,sol); new_line(file);
    end if;
  end sWrite;

  procedure sWrite ( file : in file_type;
                     sol : in QuadDobl_Complex_Solutions.Solution ) is
  begin
    if (out_code = s) or (out_code = sp) 
       or (out_code = sc) or (out_code = spc)
     then put(file,sol); new_line(file);
    end if;
  end sWrite;

  procedure sWrite ( file : in file_type;
                     sol : in Multprec_Complex_Solutions.Solution ) is
  begin
    if (out_code = s) or (out_code = sp) 
      or (out_code = sc) or (out_code = spc)
     then put(file,sol); new_line(file);
    end if;
  end sWrite;

  procedure pWrite ( file : in file_type; step : in double_float;
                     t : in Standard_Complex_Numbers.Complex_Number ) is
  begin
    if (out_code = p) or (out_code = sp) 
     or (out_code = pc) or (out_code = spc) then
      put(file,"step : "); put(file,step); put(file,"  ");
      put(file,"t : "); put(file,t); new_line(file);
    end if;
  end pWrite;

  procedure pWrite ( file : in file_type; step : in double_float;
                     t : in DoblDobl_Complex_Numbers.Complex_Number ) is
  begin
    if (out_code = p) or (out_code = sp) 
     or (out_code = pc) or (out_code = spc) then
      put(file,"step : "); put(file,step); put(file,"  ");
      put(file,"t : "); put(file,t); new_line(file);
    end if;
  end pWrite;

  procedure pWrite ( file : in file_type; step : in double_float;
                     t : in QuadDobl_Complex_Numbers.Complex_Number ) is
  begin
    if (out_code = p) or (out_code = sp) 
     or (out_code = pc) or (out_code = spc) then
      put(file,"step : "); put(file,step); put(file,"  ");
      put(file,"t : "); put(file,t); new_line(file);
    end if;
  end pWrite;

  procedure pWrite ( file : in file_type; step : in Floating_Number;
                     t : in Multprec_Complex_Numbers.Complex_Number ) is
  begin
    if (out_code = p) or (out_code = sp) 
     or (out_code = pc) or (out_code = spc) then
      put(file,"step : "); put(file,step); put(file,"  ");
      put(file,"t : "); put(file,t); new_line(file);
    end if;
  end pWrite;

  procedure pWrite ( file : in file_type; step : in double_float;
                     t : in Standard_Complex_Numbers.Complex_Number;
                     sol : in Standard_Complex_Solutions.Solution ) is
  begin
    if (out_code = p) or (out_code = sp) 
     or (out_code = pc) or (out_code = spc) then
      put(file,"step : "); put(file,step); put(file,"  ");
      put(file,"t : "); put(file,t); new_line(file);
      if (out_code = sp) or (out_code = spc) then
        put_line(file,"the predicted solution for t :");
        put_vector(file,sol);
      end if;
    end if;
  end pWrite;

  procedure pWrite ( file : in file_type; step : in double_float;
                     t : in DoblDobl_Complex_Numbers.Complex_Number;
                     sol : in DoblDobl_Complex_Solutions.Solution ) is
  begin
    if (out_code = p) or (out_code = sp) 
     or (out_code = pc) or (out_code = spc) then
      put(file,"step : "); put(file,step); put(file,"  ");
      put(file,"t : "); put(file,t); new_line(file);
      if (out_code = sp) or (out_code = spc) then
        put_line(file,"the predicted solution for t :");
        put_vector(file,sol);
      end if;
    end if;
  end pWrite;

  procedure pWrite ( file : in file_type; step : in double_float;
                     t : in QuadDobl_Complex_Numbers.Complex_Number;
                     sol : in QuadDobl_Complex_Solutions.Solution ) is
  begin
    if (out_code = p) or (out_code = sp) 
     or (out_code = pc) or (out_code = spc) then
      put(file,"step : "); put(file,step); put(file,"  ");
      put(file,"t : "); put(file,t); new_line(file);
      if (out_code = sp) or (out_code = spc) then
        put_line(file,"the predicted solution for t :");
        put_vector(file,sol);
      end if;
    end if;
  end pWrite;

  procedure pWrite ( file : in file_type; step : in Floating_Number;
                     t : in Multprec_Complex_Numbers.Complex_Number;
                     sol : in Multprec_Complex_Solutions.Solution ) is
  begin
    if (out_code = p) or (out_code = sp) 
     or (out_code = pc) or (out_code = spc) then
      put(file,"step : "); put(file,step); put(file,"  ");
      put(file,"t : "); put(file,t); new_line(file);
      if (out_code = sp) or (out_code = spc) then
        put_line(file,"the predicted solution for t :");
        put_vector(file,sol);
      end if;
    end if;
  end pWrite;

  procedure pWrite ( file : in file_type; step_number : in natural32;
                     step : in double_float;
                     t : in Standard_Complex_Numbers.Complex_Number ) is
  begin
    put(file,"step# "); put(file,step_number,1); 
    put(file,"  t = "); put(file,Standard_Complex_Numbers.REAL_PART(t),3);
    put(file,"  step = "); put(file,step,3); new_line(file);
  end pWrite;

  procedure pWrite ( file : in file_type; step_number : in natural32;
                     step : in double_double;
                     t : in DoblDobl_Complex_Numbers.Complex_Number ) is
  begin
    put(file,"step# "); put(file,step_number,1); 
    put(file,"  t = "); put(file,DoblDobl_Complex_Numbers.REAL_PART(t),3);
    put(file,"  step = "); put(file,step,3); new_line(file);
  end pWrite;

  procedure pWrite ( file : in file_type; step_number : in natural32;
                     step : in quad_double;
                     t : in QuadDobl_Complex_Numbers.Complex_Number ) is
  begin
    put(file,"step# "); put(file,step_number,1); 
    put(file,"  t = "); put(file,QuadDobl_Complex_Numbers.REAL_PART(t),3);
    put(file,"  step = "); put(file,step,3); new_line(file);
  end pWrite;

  procedure cWrite ( file : in file_type;
                     normax,normrx,normaf,normrf : in double_float ) is
  begin
    if (out_code = c) or (out_code = pc) 
      or (out_code = sc) or (out_code = spc) then
       put(file,"correction (a&r):");
       put(file,normax,3,3,3); put(file,normrx,3,3,3); put(file," ");
       put(file,"residual (a&r):");
       put(file,normaf,3,3,3); put(file,normrf,3,3,3); new_line(file);
    end if;
  end cWrite;

  procedure cWrite ( file : in file_type;
                     normax,normrx,normaf,normrf : in Floating_Number ) is
  begin
    if (out_code = c) or (out_code = pc) 
     or (out_code = sc) or (out_code = spc) then
      put(file,"correction (a&r):");
      put(file,normax,3,3,3); put(file,normrx,3,3,3); put(file," ");
      put(file,"residual (a&r):");
      put(file,normaf,3,3,3); put(file,normrf,3,3,3); new_line(file);
    end if;
  end cWrite;

  procedure cWrite ( file : in file_type; 
                     rcond : in double_float; m : in natural32 ) is
  begin
    if (out_code = c) or (out_code = sc) 
     or (out_code = pc) or (out_code = spc) then
      put(file,"rcond : "); put(file,rcond); put(file,"  ");
      put(file,"multiplicity : "); put(file,m,1); new_line(file);
    end if;
  end cWrite;

  procedure cWrite ( file : in file_type; 
                     rcond : in Floating_Number; m : in natural32 ) is
  begin
    if (out_code = c) or (out_code = sc) 
     or (out_code = pc) or (out_code = spc) then
      put(file,"rcond : "); put(file,rcond); put(file,"  ");
      put(file,"multiplicity : "); put(file,m,1); new_line(file);
    end if;
  end cWrite;

  procedure Write_Statistics ( file : in file_type;
                               nstep,nfail,niter,nsyst : in natural32 ) is
  begin
    put_line(file,"######################################################");
    put(file,"number of steps      :"); put(file,nstep); new_line(file);
    put(file,"number of failures   :"); put(file,nfail); new_line(file);
    put(file,"number of iterations :"); put(file,niter); new_line(file);
    put(file,"number of systems    :"); put(file,nsyst); new_line(file);
  end Write_Statistics;

  procedure Write_Total_Statistics 
                ( file : in file_type;
                  tnstep,tnfail,tniter,tnsyst : in natural32 ) is
  begin
    put(file,"total number of steps      :"); put(file,tnstep); new_line(file);
    put(file,"total number of failures   :"); put(file,tnfail); new_line(file);
    put(file,"total number of iterations :"); put(file,tniter); new_line(file);
    put(file,"total number of systems    :"); put(file,tnsyst); new_line(file);
  end Write_Total_Statistics;

  procedure sWrite_Solutions
              ( file : in file_type;
                sols : in Standard_Complex_Solutions.Solution_List ) is
  begin
    if out_code /= nil
     then put(file,sols);
    end if;
  end sWrite_Solutions;

  procedure sWrite_Solutions
              ( file : in file_type;
                sols : in Multprec_Complex_Solutions.Solution_List ) is
  begin
    if out_code /= nil
     then put(file,sols);
    end if;
  end sWrite_Solutions;

end Process_io;
