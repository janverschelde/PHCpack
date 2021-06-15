with Multprec_Floating64_Numbers;        use Multprec_Floating64_Numbers;
with Multprec_Floating64_Numbers_io;     use Multprec_Floating64_Numbers_io;

package body Multprec_Complex64_Numbers_io is

  procedure get ( x : in out Complex_Number ) is
  begin
    get(Standard_Input,x);
  end get;

  procedure get ( file : in file_type; x : in out Complex_Number ) is

    xre,xim : Floating_Number;

  begin
    get(file,xre);
	get(file,xim);
    x := Create(xre,xim);
    Clear(xre);
    Clear(xim);
  end get;

  procedure put ( x : in Complex_Number ) is
  begin
    put(Standard_Output,x);
  end put;

  procedure put ( file : in file_type; x : in Complex_Number ) is

    xre,xim : Floating_Number;

  begin
    xre := REAL_PART(x);
    xim := IMAG_PART(x);
    put(file,xre);
    put(file,"  ");
    put(file,xim);
    Clear(xre);
    Clear(xim);
  end put;

  procedure put ( c : in Complex_Number; fore,aft,exp : in natural32 ) is
  begin
    put(Standard_Output,c,fore,aft,exp);
  end put;

  procedure put ( file : in file_type;
                  c : in Complex_Number; fore,aft,exp : in natural32 ) is

    cre,cim : Floating_Number;

  begin
    cre := REAL_PART(c);
    cim := IMAG_PART(c);
    put(file,cre,fore,aft,exp);
    put(file,"  ");
    put(file,cim,fore,aft,exp);
    Clear(cre);
    Clear(cim);
  end put;

  procedure put ( c : in Complex_Number; dp : in natural32 ) is
  begin
    put(c,dp,dp,dp);
  end put;

  procedure put ( file : in file_type;
                  c : in Complex_Number; dp : in natural32 ) is
  begin
    put(file,c,dp,dp,dp);
  end put;                                       

end Multprec_Complex64_Numbers_io;
