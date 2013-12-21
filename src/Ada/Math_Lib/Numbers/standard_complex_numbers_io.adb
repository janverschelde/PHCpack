with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;

package body Standard_Complex_Numbers_io is

  procedure get ( c : in out Complex_Number ) is
  begin
    get(Standard_Input,c);
  end get;

  procedure get ( file : in file_type; c : in out Complex_Number ) is

    x,y : double_float := 0.0;

  begin
    get(file,x); get(file,y);
    c := Create(x,y);
  end get;

  procedure put ( c : in Complex_Number ) is
  begin
    put(Standard_Output,c);
  end put;

  procedure put ( file : in file_type; c : in Complex_Number ) is
  begin
    put(file,REAL_PART(c));
    put(file,"  ");
    put(file,IMAG_PART(c));
  end put;

  procedure put ( c : in Complex_Number; fore,aft,exp : in natural32 ) is
  begin
    put(Standard_Output,c,fore,aft,exp);
  end put;

  procedure put ( file : in file_type;
                  c : in Complex_Number; fore,aft,exp : in natural32 ) is
  begin
    put(file,REAL_PART(c),fore,aft,exp);
    put(file,"  ");
    put(file,IMAG_PART(c),fore,aft,exp);
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

end Standard_Complex_Numbers_io;
