package body Generic_Complex_Numbers_io is

  use Ring_io;  use Ring_io.Ring;

  procedure get ( x : in out Complex_Number ) is
  begin
    get(Standard_Input,x);
  end get;

  procedure get ( file : in file_type; x : in out Complex_Number ) is

    xre,xim : number;

  begin
    get(file,xre); get(file,xim);
    c := Create(xre,xim);
    Clear(xre); Clear(xim);
  end get;

  procedure put ( x : in Complex_Number ) is
  begin
    put(Standard_Output,x);
  end put;

  procedure put ( file : in file_type; x : in Complex_Number ) is
  begin
    put(file,REAL_PART(x));
    put(file,"  ");
    put(file,IMAG_PART(x));
  end put;

end Generic_Complex_Numbers_io;
