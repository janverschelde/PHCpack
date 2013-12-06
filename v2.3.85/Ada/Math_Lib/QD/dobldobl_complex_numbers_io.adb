with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;

package body DoblDobl_Complex_Numbers_io is

  procedure get ( c : in out Complex_Number ) is
  begin
    get(standard_input,c);
  end get;

  procedure get ( file : in file_type; c : in out Complex_Number ) is

    x,y : double_double := Create(0.0);

  begin
    get(file,x,y);
    c := Create(x,y);
  end get;

  procedure put ( c : in Complex_Number ) is
  begin
    put(standard_output,c);
  end put;

  procedure put ( file : in file_type; c : in Complex_Number ) is
  begin
    put(file,REAL_PART(c));
    put(file,"  ");
    put(file,IMAG_PART(c));
  end put;

  procedure put ( c : in Complex_Number; dp : in natural32 ) is
  begin
    put(standard_output,c,dp);
  end put;

  procedure put ( file : in file_type;
                  c : in Complex_Number; dp : in natural32 ) is
  begin
    put(file,REAL_PART(c),dp);
    put(file,"  ");
    put(file,IMAG_PART(c),dp);
  end put;

end DoblDobl_Complex_Numbers_io;
