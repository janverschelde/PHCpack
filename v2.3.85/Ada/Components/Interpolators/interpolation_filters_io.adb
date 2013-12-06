with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Multprec_Complex_Polynomials_io;   use Multprec_Complex_Polynomials_io;

package body Interpolation_Filters_io is

  procedure get ( f : in out Standard_Filter ) is
  begin
    get(Standard_Input,f);
  end get;

  procedure get ( file : in file_type; f : in out Standard_Filter ) is
  begin
    null;
  end get;

  procedure get ( f : in out Multprec_Filter ) is
  begin
    get(Standard_Input,f);
  end get;

  procedure get ( file : in file_type; f : in out Multprec_Filter ) is
  begin
    null;
  end get;

  procedure put ( f : in Standard_Filter ) is
  begin
    put(Standard_Output,f);
  end put;

  procedure put ( file : in file_type; f : in Standard_Filter ) is
  begin
    put(file,Dimension(f),1); put(file,"  ");
    put(file,Degree(f),1); put(file,"  ");
    put(file,Centrality(f),1); new_line(file);
    put_line(file,Interpolator(f));
  end put;

  procedure put ( f : in Multprec_Filter ) is
  begin
    put(Standard_Output,f);
  end put;

  procedure put ( file : in file_type; f : in Multprec_Filter ) is
  begin
    put(file,Dimension(f),1); put(file,"  ");
    put(file,Degree(f),1); put(file,"  ");
    put(file,Centrality(f),1); new_line(file);
    put_line(file,Interpolator(f));
  end put;

end Interpolation_Filters_io;
