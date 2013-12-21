with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;  use Multprec_Complex_Poly_Systems_io;

package body Span_of_Component_io is

  procedure get ( sp : in out Standard_Span ) is
  begin
    get(Standard_Input,sp);
  end get;

  procedure get ( file : file_type; sp : in out Standard_Span ) is

    n,d : integer32 := 0;

  begin
    get(file,n); get(file,d);
    declare
      frv : Standard_Integer_Vectors.Vector(1..d);
      equ : Standard_Complex_Poly_Systems.Poly_Sys(1..n-d);
    begin
      get(file,frv);
      get(file,equ);
      sp := Create(natural32(n),natural32(d),frv,equ);
    end;
  end get;

  procedure get ( sp : in out Multprec_Span ) is
  begin
    get(Standard_Input,sp);
  end get;

  procedure get ( file : file_type; sp : in out Multprec_Span ) is

    n,d : integer32 := 0;

  begin
    get(file,n); get(file,d);
    declare
      frv : Standard_Integer_Vectors.Vector(1..d);
      equ : Multprec_Complex_Poly_Systems.Poly_Sys(1..n-d);
    begin
      get(file,frv);
      get(file,equ);
      sp := Create(natural32(n),natural32(d),frv,equ);
    end;
  end get;

  procedure put ( sp : in Standard_Span ) is
  begin
    put(Standard_Output,sp);
  end put;

  procedure put ( file : file_type; sp : in Standard_Span ) is
  begin
    if not Empty(sp) then
      put(file,Ambient_Dimension(sp),1); put(file," ");
      put(file,Dimension(sp),1); new_line(file);
      put(file,Free_Variables(sp)); new_line(file);
      put(file,Equations(sp));
    end if;
  end put;

  procedure put ( sp : in Multprec_Span ) is
  begin
    put(Standard_Output,sp);
  end put;

  procedure put ( file : file_type; sp : in Multprec_Span ) is
  begin
    if not Empty(sp) then
      put(file,Ambient_Dimension(sp),1); put(file," ");
      put(file,Dimension(sp),1); new_line(file);
      put(file,Free_Variables(sp)); new_line(file);
      put(file,Equations(sp));
    end if;
  end put;

end Span_of_Component_io;
