with Standard_Natural_Numbers;             use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;          use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;             use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;          use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;          use Standard_Natural_Vectors_io;
with Interpolation_Filters;                use Interpolation_Filters;
with Interpolation_Filters_io;             use Interpolation_Filters_io;
with Span_of_Component_io;                 use Span_of_Component_io;

package body Irreducible_Components_io is

  procedure get ( c : out Standard_Irreducible_Component ) is
  begin
    get(Standard_Input,c);
  end get;

  procedure get ( c : out Multprec_Irreducible_Component ) is
  begin
    get(Standard_Input,c);
  end get;

  procedure get ( file : in file_type;
                  c : out Standard_Irreducible_Component ) is

    f : Standard_Filter;

  begin
    get(file,f);
    c := Create(f);
  end get;

  procedure get ( file : in file_type;
                  c : out Multprec_Irreducible_Component ) is

    f : Multprec_Filter;

  begin
    get(file,f);
    c := Create(f);
  end get;

  procedure put ( c : in Standard_Irreducible_Component ) is
  begin
    put(Standard_Output,c);
  end put;

  procedure put ( c : in Multprec_Irreducible_Component ) is
  begin
    put(Standard_Output,c);
  end put;

  procedure put ( file : in file_type;
                  c : in Standard_Irreducible_Component ) is
  begin
    put(file,Filter(c));
    put(file,Span(c));
  end put;

  procedure put ( file : in file_type; 
                  c : in Multprec_Irreducible_Component ) is
  begin
    put(file,Filter(c));
    put(file,Span(c));
  end put;

-- MINIMAL DATA :

  procedure get_labels ( c : out Standard_Irreducible_Component ) is
  begin
    get(Standard_Input,c);
  end get_labels;

  procedure get_labels ( file : in file_type;
                         c : out Standard_Irreducible_Component ) is

    d : integer32 := 0;
    ch : character := ' ';

  begin
    get(file,d);
    while ch /= ':' loop
      get(file,ch);
    end loop;
    if d /= 0 and not End_of_Line(file) then
      declare
        lab : Standard_Natural_Vectors.Vector(1..d);
      begin
        get(file,lab);
        Initialize_Labels(c,lab);
      end;
    end if;
  end get_labels;

  procedure get_labels ( c : out Multprec_Irreducible_Component ) is
  begin
    get(Standard_Input,c);
  end get_labels;

  procedure get_labels ( file : in file_type;
                         c : out Multprec_Irreducible_Component ) is

    d : integer32 := 0;
    ch : character := ' ';

  begin
    get(file,d);
    while ch /= ':' loop
      get(file,ch);
    end loop;
    if d /= 0 and not End_of_Line(file) then
      declare
        lab : Standard_Natural_Vectors.Vector(1..d);
      begin
        get(file,lab);
        Initialize_Labels(c,lab);
      end;
    end if;
  end get_labels;

  procedure put_labels ( c : in Standard_Irreducible_Component ) is
  begin
    put(Standard_Output,c);
  end put_labels;

  procedure put_labels ( file : in file_type;
                         c : in Standard_Irreducible_Component ) is

    d : constant natural32 := Degree(c);

  begin
    put(file,d,1);
    put(file," : ");
    if d /= 0
     then put(file,Labels(c)); new_line(file);
    end if;
  end put_labels;
 
end Irreducible_Components_io;
