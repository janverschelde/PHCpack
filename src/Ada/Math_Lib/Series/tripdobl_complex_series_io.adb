with Standard_Integer_Numbers_io;
with TripDobl_Complex_Vectors_io;

package body TripDobl_Complex_Series_io is

  procedure get ( s : in out Series ) is
  begin
    get(standard_input,s);
  end get;

  procedure get ( file : in file_type; s : in out Series ) is
  begin
    TripDobl_Complex_Vectors_io.get(file,s.cff(0..s.deg));
  end get;

  procedure get ( s : in out Link_to_Series; deg : in integer32 ) is
  begin
    get(standard_input,s,deg);
  end get;

  procedure get ( file : in file_type;
                  s : in out Link_to_Series; deg : in integer32 ) is
  begin
    s := Create(0,deg);
    TripDobl_Complex_Vectors_io.get(file,s.cff(0..s.deg));
  end get;

  procedure get ( s : in out Link_to_Series ) is
  begin
    get(standard_input,s);
  end get;

  procedure get ( file : in file_type; s : in out Link_to_Series ) is

    deg : integer32 := 0;

  begin
    Standard_Integer_Numbers_io.get(file,deg);
    s := Create(0,deg);
    TripDobl_Complex_Vectors_io.get(file,s.cff(0..s.deg));
  end get;

  procedure put ( s : in Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( file : in file_type; s : in Series ) is
  begin
    TripDobl_Complex_Vectors_io.put_line(file,s.cff(0..s.deg));
  end put;

  procedure put ( s : in Link_to_Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( file : in file_type; s : in Link_to_Series ) is
  begin
    if s /= null
     then TripDobl_Complex_Vectors_io.put_line(file,s.cff(0..s.deg));
    end if;
  end put;

  procedure put ( s : in Series; dp : in natural32 ) is
  begin
    put(standard_output,s,dp);
  end put;

  procedure put ( file : in file_type;
                  s : in Series; dp : in natural32 ) is
  begin
    TripDobl_Complex_Vectors_io.put_line(file,s.cff(0..s.deg),dp);
  end put;

  procedure put ( s : in Link_to_Series; dp : in natural32 ) is
  begin
    put(standard_output,s,dp);
  end put;

  procedure put ( file : in file_type;
                  s : in Link_to_Series; dp : in natural32 ) is
  begin
    if s /= null
     then TripDobl_Complex_Vectors_io.put_line(file,s.cff(0..s.deg),dp);
    end if;
  end put;

end TripDobl_Complex_Series_io;
