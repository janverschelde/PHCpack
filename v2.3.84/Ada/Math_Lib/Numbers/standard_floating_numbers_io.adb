package body Standard_Floating_Numbers_io is

  package single_float_io is new text_io.float_io(single_float);
  package double_float_io is new text_io.float_io(double_float);

  procedure get ( f : in out single_float ) is
  begin
    single_float_io.get(f);
  end get;

  procedure get ( s : in string; f : in out single_float;
                  last : out positive ) is
  begin
    single_float_io.get(s,f,last);
  end get;

  procedure get ( file : in file_type; f : in out single_float ) is
  begin
    single_float_io.get(file,f);
  end get;

  procedure get ( f : in out double_float ) is
  begin
    double_float_io.get(f);
  end get;

  procedure get ( file : in file_type; f : in out double_float ) is
  begin
    double_float_io.get(file,f);
  end get;

  procedure get ( s : in string; f : in out double_float;
                  last : out positive ) is
  begin
    double_float_io.get(s,f,last);
  end get;

  procedure put ( f : in single_float ) is
  begin
    single_float_io.put(f);
  end put;

  procedure put ( s : out string; f : in single_float ) is
  begin
    single_float_io.put(s,f);
  end put;

  procedure put ( file : in file_type; f : in single_float ) is
  begin
    single_float_io.put(file,f);
  end put;

  procedure put ( f : in double_float ) is
  begin
    double_float_io.put(f);
  end put;

  procedure put ( s : out string; f : in double_float ) is
  begin
    double_float_io.put(s,f);
  end put;

  procedure put ( file : in file_type; f : in double_float ) is
  begin
    double_float_io.put(file,f);
  end put;

  procedure put ( f : in single_float; fore,aft,exp : in natural32 ) is
  begin
    single_float_io.put(f,natural(fore),natural(aft),natural(exp));
  end put;

  procedure put ( s : out string;
                  f : in single_float; aft,exp : in natural32 ) is
  begin
    single_float_io.put(s,f,natural(aft),natural(exp));
  end put;

  procedure put ( f : in double_float; fore,aft,exp : in natural32 ) is
  begin
    double_float_io.put(f,natural(fore),natural(aft),natural(exp));
  end put;

  procedure put ( s : out string;
                  f : in double_float; aft,exp : in natural32 ) is
  begin
    double_float_io.put(s,f,natural(aft),natural(exp));
  end put;

  procedure put ( file : in file_type;
                  f : in single_float; fore,aft,exp : in natural32 ) is
  begin
    single_float_io.put(file,f,natural(fore),natural(aft),natural(exp));
  end put;

  procedure put ( file : in file_type;
                  f : in double_float; fore,aft,exp : in natural32 ) is
  begin
    double_float_io.put(file,f,natural(fore),natural(aft),natural(exp));
  end put;

  procedure put ( f : in double_float; dp : in natural32 ) is
  begin
    double_float_io.put(f,natural(dp),natural(dp),natural(dp));
  end put;

  procedure put ( s : out string; f : in double_float; dp : in natural32 ) is
  begin
    double_float_io.put(s,f,natural(dp),natural(dp));
  end put;

  procedure put ( file : in file_type;
                  f : in double_float; dp : in natural32 ) is
  begin
    double_float_io.put(file,f,natural(dp),natural(dp),natural(dp));
  end put;

end Standard_Floating_Numbers_io;
