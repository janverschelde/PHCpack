with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Floating_Vectors;
with Bits_of_Doubles;

package body Vectored_Doubles is

  procedure Balanced_Quarter_Product
              ( dim : in integer32;
                x0,x1,x2,x3 : in Standard_Floating_Vectors.Vector;
                y0,y1,y2,y3 : in Standard_Floating_Vectors.Vector;
                s0,s1,s2,s3 : out double_float ) is
  begin
    s0 := 0.0; s1 := 0.0; s2 := 0.0; s3 := 0.0;
    for i in 1..dim loop
      s0 := s0 + x0(i)*y0(i);
      s1 := s1 + x0(i)*y1(i) + x1(i)*y0(i);
      s2 := s2 + x0(i)*y2(i) + x1(i)*y1(i) + x2(i)*y0(i);
      s3 := s3 + x0(i)*y3(i) + x1(i)*y2(i) + x2(i)*y1(i) + x3(i)*y0(i);
    end loop;
  end Balanced_Quarter_Product;

  procedure Write_Subsums ( s0,s1,s2,s3 : in double_float ) is

    use Bits_of_Doubles;

  begin
    put("s0 : "); put(s0);
    put(", n0 : "); put(Last_Zero_Count(s0),1); new_line;
    put("s1 : "); put(s1);
    put(", n1 : "); put(Last_Zero_Count(s1),1); new_line;
    put("s2 : "); put(s2);
    put(", n2 : "); put(Last_Zero_Count(s2),1); new_line;
    put("s3 : "); put(s3);
    put(", n3 : "); put(Last_Zero_Count(s3),1); new_line;
  end Write_Subsums;

  function to_double ( s0,s1,s2,s3 : double_float;
                       verbose : boolean := true ) return double_float is
  begin
    if verbose
     then Write_Subsums(s0,s1,s2,s3);
    end if;
    return s3 + s2 + s1 + s0;
  end to_double;

end Vectored_Doubles; 
