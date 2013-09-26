with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Random_Matrices;          use Standard_Random_Matrices;
with Standard_Flag_Representations;     use Standard_Flag_Representations;

procedure ts_flag is

-- DESCRIPTION :
--   A flag is a nested sequence of linear spaces.

  procedure Random_Flag ( n,k : in integer32 ) is

    extrinsic_flag : constant Matrix(1..k,0..n) := Random_Matrix(1,k,0,n);
    intrinsic_flag : Matrix(1..n,0..n);

  begin
    new_line;
    put("Generating a random "); put(k,1); put("-flag in ");
    put(n,1); put_line("-space.");
    intrinsic_flag := Create_Intrinsic(extrinsic_flag);
    Test_Flag(Standard_Output,extrinsic_flag,intrinsic_flag);
  end Random_Flag;

  procedure Main is

    n,k : integer32;

  begin
    new_line;
    put_line("Extrinsic and intrinsic representations of flags");
    new_line;
    put("Give the dimension n of the ambient space : "); get(n);
    put("Give the co-dimension k of smallest space : "); get(k);
    Random_Flag(n,k);
  end Main;

begin
  Main;
end ts_flag;
