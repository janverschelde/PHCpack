with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Standard_Matrix_Inversion;          use Standard_Matrix_Inversion;

procedure ts_matinv is

-- DESCRIPTION :
--   Test on matrix inversion.

  procedure Main is
    
    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      cmpmat : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
             := Random_Matrix(natural32(n),natural32(n));
      invcmp : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
             := Inverse(cmpmat);
      fltmat : constant Standard_Floating_Matrices.Matrix(1..n,1..n)
             := Random_Matrix(natural32(n),natural32(n));
      invflt : constant Standard_Floating_Matrices.Matrix(1..n,1..n)
             := Inverse(fltmat);
    begin
      put_line("Test on complex matrix inversion : ");
      put(cmpmat*invcmp,3);
      put_line("Test on commutativity with complex matrix inversion : ");
      put(invcmp*cmpmat,3);
      put_line("Test on floating matrix inversion : ");
      put(fltmat*invflt,3);
      put_line("Test on commutativity with floating matrix inversion : ");
      put(invflt*fltmat,3);
    end;
  end Main;

begin 
  Main;
end ts_matinv;
