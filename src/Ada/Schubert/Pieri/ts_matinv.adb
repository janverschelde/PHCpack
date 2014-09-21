with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Double_Double_Matrices;
with Double_Double_Matrices_io;          use Double_Double_Matrices_io;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Matrices_io;       use DoblDobl_Complex_Matrices_io;
with Quad_Double_Matrices;
with Quad_Double_Matrices_io;            use Quad_Double_Matrices_io;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices_io;       use QuadDobl_Complex_Matrices_io;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with DoblDobl_Random_Matrices;           use DoblDobl_Random_Matrices;
with QuadDobl_Random_Matrices;           use QuadDobl_Random_Matrices;
with Standard_Matrix_Inversion;
with DoblDobl_Matrix_Inversion;
with QuadDobl_Matrix_Inversion;

procedure ts_matinv is

-- DESCRIPTION :
--   Test on matrix inversion.

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Tests matrix inversion using standard double arithmetic.

    use Standard_Floating_Matrices;
    use Standard_Complex_Matrices;
    use Standard_Matrix_Inversion;
    
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
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Tests matrix inversion using double double arithmetic.

    use Double_Double_Matrices;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Matrix_Inversion;
    
    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      cmpmat : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
             := Random_Matrix(natural32(n),natural32(n));
      invcmp : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
             := Inverse(cmpmat);
      fltmat : constant Double_Double_Matrices.Matrix(1..n,1..n)
             := Random_Matrix(natural32(n),natural32(n));
      invflt : constant Double_Double_Matrices.Matrix(1..n,1..n)
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
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Tests matrix inversion using double double arithmetic.

    use Quad_Double_Matrices;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Matrix_Inversion;
    
    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      cmpmat : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
             := Random_Matrix(natural32(n),natural32(n));
      invcmp : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
             := Inverse(cmpmat);
      fltmat : constant Quad_Double_Matrices.Matrix(1..n,1..n)
             := Random_Matrix(natural32(n),natural32(n));
      invflt : constant Quad_Double_Matrices.Matrix(1..n,1..n)
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
  end QuadDobl_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test matrix inversion :");
    put_line("  0. in standard double precision");
    put_line("  1. in double double precision");
    put_line("  2. in quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin 
  Main;
end ts_matinv;
