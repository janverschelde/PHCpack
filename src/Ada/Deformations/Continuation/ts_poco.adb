with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Multprec_Floating_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Main_Poly_Continuation;             use Main_Poly_Continuation;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;

procedure ts_poco is

-- DESCRIPTION :
--   Calls the driver for polynomial continuation.

  procedure Standard_Continuation is

  -- DESCRIPTION :
  --   This is the standard driver for polynomial continuation
  --   with standard double floats.

    use Standard_Complex_Poly_Systems;

    file : file_type;
    lp : Link_to_Poly_Sys;
    ls : constant Link_to_Array_of_Strings := null;
    sols : Standard_Complex_Solutions.Solution_List;
    ddsols : DoblDobl_Complex_Solutions.Solution_List;
    qdsols : QuadDobl_Complex_Solutions.Solution_List;
    mpsols : Multprec_Complex_Solutions.Solution_List;
    target : Complex_Number;

  begin
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Driver_for_Polynomial_Continuation
      (file,lp.all,1,ls,sols,ddsols,qdsols,mpsols,target);
  end Standard_Continuation;

  procedure DoblDobl_Continuation is

  -- DESCRIPTION :
  --   This is the standard driver for polynomial continuation
  --   with double double floats.

    use DoblDobl_Complex_Poly_Systems;

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    target : Complex_Number;

  begin
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Driver_for_Polynomial_Continuation(file,lp.all,sols,target);
  end DoblDobl_Continuation;

  procedure QuadDobl_Continuation is

  -- DESCRIPTION :
  --   This is the standard driver for polynomial continuation
  --   with quad double floats.

    use QuadDobl_Complex_Poly_Systems;

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    target : Complex_Number;

  begin
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Driver_for_Polynomial_Continuation(file,lp.all,sols,target);
  end QuadDobl_Continuation;

  procedure Multprec_Continuation is

  -- DESCRIPTION :
  --   Calls polynomial continuation drivers for arbitrary
  --   precision floating point arithmetic.

    use Multprec_Complex_Poly_Systems;

    file : file_type;
    lp : Link_to_Poly_Sys;
    dp,sz : natural32 := 0;
    sols : Multprec_Complex_Solutions.Solution_List;
    target : Complex_Number;

  begin
    new_line;
    put("Give the number of decimal places : "); get(dp);
    sz := Multprec_Floating_Numbers.Decimal_to_Size(dp);
    Multprec_Complex_Polynomials_io.Set_Working_Precision(sz);
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Driver_for_Polynomial_Continuation(file,dp,lp.all,sols,target);
  end Multprec_Continuation;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test the continuation drivers :");
    put_line("  1. use standard double float arithmetic;");
    put_line("  2. use double double float arithmetic;");
    put_line("  3. use quad double float arithmetic;");
    put_line("  4. use multiprecision arithmetic.");
    put("Type 1, 2, 3, or 4 to make your choice : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' => Standard_Continuation;
      when '2' => DoblDobl_Continuation;
      when '3' => QuadDobl_Continuation;
      when '4' => Multprec_Continuation;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_poco;
