with text_io;                            use text_io;
with Triple_Double_Numbers_io;           use Triple_Double_Numbers_io;
with TripDobl_Complex_Numbers;           use TripDobl_Complex_Numbers;
with TripDobl_Complex_Numbers_io;        use TripDobl_Complex_Numbers_io;

procedure ts_tdcmplx is

-- DESCRIPTION :
--   Tests the operations on complex numbers in triple double precision.

  procedure Test_io is

  -- DESCRIPTION :
  --   Prompts for a complex number and writes the number.

    c : Complex_Number;

  begin
    put("Give a complex number : "); get(c);
    put_line("-> the real part : "); put(REAL_PART(c)); new_line;
    put_line("-> the imaginary part : "); put(IMAG_PART(c)); new_line;
    put_line("-> your number :"); put(c); new_line;
  end Test_io;

  procedure Main is
  begin
    new_line;
    put_line("Testing triple double complex arithmetic ...");
    Test_io;
  end Main;

begin
  Main;
end ts_tdcmplx;
