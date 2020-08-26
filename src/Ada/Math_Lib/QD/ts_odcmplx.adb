with text_io;                            use text_io;
with Octo_Double_Numbers_io;             use Octo_Double_Numbers_io;
with OctoDobl_Complex_Numbers;           use OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers_io;        use OctoDobl_Complex_Numbers_io;

procedure ts_odcmplx is

-- DESCRIPTION :
--   Tests the operations on complex numbers in octo double precision.

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
    put_line("Testing octo double complex arithmetic ...");
    Test_io;
  end Main;

begin
  Main;
end ts_odcmplx;
