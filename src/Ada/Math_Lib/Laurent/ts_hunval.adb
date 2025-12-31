with Ada.Text_IO;
with Communications_with_User;
with Test_Leading_Evaluations;
with Test_Ordered_Evaluations;

procedure ts_hunval is

-- DESCRIPTION :
--   Calls the main test on the evaluation and differentiation of Laurent
--   polynomials in several variables at series with real powers for
--   a generalized Newton-Puiseux method.

  ans : character;

begin
  Ada.Text_IO.new_line;
  Ada.Text_IO.put("Test higher order evaluations ? (y/n) ");
  Communications_with_User.Ask_Yes_or_No(ans);
  if ans = 'y' then
    Test_Ordered_Evaluations.main;
  else
    Ada.Text_IO.put("Test indexed derivatives of monomial ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Leading_Evaluations.test_indexed_derivatives;
    else
      Ada.Text_IO.put("Test enumeration of numbers ? (y/n) ");
      Communications_with_User.Ask_Yes_or_No(ans);
      if ans = 'y'
       then Test_Leading_Evaluations.test_number_enumeration;
       else Test_Leading_Evaluations.main;
      end if;
    end if;
  end if;
end ts_hunval;
