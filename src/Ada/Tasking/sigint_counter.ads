-- defines an interactive counter

package Sigint_Counter is

  continue : boolean := true;
  -- to continue or not

  function Continue_Counting return boolean;
  -- prompts the user to continue or not

  task Counter is
    entry Stop;
  end Counter;

end Sigint_Counter;
