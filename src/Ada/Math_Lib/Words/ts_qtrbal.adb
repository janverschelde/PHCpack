with Ada.Text_IO;
with Communications_with_User;
with Test_Random_Balanced_Quarters;
with Test_Quarter_Balancers;

procedure ts_qtrbal is

  ans : character;

begin
  Ada.Text_IO.put("Test making of random balanced quarters ? (y/n) ");
  Communications_with_User.ask_yes_or_no(ans);
  if ans = 'y'
   then Test_Random_Balanced_Quarters.Main;
   else Test_Quarter_Balancers.Main;
  end if;
end ts_qtrbal;
