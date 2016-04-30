package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.15 released 2016-04-30";

  begin
    return res;
  end Version;

end Greeting_Banners;
