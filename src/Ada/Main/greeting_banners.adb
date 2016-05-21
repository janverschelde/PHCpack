package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.17 released 2016-05-21";

  begin
    return res;
  end Version;

end Greeting_Banners;
