package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.13 released 2016-03-18";

  begin
    return res;
  end Version;

end Greeting_Banners;
