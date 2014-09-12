package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.3.93 released 2014-09-12";

  begin
    return res;
  end Version;

end Greeting_Banners;
