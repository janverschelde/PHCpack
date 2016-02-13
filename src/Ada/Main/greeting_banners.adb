package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.11 released 2016-02-12";

  begin
    return res;
  end Version;

end Greeting_Banners;
