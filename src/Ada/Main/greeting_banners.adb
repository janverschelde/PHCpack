package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.12 released 2016-02-19";

  begin
    return res;
  end Version;

end Greeting_Banners;
