package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.18 released 2016-05-26";

  begin
    return res;
  end Version;

end Greeting_Banners;
