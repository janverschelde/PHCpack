package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.07 released 2015-12-24";

  begin
    return res;
  end Version;

end Greeting_Banners;
