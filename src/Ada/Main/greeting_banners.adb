package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.06 released 2015-12-11";

  begin
    return res;
  end Version;

end Greeting_Banners;
