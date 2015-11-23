package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.05 released 2015-11-23";

  begin
    return res;
  end Version;

end Greeting_Banners;
