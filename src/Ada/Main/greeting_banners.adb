package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.3.92 released 2014-08-31";

  begin
    return res;
  end Version;

end Greeting_Banners;
