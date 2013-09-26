package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.3.84 released 2013-09-25";

  begin
    return res;
  end Version;

end Greeting_Banners;
