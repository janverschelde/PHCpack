package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.3.94 released 2014-09-25";

  begin
    return res;
  end Version;

end Greeting_Banners;
