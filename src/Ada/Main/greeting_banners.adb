package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.3.97 released 2015-05-22";

  begin
    return res;
  end Version;

end Greeting_Banners;
