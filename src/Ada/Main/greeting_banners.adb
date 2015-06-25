package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.3.98 released 2015-06-25";

  begin
    return res;
  end Version;

end Greeting_Banners;
