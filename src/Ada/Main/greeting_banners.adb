package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.3.88 released 2014-06-06";

  begin
    return res;
  end Version;

end Greeting_Banners;
