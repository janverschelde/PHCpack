package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.3.89 released 2014-06-30";

  begin
    return res;
  end Version;

end Greeting_Banners;
