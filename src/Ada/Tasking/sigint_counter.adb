with Ada.text_io;
with Sigint_Handler;
with global_counter;

package body Sigint_Counter is

  function Continue_Counting return boolean is

    answer : character;

  begin
    Ada.text_io.put("Continue ? (y/n) ");
    Ada.text_io.get(answer);
    return (answer = 'y');
  end Continue_Counting;

  task body Counter is
  begin
    loop
      select
        accept Stop do 
          continue := Continue_Counting;
          if continue
           then requeue Sigint_Handler.Handler.Wait;
          end if;
        end Stop;
        exit when not continue;
      or delay 0.5;
      end select;
      global_counter.value := global_counter.value + 1;
      Ada.text_io.put_line(Natural'Image(global_counter.value));
    end loop;
  end Counter;

end Sigint_Counter;
