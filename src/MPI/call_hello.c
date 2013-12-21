/* IMPORTANT NOTICE:
 * The adainit and adafinal must be called only once,
 * and since call_hello will be executed many times,
 * it is better to to do these calls in the main C program */

/* extern void adainit(); */
extern int _ada_hello ( int id, const char *action, const char *message );
/* extern void adafinal(); */

void call_hello ( int id, const char *action, const char *message )
{
/* this will be done by the Ada routine : 
 *   printf("%d %s %s\n", id, action, message); */


/*   adainit(); */
  int result = _ada_hello(id,action,message);
/*   adafinal(); */
}
