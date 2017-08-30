// defines the function log2ceil, with prototype in log2ceil.h

int log2ceil ( int n )
{
   n = 2*n-1;
   int log2n = 0;
   while(n>1)
   {
      log2n++;
      n/=2;
   }
   return log2n;
}
