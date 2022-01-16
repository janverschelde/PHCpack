/* contains definitions of functions with prototypes in lists_and_strings.h */

#include<stdio.h>

int intlist2str ( int n, int *d, char *s )
{
   int bufsize = 16;
   int cnt = 0;
   int i,j;

   s[cnt++] = '[';
   for(i=0; i<n; i++)
   {
      char buf[bufsize];
      for(j=0; j<bufsize; j++) buf[j] = ' ';
      sprintf(buf,"%d",d[i]);
      if(i != 0) s[cnt++] = ' ';
      for(j=0; j<bufsize; j++)
      {
         if((buf[j] == '\0') || (buf[j] == ' ')) break;
         s[cnt++] = buf[j];
      }
      if(i < n-1) s[cnt++] = ',';
   }
   s[cnt++] = ']';
   s[cnt++] = '\0';
   return cnt;
}

int dbllist2str ( int n, double *d, char *s )
{
   int bufsize = 26;
   int cnt = 0;
   int i,j;

   s[cnt++] = '[';
   for(i=0; i<n; i++)
   {
      char buf[bufsize];
      for(j=0; j<bufsize; j++) buf[j] = ' ';
      sprintf(buf,"%.16le",d[i]);
      if(i != 0) s[cnt++] = ' ';
      for(j=0; j<bufsize; j++)
      {
         if((buf[j] == '\0') || (buf[j] == ' ')) break;
         s[cnt++] = buf[j];
      }
      if(i < n-1) s[cnt++] = ',';
   }
   s[cnt++] = ']';
   s[cnt++] = '\0';
   return cnt;
}

int itemcount ( char *s )
{
   int cnt = 1;
   int pos = 0;
   
   while(s[pos] != '\0')
      if(s[pos++] == ',') cnt++;
 
   return cnt;
}

void str2intlist ( int n, char *s, int *d )
{
   int bufsize = 16;
   char buf[bufsize];
   int spos = 0;
   int i,cnt;

   while(s[spos] != '[') spos++;
   spos++;

   for(i=0; i<n-1; i++)
   {
      cnt = 0;
      while(s[spos] == ' ') spos++;
      while(s[spos] != ',') buf[cnt++] = s[spos++];
      spos++;
      buf[cnt] = '\0';
      // printf("the buffer : %s\n",buf);
      sscanf(buf,"%d",&d[i]);
      // printf("the number : %d\n",d[i]);
   }
   cnt = 0;
   while(s[spos] == ' ') spos++;
   while(s[spos] != ']') buf[cnt++] = s[spos++];
   buf[cnt] = '\0';
   // printf("the buffer : %s\n",buf);
   sscanf(buf,"%d",&d[n-1]);
}

void str2dbllist ( int n, char *s, double *d )
{
   int bufsize = 26;
   char buf[bufsize];
   int spos = 0;
   int i,cnt;

   while(s[spos] != '[') spos++;
   spos++;

   for(i=0; i<n-1; i++)
   {
      cnt = 0;
      while(s[spos] == ' ') spos++;
      while(s[spos] != ',') buf[cnt++] = s[spos++];
      spos++;
      buf[cnt] = '\0';
      // printf("the buffer : %s\n",buf);
      sscanf(buf,"%le",&d[i]);
      // printf("the number : %.16le\n",d[i]);
   }
   cnt = 0;
   while(s[spos] == ' ') spos++;
   while(s[spos] != ']') buf[cnt++] = s[spos++];
   buf[cnt] = '\0';
   // printf("the buffer : %s\n",buf);
   sscanf(buf,"%le",&d[n-1]);
}
