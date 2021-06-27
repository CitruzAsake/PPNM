#include <stdio.h>
#include <stdlib.h>

int main()
{
   FILE *time = fopen("time.txt", "r");
   FILE *activity = fopen("activity.txt", "r");
   FILE *unceractivity = fopen("unceractivity.txt","r");

   FILE *data = fopen("data.txt", "w");
   char c;

   if (time == NULL || activity == NULL || unceractivity == NULL || data == NULL)
   {
         puts("Could not open files");
         exit(0);
   }

   // Copy contents of first file to file3.txt
   while ((c = fgetc(time)) != EOF)
      fputc(c, data);

   // Copy contents of second file to file3.txt
   while ((c = fgetc(activity)) != EOF)
      fputc(c, data);

   while ((c = fgetc(time)) != EOF)
      fputc(c,data);

   printf("Merged data for variables into data.txt");
   
   fclose(time);
   fclose(activity);
   fclose(unceractivity);
   fclose(data);
   return 0;
}
        
