/*
Title:      create nice look up table
Author:     Alexandre Schoepfer
Version:    31st August 2021, 15:45 (GMT+1)
Notes:      For Saccharomyces Cerevisiae
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BCHAR 4096
#define NCHAR 512

#define BCLEN 15

#define BARCOL 0
#define BQLCOL 1
#define MUTCOL 7



int readLine (FILE *fl, char **lin)
{
    char c;
    char *tmp;
    size_t lenChar = 0;
    size_t maxChar = BCHAR;
    
    if (!(*lin = (char *)malloc (maxChar))) return 1; 

    while ((c = fgetc(fl)) != EOF)
    {   
        (*lin)[lenChar++] = c;
        
        if (lenChar > maxChar)
        {
            if (!(tmp = (char *)realloc (*lin, maxChar += NCHAR))) return 1;
            *lin = tmp;
        }

        if (c == '\n') break;
    }
    
    (*lin)[lenChar] = '\0';

    if (c == EOF)
    {
        free (*lin);
        *lin = NULL;
        return 1;
    }
    
    return 0;
}

int allocLut (char **tok, char *str, size_t bs)
{
    char *tmp;
    size_t strl = strlen (str);
    
    if (strl > bs )
    {
        tmp = (char *)realloc (*tok, strl + 1);
        *tok = tmp;
    }
    
    strcpy(*tok, str);
    return 0;
}

int curateLut(char *blut, char **olut, char *barcd, char *barql, char *aamut, char **obarcd, char **obarql, char **oaamut, long *sbc)
{
    
    
    long phret = 0; // barcode phret score
    long ophret = 0; // old barcode phret score
    char tmp; // temporary char

    if (strlen (barcd) == BCLEN)
    {
        
        for (size_t i = 0; i < BCLEN; i++) phret += barql[i];
        for (size_t i = 0; i < BCLEN; i++)
        {
            strncpy (&tmp, *obarql+i, 1); // Did not work otherwise
            ophret += tmp; 
        }
         
        if (strcmp (barcd, *obarcd) == 0){
            if (*sbc > 0 && strcmp (aamut, *oaamut) == 0)
            {
                *sbc += 1;
            }
            else
            {
                *sbc = -1;
            }
        }
        else
        {
            if (*sbc > 0)
            {
                //printf ("%s,%s,%ld\n", *obarcd, *oaamut, *sbc);
                printf ("%s\t%ld\n", *olut, *sbc);
            }
            *sbc = 1;
        }

        if ( ! ( strcmp ( barcd, *obarcd) == 0 && phret < ophret ) )
        {
            strcpy (*obarcd, barcd); 
            strcpy (*obarql, barql);
            strcpy (*oaamut, aamut);  

            strcpy (*olut, blut);
        }

    }

    return 0;
}

int parseLut (char *lut, char *blut, char **olut, char **obarcd, char **obarql, char **oaamut, long *sbc)
{
    char *tkp; // current token of lut line
    char *barcd = (char *)malloc(BCHAR); // barcode string
    char *barql = (char *)malloc(BCHAR); // barcode quality
    char *aamut = (char *)malloc(BCHAR); // amino acid mutation
    size_t toki = 0; // token counter

    tkp = strtok (lut, "\t\n");

    while (tkp != NULL)
    {

        if (toki == BARCOL) allocLut (&barcd, tkp, BCHAR);                      
        else if (toki == BQLCOL) allocLut (&barql, tkp, BCHAR);
        else if (toki == MUTCOL) allocLut (&aamut, tkp, BCHAR);
    
        tkp = strtok (NULL, "\t\n");
        toki++;
    }
    
    curateLut(blut, olut, barcd, barql, aamut, obarcd, obarql, oaamut, sbc);

    free (barcd);
    free (barql);
    free (aamut);

    return 0;
}

int readLut (FILE *lutF, char **lut)
{
    long sbc = 0; // barocde multiplicate counter
    char *oaamut = (char *)malloc(BCHAR); // old amino acid mutation
    char *obarcd = (char *)malloc(BCHAR); // old barcode string
    char *obarql = (char *)malloc(BCHAR); // old barcode quality
    strcpy (oaamut, " ");
    strcpy (obarcd, " ");
    strcpy (obarql, "               ");

    char *olut = (char *)malloc (BCHAR*2); // old lut
    strcpy (olut, " ");

    char *blut = (char *)malloc (BCHAR*2); // unmodified lut
    strcpy (blut, " ");

    while (!(readLine (lutF, lut)))
    {          
        strcpy (blut, *lut);
        blut[strlen (blut)-2] = '\0';
        
        parseLut (*lut, blut, &olut, &obarcd, &obarql, &oaamut, &sbc);

        free (*lut);
        *lut = NULL;
    }
    //if (sbc !=0) printf ("%s,%s,%ld\n", obarcd, oaamut, sbc);
    if (sbc !=0) printf ("%s\t%ld\n", olut, sbc);
    
    free (oaamut);
    free (obarcd);
    free (obarql);

    free (olut);
    free (blut);
    
    return 0;
}

int main (int argc, char *argv[])
{

    FILE *lutFile = fopen (argv[argc-1], "r");

    char *lut = NULL;

    if (argc == 1)
    {
        printf ("Usage: program <raw_lookup_table.tsv>\n");
        printf ("Try 'mlut --help' for more information.\n");
        exit (1);
    }
    
    for (size_t i = 0; i < argc; i++)
    {
        if ( strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
        {
            printf ("Usage: program <raw_lookup_table.tsv>\n");
            printf ("Only use for .tsv files.\n");
            printf ("Barcode length default: 15\n");
            printf ("\nDocumentation: <https://github.com/Nash-Lab/ngs-ml>\n");
            exit (0);
        }
    }

    if (lutFile != NULL) readLut (lutFile, &lut);  
    fclose (lutFile);
    
    free (lut);

    return 0;
}