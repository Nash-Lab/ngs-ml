/*
Title:      Process illumina data
Author:     Alexandre Schoepfer
Version:    28th July 2021, 10:15 (GMT+1)
Notes:      For Saccharomyces Cerevisiae
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BCHAR 4096
#define NCHAR 512

#define BCLEN 15

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


int allocIll (char **tok, char *str, size_t bs)
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

size_t dlen (long num)
{
    size_t nod = 0;
    
    do
    {
        nod++;
        num /= 10;
    } while (num);
    
    return nod;
    
}

int parseCig (char *cig, char *ill, char *qal)
{
    long cigptr = 0;
    long illptr = 0;

    long shi = 0;
    char fla = ' ';

    char *barc = (char *)malloc(BCHAR);
    char *qual = (char *)malloc(BCHAR);

    while (cig[cigptr] != '\0')
    {
        sscanf (cig+cigptr, "%ld%c", &shi, &fla);
        
        if (fla == 'M' && illptr > 0)
        {
            strncpy (barc, ill + illptr, shi);
            strncpy (qual, qal + illptr, shi);
            barc[shi] = '\0';
            qual[shi] = '\0';
            printf ("@seq1\n%s\n+\n%s\n", barc, qual);
            break;
        }
        else
        {
            illptr += shi;
        }

        cigptr += dlen (shi) + 1;
    }
    
    free (barc);
    free (qual);

    return 0;
}

int parseIll (char **ill)
{
    char *tkp;
    char *cigar = (char *)malloc(BCHAR);
    char *illum = (char *)malloc(BCHAR);
    char *quali = (char *)malloc(BCHAR);
    size_t toki = 0;

    tkp = strtok (*ill, "\t\n");

    while (tkp != NULL)
    {

        if (toki == 5) allocIll (&cigar, tkp, BCHAR);
        if (toki == 9) allocIll (&illum, tkp, BCHAR);
        if (toki == 10) allocIll (&quali, tkp, BCHAR);                      

        tkp = strtok (NULL, "\t\n");
        toki++;
    }
    
    parseCig (cigar, illum, quali);

    free (cigar);
    free (illum);
    free (quali);

    return 0;
}

int readIll (FILE *illF, char **ill)
{

    char isMeta = 1;
    
    while (!(readLine (illF, ill)))
    {  
        if (*ill[0] != '@') isMeta = 0;
        if (!isMeta) parseIll (ill);

        free (*ill);
        *ill = NULL;
    }

    return 0;
}

int main (int argc, char *argv[])
{

    FILE *illFile = fopen (argv[argc-1], "r");

    char *ill = NULL;

    if (argc == 1)
    {
        printf ("Usage: program <alignment.sam>\n");
        printf ("Try 'pib --help' for more information.\n");
        exit (1);
    }
    
    for (size_t i = 0; i < argc; i++)
    {
        if ( strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
        {
            printf ("Usage: program <alignment.sam>\n");
            printf ("\nDocumentation: <https://github.com/Nash-Lab/ngs-ml>\n");
            exit (0);
        }
    }

    if (illFile != NULL) readIll (illFile, &ill);  
    fclose (illFile);
    
    free (ill);

    return 0;
}

