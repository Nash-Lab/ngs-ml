/*
Title:   Read Illumina Barcode
Author:   Alexandre Schoepfer
Version:  28rd July 2021, 10:15 (GMT+1)
Notes:
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define NOT_FOUND 0
#define FOUND 1
#define LOW_QUALITY 2
#define WRONG_SIZE 3

unsigned long lookupTableLength = 0;
unsigned long dataLength  = 0;

double qualityCutoff = 20.0;
unsigned int barcodel = 15;
unsigned int maxDataLength = 0;
char delim = '\t';

unsigned int getMaxLineLength (FILE* file, unsigned long * length)
{
    char c;
    unsigned int count = 0;
    unsigned int max = 0;

    rewind (file);

    if (file != NULL)
    {
        do
        {        
            c = fgetc (file);
            count++;
            
            if (c == '\n')
            {
                if (max < count)
                {
                    max = count;
                }
                count = 1;
                *length += 1 ;
            }

        } while (c != EOF);
    }

    return (max);
}

unsigned short getStartingIndex (FILE* file, unsigned int maxLineLength)
{

    char line[maxLineLength];
    unsigned short startingIndex = 0;

    rewind (file);

    if (file != NULL)
    {
                
        while (fgets (line, maxLineLength, file) != NULL)
        {   

            if (strcmp (line, "+\n") == 0 && startingIndex > 1)
            {
                return (startingIndex - 1);
            }
            else if (strcmp (line, "+\n") == 0 && startingIndex < 2)
            {
                printf("X Error, '+' line should not be here at line %d.\n", 
                       startingIndex+1);
                exit(1);
            }
            
            startingIndex++;   
        } 
    }

    return (0);
}

int storeLookupTable (FILE * file, unsigned int maxLineLength,
                      char * luts, char * muts)
{
    
    char line[maxLineLength];
    rewind (file);

    unsigned long i = 0;
    unsigned long j = 0;
    unsigned long k = 1;

    char sdelim[2];
    sdelim[0] = delim;
    sdelim[1] = '\0';

    if (file != NULL)
    {
                
        while (fgets (line, maxLineLength, file) != NULL)
        {   
            if ( strncmp(line + barcodel, sdelim, 1) == 0 )
            {
                strncpy (luts + i, line, barcodel);
                strncpy (muts + j, line + barcodel + 1,
                        maxLineLength - barcodel - 1);

                i += barcodel + 1;
                j += maxLineLength - barcodel -1;
                k++;
            }
            else
            {
                printf(
                    "X Error, barcode with wrong size found in lookup table.\n"
                    );
                printf(
                    "Please format the lookuptable like: '<barcode>,<rest>...'"
                    );
                printf("\nBarcode Size should be %d. Check line %ld.\n", barcodel, k);
                exit(1);
            }
            
        } 
    }

    return (0);
}

unsigned long findMatchOld (char * line, char * luts, unsigned long lutLength)
{
    int cmp;

    for (unsigned long i = 0; i < lutLength; i++)
    {
        cmp = strncmp (line, luts + (i * (barcodel + 1)), barcodel + 1);

        if (cmp == 0)
        {
            return (i);
        }
    }

    return (lutLength);
    
}

unsigned long findMatch (char * line, char * luts, unsigned long lutLength)
{
    
    long lowerB = 0;
    long upperB = lutLength-1;

    int cmp;

    do   
    {
        cmp = strncmp (line, luts + 
                       (( lowerB + upperB) / 2) * (barcodel + 1), barcodel + 1);
        
        if (cmp == 0)
        {
            return ((( lowerB + upperB ) / 2));
        }
        else if (upperB <= lowerB)
        {
            return (lutLength);
        }
        else if (cmp < 0)
        {
            upperB = (( lowerB + upperB) / 2) - 1;
        }
        else if (cmp > 0)
        {
            lowerB = (( lowerB + upperB) / 2) + 1;
        }      
    } while (lowerB <= upperB);

    return (lutLength);
}

int parseSequences (FILE* file, unsigned int maxLineLength,
                    unsigned long startIndex, char * luts, char * muts,
                    unsigned int maxLutLength, unsigned long lutLength)
{
    
    unsigned int sl;
    long sum;
    double avg;
    unsigned long matchIndex;

    char line [maxLineLength];
    char tempLine [maxLineLength];
    char * tempLineP;
    unsigned long index = startIndex;

    rewind (file);

    if (file != NULL)
    {
        //printf ("tag,barcode,aa_mutation,n_aa_substitutions\n");    
        while (fgets (line, maxLineLength, file) != NULL)
        {   
            if (index % 4 == 2)
            {
                tempLineP = strtok (line, "\n");
                strncpy (tempLine, tempLineP, maxLineLength);              
            }
            else if (index % 4 == 0)
            {
                if (barcodel != strlen (tempLine))
                {
                    printf ("%d%c%s\n", WRONG_SIZE, delim, tempLine);
                }
                else
                {
                    sl = strlen (line);
                    sum  = 0;

                    for (unsigned int i = 0; i < sl; i++)
                    {
                        sum += line [i] - 33;

                    }
                    avg = sum / sl;
                    if (avg < qualityCutoff)
                    {
                        printf ("%d%c%s\n", LOW_QUALITY, delim, tempLine);
                    }
                    else
                    {
                        matchIndex = findMatch(tempLine, luts, lutLength);
                        if (matchIndex < lutLength)
                        {
                            printf ("%d%c%s%c%s", FOUND, delim,
                                    luts + (matchIndex * (barcodel + 1)),
                                    delim,
                                    muts + (matchIndex *
                                    (maxLutLength - barcodel - 1)));
                        }
                        else
                        {
                            printf ("%d%c%s\n", NOT_FOUND, delim, tempLine);
                        }
                    }
                }            
            }
            index++;
        } 
    }

    return (0);
}

int extractSequences (FILE* file, unsigned int maxLineLength,
                    unsigned long startIndex)
{
    
    unsigned int sl;
    long sum;
    double avg;
    unsigned long matchIndex;

    char line [maxLineLength];
    char tempLine [maxLineLength];
    char * tempLineP;
    unsigned long index = startIndex;

    rewind (file);

    if (file != NULL)
    {
        printf ("barcode,quality\n");    
        while (fgets (line, maxLineLength, file) != NULL)
        {   

            if (index % 4 == 2)
            {
                tempLineP = strtok (line, "\n");
                strncpy (tempLine, tempLineP, maxLineLength);              
            }
            else if (index % 4 == 0)
            {
                sl = strlen (line);
                sum  = 0;

                for (unsigned int i = 0; i < sl; i++)
                {
                    sum += line [i] - 33;

                }
                avg = sum / sl;

                printf ("%s,%d\n", tempLine, (int) avg);        
            }
            index++;
        } 
    }

    return (0);
}

int main (int argc, char* argv[])
{
    FILE* lookupTable;
    FILE* data;
    
    short useLut = 0;

    if (argc == 1)
    {
        printf ("Usage: program [OPTION]... [lookuptable] <fastq>\n");
        printf ("Try 'rib --help' for more information.\n");
        exit (1);
    }
    
    for (unsigned int i = 0; i < argc; i++)
    {
        if ( strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
        {
            printf ("Usage: program [OPTION]... [lookuptable] <fastq>\n");
            printf ("\n");
            printf ("  -b, --barcode    Length of barcode. Default: 15\n");
            printf ("  -d, --delimiter  Tag delimiter. Default: '\\t'\n");
            printf ("  -h, --help       This prompt.\n");
            printf ("  -l, --length     Length of longest line in fastq file.");
            printf ("\n                   Default: automatic\n");
            printf ("  -q, --quality    Cutoff value for the mean Phred\n");
            printf ("                   quality score of the sequence.\n");
            printf ("                   Default: 20\n");
            printf ("  -t, --lut        Use lookup table\n");
            printf ("\nDocumentation: <https://github.com/Nash-Lab/ngs-ml>\n");
            exit (0);
        }
        else if ( strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--lut") == 0 )
        {
            useLut = i + 1;
        }
    }
    
    for (unsigned int i = 0; i < argc; i++)
    {
        if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--barcode") == 0)
        {
            sscanf (argv[i+1], "%d", & barcodel);
        }
        else if (strcmp(argv[i], "-q") == 0 ||
            strcmp(argv[i], "--quality") == 0)
        {
            sscanf (argv[i+1], "%lf", & qualityCutoff);
        }
        else if (strcmp(argv[i], "-l") == 0 || 
            strcmp(argv[i], "--length") == 0)
        {
            sscanf (argv[i+1], "%d", & maxDataLength);
        }
        else if (strcmp(argv[i], "-d") == 0 || 
            strcmp(argv[i], "--delimiter") == 0)
        {
            sscanf (argv[i+1], "%c", & delim);
        }
    }

    // Fastq file
    data = fopen (argv[argc-1], "r");
    
    if (maxDataLength == 0)
    {
        maxDataLength = getMaxLineLength (data, & dataLength);
    }
    else
    {
        maxDataLength++;
    }
    unsigned long startIndex = getStartingIndex (data, maxDataLength);
    if (startIndex == 0)
    {
        printf ("X Error, invalid fastq file, couldn't find '+' line.");
        exit(1);
    }

    // Lookup Table
    if ( useLut > 0 )
    {
        lookupTable = fopen (argv[useLut], "r");
        
        unsigned int maxLutLength = getMaxLineLength (lookupTable,
                                                    & lookupTableLength);

        char * luts = (char *)calloc (lookupTableLength, (barcodel + 1));
        if (luts == NULL)
        {
            printf ("X Error, not enough memory for lookup table barcodes.");
            exit(1);
        }
        char * muts = (char *)calloc (lookupTableLength,
        (maxLutLength - barcodel));
        if (muts == NULL)
        {
            printf ("X Error, not enough memory for lookup table mutations.");
            exit(1);
        }
        storeLookupTable(lookupTable, maxLutLength, luts, muts);

        fclose (lookupTable);

        int ps = parseSequences (data, maxDataLength, startIndex, luts, muts,
        maxLutLength, lookupTableLength);
    }
    else
    {
        int es = extractSequences (data, maxDataLength, startIndex);
    }

    fclose (data);

}
