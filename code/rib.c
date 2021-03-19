/*
Title:   Read Illumina Barcode
Author:   Alexandre Schoepfer
Version:  19th March 2021, 10:21 (GMT+1)
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

            if (strcmp (line, "+\n") == 0)
            {
                return (startingIndex - 1);
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

    if (file != NULL)
    {
                
        while (fgets (line, maxLineLength, file) != NULL)
        {   
            if ( strncmp(line + barcodel, ",", 1) == 0 )
            {
                strncpy (luts + i, line, barcodel);
                strncpy (muts + j, line + barcodel + 1,
                        maxLineLength - barcodel -1);

                i += barcodel + 1;
                j += maxLineLength - barcodel -1;
            }
            else
            {
                printf(
                    "X Error, barcode with wrong size found in lookup table.\n"
                    );
                printf(
                    "Please format the lookuptable like: '<barcode>,<rest>...'"
                    );
                printf("\nBarcode Size should be %d.\n", barcodel);
                exit(1);
            }
            
        } 
    }

    return (0);
}

unsigned long findMatch (char * line, char * luts, unsigned long lutLength)
{
    
    unsigned long lowerB = 0;
    unsigned long upperB = lutLength-1;

    int cmp;

    do   
    {
        cmp = strncmp (line, luts + 
                       (( upperB + lowerB) / 2) * (barcodel + 1), barcodel);
        
        if (cmp == 0)
        {
            return ((( upperB + lowerB) / 2));
        }
        else if (upperB - lowerB == 0)
        {
            return (lutLength);
        }
        else if (cmp < 0)
        {
            upperB = (( upperB + lowerB) / 2) - 1;
        }
        else if (cmp > 0)
        {
            lowerB = (( upperB + lowerB) / 2) + 1;
        }      
    } while (lowerB < upperB);

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
        printf ("tag,barcode,aa_mutation,n_aa_substitutions\n");    
        while (fgets (line, maxLineLength, file) != NULL)
        {   

            if (index % 4 == 2)
            {
                tempLineP = strtok (line, "\n");
                strncpy (tempLine, tempLineP, maxLineLength);              
            }
            else if (index % 4 == 0)
            {
                if (barcodel < strlen (tempLine))
                {
                    printf ("%d,%s,,\n", WRONG_SIZE, tempLine);
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
                        printf ("%d,%s,,\n", LOW_QUALITY, tempLine);
                    }
                    else
                    {
                        matchIndex = findMatch(tempLine, luts, lutLength);
                        if (matchIndex < lutLength)
                        {
                            printf ("%d,%s,%s", FOUND,
                                    luts + (matchIndex * (barcodel + 1)),
                                    muts + (matchIndex *
                                    (maxLutLength - barcodel - 1)));
                        }
                        else
                        {
                            printf ("%d,%s,,\n", NOT_FOUND, tempLine);
                        }
                    }
                }            
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
    
    if (argc == 1)
    {
        printf ("Usage: program [OPTION]... <lookuptable> <fastq>\n");
        printf ("Try 'rib --help' for more information.\n");
        exit (1);
    }
    
    for (unsigned int i = 0; i < argc; i++)
    {
        if ( strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
        {
            printf ("Usage: program [OPTION]... <lookuptable> <fastq>\n");
            printf ("\n");
            printf ("  -b, --barcode    Length of barcode. Default: 15\n");
            printf ("  -h, --help       This prompt.\n");
            printf ("  -l, --length     Length of longest line in fastq file.");
            printf ("\n                   Default: automatic\n");
            printf ("  -q, --quality    Cutoff value for the mean Phred\n");
            printf ("                   quality score of the sequence.\n");
            printf ("                   Default: 20\n");
            printf ("\nDocumentation: <https://github.com/Nash-Lab/ngs-ml>\n");
            exit (0);
            break;
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
        
    }


    lookupTable = fopen (argv[argc-2], "r");
    data = fopen (argv[argc-1], "r");

    // Lookup Table
    unsigned int maxLutLength = getMaxLineLength (lookupTable,
                                                  & lookupTableLength);
  
    char * luts = (char *)calloc (lookupTableLength, (barcodel + 1));
    if (luts == NULL)
    {
        printf ("X Not enough memory for lookup table barcodes.");
        exit(1);
    }
    char * muts = (char *)calloc (lookupTableLength, (maxLutLength - barcodel));
    if (muts == NULL)
    {
        printf ("X Not enough memory for lookup table mutations.");
        exit(1);
    }
    storeLookupTable(lookupTable, maxLutLength, luts, muts);
  
    // Fastq file
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
    
    int ps = parseSequences (data, maxDataLength, startIndex, luts, muts,
                             maxLutLength, lookupTableLength);

}