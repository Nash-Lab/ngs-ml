/*
Title:      Parse PacBio Alignement main
Author:     Alexandre Schoepfer
Version:    28th July 2021, 10:15 (GMT+1)
Notes:  - TODO: memory allocation checks
        - TODO: dynamic memory for parseCS
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "libsccodon.h"

#define BCHAR 16384
#define NCHAR 512
#define BALIN 1024

#define INDCOL 3 // Index column location
#define CIGCOL 5 // Cigar string column location
#define SEQCOL 9 // Base sequence column location in .sam
#define QALCOL 10 // Sequence quality column location in .sam
#define ALICOL 20 // Obsolete

#define START 600
#define END 2102

size_t s = 1;
size_t e = 1;

int readLine (FILE *fl, char **ref)
{
    char c;
    char *tmp;
    size_t lenChar = 0;
    size_t maxChar = BCHAR;
    
    if (!(*ref = (char *)malloc (maxChar))) return 1; 

    while ((c = fgetc(fl)) != EOF)
    {   
        
        (*ref)[lenChar++] = c;
        
        if (lenChar > maxChar)
        {
            if (!(tmp = (char *)realloc (*ref, maxChar += NCHAR))) return 1;
            *ref = tmp;
        }
        if (c == '\n') break;
    }
    
    (*ref)[lenChar] = '\0';

    if (c == EOF)
    {
        free (*ref);
        *ref = NULL;
        return 1;
    }
    //free (tmp);
    return 0;
}

int allocSam (char **tok, char *str, size_t bs)
{
    char *tmp;
    size_t strl = strlen (str);
    
    if (strl > bs )
    {
        tmp = (char *)realloc (*tok, strl + 1);
        *tok = tmp;
    }
    
    //free (tmp);

    strcpy(*tok, str);
    return 0;
}

int readRef (FILE *refF, char **ref)
{
    char refLinePointer = 0; 
    
    while (!(readLine (refF, ref)))
    {

        if (refLinePointer) break;
        if (*ref[0] == '>') refLinePointer = 1;
  
        free (*ref);
        *ref = NULL;
    }
    
    return 0;
}

int parseCS (char *ref, long ind, long cig, char *seq, char *qal, char *ali)
{   
    char mut1;
    char mut2;
    char *lts = (char *)malloc(BCHAR); // TODO dynamic size
    char *tmp = (char *)malloc(BCHAR); // TODO dynamic size
    char *fps = (char *)malloc(BCHAR); // TODO dynamic size
    char *barcodeSeq = (char *)malloc(BCHAR); // TODO dynamic size
    char *barcodeQal = (char *)malloc(BCHAR); // TODO dynamic size
    char *genMutsSeq = (char *)malloc(BCHAR); // TODO dynamic size
    char *genMutsQal = (char *)malloc(BCHAR); // TODO dynamic size
    char *genIndlSeq = (char *)malloc(BCHAR); // TODO dynamic size
    char *genIndlQal = (char *)malloc(BCHAR); // TODO dynamic size
    long aliptr = 5; // Alignement pointer TODO check if str genuine CS str
    long seqptr = cig; // Sequence pointer
    long seqind = ind; // Sequence reference index
    long seqshi = 0; // Sequence index shifter
    size_t bci = 0; // Barcode index
    size_t msi = 0; // Mutation quality index
    size_t mqi = 0; // Mutation quality index 
    size_t isi = 0; // Indels quality index
    size_t iqi = 0; // Indels quality index 
    
    size_t ii = 0;
    double seqqal = 0;
    double genqal = 0;
    double barqal = 0;


    char bcd = 0;
    char rcd[] = "000";
    char qcd[] = "000";
    long icd[] = {0,0,0};
    long iaa = 0;
    for (size_t i = 0; i < 3; i++) icd[i] = 0;

    char *cdStr = (char *)malloc(BCHAR); // TODO dynamic size
    char *aaStr = (char *)malloc(BCHAR); // TODO dynamic size
    char *nkStr = (char *)malloc(BCHAR); // TODO dynamic size
    char *frStr = (char *)malloc(BCHAR); // TODO dynamic size
    size_t cdi = 0;
    size_t aai = 0;
    size_t nki = 0;
    size_t fri = 0;

    size_t naa = 0;

    strcpy (lts, "");
    strcpy (tmp, "");
    strcpy (fps, "");
    strcpy (barcodeSeq, " ");
    strcpy (barcodeQal, " ");
    strcpy (genMutsSeq, " ");
    strcpy (genMutsQal, " ");
    strcpy (genIndlSeq, " ");
    strcpy (genIndlQal, " ");
    strcpy (cdStr, " ");
    strcpy (aaStr, " ");
    strcpy (nkStr, " ");
    strcpy (frStr, " ");

    while (ali[aliptr] != '\0')
    {

        if (strncmp (ali+aliptr, ":", 1) == 0)
        {
            sscanf (ali+aliptr+1, "%ld", &seqshi);
            seqptr += seqshi;
            seqind += seqshi;
        }
        else if (strncmp (ali+aliptr, "+", 1) == 0)
        {
            sscanf (ali+aliptr+1, "%[acgtn]", tmp);
            for (size_t i = 0; tmp[i] != '\0'; i++) tmp[i] = toupper (tmp[i]);
            sprintf (lts, "+%ld%s ",seqind,tmp);
            strcpy (genIndlSeq + isi, lts);
            isi += strlen (lts); 
            
            for (size_t i = 0; strncmp (seq+seqptr+i, tmp, strlen (tmp)) == 0; i+=strlen (tmp))
            {
                strncpy (genIndlQal+iqi,qal+seqptr+i,strlen (tmp));
                iqi += strlen (tmp);
                genIndlQal[iqi] = '\0';
            }
            genIndlQal[iqi++] = ' ';
            
            seqptr += strlen (tmp);
        }
        else if (strncmp (ali+aliptr, "-", 1) == 0)
        {
            sscanf (ali+aliptr+1, "%[acgtn]", tmp);
            for (size_t i = 0; tmp[i] != '\0'; i++) tmp[i] = toupper (tmp[i]);
            sprintf (lts, "-%ld%s ",seqind,tmp);
            strcpy (genIndlSeq + isi, lts);
            isi += strlen (lts); 
            
            for (size_t i = 0; tmp[i] != '\0'; i++) genIndlQal[iqi++] = '~';
            for (size_t i = 0; strncmp (seq+seqptr+i, tmp, strlen (tmp)) == 0; i+=strlen (tmp))
            {
                strncpy (genIndlQal+iqi,qal+seqptr+i,strlen (tmp));
                iqi += strlen (tmp);
                genIndlQal[iqi] = '\0';
            }
            genIndlQal[iqi++] = ' ';

            seqind += strlen (tmp);
        }
        else if (strncmp (ali+aliptr, "*n", 2) == 0)
        {
            barcodeSeq[bci] = ali[aliptr+2]; // seq[seqptr]
            barcodeQal[bci] = qal[seqptr];
            bci++;
            seqptr++;
            seqind++;
        }
        else if (strncmp (ali+aliptr, "*", 1) == 0)
        {
            sscanf (ali+aliptr+1, "%c%c", &mut1, &mut2);
            sprintf (lts, "%c%ld%c ", mut1-32, seqind, mut2-32);
            strcpy (genMutsSeq + msi, lts);
            msi += strlen (lts); 
    
            genMutsQal[mqi++] = qal[seqptr];
            genMutsQal[mqi++] = ' ';

            if (seqind >= START && seqind <= END)
            {
                rcd[(seqind-START)%3] = mut1-32;
                qcd[(seqind-START)%3] = mut2-32;
                icd[(seqind-START)%3] = seqind;
                bcd = 1;
                iaa = (seqind-START+3)/3;
            }

            seqptr++;
            seqind++;
        }
        else if (strncmp (ali+aliptr, "~", 1) == 0) exit (2); //There's a ~ somewhere        

        //Codon information
        if (bcd && iaa < (seqind-START+3)/3)
        {
            long ics = 0;
            
            for (size_t i = 0; i < 3; i++)
            {
                if (icd[i] != 0) ics = icd[i] - i;
            }

            for (size_t i = 0; i < 3; i++)
            {
                if (rcd[i] == '0') rcd[i] = ref[ics+i-1];
            }

            for (size_t i = 0; i < 3; i++)
            {
                if (qcd[i] == '0') qcd[i] = ref[ics+i-1];
            }
            
            sprintf (lts, "%s%ld%s ", rcd, ics, qcd);
            strcpy (cdStr + cdi, lts);
            cdi += strlen (lts); 

            sprintf (lts, "%c%ld%c ", codon2aa (rcd), iaa, codon2aa (qcd));
            strcpy (aaStr + aai, lts);
            aai += strlen (lts);             

            sprintf (lts, "%d ", codon2nnk (qcd));
            strcpy (nkStr + nki, lts);
            nki += strlen (lts);

            sprintf (lts, "%.2lf ", codon2fraction (qcd));
            strcpy (frStr + fri, lts);
            fri += strlen (lts);

            strcpy (rcd, "000");
            strcpy (qcd, "000");
            for (size_t i = 0; i < 3; i++) icd[i] = 0;
            

            bcd = 0;
        }

        aliptr++;
    }
    if (bci != 0) barcodeSeq[bci] = '\0';
    if (bci != 0) barcodeQal[bci] = '\0';
    if (msi != 0) genMutsSeq[msi] = '\0';
    if (mqi != 0) genMutsQal[mqi] = '\0';
    if (isi != 0) genIndlSeq[isi] = '\0';
    if (iqi != 0) genIndlQal[iqi] = '\0';

    if (cdi != 0) cdStr[cdi] = '\0';
    if (aai != 0) aaStr[aai] = '\0';
    if (nki != 0) nkStr[nki] = '\0';
    if (fri != 0) frStr[fri] = '\0';
    

    ii=0;
    while (barcodeQal[ii] != '\0')
    {
        barqal += barcodeQal[ii];
        ii++;
    }
    barqal /= ii;
    barqal = pow (10, (-(barqal-33)/10) );


    ii=0;
    while (qal[ii] != '\0')
    {
        seqqal += qal[ii];
        ii++;
    }
    seqqal /= ii;
    seqqal = pow (10, (-(seqqal-33)/10) );

    ii=START-1;
    while (ii < END)
    {
        genqal += qal[ii];
        ii++;
    }
    genqal /= (ii-START+1);
    genqal = pow (10, (-(genqal-33)/10) );

    for (size_t i = 0; barcodeSeq[i]!='\0'; i++) barcodeSeq[i] = toupper (barcodeSeq[i]);
    
    for (size_t i = 0; aaStr[i]!='\0'; i++)
    {
        if (aaStr[i] == ' ') naa++;
    }
    
    
    //printf ("%.4lf,%.4lf,%.4lf\n", barqal,genqal,seqqal);
    printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%ld\t%s\t%s\t%.4lf\t%.4lf\t%.4lf\n", barcodeSeq, barcodeQal, genMutsSeq, genMutsQal, genIndlSeq, genIndlQal, cdStr, aaStr, naa, nkStr, frStr, barqal, genqal, seqqal);
    //printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", barcodeSeq, barcodeQal, genMutsSeq, genMutsQal, genIndlSeq, genIndlQal, cdStr, aaStr, nkStr);
    
    free (barcodeSeq);
    free (barcodeQal);
    free (genMutsSeq);
    free (genMutsQal);
    free (genIndlSeq);
    free (genIndlQal);
    free (tmp);
    free (lts);
    free (fps);

    free (cdStr);
    free (aaStr);
    free (nkStr);
    free (frStr);

    return 0;
}

int parseSam (char *ref, char **sam)
{
    char *tkp;
    char stc;
    char *seq = (char *)malloc(BCHAR);
    char *qal = (char *)malloc(BCHAR);
    char *ali = (char *)malloc(BCHAR);
    long ind = 0;
    long cig = 0;
    size_t toki = 0;

    tkp = strtok (*sam, "\t\n");

    while (tkp != NULL)
    {

        if (toki == INDCOL)
        {
            sscanf (tkp, "%ld", &ind); //Use atoi ?
        }
        if (toki == CIGCOL)
        {
            sscanf (tkp, "%ld%c", &cig, &stc); //Use atoi ?
        }
        if (toki == SEQCOL) allocSam (&seq, tkp, BCHAR);                      
        else if (toki == QALCOL) allocSam (&qal, tkp, BCHAR);
        else if (strncmp (tkp,"cs:Z:", 5) == 0) allocSam (&ali, tkp, BCHAR);

        tkp = strtok (NULL, "\t\n");
        toki++;
    }
    
    if (stc != 'S') cig = 0;
    parseCS (ref, ind, cig, seq, qal, ali);

    free (seq); 
    free (qal);
    free (ali);

    
    return 0;
}

int readSam (FILE *samF, char *ref, char **sam)
{
    char isMeta = 1;
    
    while (!(readLine (samF, sam)))
    {  

        if (*sam[0] != '@') isMeta = 0;
        if (!isMeta) parseSam (ref, sam);

        free (*sam);
        *sam = NULL;
    }
    
    return 0;
}

int main (int argc, char *argv[])
{
    FILE *refFile = fopen (argv[argc-2], "r");
    FILE *samFile = fopen (argv[argc-1], "r");

    char *ref = NULL;
    char *sam = NULL;

    if (argc == 1)
    {
        printf ("Usage: program <reference.fa> <alignment.sam>\n");
        printf ("Try 'ppba --help' for more information.\n");
        exit (1);
    }
    
    for (size_t i = 0; i < argc; i++)
    {
        if ( strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
        {
            printf ("Usage: program <reference.fa> <alignement.sam>\n");
            printf ("\n");
            printf ("For Saccharomyces Cerevisiae.\n");
            printf ("Currently, start is set at 600 and end at 2102.\n");
            printf ("\nDocumentation: <https://github.com/Nash-Lab/ngs-ml>\n");
            exit (0);
        }
    }
    
    if (refFile != NULL) readRef (refFile, &ref);
    fclose (refFile);

    if (samFile != NULL) readSam (samFile, ref, &sam);  
    fclose (samFile);
    
    free (ref);
    free (sam);

    return 0;
}