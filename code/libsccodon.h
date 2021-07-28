/*
Title:      Parse PacBio Alignement codon library
Author:     Alexandre Schoepfer
Version:    23th July 2021, 8:30 (GMT+1)
Notes:      For Saccharomyces Cerevisiae
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char codon2aa (char *codon)
{    
    if (strcmp (codon, "TTT") == 0) return 'F';
    else if (strcmp (codon, "TTC") == 0) return 'F';
    else if (strcmp (codon, "TTA") == 0) return 'L';
    else if (strcmp (codon, "TTG") == 0) return 'L';
    else if (strcmp (codon, "TCT") == 0) return 'S';
    else if (strcmp (codon, "TCC") == 0) return 'S';
    else if (strcmp (codon, "TCA") == 0) return 'S';
    else if (strcmp (codon, "TCG") == 0) return 'S';
    else if (strcmp (codon, "TAT") == 0) return 'Y';
    else if (strcmp (codon, "TAC") == 0) return 'Y';
    else if (strcmp (codon, "TAA") == 0) return '*';
    else if (strcmp (codon, "TAG") == 0) return '*';
    else if (strcmp (codon, "TGT") == 0) return 'C';
    else if (strcmp (codon, "TGC") == 0) return 'C';
    else if (strcmp (codon, "TGA") == 0) return '*';
    else if (strcmp (codon, "TGG") == 0) return 'W';
    else if (strcmp (codon, "CTT") == 0) return 'L';
    else if (strcmp (codon, "CTC") == 0) return 'L';
    else if (strcmp (codon, "CTA") == 0) return 'L';
    else if (strcmp (codon, "CTG") == 0) return 'L';
    else if (strcmp (codon, "CCT") == 0) return 'P';
    else if (strcmp (codon, "CCC") == 0) return 'P';
    else if (strcmp (codon, "CCA") == 0) return 'P';
    else if (strcmp (codon, "CCG") == 0) return 'P';
    else if (strcmp (codon, "CAT") == 0) return 'H';
    else if (strcmp (codon, "CAC") == 0) return 'H';
    else if (strcmp (codon, "CAA") == 0) return 'Q';
    else if (strcmp (codon, "CAG") == 0) return 'Q';
    else if (strcmp (codon, "CGT") == 0) return 'R';
    else if (strcmp (codon, "CGC") == 0) return 'R';
    else if (strcmp (codon, "CGA") == 0) return 'R';
    else if (strcmp (codon, "CGG") == 0) return 'R';
    else if (strcmp (codon, "ATT") == 0) return 'I';
    else if (strcmp (codon, "ATC") == 0) return 'I';
    else if (strcmp (codon, "ATA") == 0) return 'I';
    else if (strcmp (codon, "ATG") == 0) return 'M';
    else if (strcmp (codon, "ACT") == 0) return 'T';
    else if (strcmp (codon, "ACC") == 0) return 'T';
    else if (strcmp (codon, "ACA") == 0) return 'T';
    else if (strcmp (codon, "ACG") == 0) return 'T';
    else if (strcmp (codon, "AAT") == 0) return 'N';
    else if (strcmp (codon, "AAC") == 0) return 'N';
    else if (strcmp (codon, "AAA") == 0) return 'K';
    else if (strcmp (codon, "AAG") == 0) return 'K';
    else if (strcmp (codon, "AGT") == 0) return 'S';
    else if (strcmp (codon, "AGC") == 0) return 'S';
    else if (strcmp (codon, "AGA") == 0) return 'R';
    else if (strcmp (codon, "AGG") == 0) return 'R';
    else if (strcmp (codon, "GTT") == 0) return 'V';
    else if (strcmp (codon, "GTC") == 0) return 'V';
    else if (strcmp (codon, "GTA") == 0) return 'V';
    else if (strcmp (codon, "GTG") == 0) return 'V';
    else if (strcmp (codon, "GCT") == 0) return 'A';
    else if (strcmp (codon, "GCC") == 0) return 'A';
    else if (strcmp (codon, "GCA") == 0) return 'A';
    else if (strcmp (codon, "GCG") == 0) return 'A';
    else if (strcmp (codon, "GAT") == 0) return 'D';
    else if (strcmp (codon, "GAC") == 0) return 'D';
    else if (strcmp (codon, "GAA") == 0) return 'E';
    else if (strcmp (codon, "GAG") == 0) return 'E';
    else if (strcmp (codon, "GGT") == 0) return 'G';
    else if (strcmp (codon, "GGC") == 0) return 'G';
    else if (strcmp (codon, "GGA") == 0) return 'G';
    else if (strcmp (codon, "GGG") == 0) return 'G';
    else
    {
        return 1;
        exit(1);
    }
}

double codon2fraction (char *codon) //source: https://www.genscript.com/tools/codon-frequency-table
{
    if (strcmp (codon, "TTT") == 0) return 0.59;
    else if (strcmp (codon, "TTC") == 0) return 0.41;
    else if (strcmp (codon, "TTA") == 0) return 0.28;
    else if (strcmp (codon, "TTG") == 0) return 0.29;
    else if (strcmp (codon, "TCT") == 0) return 0.26;
    else if (strcmp (codon, "TCC") == 0) return 0.16;
    else if (strcmp (codon, "TCA") == 0) return 0.21;
    else if (strcmp (codon, "TCG") == 0) return 0.10;
    else if (strcmp (codon, "TAT") == 0) return 0.56;
    else if (strcmp (codon, "TAC") == 0) return 0.44;
    else if (strcmp (codon, "TAA") == 0) return 0.48;
    else if (strcmp (codon, "TAG") == 0) return 0.22;
    else if (strcmp (codon, "TGT") == 0) return 0.63;
    else if (strcmp (codon, "TGC") == 0) return 0.37;
    else if (strcmp (codon, "TGA") == 0) return 0.30;
    else if (strcmp (codon, "TGG") == 0) return 1.0;
    else if (strcmp (codon, "CTT") == 0) return 0.13;
    else if (strcmp (codon, "CTC") == 0) return 0.06;
    else if (strcmp (codon, "CTA") == 0) return 0.14;
    else if (strcmp (codon, "CTG") == 0) return 0.11;
    else if (strcmp (codon, "CCT") == 0) return 0.31;
    else if (strcmp (codon, "CCC") == 0) return 0.15;
    else if (strcmp (codon, "CCA") == 0) return 0.42;
    else if (strcmp (codon, "CCG") == 0) return 0.12;
    else if (strcmp (codon, "CAT") == 0) return 0.64;
    else if (strcmp (codon, "CAC") == 0) return 0.36;
    else if (strcmp (codon, "CAA") == 0) return 0.69;
    else if (strcmp (codon, "CAG") == 0) return 0.31;
    else if (strcmp (codon, "CGT") == 0) return 0.14;
    else if (strcmp (codon, "CGC") == 0) return 0.06;
    else if (strcmp (codon, "CGA") == 0) return 0.07;
    else if (strcmp (codon, "CGG") == 0) return 0.04;
    else if (strcmp (codon, "ATT") == 0) return 0.46;
    else if (strcmp (codon, "ATC") == 0) return 0.26;
    else if (strcmp (codon, "ATA") == 0) return 0.27;
    else if (strcmp (codon, "ATG") == 0) return 1.0;
    else if (strcmp (codon, "ACT") == 0) return 0.35;
    else if (strcmp (codon, "ACC") == 0) return 0.22;
    else if (strcmp (codon, "ACA") == 0) return 0.30;
    else if (strcmp (codon, "ACG") == 0) return 0.14;
    else if (strcmp (codon, "AAT") == 0) return 0.59;
    else if (strcmp (codon, "AAC") == 0) return 0.41;
    else if (strcmp (codon, "AAA") == 0) return 0.58;
    else if (strcmp (codon, "AAG") == 0) return 0.42;
    else if (strcmp (codon, "AGT") == 0) return 0.16;
    else if (strcmp (codon, "AGC") == 0) return 0.11;
    else if (strcmp (codon, "AGA") == 0) return 0.48;
    else if (strcmp (codon, "AGG") == 0) return 0.21;
    else if (strcmp (codon, "GTT") == 0) return 0.39;
    else if (strcmp (codon, "GTC") == 0) return 0.21;
    else if (strcmp (codon, "GTA") == 0) return 0.21;
    else if (strcmp (codon, "GTG") == 0) return 0.19;
    else if (strcmp (codon, "GCT") == 0) return 0.38;
    else if (strcmp (codon, "GCC") == 0) return 0.22;
    else if (strcmp (codon, "GCA") == 0) return 0.29;
    else if (strcmp (codon, "GCG") == 0) return 0.11;
    else if (strcmp (codon, "GAT") == 0) return 0.65;
    else if (strcmp (codon, "GAC") == 0) return 0.35;
    else if (strcmp (codon, "GAA") == 0) return 0.70;
    else if (strcmp (codon, "GAG") == 0) return 0.30;
    else if (strcmp (codon, "GGT") == 0) return 0.47;
    else if (strcmp (codon, "GGC") == 0) return 0.19;
    else if (strcmp (codon, "GGA") == 0) return 0.22;
    else if (strcmp (codon, "GGG") == 0) return 0.12;
    else
    {
        return 1;
        exit(1);
    }
}

unsigned char codon2nnk (char *codon)
{
    if (strcmp (codon, "TTT") == 0) return 1;
    else if (strcmp (codon, "TTC") == 0) return 0;
    else if (strcmp (codon, "TTA") == 0) return 0;
    else if (strcmp (codon, "TTG") == 0) return 1;
    else if (strcmp (codon, "TCT") == 0) return 1;
    else if (strcmp (codon, "TCC") == 0) return 0;
    else if (strcmp (codon, "TCA") == 0) return 0;
    else if (strcmp (codon, "TCG") == 0) return 1;
    else if (strcmp (codon, "TAT") == 0) return 1;
    else if (strcmp (codon, "TAC") == 0) return 0;
    else if (strcmp (codon, "TAA") == 0) return 0;
    else if (strcmp (codon, "TAG") == 0) return 1;
    else if (strcmp (codon, "TGT") == 0) return 1;
    else if (strcmp (codon, "TGC") == 0) return 0;
    else if (strcmp (codon, "TGA") == 0) return 0;
    else if (strcmp (codon, "TGG") == 0) return 1;
    else if (strcmp (codon, "CTT") == 0) return 1;
    else if (strcmp (codon, "CTC") == 0) return 0;
    else if (strcmp (codon, "CTA") == 0) return 0;
    else if (strcmp (codon, "CTG") == 0) return 1;
    else if (strcmp (codon, "CCT") == 0) return 1;
    else if (strcmp (codon, "CCC") == 0) return 0;
    else if (strcmp (codon, "CCA") == 0) return 0;
    else if (strcmp (codon, "CCG") == 0) return 1;
    else if (strcmp (codon, "CAT") == 0) return 1;
    else if (strcmp (codon, "CAC") == 0) return 0;
    else if (strcmp (codon, "CAA") == 0) return 0;
    else if (strcmp (codon, "CAG") == 0) return 1;
    else if (strcmp (codon, "CGT") == 0) return 1;
    else if (strcmp (codon, "CGC") == 0) return 0;
    else if (strcmp (codon, "CGA") == 0) return 0;
    else if (strcmp (codon, "CGG") == 0) return 1;
    else if (strcmp (codon, "ATT") == 0) return 1;
    else if (strcmp (codon, "ATC") == 0) return 0;
    else if (strcmp (codon, "ATA") == 0) return 0;
    else if (strcmp (codon, "ATG") == 0) return 1;
    else if (strcmp (codon, "ACT") == 0) return 1;
    else if (strcmp (codon, "ACC") == 0) return 0;
    else if (strcmp (codon, "ACA") == 0) return 0;
    else if (strcmp (codon, "ACG") == 0) return 1;
    else if (strcmp (codon, "AAT") == 0) return 1;
    else if (strcmp (codon, "AAC") == 0) return 0;
    else if (strcmp (codon, "AAA") == 0) return 0;
    else if (strcmp (codon, "AAG") == 0) return 1;
    else if (strcmp (codon, "AGT") == 0) return 1;
    else if (strcmp (codon, "AGC") == 0) return 0;
    else if (strcmp (codon, "AGA") == 0) return 0;
    else if (strcmp (codon, "AGG") == 0) return 1;
    else if (strcmp (codon, "GTT") == 0) return 1;
    else if (strcmp (codon, "GTC") == 0) return 0;
    else if (strcmp (codon, "GTA") == 0) return 0;
    else if (strcmp (codon, "GTG") == 0) return 1;
    else if (strcmp (codon, "GCT") == 0) return 1;
    else if (strcmp (codon, "GCC") == 0) return 0;
    else if (strcmp (codon, "GCA") == 0) return 0;
    else if (strcmp (codon, "GCG") == 0) return 1;
    else if (strcmp (codon, "GAT") == 0) return 1;
    else if (strcmp (codon, "GAC") == 0) return 0;
    else if (strcmp (codon, "GAA") == 0) return 0;
    else if (strcmp (codon, "GAG") == 0) return 1;
    else if (strcmp (codon, "GGT") == 0) return 1;
    else if (strcmp (codon, "GGC") == 0) return 0;
    else if (strcmp (codon, "GGA") == 0) return 0;
    else if (strcmp (codon, "GGG") == 0) return 1;
    else
    {
        return 1;
        exit(1);
    }  
}