/* 
 * File:   main.cpp
 * Author: daniel
 *
 * Created on 20 de octubre de 2015, 19:04
 */

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
 #include <getopt.h>
#include <math.h>
#include "bowtie.h"
#include "bwa.h"
 #include "bwtSW.h"

using namespace std;
/* Execution modes */
#define MODE_BOWTIE   0
#define MODE_BWA      1
#define MODE_BWT      2

/* Global variables for options processing */
int mode            = -1;
char structure_file[100];
char  input_file[100];
char  output_file[100];
char  reads_file[100];
int mismatches      = 0;
int gap				= 0;
int seed            = 0;
int keep            = 0;
int load            = 0;

/* Stores the values for the options in the propper global variables */
int process_opts(int argc, char *const *argv);

/* Prints the program ussage */
void print_help();

int main(int argc, char** argv) {

    string id, seq, qua, aux;
    Read *read;
    BowtieFmi * bowfmi;
    BwaFmi * bwaFmi;
    BwtSW * bwt;
    if (process_opts(argc, argv) == -1) {
        printf("[ ERROR] Error parsing options.\n");
        print_help();
        return -1;
    }    
    switch(mode){
        case MODE_BOWTIE:
            if(load)
                bowfmi = loadBowtieFMI(input_file,seed);
            else
                bowfmi = createBowtieFMI(input_file,seed);
            break;
        case MODE_BWA:
            if(load)
                bwaFmi=  loadBwaFMI(input_file, gap);
            else
                bwaFmi=  createBwaFMI(input_file, gap);
                break;
        case MODE_BWT:
            if(load)
                bwt=  loadBwtSW(input_file);
            else
                bwt=  createBwtSW(input_file);
            break;
    }

    ifstream f;
    f.open(reads_file);
    while (getline(f, id)) {
        getline(f, seq);
        getline(f, aux);
        getline(f, qua);
        read = createRead(id.c_str(), seq.c_str(), qua.c_str());
        switch(mode){
            case MODE_BOWTIE:
                inexactMatchBowtie(bowfmi, *read, mismatches);
                break;
            case MODE_BWA:
                inexactMatching(*read, mismatches, *bwaFmi) ;
                break;
            case MODE_BWT:
                alignReadBwtSW(bwt, read);
                break;
            default:
                break;
        }        
        destroyRead(read);
    }
    f.close(); 
    switch(mode){
            case MODE_BOWTIE:
                if(keep == 1)
                    saveBowtieFMI(structure_file,bowfmi);
                removeBowtieFMI(bowfmi); 
                break;
            case MODE_BWA:
                if(keep == 1)
                    saveBwaFMI(structure_file,bwaFmi);
                removeBwaFmi(bwaFmi);
                break;
            case MODE_BWT:
                if(keep == 1)
                    saveBwtSW(structure_file,bwt);
                removeBwtSW(bwt);
                break;
            default:
                break;
        }
    return 0;
}

int process_opts(int argc, char *const *argv) {
    int c;

    while (1) {
        static struct option long_options[] =
            {
                {"mode",            required_argument,  0, 'm'},
                {"input-file",      required_argument,  0, 'i'},
                {"output-file",     required_argument,  0, 'o'},
                {"reads-file",     required_argument,  0, 'r'},
                {"seed",    required_argument,  0, 's'},
                {"keep",    required_argument,  0, 'k'},
                {"load",    no_argument,  0, 'l'},
				{"gap",    no_argument,  0, 'g'},
                {"mismatches",    required_argument,  0, 'n'},
                {"help",    no_argument,  0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;

        c = getopt_long (argc, argv, "m:i:o:r:n:s:hk:lg", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
        case 'm':
            printf("[ INFO ] Mode: %s\n", optarg);
            if (!strcmp(optarg, "bowtie")) {
                mode = MODE_BOWTIE;
            } else if (!strcmp(optarg, "bwa")) {
                mode = MODE_BWA;
            }else if (!strcmp(optarg, "bwt")) {
                mode = MODE_BWT;
            } else {
                printf("[ ERROR] Unrecogniced mode: %s\n", optarg);
                return -1;
            }
            break;
        case 'i':
            printf("[ INFO ] Input file: %s\n", optarg);
            strcpy(input_file,optarg);
            break;
        case 'o':
            printf("[ INFO ] Output file: %s\n", optarg);
            strcpy(output_file,optarg);
            break;
        case 'r':
            printf("[ INFO ] Read's file: %s\n", optarg);
            strcpy(reads_file,optarg);
            break;
        case 's':
            printf("[ INFO ] Size of the seed: %d\n", atoi(optarg));
            seed = atoi(optarg);
            break;
        case 'k':
            printf("[ INFO ] Structure's file: %s\n", optarg);
            strcpy(structure_file,optarg);
            keep = 1;
            break;
        case 'l':
            load = 1;
            break;
		case 'g':
            gap = 1;
            break;
        case 'n':
            printf("[ INFO ] Number of mismatches: %d\n", atoi(optarg));
            mismatches = atoi(optarg);
            break;
        case 'h':
            print_help();
            exit(1);
            break;
        default:
            return -1;
        }
    }
    if ( mode == -1) {
        return -1;
    }
    return 0;
}

void print_help() {
    printf("[ INFO ] Usage:\n");
    printf("             aligners -m MODE -i INPUT -o OUTPUT -r READS    [-n MISMATCHES -s SEED -k FILE -l ]\n");
    printf("         Options:\n");
    printf("             -m, --mode:            Aligner mode [bowtie, bwa, bwt].\n");;
    printf("             -i, --input-file:      File with the genome sequence.\n");
    printf("             -o, --ouput-file:      File to write the alignments.\n");
    printf("             -r, --reads-file:      File with the reads.\n");
    printf("             -n, --mismatches:      Number os mismatches.\n");
    printf("             -k, --keep:            Save the structure in the file specified.\n");
    printf("             -l, --load:            Reads the structure in the input file.\n");
	printf("             -g, --gap:             Enable gap mode (only with bwa mode).\n");
    printf("             -s, --seed:            Size of the seed.\n");
    printf("             -h, --help:            Shows the different options and arguments.\n");


}

