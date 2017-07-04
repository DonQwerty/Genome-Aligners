

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include "byte.h"

u_int8_t getConversion(char c, int position) {
    u_int8_t byte;
    switch (c) {
        case 'A':
            byte = 0x00;
            break;
        case 'C':
            byte = 0x01;
            break;
        case 'G':
            byte = 0x02;
            break;
        default:
            byte = 0x03;
            break;
    }
    return byte << (3 - position) * 2;
}

char getReverseConversion(u_int8_t byte, int position) {
    char c;
   
    u_int8_t b=byte;
    b=b << (position*2);
    b=b >> 6;
   
    switch (b) {
        case 0x00:
            c='A';
            break;
        case 0x01:
            c='C';
            break;
        case 0x02:
            c='G';
            break;
        default:
            c='T';
            break;
    }
    return c;
}

int getIntegerNuc(char c) {
    switch (c) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return -1;
    }
}

char getCharNuc(int n) {
    switch (n) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        default:
            return 'O';
    }
}