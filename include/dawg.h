/* 
 * File:   newfile.h
 * Author: daniel
 *
 * Created on 10 de abril de 2016, 21:51
 */


#ifndef DAWG_H
#define	DAWG_H
#include <stack> 
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <string.h>
#include "fmi.h"

 typedef struct Node_
{
	int NChilds;
    Node_ ** childs;
    long long * sa;
    char  word; 
    int vis;
    int * fathers;
    int NFathers;
	
} Node;

typedef struct Dawg_
{
  Node * root;
  Node ** arrayNodes;
  int length;
	
} Dawg;

Node * createNode(char word,long long sp, long long ep);

void freeNode(Node * node);

void setChildsNode(Node * node, Node ** childs, int n);

void destroyDawg(Dawg * g);

void destroyNode(Node * n);

Node * findNode(Node * node, long long sp, long long ep);

Dawg * createDawg(char * read);

void fillGraph(Node * node, Fmi * fmi );

void fillGraphAux(Node * node, Fmi * fmi, Node *root);

void postOrderReverse(Node * n, std::stack<Node*>& p);

void addFathers(Dawg * dawg);

#endif	/* DAWG_H */
