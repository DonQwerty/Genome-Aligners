
#include <cstdlib>
#include <stdio.h>
#include "dawg.h"



Node * createNode(char word,long long sp, long long ep){
	Node * node = (Node *)malloc(sizeof(Node));
	node->sa = (long long *)malloc(2*sizeof(long long));
	node->sa[0] = sp;
	node->sa[1] = ep;
	node->word = word;
	node->NChilds = 0;
	node->childs = NULL;
	node->NFathers = 0;
	node->fathers = NULL;
	node->vis = 0;
	
	return node;
}

void freeNode(Node * node){
	if(node){
		free(node->sa);
		if(node->childs)
			free(node->childs);
		if(node->fathers){
			free(node->fathers);
		}
		free(node);

	}
	return;
}

void setChildsNode(Node * node, Node ** childs, int n){
	node->childs = childs;
	node->NChilds = n;
}

void destroyDawg(Dawg * g){
	int i;
	for (i = 0 ; i < g->length; i++){
		freeNode(g->arrayNodes[i]);
	}
	free(g->arrayNodes);
	
	free(g);
}


Node * findNode(Node * node, long long sp, long long ep){
	if(node->sa[0]==sp && node->sa[1]==ep)
		return node;
	int i;
	for(i = 0; i< node->NChilds;i++){
		Node * n=findNode(node->childs[i],sp,ep);
		if(n)
			return n;
	}
	return NULL;
}

Dawg * createDawg(char * read){

	char * fmistring;
	int l,i;
	Node *root;
	Fmi * fmiRead;
	Dawg * dawg;

	dawg=(Dawg *)malloc(sizeof(Dawg));

	l = strlen(read);
	fmistring = (char *)malloc((l+2) * sizeof(char));
	strcpy(fmistring,read);
	strcat(fmistring,"$");
	fmiRead= createFMI(fmistring, l);
	free(fmistring);
	
	root = createNode('R',0,0);
	fillGraph(root,fmiRead);
	dawg->root = root;

	std::stack<Node*> p;
	postOrderReverse(root,p);

	dawg->arrayNodes = (Node **)malloc(p.size()*sizeof(Node*));
	i = 0;
	while (!p.empty()){  	
	    dawg->arrayNodes[i] = p.top();
	    p.pop();
	    if(i != 0){
	    	dawg->arrayNodes[i]->fathers = (int *)malloc(sizeof(int ) * dawg->arrayNodes[i]->NFathers);
	    }
	    i++;
  	}
  	dawg->length = i;
  	addFathers(dawg);
	removeFMI(fmiRead);
	return dawg;
}

void fillGraph(Node * node, Fmi * fmi ){
	int i,cont;
	long long sp, ep;
	char nuc[]="ACGT";
	cont = 0;
	
	Node ** childs = (Node**)malloc(sizeof(Node*)*4);
	cont = 0;
	for(i = 0 ; i < 4; i++){
		sp = getMinRank(nuc[i], fmi->ranking) - 1;
    	ep = getMaxRank(nuc[i], fmi->ranking) - 1;
    	if (sp < ep){
    		childs[cont] = createNode(nuc[i], sp, ep);
    		cont++;
    	}
	}
	setChildsNode(node, childs, cont);
	for(i = 0; i < cont; i++){
		fillGraphAux(node->childs[i],fmi, node);
	}
	return;

}

void fillGraphAux(Node * node, Fmi * fmi, Node *root){
	int i,cont;
	long long sp, ep;
	char nuc[]="ACGT";
	cont = 0;
	int val[]={0,0,0,0};
	for(i = 0 ; i < 4; i++){
		sp = getMinRank(nuc[i], fmi->ranking) + getCheckpoints(fmi, node->sa[0], nuc[i])-1 ;
        ep = getMinRank(nuc[i], fmi->ranking) + getCheckpoints(fmi, node->sa[1], nuc[i])-1 ;
    	if (sp < ep)
    		cont++;
	}
	if(cont > 0){
		Node ** childs = (Node**)malloc(sizeof(Node*)*cont);
		cont = 0;
		for(i = 0 ; i < 4; i++){
			sp = getMinRank(nuc[i], fmi->ranking) + getCheckpoints(fmi, node->sa[0], nuc[i]) -1;
        	ep = getMinRank(nuc[i], fmi->ranking) + getCheckpoints(fmi, node->sa[1], nuc[i]) -1;
	    	if (sp < ep){
	    		Node * n=findNode(root,sp,ep);
	    		//n=NULL;
	    		if(n)
	    			childs[cont] = n;
	    		else{
	    			val[cont] = 1;
	    			childs[cont] = createNode(nuc[i], sp, ep);
	    		}
	    		cont++;
	    	}
		}
		setChildsNode(node, childs, cont);
	}
	else
		setChildsNode(node, NULL, cont);
	for(i = 0; i < cont; i++){
		if(val[i]==1)
			fillGraphAux(node->childs[i],fmi,root);
	}
	return;

}

void postOrderReverse(Node * n, std::stack<Node*>& p){
	n->NFathers++;
	if(n->vis!=1){
		n->vis=1;
		for (int i = 0; i < n->NChilds; ++i)
		{
			postOrderReverse(n->childs[i],p);
		}
		p.push(n);
	}
	return;
}

void addFathers(Dawg * dawg){
	int i,j,k,sum;
	Node *n1,*n2;
	for(i = 0; i < dawg->length; i++){
		n1 = dawg->arrayNodes[i];
		sum = 0;
		for(j = 0; j< i; j++){
			n2 = dawg->arrayNodes[j];
			for(k = 0; k < n2->NChilds; k++){
				if(n2->childs[k] == n1){
					n1->fathers[sum] = j;
					sum++;
					break;
				}
			}
		}
	}
}


