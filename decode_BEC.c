#include <stdio.h>
#include "LDPC.h"

void printarray(int * a, int n)
{
	for(int i=0;i<n;i++)
		printf("%d ",a[i]);
	printf("\n");
}
void squarearray(int * a, int n)
{
	for(int i=0;i<n;i++)
		a[i]*=a[i];
}

/*
 * Argument descriptions :
 * bsi : start indices for list of indices of edges for each bit node
 * bei : end indices for list indices of edges for each bit node
 * bls : list of indices of edges grouped together by bit node to which
 *       they are connected
 * csi : start indices for list of indices of edges for each check node
 * cei : end indices for list indices of edges for each check node
 * n : Block length
 * nk : Rows in parity check matrix
 * edg : working memory to store messages on edges (initialized to BEC output (-1,0,1))
 * out : memory in which to store decoder output
 *
 * */
void raw_decode_bec(int* bsi, int* bei, int* bls, int* csi, int* cei,
	int n, int nk, int* edg, int* out, int max_itr)
{
	//printf("Max iterations : %d\n",max_itr);
	for(int it=0; it<max_itr ; it++)
	{
		//printf("Iteration %d\n",it);
		// Bit to check operations
		int nchanges=0; // Track the number of updates made to bit values
		for(int bn=0;bn<n;bn++)
		{
			// Operations for each bit node
			// If bit value is an erasure
			if(out[bn]==0)
			{
				//printf("BN %d\n",bn);
				int t=0;// Variable to store any non erasure value that might be received on a bit node
				for(int ei=bsi[bn]; ei<=bei[bn]; ei++)
				{
					//printf("EID %d, VAL %d\n",bls[ei],edg[bls[ei]]);
					if (edg[bls[ei]]) // If non zero, store that value and break
					{
						t=edg[bls[ei]];
						break;
					}
				}
				for(int ei=bsi[bn]; ei<=bei[bn]; ei++)// Send back whatever was stored earlier
					edg[bls[ei]]=t;
				// Update bit value as well as tracker
				out[bn]=t;
				if (t)
					nchanges++;
			}
			else
			{
				int t=out[bn];
				for(int ei=bsi[bn]; ei<=bei[bn]; ei++)// Send back known value
					edg[bls[ei]]=t;
			}
		}
		// If no erasures have been corrected, stop decoding as it will not proceed further
		if((!nchanges) && it!=0)
		{
			//printf("No changes. Breaking\n");
			break;
		}
		// Check to bit operations
		for(int cn=0;cn<nk;cn++)
		{
			// Operations for individual check nodes
			//printf("CN %d\n",cn);
			int neras=0,eras_loc=-1;
			// Count number of erasures
			for(int e=csi[cn]; e<=cei[cn]; e++)
				if (!edg[e])
				{
					neras++;
					if(neras==1) // If erasure not yet found, update
						eras_loc=e;
					else if(neras>1)
						break;
				}
			//printf("NERAS %d, ERLOC %d\n",neras,eras_loc);
			if(neras>1) // If more than one erasure, send erasures everywhere
				for(int e=csi[cn]; e<=cei[cn]; e++)
					edg[e]=0;
			if(neras==1) // If exactly 1 erasure, send appropriate messages
			{
				int xormsg=1;
				for(int e=csi[cn]; e<=cei[cn]; e++)
					if (edg[e]==-1) // Flip parity
						xormsg=-xormsg;
				// Send back erasures everywhere except the edge that gave the erasure
				for(int e=csi[cn]; e<=cei[cn]; e++)
					edg[e]=0;
				edg[eras_loc]=xormsg;
				//printf("XORMSG %d\n",xormsg);
			}
		}
	}
}
