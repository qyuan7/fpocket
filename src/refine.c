#include "../headers/refine.h"
/*
 * Copyright <2012> <Vincent Le Guilloux,Peter Schmidtke, Pierre Tuffery>
 * Copyright <2013-2018> <Peter Schmidtke, Vincent Le Guilloux>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 */
/*

## GENERAL INFORMATION
##
## FILE 					refine.c
## AUTHORS					P. Schmidtke and V. Le Guilloux
## LAST MODIFIED			30-01-12
##
## SPECIFICATIONS
##
## This file defins several routines which refine the clustering algorithm
## used by fpocket. In particular, we merge pockets too closed from each other
## we drop small pockets, and we perform reindexation after those operations.
##
## MODIFICATIONS HISTORY
##
##      30-01-12        (p)  Adding apply clustering routine for new clustering
##	09-02-09	(v)  Drop tiny pocket routine added
##	28-11-08	(v)  Comments UTD 
##	01-04-08	(v)  Added template for comments and creation of history
##	01-01-08	(vp) Created (random date...)
##	
## TODO or SUGGESTIONS
##
## (v) Improve and optimize the algorithm.
##

*/


/**
   ## FUNCTION:
	apply_clustering


   ## SPECIFICATION:
        Using prior clustering results written to vertice resid's, update and
 *      merge pockets here

   ## PARAMETRES:
	@ c_lst_pockets *pockets : The list of pockets.
	@ s_fparams *params      : Parameters

   ## RETURN:
	void

*/
c_lst_pockets *apply_clustering(c_lst_pockets *pockets, s_fparams *params, s_lst_vvertice *lvert)
{
	int npockets = 0;
	node_pocket *nextPocket;
	node_pocket *curMobilePocket;

	node_pocket *pcur = NULL ;
        s_vvertice *vcur=NULL;
        int cur_resid=-1;
        s_vvertice *vmobile=NULL;

        int cur_mobile_resid=-2;

	if(pockets){
		pcur = pockets->first;
		while(pcur){
			npockets += 1;
			pcur = pcur->next;
		}
	}
	else {
		fprintf(stderr, "! No pocket to refine! (argument NULL: %p).\n", pockets) ;
	}
	int n_resid = 0;
	int max_resid = -1;

	if(pockets){
		pcur = pockets->first;
		while(pcur){
			vcur = pcur->pocket->v_lst->first->vertice;
			cur_resid = vcur->resid;
			if(cur_resid > max_resid){
				max_resid = cur_resid;
			}
			pcur = pcur->next;
		}
	}
	int mem_table[max_resid+1];

	for(int i=0; i<=max_resid; i++){
		mem_table[i]=-1;
	}

	if(pockets){
		pcur = pockets->first;
		while(pcur){
			vcur = pcur->pocket->v_lst->first->vertice;
			cur_resid = vcur->resid;
			if(mem_table[cur_resid] == -1){
				mem_table[cur_resid] = n_resid;
				n_resid += 1;
			}
			//if(HashmapGet(map, cur_resid) == NULL){
			//	HashmapPut(map, cur_resid, (void*)(ptr+n_resid));
			//	int *tmp = (int*)HashmapGet(map, cur_resid);
			//	mem_table[cur_resid] = n_resid;
			//	//fprintf(stdout, "%d,%d\n",cur_resid, *tmp);
			//	n_resid += 1;
			//}
			pcur = pcur->next;
		}
	}
    c_lst_pockets **mergedPockets = (c_lst_pockets *)
                                               my_malloc(n_resid*sizeof(c_lst_pockets));

	//c_lst_pockets *mergedPockets[n_resid];
	for(int i=0; i<n_resid; i++){
		mergedPockets[i] = c_lst_pockets_alloc();
	}
    int mergeId;

    int i=0;
	if(pockets){
		pcur = pockets->first;
		while(pcur){

			vcur = pcur->pocket->v_lst->first->vertice;
			cur_resid = vcur->resid;
			//fprintf(stdout, "%d ", pcur->pocket->v_lst->n_vertices);
			mergeId = mem_table[cur_resid];
			//fprintf(stdout, "%d \n", mergeId);
			int cur_n_apol = pcur->pocket->nAlphaApol;
			int cur_n_pol = pcur->pocket->nAlphaPol;
			s_pocket *pock = alloc_pocket();
			pock->v_lst=c_lst_vertices_alloc();
			//pock->v_lst = pcur->pocket->v_lst;
			//c_lst_vertices_add_last(pock->v_lst, pcur->pocket->v_lst->first->vertice);
			c_lst_vertices_add_last(pock->v_lst, lvert->vertices+i);
			pock->size = pock->v_lst->n_vertices;
			c_lst_pockets_add_last(mergedPockets[mergeId], pock,cur_n_apol,cur_n_pol);
			if(mergedPockets[mergeId]->n_pockets>1){
				mergePockets(mergedPockets[mergeId]->first, mergedPockets[mergeId]->first->next,mergedPockets[mergeId]);

			}
            i += 1;
			pcur = pcur->next;
		}
	}

	if(pockets){
		pcur = pockets->first;
		while(pcur){
			my_free(pcur->pocket);
			pcur = pcur->next;
		}
	}

	c_lst_pockets *pockets1 = c_lst_pockets_alloc();
	for(int i=0; i<n_resid; i++){
	
   		int cur_n_apol = mergedPockets[i]->first->pocket->nAlphaApol;
	   	int cur_n_pol = mergedPockets[i]->first->pocket->nAlphaPol;

	 //c_lst_pockets_add_last(pockets1, pock,cur_n_apol, cur_n_pol);
	    c_lst_pockets_add_last(pockets1, mergedPockets[i]->first->pocket, cur_n_apol, cur_n_pol);
	}
	pockets1->vertices=lvert;
	//pockets1->vertices->nvert;
	//fprintf(stdout, "% verts in total\n", pockets1->vertices->nvert);
	//print_number_of_objects_in_memory();
	//for (int i=0; i<n_resid; i++){
	//	c_lst_pocket_free(mergedPockets[i]);
	//}
	//my_free(mergedPockets);
    //print_number_of_objects_in_memory();
	node_pocket *p = pockets1->first ;
	while(p) {
		p->pocket->size = p->pocket->v_lst->n_vertices ;
		p = p->next ;
	}

	if(pockets1->n_pockets > 0){
		return pockets1;
	}
	else {
		my_free(pockets1) ;
		return NULL ;
	}

}


void apply_clustering_old(c_lst_pockets *pockets, s_fparams *params)
{
	int npockets = 0;
	node_pocket *nextPocket;
	node_pocket *curMobilePocket;

	node_pocket *pcur = NULL ;
        s_vvertice *vcur=NULL;
        int cur_resid=-1;
        s_vvertice *vmobile=NULL;

        int cur_mobile_resid=-2;

	if(pockets){
		pcur = pockets->first;
		while(pcur){
			npockets += 1;
			pcur = pcur->next;
		}
	}
	else {
		fprintf(stderr, "! No pocket to refine! (argument NULL: %p).\n", pockets) ;
	}
	int n_resid = 0;
	int max_resid = -1;
	sort_pockets(pockets, cmp_pockets);
	if(pockets){
		pcur = pockets->first;
	
		while(pcur){
			vcur = pcur->pocket->v_lst->first->vertice;
			cur_resid = vcur->resid;
			curMobilePocket=pcur->next;
			if(curMobilePocket){
			    if(curMobilePocket->pocket->v_lst->first->vertice->resid == cur_resid){
				    mergePockets(pcur, curMobilePocket, pockets);
			    }
			    else{
				    pcur = curMobilePocket;
				    //cur_resid = curMobilePocket->pocket->v_lst->first->vertice->resid;
			    }
			}
			else{
				break;
			}

		}
	}
	/*if(pockets) {
            pcur = pockets->first ;
            while(pcur) {
                vcur=pcur->pocket->v_lst->first->vertice;
                curMobilePocket = pcur->next ;
                cur_resid=vcur->resid;
                while(curMobilePocket) {
                    vmobile=curMobilePocket->pocket->v_lst->first->vertice;
                    cur_mobile_resid=vmobile->resid;
                    nextPocket = curMobilePocket->next;
                    
                    if(cur_resid==cur_mobile_resid) {
                    // Merge pockets if barycentres are close to each other
                            mergePockets(pcur, curMobilePocket, pockets);
                    }
                    curMobilePocket = nextPocket ;
                }

                pcur = pcur->next ;
            }
	}*/
	
	else {
		fprintf(stderr, "! No pocket to refine! (argument NULL: %p).\n", pockets) ;
	}

}


c_lst_pockets *assign_apply_clustering( s_fparams *params, s_lst_vvertice *lvert)
{
	int i = -1,
			
    cur_n_pol=0,
    cur_n_apol=0;

	s_vvertice *vertices = lvert->vertices,
			   *vcur = NULL ;
	c_lst_pockets *pockets = c_lst_pockets_alloc();
	int npockets = lvert->nvert;
	int n_resid = 0;
	int max_resid = -1;
	int cur_resid = -1;
	for(i=0; i<lvert->nvert;i++){
		vcur = vertices + i;
		cur_resid = vcur->resid;
		if(cur_resid > max_resid){
			max_resid = cur_resid;
		}
	}
	int mem_table[max_resid+1];
	for(int i=0; i<=max_resid; i++){
		mem_table[i]=-1;
	}
	for(i=0; i<lvert->nvert; i++){
		vcur = vertices + i;
		cur_resid = vcur->resid;
		//fprintf(stdout, "%d \n", cur_resid);

		if(mem_table[cur_resid] == -1){
			mem_table[cur_resid] = n_resid;
			n_resid += 1;
		}
	}
	
	c_lst_pockets **mergedPockets = (c_lst_pockets *)
                                               my_malloc(n_resid*sizeof(c_lst_pockets));

	//c_lst_pockets *mergedPockets[n_resid];
	for(int i=0; i<n_resid; i++){
		mergedPockets[i] = c_lst_pockets_alloc();
	}
    int mergeId;
    //fprintf(stdout, "max id is %d\n", max_resid);

	//fprintf(stdout, "%d allocated.\n", n_resid);

	for(i=0;i<npockets;i++) {

       	vcur = vertices + i ;
		cur_resid = vcur->resid;
		mergeId = mem_table[cur_resid];
		//fprintf(stdout, "group %d\n", mergeId);
        cur_n_apol=0;
        cur_n_pol=0;
        //vcur->resid=i+1;
        //vcur->id=i+1;
        /* Create a new pocket */

        s_pocket *pocket = alloc_pocket();
        pocket->v_lst=c_lst_vertices_alloc();
        /* Add vertices to the pocket */

        c_lst_vertices_add_last(pocket->v_lst, vcur);

        if(vcur->type==M_APOLAR_AS) cur_n_apol++;
        else cur_n_pol++;
        //pocket->rank=i+1;

        pocket->size = pocket->v_lst->n_vertices;
        //fprintf(stdout, "%d\n", mergedPockets[mergeId]->n_pockets);
        //fprintf(stdout, "before.\n");
        //fprintf(stdout, "%d pockets\n", mergedPockets[mergeId]->n_pockets);
        c_lst_pockets_add_last(mergedPockets[mergeId], pocket,cur_n_apol,cur_n_pol);
		//fprintf(stdout, "added.\n");
		if(mergedPockets[mergeId]->n_pockets>1){
			mergePockets(mergedPockets[mergeId]->first, mergedPockets[mergeId]->first->next,mergedPockets[mergeId]);
		}

	}

	for(int i=0; i<n_resid; i++){
	
   		int cur_n_apol = mergedPockets[i]->first->pocket->nAlphaApol;
	   	int cur_n_pol = mergedPockets[i]->first->pocket->nAlphaPol;

	    c_lst_pockets_add_last(pockets, mergedPockets[i]->first->pocket, cur_n_apol, cur_n_pol);
	}
	pockets->vertices=lvert;
	node_pocket *p = pockets->first ;
	while(p) {
		p->pocket->size = p->pocket->v_lst->n_vertices ;
		p = p->next ;
	}
    for(int i=0; i<n_resid; i++){
        my_free(mergedPockets[i]);
	}
	my_free(mergedPockets);
	if(pockets->n_pockets > 0){
		return pockets;
	}
	else {
		fprintf(stderr, "! No pocket to refine! (argument NULL: %p).\n", pockets) ;
		my_free(pockets) ;
		return NULL ;
	}

}

/**
   ## FUNCTION: 
	dropSmallNpolarPockets
	--
   ## SPECIFICATION:
	Refine algorithm: will remove small pockets (depends on the corresponding
	parameters in params), pockets containing less than NB apolar alpha spheres
	(given in params)..
  
   ## PARAMETRES:
	@ c_lst_pockets *pockets : The list of pockets.
	@ s_fparams *params      : Parameters
  
   ## RETURN:
	void
  
*/
void dropSmallNpolarPockets(c_lst_pockets *pockets, s_fparams *params)
{
	double pasph = 0.0 ;
	node_pocket *npcur = NULL,
				*nextPocket1 = NULL;
	s_pocket *pcur = NULL ;
	//fprintf(stdout, "%d pockets processed.\n", pockets->n_pockets);
	if(pockets) {
		npcur = pockets->first ;
		while(npcur) {
			pcur = npcur->pocket ;
			nextPocket1 = npcur->next ;
			pasph = (float)((float)pcur->nAlphaApol/(float)pcur->v_lst->n_vertices) ;


			if(pcur->v_lst->n_vertices < (size_t) params->min_pock_nb_asph 
				||  pasph <  (params->refine_min_apolar_asphere_prop) || pcur->pdesc->as_density<(params->min_as_density)){
			/* If the pocket is too small or has not enough apolar alpha
			 * spheres, drop it: QY:correct code to work as expected*/
			    //fprintf(stdout, "size is %d\n", pcur->v_lst->n_vertices);
				dropPocket(pockets, npcur);		
			}

			if(pockets->n_pockets <= 0) fprintf(stderr, "! No Pockets Found while refining\n");
			npcur = nextPocket1 ;
		}
	}
	else {
		fprintf(stderr, "! No pockets to drop from (argument NULL: %p).\n", pockets) ;
	}

}

void dropSmallPockets(c_lst_pockets *pockets, s_fparams *params){
	node_pocket *npcur = NULL,
				*nextPocket1 = NULL;
	s_pocket *pcur = NULL ;
	
	if(pockets) {
		npcur = pockets->first ;
		while(npcur) {
			pcur = npcur->pocket ;
			nextPocket1 = npcur->next ;

			if(pcur->v_lst->n_vertices < (size_t) params->min_pock_nb_asph){
			/* If the pocket is too small or has not enough apolar alpha
			 * spheres, drop it */
				dropPocket(pockets, npcur);		
			}

			if(pockets->n_pockets <= 0) fprintf(stderr, "! No Pockets Found while refining\n");
			npcur = nextPocket1 ;
		}
	}
	else {
		fprintf(stderr, "! No pockets to drop from (argument NULL: %p).\n", pockets) ;
	}
}

/**
   ## FUNCTION:
	drop_tiny
	--
   ## SPECIFICATION:
	Drop really tiny pockets (< 5 alpha spheres)
  
   ## PARAMETRES:
	@ c_lst_pockets *pockets : The list of pockets.
  
   ## RETURN:
	void
  
*/
void drop_tiny(c_lst_pockets *pockets, s_fparams *params)
{
	node_pocket *npcur = pockets->first,
				*npnext = NULL ;
	while(npcur) {
		npnext = npcur->next ;

		if(npcur->pocket->v_lst->n_vertices < params->min_pock_nb_asph){
		/* If the pocket is really small, drop it */
			dropPocket(pockets, npcur);
		}

		if(pockets->n_pockets <= 0) fprintf(stderr, "! No Pockets Found while refining\n");
		npcur = npnext ;
	}
}

/**
   ## FUNCTION: 
	reIndexPockets
  
   ## SPECIFICATION:
	Reindex pockets, after dropping and merging operations on pockets and
	recalculate barycentres in the same loop
  
   ## PARAMETRES:
	@ c_lst_pockets *pockets: The list of pockets.
  
   ## RETURN:
	void
  
*/
void reIndexPockets(c_lst_pockets *pockets)
{
	node_vertice *vcur = NULL ;
	node_pocket *pcur = NULL ;
	s_pocket *pock_cur = NULL ;

	int curPocket = 0,
		n_vert ;

	float posSum[3];

	if(pockets && pockets->n_pockets > 0) {
		pcur = pockets->first ;
		while(pcur) {
			curPocket++ ;						//new index counter
			n_vert = 0 ;

			pock_cur = pcur->pocket ;
			pock_cur->bary[0]=0 ;
			pock_cur->bary[1]=0 ;
			pock_cur->bary[2]=0 ;
				
			posSum[0]=0; posSum[1]=0; posSum[2]=0;
			if(pock_cur->v_lst){
				vcur = pock_cur->v_lst->first;
				while(vcur){
					posSum[0] += vcur->vertice->x;
					posSum[1] += vcur->vertice->y;
					posSum[2] += vcur->vertice->z;
					n_vert++;
	
					vcur->vertice->resid = curPocket;	//set new index
					vcur = vcur->next ;
				}
				//fprintf(stdout, "%d ", n_vert);
			}
			else {
				fprintf(stderr, "! Empty Pocket...something might be wrong over here ;).\n") ; 
			}
			//set new barycentre
			pock_cur->bary[0] = posSum[0]/(float)n_vert;
			pock_cur->bary[1] = posSum[1]/(float)n_vert;
			pock_cur->bary[2] = posSum[2]/(float)n_vert;
			pcur = pcur->next ;
			//fprintf(stdout, "%d\n", curPocket);
		}
	}
	else {
		fprintf(stderr, "! No pocket to reindex.\n") ;
	}
}

int cmp_pockets(const node_pocket *p1, const node_pocket *p2){
    return ((p1->pocket->v_lst->first->vertice->resid < p2->pocket->v_lst->first->vertice->resid )? -1 : 1);
}
