#ifdef NDEBUG  
#undef NDEBUG /* re-enable asserts! */
#endif

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "alternatingdirection_c.h"


void delete_adi_step(struct adi_step *step)
{
  int Cnt;

  if (!step) return;
  
  if (step->Dvec) free(step->Dvec);
  step->Dvec=NULL;

  for (Cnt=0;step->CmatsData && step->CmatsIndices && Cnt < step->NumCmats;Cnt++) {
    if (step->CmatsData[Cnt]) free(step->CmatsData[Cnt]);
    step->CmatsData[Cnt]=NULL;

    if (step->CmatsIndices[Cnt]) free(step->CmatsIndices[Cnt]);
    step->CmatsIndices[Cnt]=NULL;
  }
  



  if (step->CmatsEntriesUsed) free(step->CmatsEntriesUsed);
  step->CmatsEntriesUsed=NULL;
  
  if (step->CmatsData) free(step->CmatsData);
  step->CmatsData=NULL;
  
  if (step->CmatsIndices) free(step->CmatsIndices);
  step->CmatsIndices=NULL;

  if (step->BmatData) free(step->BmatData);
  step->BmatData=NULL;
  
  if (step->BmatIndices) free(step->BmatIndices);
  step->BmatIndices=NULL;
  
  if (step->Amat) free(step->Amat);
  step->Amat=NULL;

  free(step);

}

struct adi_step *create_adi_step(uint64_t shape[3],uint64_t permutedshape[3],int permuteorder[3],int invpermuteorder[3],int stepnum, int NumCmats)
{
  struct adi_step *step;
  int Cnt;
  
  step=calloc(sizeof(*step),1);

  assert(sizeof(*shape)==sizeof(*step->shape));
  memcpy(step->shape,shape,sizeof(step->shape));
  step->n=step->shape[0]*step->shape[1]*step->shape[2];
  //fprintf(stderr,"n=%d\n",(int)step->n);

  assert(sizeof(*permutedshape)==sizeof(*step->permutedshape));  
  memcpy(step->permutedshape,permutedshape,sizeof(step->permutedshape));

  assert(sizeof(*permuteorder)==sizeof(*step->permuteorder));  
  memcpy(step->permuteorder,permuteorder,sizeof(step->permuteorder));

  assert(sizeof(*invpermuteorder)==sizeof(*step->invpermuteorder));  
  memcpy(step->invpermuteorder,invpermuteorder,sizeof(step->invpermuteorder));


  step->stepnum=stepnum;
  step->NumCmats=NumCmats;

  /* allocate arrays */
  step->Amat=calloc(sizeof(*step->Amat),step->n*3);
  step->BmatIndices=calloc(sizeof(*step->BmatIndices),BmatNRows(step)*2);
  step->BmatData=calloc(sizeof(*step->BmatData),BmatNRows(step));
  step->CmatsIndices=calloc(sizeof(*step->CmatsIndices),step->NumCmats);
  step->CmatsData=calloc(sizeof(*step->CmatsData),step->NumCmats);
  step->CmatsEntriesUsed=calloc(sizeof(*step->CmatsEntriesUsed),step->NumCmats);
  for (Cnt=0;Cnt < step->NumCmats;Cnt++) {
    step->CmatsData[Cnt]=calloc(sizeof(*step->CmatsData[Cnt]),BmatNRows(step));
    step->CmatsIndices[Cnt]=calloc(sizeof(*step->CmatsIndices[Cnt]),BmatNRows(step)*2);
    step->CmatsEntriesUsed[Cnt]=0;
  }
  
  step->Dvec=calloc(sizeof(*step->Dvec),step->n);
  //step->Lmat=calloc(sizeof(*step->Lmat),step->n*2);
  //step->Umat=calloc(sizeof(*step->Umat),step->n);

  return step;
}


void add_equation(struct adi_step *step,uint64_t posindex[3],char **eqvarnames,double *eqvalues,int numeqvars)
{

  uint64_t rownum;
  int eqvarcnt;
  int64_t shift[3];
  long solnum;
  uint64_t Bmatentry;
  uint64_t Cmatentry;
  char tshift;
  int64_t shiftedcolnum;
  
  rownum=step->permutedshape[2]*step->permutedshape[1]*posindex[step->permuteorder[0]] + step->permutedshape[2]*posindex[step->permuteorder[1]] + posindex[step->permuteorder[2]];

  assert(rownum < step->n);
  
  for (eqvarcnt=0;eqvarcnt < numeqvars;eqvarcnt++) {
    if (!strcmp(eqvarnames[eqvarcnt],"volumetric_source")) {
      /* constant term... drop in to Dvec */
      step->Dvec[rownum]=eqvalues[eqvarcnt];
    } else if (eqvarnames[eqvarcnt][0]==0) {
      /* Blank... ignore constant term if it has zero value */
      assert(eqvalues[eqvarcnt]==0.0);
    } else {
      assert(eqvarnames[eqvarcnt][0]=='T'); /* all other vars should be temperature */
      
      /* Evaluate shift from variable name */
      shift[0]=((int)eqvarnames[eqvarcnt][1])-'5';
      shift[1]=((int)eqvarnames[eqvarcnt][2])-'5';
      shift[2]=((int)eqvarnames[eqvarcnt][3])-'5';

      
      /* check tshift... positive? ... if so may refer to explicit solnum */
      tshift=eqvarnames[eqvarcnt][4];
      if (tshift=='p') {
        solnum=strtol(eqvarnames[eqvarcnt]+5,NULL,10);
      } else {
        solnum=0;
      }

      /* Check if this should go on the left side of the equals sign */
      if (tshift=='p' && solnum==step->stepnum) {
        /* it is to be solved in this step */
	//fprintf(stderr,"%s\n",eqvarnames[eqvarcnt]);
        assert(shift[step->permuteorder[0]]==0 && shift[step->permuteorder[1]]==0);
	step->Amat[rownum*3+(shift[step->permuteorder[2]]+1)]=-eqvalues[eqvarcnt];
	//fprintf(stderr,"A row %d/%d\n",(int)rownum,(int)step->n);
	
      } else {


        /* add shift to posindex to identify absolute location */
        shift[0]+=posindex[0];
        shift[1]+=posindex[1];
        shift[2]+=posindex[2];
      
	shiftedcolnum=step->permutedshape[2]*step->permutedshape[1]*shift[step->permuteorder[0]] + step->permutedshape[2]*shift[step->permuteorder[1]] + shift[step->permuteorder[2]];

	if (shiftedcolnum < 0 || shiftedcolnum >= step->n) {
	  fprintf(stderr,"heatsim2/alternatingdirection_c/add_equation: Equation exceeds bounds of domain. Are external boundaries set correctly?\n");
	  exit(1);
	}

	if (tshift=='m' || tshift==0) {
          /* Goes into 'B' matrix */
          Bmatentry=step->BmatEntriesUsed;
	  assert(Bmatentry < BmatNRows(step));

	  //fprintf(stderr,"Bmatentry=%d;BmatNRows(step)=%d\n",(int)Bmatentry,(int)BmatNRows(step));
	  //fprintf(stderr,"rownum=%lld;shiftedcolunum=%lld\n",(long long)rownum,(long long)shiftedcolnum);

	  
	  step->BmatIndices[Bmatentry*2]=rownum;
	  step->BmatIndices[Bmatentry*2+1]=shiftedcolnum;
	  
	  step->BmatData[Bmatentry]=eqvalues[eqvarcnt];
	  step->BmatEntriesUsed++;
	  
        } else {
          /* Put in the appropriate 'C' matrix */
          assert(solnum < step->NumCmats);
          Cmatentry=step->CmatsEntriesUsed[solnum];
	  assert(Cmatentry < BmatNRows(step));

	  step->CmatsIndices[solnum][Cmatentry*2]=rownum;
	  step->CmatsIndices[solnum][Cmatentry*2+1]=shiftedcolnum;
	  
	  step->CmatsData[solnum][Cmatentry]=eqvalues[eqvarcnt];
	  step->CmatsEntriesUsed[solnum]++;
          
        }
	
      }
    }
  }
  
  
}
