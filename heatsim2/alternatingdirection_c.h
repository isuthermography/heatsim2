#include <stdint.h>


typedef struct adi_step {
  uint64_t shape[3];
  uint64_t n; /* prod(shape) */

  int permuteorder[3];

  int invpermuteorder[3];

  int stepnum;
  
  uint64_t permutedshape[3];

  double *Amat; /*  tridiagonal... each column represents a diagonal
		    each row represents an equation. First element of
                    first row and last element of last row must be 0
		    dimensions: n x 3
                    Each row represents an equation, indexed
		    according to unwrapped permuted indices. 
		*/
  uint64_t *BmatIndices; /* n*6 x 2 array of 
			    indices for Bmat */
  double *BmatData; /* n*6 data array for Bmat */
  uint64_t BmatEntriesUsed; /* How many of the rows in BmatIndices/BmatData are used? */

  int NumCmats; 
  uint64_t **CmatsIndices;  /* like Bmat, but for C matrices */
  double **CmatsData; 
  uint64_t *CmatsEntriesUsed;
  
  double *Dvec; /* length n vector for Dvec */

  //double *Lmat; /* Tridiagonal L for LU factorization... first subdiagonal and diagonal */
  //double *Umat; /* Tridiagonal U for LU factorization: diagonal assumed to be ones, contains first superdiagonal */
  
  

} adi_step;


static inline uint64_t BmatNRows(struct adi_step *step)
{

  return step->n*12;
}


void delete_adi_step(struct adi_step *step);

struct adi_step *create_adi_step(uint64_t shape[3],uint64_t permutedshape[3],int permuteorder[3],int invpermuteorder[3],int stepnum, int NumCmats);

void add_equation(struct adi_step *step,uint64_t posindex[3],char **eqvarnames,double *eqvalues,int numeqvars);
