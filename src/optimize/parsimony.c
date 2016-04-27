#include "pll_tree.h"

static void newviewParsimonyIterativeFast (pll_partition_t * partition,
                                           const pll_operation_t *ops,
                                           unsigned int count,
                                           double *pars_vec,
                                           double *pars_scores)
{
  unsigned int index, i, k;

  for (index = 0; index < count; ++index)
  {
    unsigned int totalScore = 0;

    unsigned int
        pNumber = ops[index]->child1_clv_index,
        qNumber = ops[index]->child2_clv_index,
        rNumber = ops[index]->parent_clv_index;

    unsigned int states = partition->states, width = partition->sites; // parsimony length??

    unsigned int **left = (unsigned int **) alloca(
        states * sizeof(unsigned int *)), **right = (unsigned int **) alloca(
        states * sizeof(unsigned int *)), **this = (unsigned int **) alloca(
        states * sizeof(unsigned int *)), *o_A = (unsigned int *) alloca(
        states * sizeof(unsigned int)), *t_A = (unsigned int *) alloca(
        states * sizeof(unsigned int)), t_N;

    for (i = 0; i < states; i++)
    {
      left[i]  = &(pars_vec[(width * states * qNumber) + width * i]);
      right[i] = &(pars_vec[(width * states * rNumber) + width * i]);
      this[i]  = &(pars_vec[(width * states * pNumber) + width * i]);
    }

    for (i = 0; i < width; i++)
    {
      t_N = 0;

      for (k = 0; k < states; k++)
      {
        t_A[k] = left[k][i] & right[k][i];
        o_A[k] = left[k][i] | right[k][i];
        t_N = t_N | t_A[k];
      }

      t_N = ~t_N;

      for (k = 0; k < states; k++)
        this[k][i] = t_A[k] | (t_N & o_A[k]);

      totalScore += ((unsigned int) __builtin_popcount (t_N));
    }

    pars_scores[pNumber] = totalScore + pars_scores[rNumber]
        + pars_scores[qNumber];
  }
}

static double randum (long  *seed)
{
  long  sum, mult0, mult1, seed0, seed1, seed2, newseed0, newseed1, newseed2;
  double res;

  mult0 = 1549;
  seed0 = *seed & 4095;
  sum  = mult0 * seed0;
  newseed0 = sum & 4095;
  sum >>= 12;
  seed1 = (*seed >> 12) & 4095;
  mult1 =  406;
  sum += mult0 * seed1 + mult1 * seed0;
  newseed1 = sum & 4095;
  sum >>= 12;
  seed2 = (*seed >> 24) & 255;
  sum += mult0 * seed2 + mult1 * seed1;
  newseed2 = sum & 255;

  *seed = newseed2 << 24 | newseed1 << 12 | newseed0;
  res = 0.00390625 * (newseed2 + 0.000244140625 * (newseed1 + 0.000244140625 * newseed0));

  return res;
}

static void make_permutation_fast(int *perm, int n, long rng_seed)
{
  int
    i,
    j,
    k;

  for (i = 1; i <= n; i++)
    perm[i] = i;

  for (i = 1; i <= n; i++)
    {
      double d =  randum(&rng_seed);

      k =  (int)((double)(n + 1 - i) * d);

      j        = perm[i];

      perm[i]     = perm[i + k];
      perm[i + k] = j;
    }
}

static pll_utree_t * build_simple_tree (pll_partition_t *partition, int *ntips)
{
  pll_utree_t  *tree, *newinner, *newtip;
  int  i;

  *ntips = 3;
  tree = (pll_utree_t *) malloc (sizeof(pll_utree_t));
  newinner = pll_utree_create_node(0, 0, 0, 0);

  tree->next = NULL;
  pll_utree_conect_nodes(tree, newinner, 0);
  newtip = (pll_utree_t *) malloc (sizeof(pll_utree_t));
  pll_utree_conect_nodes(newinner->next, newtip, 0);
  newtip = (pll_utree_t *) malloc (sizeof(pll_utree_t));
  pll_utree_conect_nodes(newinner->next->next, newtip, 0);

  insertParsimony(tr, pr, s, p);
  return tree;
}

void pll_utree_make_parsimony_fast(pll_utree_t *tree, pll_partition_t *partition, long rng_seed)
{
  pll_utree_t
    *p,
    *f;

  int
    i,
    nextsp,
    *perm        = (int *)malloc((size_t)(partition->tips + 1) * sizeof(int));

  unsigned int
    randomMP,
    startMP;

  make_permutation_fast(perm, partition->tips, rng_seed);

  tr->ntips = 0;

  tr->nextnode = tr->mxtips + 1;

  buildSimpleTree(tr, pr, perm[1], perm[2], perm[3]);

  f = tr->start;

  while(tr->ntips < tr->mxtips)
    {
      nodeptr q;

      tr->bestParsimony = INT_MAX;
      nextsp = ++(tr->ntips);
      p = tr->nodep[perm[nextsp]];
      q = tr->nodep[(tr->nextnode)++];
      p->back = q;
      q->back = p;

      if(tr->grouped)
        {
          int
            number = p->back->number;

          tr->constraintVector[number] = -9;
        }

      stepwiseAddition(tr, pr, q, f->back);

      {
        nodeptr
          r = tr->insertNode->back;

        int counter = 4;

        hookupDefault(q->next,       tr->insertNode);
        hookupDefault(q->next->next, r);

        computeTraversalInfoParsimony(q, tr->ti, &counter, tr->mxtips, PLL_FALSE);
        tr->ti[0] = counter;

        newviewParsimonyIterativeFast(tr, pr);
      }
    }

  nodeRectifierPars(tr);

  randomMP = tr->bestParsimony;

  do
    {
      startMP = randomMP;
      nodeRectifierPars(tr);
      for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
        {
          rearrangeParsimony(tr, pr, tr->nodep[i], 1, 20, PLL_FALSE);
          if(tr->bestParsimony < randomMP)
            {
              restoreTreeRearrangeParsimony(tr, pr);
              randomMP = tr->bestParsimony;
            }
        }
    }
  while(randomMP < startMP);

  rax_free(perm);
}
