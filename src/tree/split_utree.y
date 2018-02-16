/*
    Copyright (C) 2017 Tomas Flouri, Diego Darriba

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

%{
#include "../pllmod_common.h"
#include "tree_hashtable.h"

extern int pllmod_utree_lex();
extern FILE * pllmod_utree_in;
extern void pllmod_utree_lex_destroy();
extern int pllmod_utree_lineno;
extern int pllmod_utree_colstart;
extern int pllmod_utree_colend;

struct parse_params_t
{
  struct split_system_t * split_stack;
  void ** splits;
  unsigned int * split_count;
  unsigned int split_size;
  unsigned int split_len;
  void * tipnames_hash;
};

static void merge_split(pll_split_t to,
                        const pll_split_t from,
                        unsigned int split_len)
{
  unsigned int i;
  for (i=0;i<split_len;++i)
    to[i] |= from[i];
}

static void pllmod_utree_error(struct parse_params_t * paramas,
                               const char * s)
{
  pll_errno = PLL_ERROR_NEWICK_SYNTAX;
  if (pllmod_utree_colstart == pllmod_utree_colend)
    snprintf(pll_errmsg, 200, "%s. (line %d column %d)\n",
             s, pllmod_utree_lineno, pllmod_utree_colstart);
  else
    snprintf(pll_errmsg, 200, "%s. (line %d column %d-%d)\n", s,
             pllmod_utree_lineno, pllmod_utree_colstart, pllmod_utree_colend);
}

%}

%union
{
  char * s;
  char * d;
  void * nulval;
}

%error-verbose

%parse-param {void * params_ptr}
%destructor { } subtree


%token<s> STRING
%token<d> NUMBER
%type<s> label optional_label
%type<d> number optional_length
%type<nulval> inputstr
%type<nulval> subtree

%start inputstr
%%

inputstr: '(' subtree ',' subtree ',' subtree ')' optional_label optional_length ';'
{
  PLLMOD_UNUSED($$);PLLMOD_UNUSED($2);
  PLLMOD_UNUSED($4);PLLMOD_UNUSED($6); /* ignore result */
};

subtree: '(' subtree ',' subtree ')' optional_label optional_length
{
  PLLMOD_UNUSED($$);PLLMOD_UNUSED($2);PLLMOD_UNUSED($4);
  struct parse_params_t * params = (struct parse_params_t *) params_ptr;

  assert(params->split_stack->split_count >= 2);
  pll_split_t split1;
  split1 = params->split_stack->splits[--params->split_stack->split_count];

  // merge splits and keep at the top
  merge_split(params->split_stack->splits[params->split_stack->split_count-1],
              split1,
              params->split_len);

  // append the new split to the final splits list
  memcpy(params->splits[(*params->split_count)++],
         params->split_stack->splits[params->split_stack->split_count-1],
         params->split_len * sizeof(pll_split_base_t));
}
       | label optional_length
{
  PLLMOD_UNUSED($$);
  struct parse_params_t * params = (struct parse_params_t *) params_ptr;

  // find tip index
  int tip_id = string_hash_lookup($1,
                                  (string_hashtable_t *) params->tipnames_hash);
  if (tip_id == -1)
    return PLL_FAILURE;

  params->split_stack->splits[params->split_stack->split_count][0] = 0;

  unsigned int split_id   = tip_id / params->split_size;
  tip_id   %= params->split_size;
  memset(params->split_stack->splits[params->split_stack->split_count], 0,
         params->split_len * sizeof(pll_split_base_t));
  params->split_stack->splits[params->split_stack->split_count][split_id] =
         (1 << tip_id);
  params->split_stack->split_count++;

  free($1);
};


optional_label:  { $$ = NULL;} | label  {$$ = $1;};
optional_length: { $$ = NULL;} | ':' number {$$ = $2;};
label: STRING    { $$=$1;} | NUMBER {$$=$1;};
number: NUMBER   { $$=$1;};

%%

#ifdef __linux__
PLL_EXPORT pll_split_t * pll_utree_split_newick_string(char * s,
                                               unsigned int tip_count,
                                               string_hashtable_t * names_hash)
{
  unsigned int i;
  pll_split_system_t * split_stack;
  pll_split_t * splits;
  unsigned int max_stack_size = 2*tip_count - 3;
  unsigned int max_splits = tip_count - 3;

  unsigned int split_count = 0;
  unsigned int split_size = (sizeof(pll_split_base_t) * 8);
  unsigned int split_len  = (tip_count / split_size) +
                            (tip_count % (sizeof(pll_split_base_t) * 8) > 0);

  pll_split_t splitschunk = (pll_split_t) calloc(max_splits * split_len,
                                                 sizeof(pll_split_base_t));
  pll_split_t stackchunk = (pll_split_t) calloc(max_stack_size * split_len,
                                                  sizeof(pll_split_base_t));

  splits = (pll_split_t *) calloc(tip_count-3, sizeof(pll_split_t));

  split_stack = (pll_split_system_t *) calloc(1, sizeof(pll_split_system_t));
  split_stack->splits = (pll_split_t *) calloc(max_stack_size,
                                               sizeof(pll_split_t));
  split_stack->support = 0;
  split_stack->split_count = 0;

  for (i=0; i<max_splits; ++i)
  {
    splits[i] = splitschunk + i*split_len;
  }
  for (i=0; i<max_stack_size; ++i)
    split_stack->splits[i] = stackchunk + i*split_len;

  pllmod_utree_in = fmemopen(s, strlen(s), "r");

  if (!pllmod_utree_in)
  {
    free(split_stack->splits);
    free(split_stack);
    free(stackchunk);
    pll_errno = PLL_ERROR_FILE_OPEN;
    snprintf(pll_errmsg, 200, "Unable to map string (%s)", s);
    return PLL_FAILURE;
  }
  else
  {
    struct parse_params_t parse_params;
    parse_params.split_stack = split_stack;
    parse_params.splits = (void **) splits;
    parse_params.split_count = &split_count;
    parse_params.split_size = split_size;
    parse_params.split_len = split_len;
    parse_params.tipnames_hash = (void *) names_hash;

    if (pllmod_utree_parse(&parse_params))
    {
      free(split_stack->splits);
      free(split_stack);
      free(stackchunk);
      splits = NULL;
      fclose(pllmod_utree_in);
      pllmod_utree_lex_destroy();
      return PLL_FAILURE;
    }
  }

  if (pllmod_utree_in) fclose(pllmod_utree_in);

  free(split_stack->splits);
  free(split_stack);
  free(stackchunk);
  pllmod_utree_lex_destroy();

  assert(split_count == max_splits);

  return splits;
}
#endif
