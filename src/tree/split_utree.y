/*
    Copyright (C) 2015 Tomas Flouri

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
#include "tree_hashtable.h"

extern int pllmod_utree_lex();
extern FILE * pllmod_utree_in;
extern void pllmod_utree_lex_destroy();
extern int pllmod_utree_lineno;
extern int pllmod_utree_colstart;
extern int pllmod_utree_colend;

static void merge_split(pll_split_t to,
                        const pll_split_t from,
                        unsigned int split_len)
{
  unsigned int i;
  for (i=0;i<split_len;++i)
    to[i] |= from[i];
}

static void pllmod_utree_error(pll_split_system_t * split_stack,
                               void ** splits_ptr,
                               unsigned int * split_count,
                               unsigned int split_size,
                               unsigned int split_len,
                               void * string_hashtable_ptr,
                               const char * s)
{
  pll_errno = PLL_ERROR_NEWICK_SYNTAX;
  if (pllmod_utree_colstart == pllmod_utree_colend)
    snprintf(pll_errmsg, 200, "%s. (line %d column %d)\n",
             s, pllmod_utree_lineno, pllmod_utree_colstart);
  else
    snprintf(pll_errmsg, 200, "%s. (line %d column %d-%d)\n",
             s, pllmod_utree_lineno, pllmod_utree_colstart, pllmod_utree_colend);
}

%}

%union
{
  char * s;
  char * d;
}

%error-verbose
%parse-param {struct split_system_t * split_stack}
             {void ** splits}
             {unsigned int * split_count}
             {unsigned int split_size}
             {unsigned int split_len}
             {void * tipnames_hash}
%destructor { } subtree

%token OPAR
%token CPAR
%token COMMA
%token COLON SEMICOLON
%token<s> STRING
%token<d> NUMBER
%type<s> label optional_label
%type<d> number optional_length

%start inputstr
%%

inputstr: OPAR subtree COMMA subtree COMMA subtree CPAR optional_label optional_length SEMICOLON
{

};

subtree: OPAR subtree COMMA subtree CPAR optional_label optional_length
{
  assert(split_stack->split_count >= 2);
  pll_split_t split1;
  split1 = split_stack->splits[--split_stack->split_count];

  // merge splits and keep at the top
  merge_split(split_stack->splits[split_stack->split_count-1], split1, split_len);

  // append the new split to the final splits list
  memcpy(splits[(*split_count)++], split_stack->splits[split_stack->split_count-1], split_len * sizeof(pll_split_base_t));
}
       | label optional_length
{

  // find tip index
  int tip_id = string_hash_lookup($1, (string_hashtable_t *) tipnames_hash);
  if (tip_id == -1)
    return PLL_FAILURE;

  split_stack->splits[split_stack->split_count][0] = 0;

  unsigned int split_id   = tip_id / split_size;
  tip_id   %= split_size;
  memset(split_stack->splits[split_stack->split_count], 0, split_len * sizeof(pll_split_base_t));
  split_stack->splits[split_stack->split_count][split_id] = (1 << tip_id);
  split_stack->split_count++;

  free($1);
};


optional_label:  { $$ = NULL;} | label  {$$ = $1;};
optional_length: { $$ = NULL;} | COLON number {$$ = $2;};
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
  unsigned int split_len  = (tip_count / split_size) + (tip_count % (sizeof(pll_split_base_t) * 8) > 0);

  pll_split_t splitschunk = (pll_split_t) calloc(max_splits * split_len,
                                                 sizeof(pll_split_base_t));
  pll_split_t stackchunk = (pll_split_t) calloc(max_stack_size * split_len,
                                                  sizeof(pll_split_base_t));

  splits = (pll_split_t *) calloc(tip_count-3, sizeof(pll_split_t));

  split_stack = (pll_split_system_t *) calloc(1, sizeof(pll_split_system_t));
  split_stack->splits = (pll_split_t *) calloc(max_stack_size, sizeof(pll_split_t));
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
  else if (pllmod_utree_parse(split_stack,
                              (void **) splits,
                              &split_count,
                              split_size, split_len,
                              (void *) names_hash))
  {
    free(split_stack->splits);
    free(split_stack);
    free(stackchunk);
    splits = NULL;
    fclose(pllmod_utree_in);
    pllmod_utree_lex_destroy();
    return PLL_FAILURE;
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
