/*
 Copyright (C) 2018 Alexey Kozlov

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

 Contact: Alexey Kozlov <Alexey.Kozlov@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */

#include "pll_tree.h"
#include "../common.h"

#include <assert.h>

#define TREEFILE  "testdata/medium.tree"

void run_ham_test(unsigned int test_num, pll_split_t s1, pll_split_t s2,
              unsigned int num_tips)
{
  pllmod_utree_split_show(s1, num_tips);
  printf("\n");
  pllmod_utree_split_show(s2, num_tips);
  printf("\n");

  unsigned int p1 =  pllmod_utree_split_lightside(s1, num_tips);
  unsigned int p2 =  pllmod_utree_split_lightside(s2, num_tips);

  unsigned int hdist = pllmod_utree_split_hamming_distance(s1, s2, num_tips);

  printf("TEST #%u: LIGHT SIDE: %3u %3u, HAMMING_DIST: %u\n", test_num, p1, p2, hdist);
}

void test_hamming()
{
  pll_split_base_t sb1[2] = {0x1, 0xFFFFFFFF};
  pll_split_base_t sb2[2] = {0x2, 0xFFFFFFF0};

  pll_split_t s1 = (pll_split_t) &sb1;
  pll_split_t s2 = (pll_split_t) &sb2;

  run_ham_test(1, s1, s2, 64);

  printf("\n");

  run_ham_test(1, s1, s2, 23);
}

void run_tbe_test(char* tree1_str, char* tree2_str)
{
  pll_utree_t * tree1 = pll_utree_parse_newick_string (tree1_str);
  pll_utree_t * tree2 = pll_utree_parse_newick_string (tree2_str);

  unsigned int split_count = tree1->tip_count - 3;
  pll_unode_t ** node_split_map = calloc(split_count, sizeof(pll_unode_t *));
  double * tbe = calloc(split_count, sizeof(double));

  pllmod_utree_consistency_set(tree1, tree2);

  pll_split_t * splits1 = pllmod_utree_split_create(tree1->vroot,
                                                    tree1->tip_count,
                                                    node_split_map);

  pll_split_t * splits2 = pllmod_utree_split_create(tree2->vroot,
                                                    tree2->tip_count,
                                                    NULL);


  assert(tree1->tip_count == tree2->tip_count);

  pllmod_utree_tbe_naive(splits1, splits2, tree1->tip_count, tbe);

  pllmod_utree_draw_support(tree1, tbe, node_split_map, NULL);

  printf("TBE: ");

  for (unsigned int i = 0; i < split_count; ++i)
    printf("%.6lf ", tbe[i]);

  printf("\n\n");

  char * newick = pll_utree_export_newick(tree1->vroot, NULL);
  printf("TBE tree: %s\n", newick);

  free(newick);
  free(node_split_map);
  free(tbe);

  pllmod_utree_split_destroy(splits1);
  pllmod_utree_split_destroy(splits2);

  pll_utree_destroy (tree1, NULL);
  pll_utree_destroy (tree2, NULL);
}

void test_tbe()
{
  char *ref_tree = "(Woolly:0.02000173,Spider:0.01195957,(Howler:0.03921588,"
      "(((Squirrel:0.04951841,(Tamarin:0.01882103,PMarmoset:0.01872779)1000:0.01620522)432:0.00209062,"
      "(Titi:0.01974091,Saki:0.02183432)999:0.01197670)385:0.00073575,(((Gorilla:0.00549912,"
      "(Human:0.00667950,Chimp:0.00208720)792:0.00128616)986:0.00708195,"
      "(Gibbon:0.02407730,Orangutan:0.01258485)738:0.00147021)937:0.01302782,"
      "(Colobus:0.00276602,(DLangur:0.00477650,(Patas:0.01102645,"
      "((Tant_cDNA:0.00133132,AGM_cDNA:0.00133913)998:0.00516221,"
      "(Rhes_cDNA:0.00595363,Baboon:0.00312241)969:0.00413146)657:0.00250131)1000:0.01235639"
      ")505:0.00123650)1000:0.03064698)1000:0.13115789)998:0.01474962)1000:0.00860350);";
  char *boot1_tree = "((Squirrel:0.04749782,((Saki:0.02577556,Titi:0.02534069):0.01417705,"
      "(Tamarin:0.01830913,PMarmoset:0.01752493):0.01595714):0.00164378):0.00319885,"
      "(Howler:0.03662786,(Spider:0.01128245,Woolly:0.02588956):0.00481877):0.01827684,"
      "(((Gorilla:0.00609643,(Chimp:0.00068926,Human:0.01011787):0.00064788):0.00456013,"
      "(Gibbon:0.02515313,Orangutan:0.00762452):0.00213596):0.01362313,"
      "((DLangur:0.00941860,Colobus:0.00415358):0.00389312,(Patas:0.01861160,"
      "((Baboon:0.00583652,Rhes_cDNA:0.00860553):0.00375633,(Tant_cDNA:0.00133482,"
      "AGM_cDNA:0.00001389):0.00461931):0.00341803):0.01152701):0.03383894):0.15261034);";
  char *boot2_tree = "((Baboon:0.100000,(Colobus:0.100000,(Gibbon:0.100000,"
      "(Tamarin:0.100000,Human:0.100000):0.100000):0.100000):0.100000):0.100000,"
      "(DLangur:0.100000,(AGM_cDNA:0.100000,(Saki:0.100000,((Woolly:0.100000,"
      "Rhes_cDNA:0.100000):0.100000,Chimp:0.100000):0.100000):0.100000):0.100000):0.100000,"
      "(Squirrel:0.100000,((PMarmoset:0.100000,((Patas:0.100000,Tant_cDNA:0.100000):0.100000,"
      "(Spider:0.100000,(Titi:0.100000,"
      "(Howler:0.100000,Orangutan:0.100000):0.100000):0.100000):0.100000):0.100000):0.100000,"
      "Gorilla:0.100000):0.100000):0.100000):0.0;";

  printf("TEST #1:\n");

  run_tbe_test(ref_tree, boot1_tree);

  printf("\nTEST #2:\n");

  run_tbe_test(ref_tree, boot2_tree);
}

int main (int argc, char * argv[])
{
  unsigned int attributes = get_attributes(argc, argv);

  if (attributes != PLL_ATTRIB_ARCH_CPU)
  {
    skip_test();
  }

  printf("Testing Hamming distance:\n\n");

  test_hamming();

  printf("\nTesting Transfer Boostrap Estimate (TBE):\n\n");

  test_tbe();

  return 0;
}
