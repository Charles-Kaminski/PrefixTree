# PrefixTree
CFK 02/23/2016 charles.kaminski@lexisnexis.com

This project is a code bundle for the HPCCSystem big data platform.  The code showcases prefix trees on a big data platform and shows how prefix trees can be used to significantly improve the performance of certain algorithms such as the Levenshtein edit distance algorithm.  The code is a set of three function macros that do the heavy lifting for you.

They are as follows:  
 1. Create - Efficiently creates a prefix tree from a dataset
 2. QueryThorLevenshtein - Uses a dataset to query a prefix tree using Levenshtein
 3. QueryRoxieLevenshtein - Uses a string to query a prefix tree using Levenshtein

Usage examples are at the end of the file PrefixTree.ecl

Read the original code walk-through on the HPCC Systems blog:
https://hpccsystems.com/resources/blog?uid=225 

Read how to install code bundles onto the HPCC platform:
https://github.com/hpcc-systems/HPCC-Platform/blob/master/ecl/ecl-bundle/BUNDLES.rst#installing-a-bundle

Read about other bundles:
https://github.com/hpcc-systems/ecl-bundles
