// CFK 02/23/2016 charles.kaminski@lexisnexis.com
// The below is a code bundle for the HPCCSystem big data platform.
// The code showcases prefix trees on a big data platform
//   and shows how prefix trees can be used to significantly improve
//   the performance of certain algorithms.  The below code 
//   significantly improves the performance of the Levenshtein
//   edit distance algorithm over the standard approach. 
//     
// It has three parts:  
//   1. Create - efficiently creates a prefix tree on a
//               horizontally scalable big data platform
//               from a dataset
//   2. QueryThorLevenshtein - Uses a dataset to query a
//                             prefix tree using Levenshtein
//   3. QueryRoxieLevenshtein - Uses a string to query a
//                              prefix tree using Levenshtein
//
// Usage examples are at the end of this file.
//
// Read the original code walk-through on the HPCC Systems blog:
// https://hpccsystems.com/resources/blog?uid=225 
//
// Read how to install code bundles here:
// https://github.com/hpcc-systems/HPCC-Platform/blob/master/ecl/ecl-bundle/BUNDLES.rst#installing-a-bundle
//
// Read about other bundles:
// https://github.com/hpcc-systems/ecl-bundles

EXPORT PrefixTree := MODULE

	IMPORT STD;

	EXPORT Bundle := MODULE(Std.BundleBase)
		EXPORT Name := 'PrefixTree';
		EXPORT Description := 'Function Macros to create a prefix tree and to search using Edit Distance.';
		EXPORT Authors := ['Charles Kaminski'];
		EXPORT License := 'http://www.apache.org/licenses/LICENSE-2.0';
		EXPORT Copyright := 'Copyright 2016-2017';
		EXPORT DependsOn := [];
		EXPORT Version := '0.1.0';
		EXPORT PlatformVersion := '6.0.8';
	END;
								
//================================================================================================================================================================
//=============================== Create - Create Prefix Tree ====================================================================================================
//================================================================================================================================================================

// The "Create" function macro creates a prefix.

	EXPORT Create(infile, infield, dist_length) := FUNCTIONMACRO
		 // This Function Macro creates a prefix tree 
		 // The prefix tree can be used to significantly speed up
		 // Edit-Distance searches and  also enable those searches at query time
		 // Inputs:
		 // infile  - The Thor file holding your data
		 // infield - The name of the field in the infile holding your data
		 // dist_length - The number of characters in the infield
		 //     used to distribute the data across the Thor cluster.
		 //     You want this number to be small enough so that like
		 //       records are on the same Thor.  But large enough so that
		 //       all of your Thor nodes are being used.
		 //       1 is a good value for small datasets.
		 //       2 is a good value for large datasets and clusters.
		 
		 // The code below assumes a maximum number of Thor compute nodes  
		 //  and a maximum number of prefix tree nodes on each Thor node.
		 //  These maximums are suffient for any of our datasets and can 
		 //  be adjusted to fit any future needs.  The current maximum
		 //  number of Thor nodes assumed below is 92,233 compute nodes.  
		 //  The code below assumes no more than 72 Trillion prefix tree nodes
		 //  on <<EACH>> of the 92,233 Thor nodes.
		 
		 WordLayout := RECORD
				STRING   word;
		 END;
		 WordLayout FormatWordsTransform(RECORDOF(infile) l) := TRANSFORM
				SELF.word := TRIM(l.infield);
		 END;
		 
		 // Trim the infield and get rid of any other fields
		 cleaned_ds := PROJECT(infile, FormatWordsTransform(LEFT));
		 
		 // Distribute the files across the Thor based on the dist_length parameter
		 distributed_ds := DISTRIBUTE(cleaned_ds, HASH(cleaned_ds.word[dist_length]));
		 
		 // Get rid of any empty entries
		 input_ds := distributed_ds(word <> '');
		 
		 NodesLayout := RECORD
				input_ds;
				// The nodes field is where we keep track of when nodes begin.  We're initializing the field to 
				// zero fill with the same number of zeros as characters in each word.
				// The extra '1' at the end is a place holder for an end-cap node representing the word itself.
				STRING    nodes        := REGEXREPLACE('.', input_ds.Word, '0') + '1'; 
				UNSIGNED4 node_cnt     := 0;        // The number of nodes in the word
				UNSIGNED4 compute_node := 0;        // The compute_node this record is on			
				DATA      node_ids     := (DATA)''; // Where we pack up the node IDs right before we normalize the prefix tree records
		 END;
		 
		 nodes_format_ds := TABLE(input_ds, NodesLayout, UNSORTED, LOCAL);
		 STRING MarkNodes(STRING l, STRING l_nodes, STRING n, STRING n_nodes) := BEGINC++
				/// This C++ function is designed for use with ECL's Iterate to take a group of 
				/// strings and return the prefix tree node structures in a separate field.
				/// The function requires a sorted list and two passes (one in each direction) 
				/// l is for the last field, n is for the next field, out is for output
				/// nodes contain the indicator of if we should start a node
				/// An initial value of the node for each field is a string of 0 [zeros]
				/// the same length as the string in question -- an optional 1 can be
				/// added at the end as a word end-cap.
				/// CFK 08/05/2015   
				#option pure
				#body
				unsigned long i = 0;
				char * out = (char*)rtlMalloc(lenN_nodes);
				memcpy(out, n_nodes, lenN_nodes);
				while ( i < lenN){
					 /// If the characters don't match or
					 ///   if we are at the end of the last string then
					 ///   flip the 0 to a 1 and break out of the loop
					 if (l[i] != n[i] || i >= lenL){
							out[i] = '1';
							break;
					 }
					 if (l_nodes[i] != n_nodes[i]){
							out[i] = '1';
					 }
					 i += 1;
				}
				__lenResult = lenN_nodes;
				__result = out;
		 ENDC++;
		 
		 NodesLayout GetPTNodesTransform(NodesLayout L, NodesLayout R) := TRANSFORM
				SELF.nodes        := MarkNodes(L.word, L.nodes, R.word, R.nodes);
				SELF.node_cnt     := stringlib.StringFindCount(SELF.nodes,'1');
				SELF.compute_node := Thorlib.Node();			
				SELF := R;
		 END;
		 // Identify where the nodes for the prefix tree will be in your word field
		 pt_nodes_1_ds := ITERATE(SORT(nodes_format_ds, -word, LOCAL), GetPTNodesTransform(LEFT, RIGHT), LOCAL);
		 pt_nodes_2_ds := ITERATE(SORT(pt_nodes_1_ds,    word, LOCAL), GetPTNodesTransform(LEFT, RIGHT), LOCAL);
		 
		 STRING AssignNodeIDs(STRING l, STRING l_nodes, UNSIGNED4 l_node_cnt, STRING l_node_ids,
													STRING n, STRING n_nodes, UNSIGNED4 n_node_cnt, UNSIGNED4 compute_node) := BEGINC++
				/// This C++ function is designed for use with ECL's Iterate
				/// It assigns ID's to prefix tree nodes so that the nodes can later be broken out into 
				/// individual records with ECL's Normalize.
				/// CFK 08/05/2015   
				#option pure
				#body      
				unsigned long long last_node_id = 0;
				unsigned long i = 0;
				unsigned long node_cnt = 0;
				// The new value is empty.  
				// Could happen if there are blanks in your word list.
				if (n_node_cnt == 0){
					 __lenResult = 0;
					 return;
				}
				
				// Initialize "out" with node ID's from the last iteration
				//  But don't go over if there are too many node ID's from the last iteration
				char * out = (char*)rtlMalloc(n_node_cnt * 8);   
				i = (l_node_cnt < n_node_cnt) ? l_node_cnt : n_node_cnt;
				memcpy(out, l_node_ids, i * 8);
				// Find the ID assigned to the last node from the last word.
				//  This will be the starting point to assign new node IDs
				if (l_node_cnt != 0){
					 memcpy(&last_node_id, &l_node_ids[(l_node_cnt - 1) * 8], 8);
				}
				else{
					 // This is the starting point for node id on each compute node.
					 // Setting last_node_id to a value based on the compute node's
					 // unique id ensures that each prefix tree built locally on each 
					 // compute node has unique non-overlapping ID's.  An unsigned long long
					 // has a maximum value of 9,223,372,036,854,775,807 or 2^64 - 1.   
					 // Multiplying compute_node by 10^14 allows for 92,233 compute nodes  
					 // and 72 Trillion prefix tree nodes on each compute node.  These maximums
					 // can be adjusted to fit your needs.
					 last_node_id = (compute_node * 100000000000000);
				}	  
				// Find the first location where the two words don't match
				//  And keep track of how many nodes we've seen
				i = 0;
				while (i < lenL && i < lenN){
					 if (l[i] != n[i]){
							break;
					 }
					 if (n_nodes[i] == '1'){
							node_cnt += 1;
					 }
					 i += 1;
				}
				
				// Start adding new node ids to  from where we left off
				while (i < lenN){
					 if (n_nodes[i] == '1'){
							last_node_id += 1;
							node_cnt += 1;
							memcpy(&out[(node_cnt-1)*8], &last_node_id, 8);
					 }
					 i += 1;
				}
				
				// Add the final node id for the end-cap (referenced in comments at the top)
				last_node_id += 1;
				node_cnt +=1;
				memcpy(&out[(node_cnt-1) * 8], &last_node_id, 8);
				__lenResult = n_node_cnt * 8;
				__result    = out;
		 ENDC++;
		 NodesLayout GetNodeIdsTransform(NodesLayout L, NodesLayout R) := TRANSFORM
				SELF.compute_node := Thorlib.Node();
				SELF.node_ids     := (DATA) AssignNodeIDs(L.word, L.nodes, L.node_cnt, (STRING) L.node_ids,
																									R.word, R.nodes, R.node_cnt, SELF.compute_node);
				SELF := R;
		 END;
		 node_ids_ds := ITERATE(SORT(pt_nodes_2_ds,  word, LOCAL), GetNodeIdsTransform(LEFT, RIGHT), LOCAL);
		 UNSIGNED GetID(STRING node_ids, UNSIGNED4 node_to_retrieve) := BEGINC++
				/// Function returns a numeric node id from a group of node ids
				/// CFK 08/05/2015   
				#option pure
				#body   
				unsigned long long node_id = 0;
				if (node_to_retrieve != 0){
					 memcpy(&node_id, &node_ids[(node_to_retrieve - 1) * 8], 8);
				}
				return node_id;
		 ENDC++;
		 STRING GetNode(STRING word, STRING nodes, UNSIGNED4 node_to_retrieve):= BEGINC++
				/// Function returns the substring associated with a node
				/// CFK 08/05/2015
				#option pure
				#body   
				unsigned int i = 0;
				unsigned int j = 0;
				unsigned int k = 0;
		 
				// find the start (i) of the node in question
				while(i < lenNodes){
					 if(nodes[i] == '1'){
							j +=1;
					 }
					 if (j==node_to_retrieve){
							break;
					 }
					 i += 1 ;
				}
		 
				j = i + 1;
		 
				// find the end (j) of the node in question
				while(j < lenNodes){
					 if (nodes[j] == '1'){
							break;
					 }
					 j += 1;
				}
		 
				k = j - i;
				char * out = (char*)rtlMalloc(k);
				memcpy(&out[0], &word[i], k);
				__lenResult = k;
				__result   = out;
		 
		 ENDC++;
		 PTLayout := RECORD
				UNSIGNED  id;            // Primary Key
				UNSIGNED  parent_id;     // The parent for this node.  Zero is a root node
				UNSIGNED1 _max;          // The max length of any child word for this node
				UNSIGNED1 _min;          // The min length of any child word for this node	 
				BOOLEAN   is_word;       // Indicates if this node is a word (an end-cap with no children)
				UNSIGNED4 compute_node;  // The compute-node ID this record was processed on			
				STRING    node;          // Contains the payload for this node			
				//UNSIGNED8 RecPtr {virtual(fileposition)};			
		 END;
		 PTLayout NormalizePTTransform(NodesLayout L, UNSIGNED4 C):= TRANSFORM
				SELF.id           := GetID((STRING)L.node_ids, C);
				SELF.parent_id    := GetID((STRING)L.node_ids, C-1);
				SELF.node         := IF(C=L.node_cnt, L.word, GetNode(L.word, L.nodes, C));
				SELF._max         := LENGTH(L.word);
				SELF._min         := LENGTH(L.word);			
				SELF.is_word      := IF(C=L.node_cnt, True, False);
				SELF.compute_node := Thorlib.Node();			
		 END; 
		 // Break each word record into multiple node records
		 normalized_pt_ds := NORMALIZE(node_ids_ds, LEFT.node_cnt, NormalizePTTransform(LEFT,COUNTER));
		 
		 // Sort the dataset for a rollup
		 sorted_normalized_pt_ds := SORT(normalized_pt_ds, id, LOCAL);
		 PTLayout RollupPTTransform(PTLayout L, PTLayout R) := TRANSFORM
				SELF._max := MAX(L._max, R._max);
				SELF._min := MIN(L._min, R._min);
				SELF      := R;
		 END;
		 //  Remove duplicates on each Thor Node and calculate Max and Min lengths
		 //  https://hpccsystems.com/resources/blog/ckaminski/accelerating-prefix-trees
		 deduped_pt_ds := ROLLUP(sorted_normalized_pt_ds, LEFT.id = RIGHT.id, 
													 RollupPTTransform(LEFT, RIGHT), LOCAL);
		 RETURN deduped_pt_ds;
	ENDMACRO;
	
	//================================================================================================================================================================
	//================ QueryThorLevenshtein - Query a Prefix Tree on a Thor big data platform using a Dataset and the Levenshtein edit distance algorithm ============
	//================================================================================================================================================================
	
	// CFK 02/23/2016 charles.kaminski@lexisnexis.com
	// This code will use a prefix tree and a query dataset
	//   in conjunction with the Levenshtein string matching 
	//   algorithm to significantly improve the performance of the 
	//   edit distance algorythm.  The code is written to run on the Thor
	//   big data platform.

	Export QueryThorLevenshtein(query_ds, query_field, pt_index, max_edit_distance) := FUNCTIONMACRO
		 STRING CalculateLevenshteinVector(STRING word, STRING node, STRING state) := BEGINC++
				/// CFK 08/20/2015 charles.kaminski@lexisnexis.com                 
				/// This C++ returns the vector used to calculate
				///  the Levenshtein distance in an interative process so that
				///  an efficient trie structure can be traversed when finding candidates
				/// Node values can be passed in and the current state can be saved 
				///   in the return vector so that other nodes can be further processed in series
				/// We're using char arrays to keep the size down
				/// That also means that all our words (both in the trie and our query words)
				///  must be less than 255 characters in length.  If that limitation becomes a problem
				///  we can easily use an array of integers or another other structure to store
				///  intermediate values.
				#option pure
				#include <algorithm>
				#body
				//  The last element in the array will hold
				//   the min value.  The min value will be used later.
				int v_size = lenWord + 2;
		
				unsigned char min_value = 255;
				unsigned char new_value = 0;
		
				// Since v0 is not returned, we should not use rtlMalloc
				//unsigned char * v0 = (unsigned char*)rtlMalloc(v_size);	
				unsigned char v0[v_size];
				// We use rtlMalloc helper function when a variable is returned by the function
				unsigned char * v1 = (unsigned char*)rtlMalloc(v_size);
		
				if (lenState < 1){
					 for (int i = 0; i < v_size; i++){
							v0[i] = i;
					 }
				}
				else{
					 memcpy(&v0[0], &state[0], v_size);
				}
		
				int cost = 0;
				int k = v0[0];
						
				for (int i=k; i<k+lenNode; i++)
				{
					 min_value = 255;
					 v1[0] = i + 1;
					 for (int j=0; j<lenWord; j++)
					 {
							cost = (node[i-k] == word[j]) ? 0 : 1;
							new_value = std::min(v1[j] + 1, std::min(v0[j+1] + 1, v0[j] + cost));
							v1[j+1] = new_value;
							if (new_value < min_value){
								 min_value=new_value;
							}
					 }
					 memcpy(&v0[0], &v1[0], lenState);
				}
		
				// Store the min_value;
				v1[v_size-1] = min_value;
		
				__lenResult = v_size;
				__result    = (char *) v1;
																																	
		 ENDC++;
	 
		 UNSIGNED1 GetMinDistance(STRING state) := BEGINC++
				/// CFK 08/20/2015 charles.kaminski@lexisnexis.com	
				///  Get the Minimum Edit Distance
				#option pure
				#body		
				//return (unsigned char) state[(unsigned int) position];
				return (unsigned char) state[lenState-1];
		 ENDC++;
		 UNSIGNED1 GetFinalDistance(STRING state) := BEGINC++
				/// CFK 08/20/2015 charles.kaminski@lexisnexis.com	
				///  Get the Final Edit Distance
				#option pure
				#body		
				//return (unsigned char) state[(unsigned int) position];
				return (unsigned char) state[lenState-2];
		 ENDC++;
		 
		 QueryPTLayout := RECORD
				STRING    query_string;
				STRING    state                := '';
				DATA      state_data           := (DATA)'';
				UNSIGNED  node_id              := 0;
				STRING    node                 := '';
				UNSIGNED1 _max                 := 255;
				UNSIGNED1 _min                 := 0;
				BOOLEAN   is_word              := False;
				STRING    cumulative_nodes     := '';
				UNSIGNED1 cumulative_node_size := 0;
				UNSIGNED1 current_distance     := 0;
				UNSIGNED1 final_distance       := 0;
		 END;
			
		 QueryPTLayout InitialQueryDSTransform(RECORDOF(query_ds) L) := TRANSFORM
				SELF.query_string := L.query_field;
		 END;
		 preped_query_ds := PROJECT(query_ds, InitialQueryDSTransform(LEFT));
		 QueryPTLayout QueryPTTransform(QueryPTLayout L, RECORDOF(pt_index) R) := TRANSFORM
				SELF.query_string         := L.query_string;
				SELF.state                := IF(R.is_word, L.state, CalculateLevenshteinVector(L.query_string, R.node, L.state));
				SELF.state_data           := (DATA)SELF.state;
				SELF.node_id              := R.id;
				SELF.node                 := R.node;
				SELF._max                 := R._max;    // The minimum word size for all the children of this node
				SELF._min                 := R._min;    // The maximum word size for all the children of this node
				SELF.is_word              := R.is_word;
				SELF.cumulative_node_size := IF(R.is_word, LENGTH(R.node), LENGTH(R.node) + L.cumulative_node_size);
				SELF.cumulative_nodes     := IF(R.is_word, R.node, L.cumulative_nodes + R.node);
				SELF.current_distance     := IF(R.is_word, L.current_distance, GetMinDistance(SELF.state));
				SELF.final_distance       := IF(R.is_word, GetFinalDistance(SELF.state), L.final_distance);	
		 END;
	 
		 looped_ds := LOOP(preped_query_ds, 
				LEFT.is_word = False,
				EXISTS(ROWS(LEFT)) = True,
				JOIN(ROWS(LEFT), pt_index, LEFT.node_id = RIGHT.parent_id AND
																	 LEFT.current_distance <= max_edit_distance AND
																	 // Current_distance alone cannot be incorporated into the next two expressions
																	 // To do so would "double count" certain edit-distance situations;
																	 // One example would be comparing "dog" and "drop" as seen in the example 1_ThorInlineDataset.ecl from
																	 // https://github.com/Charles-Kaminski/FastEditDistanceUsingBigDataPrefixTrees
																	 // Don't do this -> (LENGTH(LEFT.word)) >= (RIGHT._min - (max_edit_distance - LEFT.current_distance))
																	 LENGTH(LEFT.query_string) <= (RIGHT._max + max_edit_distance) AND
																	 LENGTH(LEFT.query_string) >= (RIGHT._min - max_edit_distance),
					 QueryPTTransform(LEFT,RIGHT), INNER)); 														 
		 interim_results_ds := looped_ds(looped_ds.final_distance <= max_edit_distance);
		 
		 RETURN PROJECT(interim_results_ds, {QueryPTLayout.query_string, QueryPTLayout.node, QueryPTLayout.final_distance});
	ENDMACRO;
	
  //=================================================================================================================================================================
	//================ QueryRoxieLevenshtein - Query a prefix tree with a string using a Roxie data delivery engine and the Levenshtein edit distance algorithm ========
	//=================================================================================================================================================================	
	
	// CFK 02/23/2016 charles.kaminski@lexisnexis.com
  // This code uses a prefix tree and a query string
	//   in conjunction with the Levenshtein string matching 
	//   algorithm to significantly improve the performance of the 
	//   edit distance algorythm.  The code is written to run on the Roxie
	//   big data delivery platform.
	
	Export QueryRoxieLevenshtein(input_query_string, pt_index, max_edit_distance) := FUNCTIONMACRO
		 STRING CalculateLevenshteinVector(STRING word, STRING node, STRING state) := BEGINC++
				/// CFK 08/20/2015 charles.kaminski@lexisnexis.com                 
				/// This C++ returns the vector used to calculate
				///  the Levenshtein distance in an interative process so that
				///  an efficient trie structure can be traversed when finding candidates
				/// Node values can be passed in and the current state can be saved 
				///   in the return vector so that other nodes can be further processed in series
				/// We're using char arrays to keep the size down
				/// That also means that all our words (both in the trie and our query words)
				///  must be less than 255 characters in length.  If that limitation becomes a problem
				///  we can easily use an array of integers or another other structure to store
				///  intermediate values.
				#option pure
				#include <algorithm>
				#body
				//  The last element in the array will hold
				//   the min value.  The min value will be used later.
				int v_size = lenWord + 2;
		
				unsigned char min_value = 255;
				unsigned char new_value = 0;
		
				// Since v0 is not returned, we should not use rtlMalloc
				//unsigned char * v0 = (unsigned char*)rtlMalloc(v_size);	
				unsigned char v0[v_size];
				// We use rtlMalloc helper function when a variable is returned by the function
				unsigned char * v1 = (unsigned char*)rtlMalloc(v_size);
		
				if (lenState < 1){
					 for (int i = 0; i < v_size; i++){
							v0[i] = i;
					 }
				}
				else{
					 memcpy(&v0[0], &state[0], v_size);
				}
		
				int cost = 0;
				int k = v0[0];
						
				for (int i=k; i<k+lenNode; i++)
				{
					 min_value = 255;
					 v1[0] = i + 1;
					 for (int j=0; j<lenWord; j++)
					 {
							cost = (node[i-k] == word[j]) ? 0 : 1;
							new_value = std::min(v1[j] + 1, std::min(v0[j+1] + 1, v0[j] + cost));
							v1[j+1] = new_value;
							if (new_value < min_value){
								 min_value=new_value;
							}
					 }
					 memcpy(&v0[0], &v1[0], lenState);
				}
		
				// Store the min_value;
				v1[v_size-1] = min_value;
		
				__lenResult = v_size;
				__result    = (char *) v1;
																																	
		 ENDC++;
	 
		 UNSIGNED1 GetMinDistance(STRING state) := BEGINC++
				/// CFK 08/20/2015 charles.kaminski@lexisnexis.com	
				///  Get the Minimum Edit Distance
				#option pure
				#body		
				//return (unsigned char) state[(unsigned int) position];
				return (unsigned char) state[lenState-1];
		 ENDC++;
		 UNSIGNED1 GetFinalDistance(STRING state) := BEGINC++
				/// CFK 08/20/2015 charles.kaminski@lexisnexis.com	
				///  Get the Final Edit Distance
				#option pure
				#body		
				//return (unsigned char) state[(unsigned int) position];
				return (unsigned char) state[lenState-2];
		 ENDC++;
		 
		 QueryPTLayout := RECORD
				STRING    query_string;
				STRING    state                := '';
				DATA      state_data           := (DATA)'';
				UNSIGNED  node_id              := 0;
				STRING    node                 := '';
				UNSIGNED1 _max                 := 255;
				UNSIGNED1 _min                 := 0;
				BOOLEAN   is_word              := False;
				STRING    cumulative_nodes     := '';
				UNSIGNED1 cumulative_node_size := 0;
				UNSIGNED1 current_distance     := 0;
				UNSIGNED1 final_distance       := 0;
		 END;
			
		 preped_query_ds := DATASET([{input_query_string}], QueryPTLayout);
		 QueryPTLayout QueryPTTransform(QueryPTLayout L, RECORDOF(pt_index) R) := TRANSFORM
				SELF.query_string         := L.query_string;
				SELF.state                := IF(R.is_word, L.state, CalculateLevenshteinVector(L.query_string, R.node, L.state));
				SELF.state_data           := (DATA)SELF.state;
				SELF.node_id              := R.id;
				SELF.node                 := R.node;
				SELF._max                 := R._max;    // The minimum word size for all the children of this node
				SELF._min                 := R._min;    // The maximum word size for all the children of this node
				SELF.is_word              := R.is_word;
				SELF.cumulative_node_size := IF(R.is_word, LENGTH(R.node), LENGTH(R.node) + L.cumulative_node_size);
				SELF.cumulative_nodes     := IF(R.is_word, R.node, L.cumulative_nodes + R.node);
				SELF.current_distance     := IF(R.is_word, L.current_distance, GetMinDistance(SELF.state));
				SELF.final_distance       := IF(R.is_word, GetFinalDistance(SELF.state), L.final_distance);	
		 END;
	 
		 looped_ds := LOOP(preped_query_ds, 
				LEFT.is_word = False,
				EXISTS(ROWS(LEFT)) = True,
				JOIN(ROWS(LEFT), pt_index, LEFT.node_id = RIGHT.parent_id AND
																	 LEFT.current_distance <= max_edit_distance AND
																	 // Current_distance alone cannot be incorporated into the next two expressions
																	 // To do so would "double count" certain edit-distance situations;
																	 // One example would be comparing "dog" and "drop" as seen in the example 1_ThorInlineDataset.ecl from
																	 // https://github.com/Charles-Kaminski/FastEditDistanceUsingBigDataPrefixTrees
																	 // Don't do this -> (LENGTH(LEFT.word)) >= (RIGHT._min - (max_edit_distance - LEFT.current_distance))
																	 LENGTH(LEFT.query_string) <= (RIGHT._max + max_edit_distance) AND
																	 LENGTH(LEFT.query_string) >= (RIGHT._min - max_edit_distance),
					 QueryPTTransform(LEFT,RIGHT), INNER)); 														 
		 interim_results_ds := looped_ds(looped_ds.final_distance <= max_edit_distance);
		 
		 RETURN PROJECT(interim_results_ds, {QueryPTLayout.node, QueryPTLayout.final_distance});
	ENDMACRO;	

END;

/*
//======================================================================================================================================
//=============== Usage Example - Create Prefix Tree ===================================================================================
//======================================================================================================================================

InFileLayout := RECORD
   STRING   word;
END;

in_ds := DATASET('~OUT::CFK::LNames_Atl_GA_Words', InFileLayout, THOR);
SMALL_DATASET_OR_COMPUTE_CLUSTER := 1;
LARGE_DATASET_AND_CLUSTER := 2;

// This dataset will be distributed based on the first 2 letters in each word
//  Query performance can be improved by using a 1 for small datasets or small compute clusters.  
//  2 is appropriate for large datasets and large compute clusters

pt_ds := <<MODULE LOCATION>>.PrefixTree.Create(in_ds, word, LARGE_DATASET_AND_CLUSTER);

OUTPUT(pt_ds, , '~OUT::CFK::LNames_Atl_GA_PT', OVERWRITE);
pt_index := INDEX(pt_ds, {parent_id}, {id, _max, _min, is_word, node}, '~OUT::CFK::LNames_Atl_GA_PT_Key');
BUILDINDEX(pt_index, OVERWRITE);
OUTPUT(pt_index);

//======================================================================================================================================
//== Usage Example - Query Prefix Tree on Thor using a query dataset ===================================================================
//======================================================================================================================================

PTTIndexLayoutThor := RECORD
   unsigned8 parent_id;
   unsigned8 id;
   unsigned1 _max;
   unsigned1 _min;
   boolean is_word;
   string node;
END;

pt_index_thor := INDEX(DATASET([], PTTIndexLayoutThor), {parent_id}, {id, _max, _min, is_word, node}, '~OUT::CFK::LNames_Atl_GA_PT_Key'); 

WordLayout := RECORD
   STRING   word;
END;

// An inline dataset is used here for simplicity sake.  An existing on-disk dataset will work the same.
query_ds := DATASET([{'KAMINSKI'}, {'BAYLISS'}, {'MUHAREMAGIC'}], WordLayout);

results_thor := <<MODULE LOCATION>>.PrefixTree.QueryThorLevenshtein(query_ds, word, pt_index_thor, 2);                                      

OUTPUT(results_thor);

//======================================================================================================================================
//== Usage Example - Query Prefix Tree on Roxie using a string =========================================================================
//======================================================================================================================================

PTIndexLayoutRoxie := RECORD
   unsigned8 parent_id;
   unsigned8 id;
   unsigned1 _max;
   unsigned1 _min;
   boolean is_word;
   string node;
END;

pt_index_roxie := INDEX(DATASET([], PTIndexLayoutRoxie), {parent_id}, {id, _max, _min, is_word, node}, '~OUT::CFK::LNames_Atl_GA_PT_Key'); 

results_roxie := <<MODULE LOCATION>>.PrefixTree.QueryRoxieLevenshtein('KAMINSKI', pt_index_roxie, 1);

OUTPUT(results_roxie);
*/
