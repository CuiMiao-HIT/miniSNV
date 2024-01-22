#ifndef DEBGA_INDEX_HPP_
#define DEBGA_INDEX_HPP_

#include <stdint.h>
#include <stdlib.h>
//#include <string>
#include <string.h>
#include <vector>

extern "C"
{
#include "../clib/utils.h"
// #include "../abPOA-v1.4.1/include/abpoa.h"
// #include "../abPOA-v1.4.1/include/simd_instruction.h"
}

#include "math_func.hpp"

//#define LEN_KMER 20 //kmer length used for search
#define ROUTE_LENGTH_MAX 1024
#define MAX_CHR_NAME_LENGTH 200
#define MAX_CHR_NUM 6000000
#define START_POS_REF 0

struct deBGA_INDEX{
private:
	//uint64_t* buffer_ref_seq = NULL;//original reference
	//uint64_t* buffer_seq = NULL;// UNITIG sequence
	uint64_t* buffer_seqf = NULL;//start offset of [UINTIG] in [buffer_seq]
	//uint64_t* buffer_off_g = NULL;//start offset of [KMER] in [buffer_seq]
	uint64_t* buffer_p = NULL;//start offset of [UINTIG] in [buffer_ref_seq]
	uint64_t* buffer_pp = NULL;//start index of [UINTIG] in [buffer_p], used to search offset of [UINTIG] in [buffer_ref_seq]
	//uint64_t* buffer_hash_g = NULL;//index for the first 14 base pair, used to search kmer in the hash table
	//uint32_t* buffer_kmer_g = NULL;//used to store other part of a kmer(except the first 14 byte)
	uint32_t chr_search_index_size = 0;
	uint32_t * chr_search_index = NULL;

	uint64_t reference_len = 0;
	uint64_t *chr_end_n;
	//char chr_names[MAX_CHR_NUM][MAX_CHR_NAME_LENGTH];
	char chr_line_content[MAX_CHR_NAME_LENGTH];
	int chr_file_n = 1;

//	bam_hdr_t * header = NULL;

//	uint64_t result_ref_seq = 0;
//	uint64_t result_seq = 0;
	uint64_t result_seqf = 0;
	uint64_t result_p = 0;
	uint64_t result_pp = 0;
//	uint64_t result_pu = 0;
//	uint64_t result_hash_g = 0;
//	uint64_t result_kmer_g = 0;
//	uint64_t result_off_g = 0;
//	uint64_t result_ref_g = 0;

public:

	//load from file, return the true load data size(in byte)
	uint64_t load_data_from_file(const char * path_name, const char * fn, void ** data, bool useCalloc, int additional_byte){
		fprintf(stderr, "END loading %s\n", fn);
		char full_fn[ROUTE_LENGTH_MAX] = {0};
		strcpy(full_fn, path_name);
		if(full_fn[strlen(full_fn) - 1] != '/')	 strcat(full_fn, "/");
		strcat(full_fn, fn);
		FILE *fp_us_b = xopen (full_fn, "rb" );
		fseek(fp_us_b, 0, SEEK_END);// non-portable
		uint64_t file_size = ftell(fp_us_b);
		rewind(fp_us_b);
		if(useCalloc)	(*data) = (uint64_t* ) xcalloc (file_size + additional_byte, 1);
		else			(*data) = (uint64_t* ) xmalloc (file_size + additional_byte);
		xread ((*data), 1, file_size, fp_us_b);
		fclose(fp_us_b);
		return file_size;
	}

	int load_index_file(const char *unitig_index_dir)
	{
//	    //        buffer_ref_seq:store reference seq;load from dm file ref.seq
//	    result_ref_seq = (load_data_from_file(index_dir, "ref.seq", ((void ** )&buffer_ref_seq), true, 536) >> 3);
//	    //        buffer_seq:store unipath seq(concatenate all unipaths's seq);load from dm file unipath.seqb
//	    result_seq = (load_data_from_file(index_dir, "unipath.seqb", ((void ** )&buffer_seq), false, 0) >> 3);
//      buffer_seqf:store unipath's offset on unipath seq;load from dm file unipath.seqfb
	    result_seqf = (load_data_from_file(unitig_index_dir, "unipath.seqfb", ((void ** )&buffer_seqf), false, 0) >> 3);
	    //buffer_p: store every unipath's set of positions on reference;load from dm file unipath.pos
	    result_p = (load_data_from_file(unitig_index_dir,  "unipath.pos", ((void ** )&buffer_p), false, 0) >> 3);
	    //buffer_pp: unipath's pointer to array buffer_p;load from dm file unipath.posp
	    result_pp = (load_data_from_file(unitig_index_dir,   "unipath.posp", ((void ** )&buffer_pp), false, 0) >> 3);

//	    for(int u = 0; u < 100; u++){
//			uint64_t ref_pos_n = buffer_pp[u + 1] - buffer_pp[u];
//			ref_pos_n = MIN(100, ref_pos_n);
//			bool find_match = false;
//			for(uint64_t i = 0; i <ref_pos_n; i++){
//				uint64_t gloabl_offset = buffer_p[i + buffer_pp[u]];
//				fprintf(stderr, "u %d, i %d , gloabl_offset %ld \n", u, i, gloabl_offset);
//			}
//	    }
//	    exit(-1);

//	    //buffer_hash_g: store kmer's hash part; load from dm file unipath_g.hash
//	    result_hash_g = (load_data_from_file(index_dir,  "unipath_g.hash", ((void ** )&buffer_hash_g), false, 0) >> 3);
//	// buffer_kmer_g:store kmer's kmer part;load from dm file unipath_g.kmer
//	    result_kmer_g = (load_data_from_file(index_dir,  "unipath_g.kmer", ((void ** )&buffer_kmer_g), false, 0) >> 2);
//	//buffer_off_g: store kmer's offset on unipath seq; load from dm file unipath_g.offset
//	    result_off_g = (load_data_from_file(index_dir,  "unipath_g.offset", ((void ** )&buffer_off_g), false, 0) >> 2);

	    //********************************************************************************************
	    //read chr names and length from unipath.chr
	    fprintf(stderr, "LOADING unipath.chr\n");

	    char unichr[ROUTE_LENGTH_MAX] = {0};
	    strcpy(unichr, unitig_index_dir);
	    strcat(unichr, "unipath.chr");
	    //fprintf(stderr,"fn: %s",unichr);
	    std::vector<std::string> unichr_str_v;
	    load_strings_from_file(unichr, unichr_str_v);
	    chr_end_n = (uint64_t *)xcalloc((unichr_str_v.size()/2 + 10), sizeof(uint64_t));
		for(uint i = 0; i < unichr_str_v.size();){
			//sscanf(unichr_str_v[i].c_str(),"%s",chr_names[chr_file_n]);
			i++;
			chr_end_n[chr_file_n++] = atol(unichr_str_v[i].c_str()); i++;
		}
	    if(false){
	    	//for(uint i = 0; i < chr_file_n;i++){
				//fprintf(stderr,"%s %ld\n",chr_names[i], chr_end_n[i]);
			//}
	    }

	    chr_end_n[0] = START_POS_REF + 1; //START_POS_REF = 0 record the start position of reference
	    //strcpy(chr_names[chr_file_n], "*");
	    reference_len = chr_end_n[chr_file_n - 1];
	    //building search index for CHR
	    building_chr_index();
	    return 0;
	}

//	//function: given an index of kmer, search the MEM within a UNITIG
//	//return 0 when successfully stored data, -1 when failed
//	int UNITIG_MEM_search(int32_t unitig_ID, int32_t unitig_offset){
//		uint64_t ref_pos_n = buffer_pp[unitig_ID + 1] - buffer_pp[unitig_ID];
//		for(uint64_t i = 0; i <ref_pos_n; i++){
//			uint64_t gloabl_offset =  buffer_p[i + buffer_pp[unitig_ID]] + unitig_offset - 1;
//			  int chr_ID = get_chromosome_ID(gloabl_offset);
//			  int POS = gloabl_offset - chr_end_n[chr_ID - 1];//end of last chr is the begin of this chr
//		}
//		return 0;
//	}

	void building_chr_index(){
		chr_search_index_size = (reference_len >> 14) + 2;
		chr_search_index = (uint32_t *)xcalloc(chr_search_index_size, sizeof(uint32_t));
		uint32_t pos_index_size = 0;
		for(int i = 0; i < chr_file_n; i++){
			int pos_index = chr_end_n[i] / 0x4000;
			while(pos_index >= pos_index_size){
				chr_search_index[pos_index_size++] = i;
			}
		}
		chr_search_index[pos_index_size] = chr_file_n;// store final block
	}

	inline int get_chromosome_ID(uint64_t position){
		int chr_n = 0;
		int pos_index = position / 0x4000;
		int low = chr_search_index[pos_index];
		int high = chr_search_index[pos_index + 1];
		int mid;
		uint64_t pos = position + 1;

		while ( low <= high ){
			mid = (low + high) >> 1;
			if(pos < (chr_end_n[mid] - 1))		high = mid - 1;
			else if(pos > (chr_end_n[mid] - 1))	low = mid + 1;
			else								return mid;
			chr_n = low;
		}
		return chr_n;
	}

	uint64_t get_ref_N(uint32_t unitig_ID){ return buffer_pp[unitig_ID + 1] - buffer_pp[unitig_ID]; }

	uint64_t get_global_offset(uint32_t unitig_ID, uint32_t unitig_offset, uint32_t ref_idx){
		return buffer_p[ref_idx + buffer_pp[unitig_ID]] + unitig_offset;
	}
	inline void get_chr_ID_POS(uint64_t global_offset, int &chr_ID, int &POS){
		  chr_ID = get_chromosome_ID(global_offset);
		  POS = global_offset - chr_end_n[chr_ID - 1];//end of last chr is the begin of this chr
		  if(chr_ID == 1) POS -= 2048;
		  chr_ID -=1;
	}

	uint32_t get_UNITIG_length(uint32_t unitig_ID)
	{
		return buffer_seqf[unitig_ID + 1] - buffer_seqf[unitig_ID];
	}

	void free_memory(){
//		if(buffer_ref_seq 	!= NULL){ free(buffer_ref_seq); buffer_ref_seq = NULL; 	} //original reference
//		if(buffer_seq 		!= NULL){ free(buffer_seq); 	buffer_seq = NULL; 		}// UNITIG sequence
		if(buffer_seqf != NULL)		{ free(buffer_seqf); 	buffer_seqf = NULL; }//start offset of [UINTIG] in [buffer_seq]
//		if(buffer_off_g != NULL)	{ free(buffer_off_g); 	buffer_off_g = NULL; }//start offset of [KMER] in [buffer_seq]
		if(buffer_p != NULL)		{ free(buffer_p); 		buffer_p = NULL; }//start offset of [UINTIG] in [buffer_ref_seq]
		if(buffer_pp != NULL)		{ free(buffer_pp); 		buffer_pp = NULL; }//start index of [UINTIG] in [buffer_p], used to search offset of [UINTIG] in [buffer_ref_seq]
//		if(buffer_hash_g != NULL)	{ free(buffer_hash_g); 	buffer_hash_g = NULL; }//index for the first 14 base pair, used to search kmer in the hash table
//		if(buffer_kmer_g != NULL)	{ free(buffer_kmer_g); 	buffer_kmer_g = NULL; }//used to store other part of a kmer(except the first 14 byte)
	}

};
// uint64_t load_data_from_file(const char * path_name, const char * fn, void ** data, bool useCalloc, int additional_byte);
// int build_deBGA_index(int argc, char *argv[]);

#endif /* DEBGA_INDEX_HPP_ */
