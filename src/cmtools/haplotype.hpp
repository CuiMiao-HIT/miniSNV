/*
 * haplotype.hpp
 *
 *  Created on: 2022年4月11日
 *      Author: fenghe
 */

#ifndef HAPLOTYPE_HPP_
#define HAPLOTYPE_HPP_
#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <map>
#include <set>
#include <zlib.h>
#include <fstream>
extern "C"{
	#include "../htslib/htslib/faidx.h"
}

#include "deBGA_index.hpp"
#include "math_func.hpp"

struct window_block_info{
	window_block_info(uint8_t chrID_, uint32_t region_st_, uint32_t region_ed_, uint64_t global_offset_){
		global_offset = global_offset_;
		region_st = region_st_;
		region_length = region_ed_ - region_st_;
		chrID = chrID_;
		flag = 0;
	}

	void set(uint8_t chrID_, uint32_t region_st_, uint32_t region_ed_, uint64_t global_offset_){
		global_offset = global_offset_;
		region_st = region_st_;
		region_length = region_ed_ - region_st_;
		chrID = chrID_;
		flag = 0;
	}
	window_block_info(const window_block_info & from){
		global_offset = from.global_offset;
		region_st = from.region_st;
		region_length = from.region_length;
		chrID = from.chrID;
		flag = from.flag;
	}
	window_block_info(){
		global_offset = 0;
		region_st = 0;
		region_length = 0;
		chrID = 0;
		flag = 0;
	}

	void show_wb(FILE * out){
		fprintf(out, "@ window block, global_offset %ld, region_st %d, region_length %d, chrID %d, flag %d\n", global_offset, region_st,region_length, chrID, flag);
	}

	bool is_SV(){
		return ((flag & 0x1) == 1);
	}

	uint64_t global_offset; //window block DATA in hap_database
	uint32_t region_st; //window block POS in reference
	uint16_t region_length;//window block length
	uint8_t chrID;
	uint8_t flag; //the bit7:( == 0) when it is from SNP/INDEL; ( == 1) when it is from SV;
};

struct Hap_seq_address{
	uint8_t *p;
	uint64_t len;
	void clear(){
		p = NULL;
		len = 0;
	}
};

struct Hap_seq_modify_info{
	Hap_seq_modify_info(int original_pos_){
		original_pos = original_pos_;
		is_modified = false;
	}
	Hap_seq_modify_info(){
		original_pos = -1;
		is_modified = false;
	}
	char is_modified;
	int original_pos;
	void print(){
		fprintf(stderr, " %s @ pos %d\n", is_modified?"TRUE":"False", original_pos);
	}
};

// X：[30 bit hash value] + [1 bit REF or ALT] + [1 bit Reverse or Forward]
// Y + Z：[25bit window ID @ MAX 32M space] + [19bit: hap ID @ MAX 0.5M space] + [20bit: hap offset @ MAX 1M]
struct mm_96_t{
	mm_96_t(uint32_t x_, uint32_t y_,uint32_t z_){
		x = x_;
		y = y_;
		z = z_;
	}
	uint32_t x, y, z;
};

struct HAP_VAR_ITEM{
	HAP_VAR_ITEM(uint32_t ref_pos_, std::string &REF_, std::string &ALT_){
		ref_pos = ref_pos_;
		std::swap(REF,REF_);
		std::swap(ALT,ALT_);
	}
	uint32_t ref_pos = 0;
	std::string REF;
	std::string ALT;

	void print(FILE * out){
		fprintf(out, "%s_%s_%d\t", REF.c_str(), ALT.c_str(), ref_pos);
	}
};

struct Window_block_handler{
private:
	//
	std::vector<char> rst_haplotype_BUFF;
	/***
	 * unfold a hap seq and generate a haplotype string
	 * input:: (1) a pointer to the bin hap seq
	 * input:: (2) length of the bin hap seq
	 * input:: (3) the string of reference
	 *
	 * Output: the haplotype sequence
	 */

	std::vector<Hap_seq_modify_info> modify_info;
public:
	std::vector<char> & get_string(){ return rst_haplotype_BUFF; }
	std::vector<Hap_seq_modify_info> & get_modify_info(){ return modify_info; }
	//debug code:: print modify inf0
	void printf_modify_info(FILE* out){
		for(auto & h : rst_haplotype_BUFF){ fprintf(out, "%c", h); }
		fprintf(out, "\n");
		for(auto & mi : modify_info){ fprintf(out, "%c", (mi.is_modified == true)?'1':'0'); }
		fprintf(out, "\n");
	}

	std::vector<HAP_VAR_ITEM> var_l;
	std::vector<char> & unfold_hap_seq(uint8_t * bin_hap_seq, uint64_t hap_seq_len, std::string & ref, bool generate_modify_string, bool generate_var_list, int kmer_size_alt_string);

	std::vector<char> & unfold_hap_seq_SV(uint8_t * bin_hap_seq, uint64_t hap_seq_len, std::string & ref, window_block_info & c_wb, bool generate_modify_string, uint32_t kmer_size_alt_string);

	uint8_t * window_block_p_old = NULL;
	std::vector<Hap_seq_address> block_all_hap_list;
	void unfold_window_block_all(uint8_t * window_block_p);
	/***
	 * unfold one hap seq in a window block
	 * input:: (1) a pointer to the window block
	 * input:: (2) ID of the hap seq
	 *
	 * Output: block_i_rst_address of the haplotype sequence
	 */
	Hap_seq_address block_i_rst_address;
	void unfold_window_block_i(uint8_t * window_block_p, uint64_t hap_i);
};


struct Simple_ref_handler{
private:
	//flag:
	bool ref_store_in_binary;
	//reference stored in ACGT
	std::vector<std::string> all_ref_data_list_ACGT;

	//reference store in binary
	struct ref_info{
		ref_info(){
			ref_length = 0;
			ref_offset = 0;
		}
		ref_info(std::string &s){
			std::vector<std::string> item_value;
			char * temp = (char *)xcalloc(1024,1);
			split_string(item_value, temp, s.c_str(), " ");
			name = item_value[0];
			ref_length = atol(item_value[1].c_str());
			ref_offset = atol(item_value[2].c_str());
			if(temp != NULL){
				free(temp); temp = NULL;
			}
		}
		uint32_t ref_length;
		uint64_t ref_offset;
		std::string name;
		void dump(FILE * out){
			fprintf(out, "%s %d %ld\n", name.c_str(), ref_length, ref_offset);
		}
	};
	std::vector<ref_info> all_ref_data_list;
	uint8_t *bin_ref_data;

	inline uint8_t char2bin(char c){
		switch(c){
		case 'A': case 'a': return 0;  break;
		case 'C': case 'c': return 1;  break;
		case 'G': case 'g': return 2;  break;
		case 'T': case 't': return 3;  break;
		case 'N': case 'n': return 2;  break;
		}
		return 0;
	}
public:
	//dump the reference info into bin_ref_info_f, dump the bin ref into bin_ref_f
	void dump_bin_ref(FILE * bin_ref_f, FILE * bin_ref_info_f){
		xassert(ref_store_in_binary == false, "");
		uint64_t total_uint8_t_size = all_ref_data_list.back().ref_offset + (all_ref_data_list.back().ref_length + 3)/4;
		bin_ref_data = (uint8_t *)xcalloc(total_uint8_t_size, 1);
		for(uint32_t chr_ID = 0;chr_ID < all_ref_data_list_ACGT.size(); chr_ID ++){
			uint64_t ref_global_offset = all_ref_data_list[chr_ID].ref_offset;
			std::string& ref_string = all_ref_data_list_ACGT[chr_ID];
			for(uint32_t offset = 0; offset < ref_string.size(); offset++){
				uint8_t bin_char = char2bin(ref_string[offset]);
				bin_ref_data[((offset) >> 2) + ref_global_offset] |= bin_char << ((3 - ((offset) & 0X3)) << 1);
			}
		}

		//dump files
		for(uint32_t chr_ID = 0;chr_ID < all_ref_data_list_ACGT.size(); chr_ID ++){
			all_ref_data_list[chr_ID].dump(bin_ref_info_f);
		}
		fwrite(bin_ref_data, 1, total_uint8_t_size, bin_ref_f);
		if(bin_ref_data != NULL){
			free(bin_ref_data); bin_ref_data = NULL;
		}
	}

	void load_bin_ref(const char * path_name){
		ref_store_in_binary = true;
		std::vector<std::string> chr_if_str;
		//load info file
		char full_fn[1024] = {0};
		strcpy(full_fn, path_name);
		if(full_fn[strlen(full_fn) - 1] != '/')	 strcat(full_fn, "/");
		strcat(full_fn, "/ref_info.bin");
		load_strings_from_file(full_fn, chr_if_str);
		all_ref_data_list.clear();
		for(auto & s: chr_if_str){
			all_ref_data_list.emplace_back(s);
		}
		//load data file:
		strcpy(full_fn, path_name);
		if(full_fn[strlen(full_fn) - 1] != '/')	 strcat(full_fn, "/");
		strcat(full_fn, "/ref.bin");
		FILE * ref_data_f = xopen(full_fn, "rb");
		uint64_t total_uint8_t_size = all_ref_data_list.back().ref_offset + (all_ref_data_list.back().ref_length + 3)/4;
		bin_ref_data = (uint8_t *) xmalloc(total_uint8_t_size);
		xread(bin_ref_data, 1, total_uint8_t_size, ref_data_f);
		fclose(ref_data_f);
	}

	//window block part loader, simple load reference from ref.fa
	void load_reference_from_file(const char * ref_fn){
		ref_store_in_binary = false;
		faidx_t * fai = fai_load(ref_fn);
		uint32_t ref_n = faidx_nseq(fai);
		uint64_t total_uint8_t_size = 0;
//		if(true){ debug_load_ref = 9; }//debug code
		for(uint32_t chr_ID = 0; chr_ID < ref_n; chr_ID++){
			int true_region_load_len = 0;
			char * ref_pointer = fai_fetch(fai, faidx_iseq(fai, chr_ID), &true_region_load_len);
			all_ref_data_list_ACGT.emplace_back(ref_pointer);

			std::string& ref_string = all_ref_data_list_ACGT[chr_ID];
			all_ref_data_list.emplace_back();
			all_ref_data_list.back().name.append(faidx_iseq(fai, chr_ID));
			all_ref_data_list.back().ref_length = ref_string.size();
			all_ref_data_list.back().ref_offset = total_uint8_t_size;
			total_uint8_t_size += (ref_string.size() + 3)/4;

			if(ref_pointer != NULL) free(ref_pointer);
		}
		fai_destroy(fai);
	}

	void load_ref_from_buff_global_offset(uint64_t global_offset, int length, std::string & rst, deBGA_INDEX * dbi){
		int chr_ID, POS;
		dbi->get_chr_ID_POS(global_offset, chr_ID, POS);
		load_ref_from_buff(chr_ID, POS, length, rst);
	}

	//window block part loader
	void load_ref_from_buff(int chr_ID, uint32_t st_pos, int length, std::string & rst){
		rst.clear();
		if(chr_ID < 0 || chr_ID >= (int)all_ref_data_list.size())	return;
		if(ref_store_in_binary){
			int32_t max_load = (all_ref_data_list[chr_ID].ref_length - st_pos + 1);
			max_load = MIN(max_load, length);
			if(max_load < 0)
				return;
			rst.resize(max_load);
			uint64_t c_ref_global_offset = all_ref_data_list[chr_ID].ref_offset;
			for(int32_t i = 0; i < max_load; i++){
				uint8_t b_c = ((bin_ref_data[((st_pos + i - 1) >> 2) + c_ref_global_offset]) >> (((3 - ((st_pos + i - 1) & 0X3)) << 1))) & 0x3;
				rst[i] = "ACGT"[b_c];
			}
		}else{
			rst.insert(rst.begin(), all_ref_data_list_ACGT[chr_ID].begin() + st_pos - 1, all_ref_data_list_ACGT[chr_ID].begin() + st_pos - 1 + length);
		}
	}

	std::string & get_chr_name(int chr_ID){
		return all_ref_data_list[chr_ID].name;
	}
	uint64_t get_chr_length(int chr_ID){
		return all_ref_data_list[chr_ID].ref_length;
	}
	uint64_t get_chr_N(){
		return all_ref_data_list.size();
	}
};


#define SV_IDX_bucket_size 0x1fff //8K, 3G / 8k = 0.375M(blocks)
struct MM_idx_loader{
	//PART1： deBGA index
	//(2)UNITIG IDs
	deBGA_INDEX unitig_idx;
	//PART2:UNITIG MM index
	mm_96_t* ref_mm;		uint64_t ref_mm_size;
	uint32_t* ref_mm_idx;	uint64_t ref_mm_idx_size;
	//PART3: hap strings
	window_block_info* wb_info; uint64_t wb_info_size;
	uint8_t * wb_data; uint64_t wb_data_size;
	//PART3.1: hap search
	std::vector<int> chr_bg_wb_ID;
	int total_chr_number = 0;//total number of chr
	int step_length = 0;
	//uint32_t* sv_wb_idx;	uint64_t sv_wb_idx_size;
	bool idx_with_SV = false;
	std::vector<std::vector<uint32_t>> sv_search_idx;

	//part4: ALT SEQ MM IDX
	mm_96_t* alt_mm; uint64_t alt_mm_size;
	uint32_t* alt_mm_idx; uint64_t alt_mm_idx_size;
	//part5: reference
	Simple_ref_handler ref;

	//global align buff
	// ALIGN_BUFF_global g;

private:
	//load from file, return the true load data size(in byte)
	uint64_t load_data_from_file(const char * path_name, const char * fn, void ** data, int skip_byte, bool must_load){
		// fprintf(stderr, "END loading %s\n", fn);
		char full_fn[1024] = {0};
		strcpy(full_fn, path_name);
		if(full_fn[strlen(full_fn) - 1] != '/')	 strcat(full_fn, "/");
		strcat(full_fn, fn);
		FILE *fp_us_b = NULL;
		if(must_load){
			fp_us_b = xopen (full_fn, "rb" );
		}else{
			fp_us_b = fopen (full_fn, "rb" );
			if(fp_us_b == NULL){
				(*data) = NULL;
				return 0;
			}
		}
		fseek(fp_us_b, 0, SEEK_END);// non-portable
		uint64_t file_size = ftell(fp_us_b);
		rewind(fp_us_b);
		fseek(fp_us_b, skip_byte, SEEK_SET);// non-portable
		(*data) = (uint64_t* ) xmalloc (file_size - skip_byte);
		xread ((*data), 1, file_size - skip_byte, fp_us_b);
		fclose(fp_us_b);
		return file_size;
	}

	uint64_t load_data_from_file_part(const char * path_name, const char * fn, void ** data, uint64_t bg, uint64_t ed){
		char full_fn[1024] = {0};
		strcpy(full_fn, path_name);
		if(full_fn[strlen(full_fn) - 1] != '/')	 strcat(full_fn, "/");
		strcat(full_fn, fn);
		FILE *fp_us_b = xopen (full_fn, "rb" );
		fseek(fp_us_b, bg, SEEK_SET);// non-portable
		(*data) = (uint64_t* ) xmalloc (ed - bg);
		xread ((*data), 1, ed - bg, fp_us_b);
		fclose(fp_us_b);
		return ed - bg;
	}

	void* load_chr_bg_wb_ID(const char * idx_path_name){
		chr_bg_wb_ID.clear();
		char full_fn[1024] = {0};
		strcpy(full_fn, idx_path_name);
		if(full_fn[strlen(full_fn) - 1] != '/')	 strcat(full_fn, "/");
		strcat(full_fn, "chr_bg_wb_ID.bin");
		load_int_from_file(full_fn, chr_bg_wb_ID);
		return (void *)this;
	}
public:
	//bool skip_unitig_mm_idx, when set to be true, only load wb index, when set false, load all index
	//set unitig_path_name to NULL when only load wb data
	void * load_all_index(const char * idx_path_name, const char * unitig_path_name, bool skip_unitig_mm_idx){
		//PART3: hap strings
		wb_info_size = load_data_from_file(idx_path_name, "wb_info.bin", (void **)&wb_info, 8, true) / sizeof(window_block_info);
		wb_data_size = load_data_from_file(idx_path_name, "wb_data.bin", (void **)&wb_data, 0, true);
		//PART3.1: wb data
		step_length = wb_info[0].region_length/2;
		load_chr_bg_wb_ID(idx_path_name);

		if(wb_info[wb_info_size - 1].is_SV()){
			idx_with_SV = true;
			total_chr_number = chr_bg_wb_ID.size() - 2;
		}
		else
			total_chr_number = chr_bg_wb_ID.size() - 1;

		if(!skip_unitig_mm_idx){
			//PART1： deBGA index
			unitig_idx.load_index_file(unitig_path_name);
			//part5: reference
			ref.load_bin_ref(idx_path_name);
			//PART2:UNITIG MM index
			ref_mm_size = load_data_from_file(idx_path_name, "ref_mm.bin", (void **)&ref_mm, 0, true) / sizeof(mm_96_t);
			ref_mm_idx_size = load_data_from_file(idx_path_name, "ref_mm_idx.bin", (void **)&ref_mm_idx, 0, true) / sizeof(uint32_t);
			//part4: ALT SEQ MM IDX
			alt_mm_size = load_data_from_file(idx_path_name, "alt_mm.bin", (void **)&alt_mm, 0, true) / sizeof(mm_96_t);
			alt_mm_idx_size = load_data_from_file(idx_path_name, "alt_mm_idx.bin", (void **)&alt_mm_idx, 0, true) / sizeof(uint32_t);
			//part 3.1: sv idx
			//sv_wb_idx_size = load_data_from_file(idx_path_name, "sv_wb_idx.bin", (void **)&sv_wb_idx, 0, false);
			if(idx_with_SV){//with SV, and build index for it
				uint32_t SV_wb_bg = chr_bg_wb_ID[total_chr_number];
				uint32_t SV_wb_ed = chr_bg_wb_ID[total_chr_number + 1];
				sv_search_idx.resize(ref.get_chr_N());
				for(uint chr_ID = 0; chr_ID < sv_search_idx.size(); chr_ID++){
					uint32_t block_number = ref.get_chr_length(chr_ID)/(SV_IDX_bucket_size) + 2;
					sv_search_idx[chr_ID].resize(block_number);
				}
				for(uint wb_id = SV_wb_bg; wb_id < SV_wb_ed; wb_id++){
					uint32_t blockID = wb_info[wb_id].region_st/(SV_IDX_bucket_size);
					sv_search_idx[wb_info[wb_id].chrID][blockID + 1] = wb_id + 1;
				}
				uint32_t old_idx = SV_wb_bg;
				for(std::vector<uint32_t> & chr :sv_search_idx){
					for(uint32_t & idx : chr){
						if(idx == 0){	idx = old_idx;}
						else{			old_idx = idx;}
					}
				}
			}
		}
		return (void *)this;
	}

	void load_window_ID_idx(const char * idx_path_name){
		//try load wb_info
		wb_info_size = load_data_from_file_part(idx_path_name, "wb_info.bin", (void **)&wb_info, sizeof(window_block_info)*0 + 8, sizeof(window_block_info)*1 + 8) / sizeof(window_block_info);
		step_length = wb_info[0].region_length/2;
		if(wb_info != NULL) {free (wb_info); wb_info = NULL;}
		load_chr_bg_wb_ID(idx_path_name);
	}

	void load_part_index(const char * idx_path_name, uint32_t wb_bg, uint32_t wb_ed){
		wb_info_size = load_data_from_file_part(idx_path_name, "wb_info.bin", (void **)&wb_info, sizeof(window_block_info)*wb_bg + 8, sizeof(window_block_info)*(wb_ed+1) + 8) / sizeof(window_block_info);
		wb_data_size = load_data_from_file_part(idx_path_name, "wb_data.bin", (void **)&wb_data, wb_info[0].global_offset, wb_info[wb_ed - wb_bg].global_offset);
		wb_info_size -= 1;
		uint64_t wb_data_bg_offset = wb_info[0].global_offset;
		for(int i = 0; i < wb_info_size; i++){
			wb_info[i].global_offset -= wb_data_bg_offset;
		}
	}

};


struct String_list_and_var_list{
	std::vector<std::string> hap_string_l;
	std::vector<std::vector<HAP_VAR_ITEM>> var_l;
	void print(FILE * out){
		for(int i = 0; i < hap_string_l.size(); i++){
			for(auto & v: var_l[i]){
				v.print(out);
			}
			fprintf(out, "\n: %s\n", hap_string_l[i].c_str());
		}
		fprintf(out, "\n");
	}

};

struct hap_string_loader_single_thread{
private:
	//PART3: hap strings
	Window_block_handler wbh;
	//load from file, return the true load data size(in byte)
	std::vector<std::string> hap_string_buff;
	String_list_and_var_list slvl;
	std::vector<std::vector<HAP_VAR_ITEM>> var_l;
	std::vector<std::string> & get_string_list_core(uint32_t window_ID, std::string & window_ref, window_block_info* wb_info, uint8_t* wb_data, bool generate_var_list){
		hap_string_buff.clear();
		window_block_info & c_wb = wb_info[window_ID];
		wbh.unfold_window_block_all(wb_data + c_wb.global_offset);
		//load reference:
		std::vector<Hap_seq_address> & hap_l = wbh.block_all_hap_list;
		//clear minimizer builder:
		if(generate_var_list){
			var_l.clear();
		}
		if(c_wb.is_SV()){
			std::vector<char> & hap_string = wbh.unfold_hap_seq_SV(hap_l[0].p, hap_l[0].len, window_ref, c_wb, false, 0);
			hap_string_buff.emplace_back(&(hap_string[0]));
		}else{
			for(auto & hap:hap_l){//load all haplotypes
				std::vector<char> & hap_string = wbh.unfold_hap_seq(hap.p, hap.len, window_ref, false, generate_var_list, 0);
				if(generate_var_list){
					var_l.emplace_back(wbh.var_l);
				}
				hap_string.emplace_back(0);
				hap_string_buff.emplace_back(&(hap_string[0]));
			}
		}
		if(generate_var_list){
			std::swap(slvl.hap_string_l, hap_string_buff);
			std::swap(slvl.var_l, var_l);
		}
		return hap_string_buff;
	}

	std::vector<char> & get_string_core(uint32_t window_ID, uint32_t hap_ID, std::string & window_ref, window_block_info* wb_info, uint8_t* wb_data){
		window_block_info & c_wb = wb_info[window_ID];
		wbh.unfold_window_block_all(wb_data + c_wb.global_offset);
		//load reference:
		std::vector<Hap_seq_address> & hap_l = wbh.block_all_hap_list;
		if(c_wb.is_SV()){
			xassert(hap_ID == 0, "SV has only one haplotype in one window block");
			return wbh.unfold_hap_seq_SV(hap_l[0].p, hap_l[0].len, window_ref, c_wb, false, 0);
		}else{
			return wbh.unfold_hap_seq(hap_l[hap_ID].p, hap_l[hap_ID].len, window_ref, false, false, 0);
		}
	}


	//step_length is 150 normally
	//return the corresponding wb ID
	//
	int get_windows_ID_core(uint8_t chr_ID, int POS, std::vector<int> &chr_bg_wb_ID, int step_length){
		return chr_bg_wb_ID[chr_ID] + (POS - 1)/step_length;
	}

	std::vector<uint32_t> SV_wb_list_buff;
	//search SVs in a region:
	//search all SVs in region chr_ID:POS-END;
	//the result (WB ID) stored in wb_bg(include) and wb_ed (not include)
	void search_SV_IN_region_core(uint8_t chr_ID, uint32_t POS, uint32_t END,  window_block_info* wb_info, std::vector<std::vector<uint32_t>> &sv_search_idx, uint32_t &wb_bg, uint32_t &wb_ed){
		uint32_t search_bg = sv_search_idx[chr_ID][POS/SV_IDX_bucket_size];
		uint32_t search_ed = sv_search_idx[chr_ID][END/SV_IDX_bucket_size + 1];
		for(uint32_t i = search_bg; i < search_ed;i++){
			if(wb_info[i].region_st >= POS){
				wb_bg = i; break;
			}
		}
		for(uint32_t i = search_ed - 1; i >= search_bg; i--){
			if(wb_info[i].region_st < END){
				wb_ed = i; break;
			}
		}
		wb_ed += 1;
	}
public:
	//function used in python : 1
	int get_windows_ID(uint8_t chr_ID, int POS, void* i_){
		MM_idx_loader * i = (MM_idx_loader *)i_;
		return get_windows_ID_core(chr_ID, POS, i->chr_bg_wb_ID, i->step_length);
	}
	//search SVs in a region:
	//search all SVs in region chr_ID:POS-END;
	//the result (WB ID) stored in wb_bg(include) and wb_ed (not include)
	//return the SV number
	std::vector<uint32_t> &search_SV_IN_region(uint8_t chr_ID, uint32_t POS, uint32_t END, void* i_){
		MM_idx_loader * i = (MM_idx_loader *)i_;
		uint32_t wb_bg = UINT32_MAX; uint32_t wb_ed = 0;
		search_SV_IN_region_core(chr_ID, POS, END, i->wb_info, i->sv_search_idx, wb_bg, wb_ed);
		SV_wb_list_buff.clear();
		for(uint32_t i = wb_bg; i < wb_ed; i++){
			SV_wb_list_buff.emplace_back(i);
		}
		return SV_wb_list_buff;
	}
	std::vector<char> & get_string(uint32_t window_ID, uint32_t hap_ID, std::string & window_ref, void* i_){
		MM_idx_loader * i = (MM_idx_loader *)i_;
		return get_string_core(window_ID, hap_ID, window_ref, i->wb_info, i->wb_data);
	}
	//function used in python : 1
	String_list_and_var_list & get_string_list_and_var_list(uint32_t window_ID, std::string & window_ref, void* i_){
		MM_idx_loader * i = (MM_idx_loader *)i_;
		get_string_list_core(window_ID, window_ref, i->wb_info, i->wb_data, true);
		return slvl;
	}

	//function used in python : 1
	std::vector<std::string> & get_string_list(uint32_t window_ID, std::string & window_ref, void* i_){
		MM_idx_loader * i = (MM_idx_loader *)i_;
		return get_string_list_core(window_ID, window_ref, i->wb_info, i->wb_data, false);
	}

	//function used in python : 2
	// std::vector<std::vector<uint8_t>> & get_bin_data_list(uint32_t window_ID,  void* i_){
	// 	MM_idx_loader * i = (MM_idx_loader *)i_;
	// 	return get_bin_data_core(window_ID,  i->wb_info, i->wb_data);
	// }

	uint64_t get_hap_N(uint32_t window_ID,  void* i_){
		MM_idx_loader * i = (MM_idx_loader *)i_;
		hap_string_buff.clear();
		window_block_info & c_wb = i->wb_info[window_ID];
		wbh.unfold_window_block_all(i->wb_data + c_wb.global_offset);
		return wbh.block_all_hap_list.size();
	}

public:
	//debug code, not use
	void printf_modify_info_after_get_string(){
		wbh.printf_modify_info(stderr);
	}

};

#endif /* HAPLOTYPE_HPP_ */
