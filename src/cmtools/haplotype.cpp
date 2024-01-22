/*
 * haplotype.cpp
 *
 *  Created on: 2022年4月14日
 *      Author: fenghe
 */

#include "haplotype.hpp"
#include "SV_handler.hpp"

#define MIN_sweet_region 0
#define MAX_sweet_region 300

std::vector<char> & Window_block_handler::unfold_hap_seq(uint8_t * bin_hap_seq, uint64_t hap_seq_len, std::string & ref, bool generate_modify_string, bool generate_var_list, int kmer_size_alt_string){
	//get hap seq:
	std::vector<char> & rst_haplotype = rst_haplotype_BUFF;
	rst_haplotype.clear();
	//modify string init:
	if(generate_modify_string){
		kmer_size_alt_string = MAX(kmer_size_alt_string, 22);
		modify_info.clear();
		int modify_info_size = ref.size();
		for(int i = 0; i < modify_info_size; i++){
			modify_info.emplace_back(i);
		}
	}
	uint32_t c_ref_pos = 0;
	uint32_t c_alt_pos = 0;
	if(generate_var_list){
		var_l.clear();
	}
	//generate hap strings
	for(uint32_t i = 0; i < hap_seq_len;){
		if((bin_hap_seq[i] & 1) == 1){//match
			int match_len = 0;
			for(;i < hap_seq_len && (bin_hap_seq[i] & 1) == 1; i++){ match_len += (bin_hap_seq[i] >> 1); }
			//debug code:
			match_len = (c_ref_pos + match_len > ref.size())?(ref.size() - c_ref_pos):match_len;
			uint32_t max_load_pos = MIN(c_ref_pos + match_len, ref.size());
			if(c_ref_pos < max_load_pos)
				rst_haplotype.insert(rst_haplotype.end(), ref.begin() + c_ref_pos, ref.begin() + max_load_pos);
			c_ref_pos += match_len;
			c_alt_pos += match_len;
		}
		else if((bin_hap_seq[i] & 2) == 0){//SNP or DEL
			uint8_t real_data = (bin_hap_seq[i] >> 2);
			if(real_data <= 63 && real_data >= 60){//SNP
				char SNV_char = "ACGT"[(real_data) - 60];
				i++;
				rst_haplotype.emplace_back(SNV_char);
				if(generate_modify_string){
					for(int mod_kmer_i = 0; mod_kmer_i < kmer_size_alt_string; mod_kmer_i++){
						if(c_alt_pos - mod_kmer_i >= 0 && c_alt_pos - mod_kmer_i < (long int)modify_info.size()){
							modify_info[c_alt_pos - mod_kmer_i].is_modified =  true;
						}
					}
				}
				if(generate_var_list){
					std::string REF; REF.insert(REF.begin(),ref.begin() + c_ref_pos, ref.begin() + c_ref_pos + 1);
					std::string ALT; ALT.insert(ALT.begin(),rst_haplotype.begin() + c_alt_pos, rst_haplotype.begin() + c_alt_pos + 1);
					var_l.emplace_back(c_ref_pos, REF, ALT);
					if(c_ref_pos < MIN_sweet_region || c_ref_pos > MAX_sweet_region ){
						rst_haplotype.erase(rst_haplotype.begin() + c_alt_pos, rst_haplotype.begin() + c_alt_pos + var_l.back().ALT.size());
						rst_haplotype.insert(rst_haplotype.begin() + c_alt_pos, var_l.back().REF.begin(), var_l.back().REF.end());
						var_l.erase(var_l.end() - 1);
					}
				}
				c_ref_pos += 1;
				c_alt_pos += 1;
			}else{//short del
				uint32_t del_length = real_data;
				xassert(del_length > 0 && del_length < 50,"");
				i++;
				if(generate_modify_string){
					for(int mod_kmer_i = 1; mod_kmer_i < kmer_size_alt_string; mod_kmer_i++){
						if(c_alt_pos - mod_kmer_i >= 0 && c_alt_pos - mod_kmer_i < (long int)modify_info.size()){
							modify_info[c_alt_pos - mod_kmer_i].is_modified =  true;
						}
					}
				}
				if(generate_var_list){
					std::string REF;  std::string ALT;
					if(c_ref_pos >= 1){
						REF.insert(REF.begin(),ref.begin() + c_ref_pos - 1, ref.begin() + c_ref_pos + del_length);
						ALT.insert(ALT.begin(),ref.begin() + c_ref_pos - 1, ref.begin() + c_ref_pos);
					}else{//deletion reach the edge of windows block
						REF.insert(REF.begin(),ref.begin() + c_ref_pos, ref.begin() + c_ref_pos + del_length);
						ALT.append("_");
					}
					var_l.emplace_back(c_ref_pos - 1, REF, ALT);

					if(c_ref_pos < MIN_sweet_region || c_ref_pos > MAX_sweet_region ){
						//rst_haplotype.erase(rst_haplotype.begin() + c_alt_pos, rst_haplotype.begin() + c_alt_pos + var_l.back().ALT.size());
						rst_haplotype.insert(rst_haplotype.begin() + c_alt_pos, var_l.back().REF.begin() + 1, var_l.back().REF.end());
						var_l.erase(var_l.end() - 1);
						//c_ref_pos -= del_length;
						c_alt_pos += del_length;//
					}
				}
				c_ref_pos += del_length;
			}
		}else{//short INS or SV
			uint8_t real_data = (bin_hap_seq[i] >> 2);
			if((real_data) == 63){//long DEL
				uint32_t del_length = ((bin_hap_seq[i+1]) << 16) + (bin_hap_seq[i + 2] << 8) + (bin_hap_seq[i+3]);
				i+= 4;
				//skip the long DEL in this function
				if(false){
					if(generate_modify_string){
						for(int mod_kmer_i = 1; mod_kmer_i < kmer_size_alt_string; mod_kmer_i++){
							if(c_alt_pos - mod_kmer_i >= 0 && c_alt_pos - mod_kmer_i < (long int)modify_info.size()){
								modify_info[c_alt_pos - mod_kmer_i].is_modified =  true;
							}
						}
					}
					c_ref_pos += del_length;//todo:: consider the condition when the deletion is too long
				}
			}else{//INS
				bool is_long_ins = ((real_data) == 62);
				uint32_t ins_length = 0;
				if(is_long_ins){
					ins_length = ((bin_hap_seq[i+1]) << 16) + (bin_hap_seq[i + 2] << 8) + (bin_hap_seq[i+3]);
					i += 4 + ins_length;
				}
				else{
					ins_length = real_data;
					i += 1 + ins_length;
				}
				//skip the long INS in this function
				if(!is_long_ins){
					uint64_t be_ins_pos = rst_haplotype.size();
					rst_haplotype.resize(rst_haplotype.size() + ins_length);
					for(uint32_t i_char_idx = 0; i_char_idx < ins_length; i_char_idx++){
						rst_haplotype[be_ins_pos + i_char_idx] = bin_hap_seq[i - ins_length + i_char_idx];
					}
					if(generate_modify_string){
						modify_info.resize(rst_haplotype.size());
						for(uint32_t mod_kmer_i = 0; mod_kmer_i < ins_length; mod_kmer_i++){
							if(c_alt_pos + mod_kmer_i < (long int)modify_info.size()){
								modify_info[c_alt_pos + mod_kmer_i].is_modified =  true;
							}
						}
						for(int mod_kmer_i = 1; mod_kmer_i < kmer_size_alt_string; mod_kmer_i++){
							if(c_alt_pos - mod_kmer_i >= 0 && c_alt_pos - mod_kmer_i < (long int)modify_info.size()){
								modify_info[c_alt_pos - mod_kmer_i].is_modified =  true;
							}
						}
					}

					if(generate_var_list){
						std::string REF;  std::string ALT;
						if(c_ref_pos >= 1){
							REF.insert(REF.begin(),rst_haplotype.begin() + c_alt_pos - 1, rst_haplotype.begin() + c_alt_pos);
							ALT.insert(ALT.begin(),rst_haplotype.begin() + c_alt_pos - 1, rst_haplotype.begin() + c_alt_pos + ins_length);
						}else{//deletion reach the edge of windows block
							REF.append("_");
							ALT.insert(ALT.begin(),rst_haplotype.begin() + c_alt_pos, rst_haplotype.begin() + c_alt_pos + ins_length);
						}
						var_l.emplace_back(c_ref_pos - 1, REF, ALT);

						if(c_ref_pos < MIN_sweet_region || c_ref_pos > MAX_sweet_region ){
							rst_haplotype.erase(rst_haplotype.begin() + c_alt_pos, rst_haplotype.begin() + c_alt_pos + var_l.back().ALT.size() - 1);
							//rst_haplotype.insert(rst_haplotype.begin() + c_alt_pos, var_l.back().REF.begin(), var_l.back().REF.end());
							var_l.erase(var_l.end() - 1);
							c_alt_pos -= ins_length;
						}
					}
					c_alt_pos += ins_length;
				}
			}
		}
	}
	if(generate_modify_string){
		modify_info.resize(rst_haplotype.size());
	}
	return rst_haplotype;
}

std::vector<char> & Window_block_handler::unfold_hap_seq_SV(uint8_t * bin_hap_seq, uint64_t hap_seq_len, std::string & ref, window_block_info & c_wb, bool generate_modify_string, uint32_t kmer_size_alt_string){
	//get hap seq:
	std::vector<char> & rst_haplotype = rst_haplotype_BUFF;
	rst_haplotype.clear();

	int bin_hap_se_idx = 0;
	int sv_type = SV_TYPE::flag_2_type(bin_hap_seq[bin_hap_se_idx]); bin_hap_se_idx += 1;//SV type
	bin_hap_se_idx += 2;//

	int sv_length = SV_STORE_BASIC::get_32bit(bin_hap_seq + bin_hap_se_idx); bin_hap_se_idx += 4;//SV length
	bin_hap_se_idx += 6;//pos of BP 1
	bin_hap_se_idx += 6;//pos of BP 2
	uint32_t hap_seq_length = SV_STORE_BASIC::get_32bit(bin_hap_seq + bin_hap_se_idx);

	bin_hap_se_idx += 4;//hap seq_length
	rst_haplotype.resize(hap_seq_length + 1);
	for(uint32_t i = 0; i < hap_seq_length; i++){
		rst_haplotype[i] = (char)bin_hap_seq[bin_hap_se_idx + i];
	}
	rst_haplotype[hap_seq_length] = 0;
	bin_hap_se_idx += hap_seq_length;

	if(generate_modify_string){
		kmer_size_alt_string = MAX(kmer_size_alt_string, 22);
		modify_info.clear();
		uint32_t modify_info_size = hap_seq_length;
		for(uint32_t i = 0; i < modify_info_size; i++){
			modify_info.emplace_back(i);
		}
		//reset
		bin_hap_se_idx = 0;
		bin_hap_se_idx += 1;//SV type
		uint32_t alt_base_bg = SV_STORE_BASIC::get_16bit(bin_hap_seq + bin_hap_se_idx);	bin_hap_se_idx += 2;//SV type
		bin_hap_se_idx += 4;//SV length

		uint16_t chr_ID1, chr_ID2;
		uint32_t POS1, POS2;
		SV_STORE_BASIC::get_chrID_POS_pair(bin_hap_seq + bin_hap_se_idx, chr_ID1, POS1); bin_hap_se_idx += 6;//pos of BP 1
		SV_STORE_BASIC::get_chrID_POS_pair(bin_hap_seq + bin_hap_se_idx, chr_ID2, POS2); bin_hap_se_idx += 6;//pos of BP 1
		if(sv_type == SV_TYPE::INS || sv_type == SV_TYPE::DEL){
			uint32_t modified_base_bg = alt_base_bg;
			uint32_t modified_base_ed = alt_base_bg + ((sv_type == SV_TYPE::DEL)?(1):(sv_length + 1));

			//fprintf(stderr, "modify_info.size() %d\n", modify_info.size());
			for(uint32_t mod_kmer_i = modified_base_bg; mod_kmer_i < modified_base_ed; mod_kmer_i++){
				if(mod_kmer_i < modify_info.size()){
					//fprintf(stderr, "%d\n", mod_kmer_i);
					modify_info[mod_kmer_i].is_modified =  true;
				}
			}
			for(uint32_t mod_kmer_i = 1; mod_kmer_i < kmer_size_alt_string; mod_kmer_i++){
				if(alt_base_bg >= mod_kmer_i && alt_base_bg - mod_kmer_i < modify_info.size()){
					//fprintf(stderr, "%d\n", alt_base_bg - mod_kmer_i);
					modify_info[alt_base_bg - mod_kmer_i].is_modified =  true;
				}
			}
		}
	}
	return rst_haplotype;
}

void Window_block_handler::unfold_window_block_all(uint8_t * window_block_p){
	if(window_block_p_old == window_block_p)	return;
	else										window_block_p_old = window_block_p;

	block_all_hap_list.clear();
	uint32_t block_count = (window_block_p[0] << 8) + window_block_p[1];
	int c_hap_id = 0;
	uint8_t * hap_seq_load_p = window_block_p + 2;
	uint64_t hap_offset = 0;
	uint64_t hap_seq_len = 0;
	//get hap_seq_len and hap_offset
	while(true){
		if(hap_seq_load_p - window_block_p >= block_count + 2){ hap_offset = -1; break; }//out of boundary
		if(((hap_seq_load_p[0]) & 0x1) == 0){ 	hap_seq_len = ( hap_seq_load_p[0] >> 1); 															hap_seq_load_p ++;	} //one byte
		else{									hap_seq_len = ((hap_seq_load_p[0] >> 1) << 16) + (hap_seq_load_p[1] << 8) + (hap_seq_load_p[2]);	hap_seq_load_p += 3;} //three bytes
		block_all_hap_list.emplace_back();
		block_all_hap_list.back().p = window_block_p + 2 + block_count + hap_offset;
		block_all_hap_list.back().len = hap_seq_len;
		hap_offset += hap_seq_len;
		c_hap_id++;
	}
}
