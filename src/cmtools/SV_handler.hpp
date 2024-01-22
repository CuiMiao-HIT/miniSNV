#ifndef SV_HANDLER_HPP_
#define SV_HANDLER_HPP_

#include "haplotype.hpp"

//------------------------------MAIN--------------------------------------//
#define DUP_MAX_LEN 2000

//'ALL' 'DEL' 'INS' 'DUP' 'TRA' 'INV' 'BND'");
struct SV_TYPE{
	static const int ALL = 0;
	static const int DEL = 1;
	static const int INS = 2;
	static const int DUP = 3;
	static const int TRA = 4;
	static const int INV = 5;
	static const int BND = 6;
	static const int UNKNOWN = 7;

	static int get_sv_type_int(char * sv_type_str){
			 if(strcmp(sv_type_str, "DEL") == 0) return SV_TYPE::DEL;
		else if(strcmp(sv_type_str, "INS") == 0) return SV_TYPE::INS;
		else if(strcmp(sv_type_str, "DUP") == 0) return SV_TYPE::DUP;
		else if(strcmp(sv_type_str, "TRA") == 0) return SV_TYPE::TRA;
		else if(strcmp(sv_type_str, "INV") == 0) return SV_TYPE::INV;
		else if(strcmp(sv_type_str, "BND") == 0) return SV_TYPE::BND;
		else if(strcmp(sv_type_str, "ALL") == 0) return SV_TYPE::ALL;
		else if(strcmp(sv_type_str, "all") == 0) return SV_TYPE::ALL;
		return SV_TYPE::UNKNOWN;
	}

	static uint8_t type_2_flag(int SV_type){
		uint8_t SV_flag_number = 0;
		if(SV_type == SV_TYPE::INV){	SV_flag_number = 59;	}
		if(SV_type == SV_TYPE::DUP){	SV_flag_number = 60;	}
		if(SV_type == SV_TYPE::TRA){	SV_flag_number = 61;	}
		if(SV_type == SV_TYPE::INS){	SV_flag_number = 62;	}
		if(SV_type == SV_TYPE::DEL){	SV_flag_number = 63;	}
		return (SV_flag_number << 2) + 2;
	}

	static int flag_2_type(int SV_flag){
		SV_flag = (SV_flag >> 2) & 0x3f;
		if(SV_flag == 59) return SV_TYPE::INV;
		if(SV_flag == 60) return SV_TYPE::DUP;
		if(SV_flag == 61) return SV_TYPE::TRA;
		if(SV_flag == 62) return SV_TYPE::INS;
		if(SV_flag == 63) return SV_TYPE::DEL;

		return SV_TYPE::UNKNOWN;
	}


};

struct SV_STORE_BASIC{

	static void store_16bit(std::vector<uint8_t> &hap_seq, uint16_t val){
		hap_seq.emplace_back((val >> 8) & 0xffff);
		hap_seq.emplace_back((val) & 0xff);
	}

	static void store_32bit(std::vector<uint8_t> &hap_seq, uint32_t val){
		hap_seq.emplace_back(val >> 24);//store the SV length, 24 bit
		hap_seq.emplace_back(val >> 16 & 0xffffff);//store the SV length, 24 bit
		hap_seq.emplace_back((val >> 8) & 0xffff);
		hap_seq.emplace_back((val) & 0xff);
	}

	static uint16_t get_16bit(uint8_t * hap_seq){
		return ((uint16_t)hap_seq[0] << 8) + ((uint16_t)hap_seq[1]);
	}

	static uint32_t get_32bit(uint8_t * hap_seq){
		return (((uint32_t)hap_seq[0]) << 24) + (((uint32_t)hap_seq[1]) << 16) + (((uint32_t)hap_seq[2]) << 8)  + (((uint32_t)hap_seq[3]));
	}

	static void get_chrID_POS_pair(uint8_t * hap_seq, uint16_t &chrID, uint32_t &POS){
		chrID = get_16bit(hap_seq);
		POS = get_32bit(hap_seq + 2);
	}

	static void store_chrID_POS_pair(std::vector<uint8_t> &hap_seq, uint16_t chrID, uint32_t POS){
		store_16bit(hap_seq, chrID);
		store_32bit(hap_seq, POS);
	}

};

#endif /* SV_HANDLER_HPP_ */
