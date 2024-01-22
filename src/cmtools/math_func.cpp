/*
 * math_func.cpp
 *
 *  Created on: 2023年02月15日
 *      Author: miaocui
 */

#include <string>
#include <stdio.h>
#include <string.h>
#include <fstream>

extern "C"{
	#include "../clib/utils.h"
}

#include "math_func.hpp"

#define LOG_10 log(10.0);

void read2higher_stand(string seq, std::string & sequence){
    reverse(seq.begin(), seq.end());
    for (int i = 0; i < seq.length(); i++){
        switch (seq[i]) {
        case 'A': case 'a': seq[i] = 'T'; break;
        case 'C': case 'c': seq[i] = 'G'; break;
        case 'G': case 'g': seq[i] = 'C'; break;
        case 'T': case 't': seq[i] = 'A'; break;
        default:   seq[i] = 'N'; break;
        }
    }
    sequence = seq.substr(0);    
}

void rescale_read_counts(int &n_alts, int &n_total){

    int max_allowed_reads=100;
    if(n_total > max_allowed_reads){
        float ratio = (float)n_alts / n_total;
        n_alts = int(floor(ratio * max_allowed_reads));
        n_total = max_allowed_reads;
    }
}

float log10sumexp(std :: vector<float> log10_probs){
    float maxvalue = *max_element(log10_probs.begin(), log10_probs.end());
    float sumexp = 0;
    for(auto a : log10_probs){
        sumexp = sumexp + pow(10.0, a-maxvalue);
    }
    sumexp = log10(sumexp);
    return maxvalue + sumexp;
}

std :: vector<float> normalize_log10_probs(std :: vector<float> log10_probs){
    float maxvalue = *max_element(log10_probs.begin(), log10_probs.end());
    if (maxvalue > 0)
    {
        throw std::invalid_argument( "log10_probs all must be <=0" );
    }
    
    float lse = log10sumexp(log10_probs);
    std :: vector<float> result;
    float zero = 0.0;
    for(auto a : log10_probs){
        float m = min(a-lse,zero);
        result.push_back(m);
    }
    return result;
}

std :: vector<float> toRealSpace(std :: vector<float> log10_probs){
    std :: vector<float> result;
    for(auto a : log10_probs){
        float cl = pow(10.0, a);
        result.push_back(cl);
    }
    return result;
}

std :: vector<float> cal_GL(int n_alts, int n_total, float err){
    int n_ref = n_total - n_alts;
    float logp = log(err) / LOG_10;
    float log1p = log(1-err) / LOG_10;
    float ori_GL00 = n_ref * log1p + n_alts * logp;
    float ori_GL01 = (-n_total) * log(2) / LOG_10;
    float ori_GL11 = n_ref *logp + n_alts * log1p;
    std :: vector<float> log10_;
    log10_.push_back(ori_GL00);
    log10_.push_back(ori_GL01);
    log10_.push_back(ori_GL11);
    std :: vector<float> log10_probs = normalize_log10_probs(log10_);
    return toRealSpace(log10_probs);
}

FormatInfo computer_quals(std :: vector<float> GL_P){
    FormatInfo FMinfo;
    std :: string Genotype[6] = {"0/0", "0/1", "1/1", "0/2", "1/2", "2/2"};
    int index_of_max = max_element(GL_P.begin(),GL_P.end())-GL_P.begin();
    float sum = accumulate(GL_P.begin(),GL_P.end(),0.0) - GL_P[index_of_max];
    if (sum < pow(10,-38))
    {
        sum = pow(10,-38);
    }
    int GQ = (int)((-10)*log10(sum));
    float clp0 = GL_P[0];
    if (clp0 < pow(10,-38))
    {
        clp0 = pow(10,-38);
    }
    float QUAL = abs(round((-10)*log10(clp0)));//??????????????????????????????????????
    FMinfo.GQ = GQ;
    FMinfo.QUAL = QUAL;
    FMinfo.GT = Genotype[index_of_max];

    return FMinfo;
}

//some tool function-----------------------------------------------------------------------------------------------------------------------------------

vector<string> split(const string &str, const string &pattern)
{
    vector<string> res;
    if(str == "")
        return res;
    //在字符串末尾也加入分隔符，方便截取最后一段
    string strs = str + pattern;
    size_t pos = strs.find(pattern);

    while(pos != strs.npos)
    {
        string temp = strs.substr(0, pos);
        res.push_back(temp);
        //去掉已分割的字符串,在剩下的字符串中进行分割
        strs = strs.substr(pos+1, strs.size());
        pos = strs.find(pattern);
    }

    return res;
}

void divide(string s,vector<string> &v,regex reg)
{
	smatch mResult;
	if (regex_search(s, mResult, reg)) {
        string pushs = "";
		for (auto it = mResult.prefix().first; it != mResult.prefix().second; it++) {
			if (*it == ',')
				continue;
			// v.push_back(string(1, *it));
            pushs = pushs + *it;
		}
		v.push_back(pushs);
		divide(mResult.suffix(), v, reg);
	}
	else {
        v.push_back(s);
	}
}

// eg : chr1->0
// int chrID2int(string chrID_str){ 
//     string ID = chrID_str.substr(3);
//     int chrID = stoi(ID)-1;
//     return chrID;
// }

int base2int(std :: string base_){ 
    int b2i;
    char base_char[1];
    base_.copy(base_char,1,0);
    switch (base_char[0])
    {
    case 'a':
    case 'A':
        b2i = 0;
        break;
    case 'c':
    case 'C':
        b2i = 1;
        break;
    case 'g':
    case 'G':
        b2i = 2;
        break;
    case 't':
    case 'T':
        b2i = 3;
        break;
    default:
        break;
    }
    return b2i;
}

string int2base(int num){ 
    string i2b;
    switch (num)
    {
    case 0:
        i2b = "A";
        break;
    case 1:
        i2b = "C";
        break;
    case 2:
        i2b = "G";
        break;
    case 3:
        i2b = "T";
        break;
    default:
        break;
    }
    return i2b;
}

void split_string_append(std::vector<std::string> &item_value, char * temp, const char * split_line, const char *split_str){
	if(strlen(split_line) == 0) return;
	strcpy(temp, split_line);
	char * token_value = NULL;
	char *save_ptr = NULL;
	token_value = strtok_r(temp, split_str, &save_ptr); item_value.emplace_back(token_value);
	for(int item_idx = 1; item_idx < MAX_LINE_ITEM_NUM; item_idx++){
		token_value = strtok_r(NULL, split_str, &save_ptr);
		if(token_value == NULL)
			break;
		item_value.emplace_back(token_value);
	}
}

//fenghe

void split_string(std::vector<std::string> &item_value, char * temp, const char * split_line, const char *split_str){
	item_value.clear();
	split_string_append(item_value, temp, split_line, split_str);
}

vector<string> split2(const string &str, const string &pattern)
{
    char * strc = new char[strlen(str.c_str())+1];
    strcpy(strc, str.c_str());   //string转换成C-string
    vector<string> res;
    char* temp = strtok(strc, pattern.c_str());
    while(temp != NULL)
    {
        res.push_back(string(temp));
        temp = strtok(NULL, pattern.c_str());
    }
    delete[] strc;
    return res;
}

void load_strings_from_file(char * string_fn, std::vector<std::string> &v){
	//get file names
	FILE * try_f = xopen(string_fn, "rb");
	fclose(try_f);
	char *temp = new char[MAX_LINE_LENGTH];//10M
	std::ifstream name_list_File(string_fn);
	while(true){
		name_list_File.getline(temp, MAX_LINE_LENGTH);
		if(name_list_File.eof() && (temp[0]) == 0) break;
		v.emplace_back(temp);
		if(name_list_File.eof())	break;
	}
	name_list_File.close();
}

void load_int_from_file(char * int_fn, std::vector<int> &v){
	FILE * try_f = xopen(int_fn, "rb");
	fclose(try_f);
	//get file names
	char *temp = new char[MAX_LINE_LENGTH];//10M
	std::ifstream name_list_File(int_fn);
	while(true){
		name_list_File.getline(temp, MAX_LINE_LENGTH);
		if(*temp == 0)	break;
		v.emplace_back(atoi(temp));
	}
	name_list_File.close();
}

int string_bias(int line0, int line16){//1:delate 0:use
    int del = 1;
    float bias;
    int all = line0 + line16;
    int diff = abs(line0 - line16);
    bias = (float)diff/all;
    if(bias < 0.75) del = 0;
    return del;
}