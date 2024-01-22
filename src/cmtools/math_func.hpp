/*
 * math_func.hpp
 *
 *  Created on: 2023年02月15日
 *      Author: miaocui
 */

#ifndef MATH_FUNC_HPP_
#define MATH_FUNC_HPP_

#include <vector>
#include <cmath>
#include <regex>
#include <numeric>

using namespace std;

struct FormatInfo
{
    int GQ;
    int QUAL;
    string GT;
};

void read2higher_stand(string seq, std::string & sequence);
void rescale_read_counts(int &n_alts, int &n_total);
float log10sumexp(std :: vector<float> log10_probs);
std :: vector<float> normalize_log10_probs(std :: vector<float> log10_probs);
std :: vector<float> toRealSpace(std :: vector<float> log10_probs);
vector<string> split(const string &str, const string &pattern);
void divide(string s,vector<string> &v,regex reg);
string int2base(int num);
int base2int(std :: string base_);
// int chrID2int(string chrID_str);

std :: vector<float> cal_GL(int n_alts, int n_total, float err);
FormatInfo computer_quals(std :: vector<float> GL_P);

//fenghe
#define MAX_LINE_LENGTH 10000000
#define MAX_LINE_ITEM_NUM 2000000

void split_string_append(std::vector<std::string> &item_value, char * temp, const char * split_line, const char *split_str);
void split_string(std::vector<std::string> &item_value, char * temp, const char * split_line, const char *split_str);
vector<string> split2(const string &str, const string &pattern);

void load_strings_from_file(char * string_fn, std::vector<std::string> &v);
void load_int_from_file(char * int_fn, std::vector<int> &v);

int string_bias(int line0, int line16);

#endif /* MATH_FUNC_HPP_ */
