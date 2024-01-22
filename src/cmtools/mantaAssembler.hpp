/*
 * assembly.hpp
 *
 *  Created on: 2020年4月24日
 *      Author: fenghe
 */

#pragma once

#include <iosfwd>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <cassert>
#include <cstring>
#include <algorithm>
#include <unordered_map>

typedef int32_t pos_t;

struct AssembledContig {
  std::string  seq;  ///< contigsequence
  unsigned seedReadCount = 0;  ///< no of reads containing the seeding kmer
  std::set<unsigned> supportReads;
  std::set<unsigned> rejectReads;
  pos_t conservativeRange_bgn = 0;///< subsection of the contig with conservative coverage
  pos_t conservativeRange_end = 0;

  //friend std::ostream& operator<<(std::ostream& os, const AssembledContig& contig)
  //{
  //  os << "CONTIG size: " << contig.seq.size() << " seedCount: " << contig.seedReadCount
  //     << " supportReads: " << contig.supportReads.size() << " seq:\n";
  //  {
  //  	const char* seq = contig.seq.c_str();
  //  	static const unsigned rowSize(100);
  //  	static const unsigned sectionSize(10);
  //  	assert(nullptr != seq);
  //  	const unsigned seqLen(strlen(seq));
  //  	for (unsigned i(0); i < seqLen; ++i) {
  //  		//if (i) {
  //  		//  if (0 == (i % rowSize))
  //  		//	os << '\n';
  //  		//  else if (0 == (i % sectionSize))
  //  		//	os << ' ';
  //  		//}
  //  		os << seq[i];
  //  	}
  //  }
  //  os << "\n";
  //  return os;
  //}

	//void debug_print(std::ostream &os) const{
	//	os << "CONTIG size: " << seq.size() << " seedCount: " << seedReadCount
	//			<< " supportReads: " << supportReads.size() << " seq:\n";
	//	{
	//		const char *seq_ = seq.c_str();
	//		assert(nullptr != seq_);
	//		const unsigned seqLen(strlen(seq_));
	//		for (unsigned i(0); i < seqLen; ++i) {
	//			os << seq_[i];
	//		}
	//	}
    //    os << "\n";
	//}

    void debug_print(std::ostream &os) const{
		os << "CONTIG size: " << seq.size() << " seedCount: " << seedReadCount
				<< " supportReads: " << supportReads.size() << " seq:\n";
		//{
		//	const char *seq_ = seq.c_str();
		//	static const unsigned rowSize(100);
		//	static const unsigned sectionSize(10);
		//	assert(nullptr != seq_);
		//	const unsigned seqLen(strlen(seq_));
		//	for (unsigned i(0); i < seqLen; ++i) {
		//		if (i) {
		//			if (0 == (i % rowSize))
		//				os << '\n';
		//			else if (0 == (i % sectionSize))
		//				os << ' ';
		//		}
		//		os << seq_[i];
		//	}
		//}
		os << "\n";
		os << "supportReads: ";
		for(auto & r: supportReads)
			os << r << "\t";
		os << "\n";

		os << "rejectReads: ";
		for(auto & r: rejectReads)
			os << r << "\t";
		os << "\n";
	}

    const char *getSeq() 
    {
        return seq.c_str();
    }

    unsigned getSupportReadCnt()
    {
        return supportReads.size();
    }
};

/**************************************************************/
/// Information added to each read in the process of assembly
struct AssemblyReadInfo {
	AssemblyReadInfo(bool isPseudo_ = false):isPseudo(isPseudo_){}
	bool isUsed = false;
	/// If true, the read was 'used' but filtered out, so there is no meaningful contig id association
	bool isFiltered = false;
	/// If true, the read was an assembled contig
	bool isPseudo = false;
	/// Index of the contigs that this read is used in
	std::vector<unsigned> contigIds;
};

typedef std::vector<AssembledContig> Assembly;
typedef std::vector<std::string>      AssemblyReadInput;
typedef std::vector<AssemblyReadInfo> AssemblyReadOutput;

/// Input parameters for IterativeAssembler
///
struct IterativeAssemblerOptions {
  IterativeAssemblerOptions() {}

  /// the symbol set used during assembly
  std::string alphabet = "ACGT";
  /// minimum basecall quality for assembly input
  int minQval = 5;
  /// initial word (kmer) length
  unsigned minWordLength = 41; //41;
  unsigned maxWordLength   =  76; //76;
  unsigned wordStepSize    = 5;
  unsigned minContigLength = 15;
  /// min. coverage required for contig extension
  unsigned minCoverage = 2; //1;
  /// coverage required for conservative contig sub-range
  unsigned minConservativeCoverage = 2;
  /// max error rates allowed during contig extension
  double maxError = 0.35; //0.35;
  /// min. number of unused reads to enable search for more contigs
  unsigned minUnusedReads = 3; //3
  /// min. number of reads required to start assembly
  unsigned minSupportReads = 2; //2
  /// Max. number of assembly returned for a given set of reads
  unsigned maxAssemblyCount = 1;
};

void runIterativeAssembler(
    const IterativeAssemblerOptions& opt,
    AssemblyReadInput&               reads,
    AssemblyReadOutput&              assembledReadInfo,
    Assembly&                        contigs);

typedef std::unordered_map<std::string, unsigned> str_uint_map_t;
// maps kmers to supporting reads
typedef std::unordered_map<std::string, std::set<unsigned>> str_set_uint_map_t;
typedef std::unordered_map<std::string, std::pair<unsigned, unsigned>> str_pair_uint_map_t;

struct AssemblyManager{

	std::vector<std::string>      reads;
	void addRead(const char * seq){ reads.push_back(std::string(seq));}
	void clear(){  reads.clear(); readInfo.clear(); }
	std::vector<AssembledContig> & getResults(){return contigs;}
	void assembley();

	IterativeAssemblerOptions o;
private:
	//************************get kmer count***************************************/
	void getKmerCounts(const unsigned wordLength, str_uint_map_t &wordCount, str_set_uint_map_t &wordSupportReads);
	//buffs for getKmerCounts
	std::set<std::string> getKmerCounts_readWords_BUFF;
	//****************************getRepeatKmers*******************************************/
	void getRepeatKmers(const str_uint_map_t &wordCount, std::set<std::string> &repeatWords);
	//buffs
	str_pair_uint_map_t getRepeatKmers_wordIndices;
	std::vector<std::string> getRepeatKmers_wordStack;

	//**************************buildContigs********************************/
	bool buildContigs(const unsigned wordLength);
	//buffs
	str_uint_map_t buildContigs_wordCount;
	// records the supporting reads for each kmer
	str_set_uint_map_t buildContigs_wordSupportReads;
	// identify repeat kmers (i.e. circles from the de bruijn graph)
	std::set<std::string> buildContigs_repeatWords;
	std::set<std::string> buildContigs_unusedWords;
	/************************selectContigs*******************************/
	void selectContigs(const unsigned normalReadCount);
	//buffs
	std::set<unsigned> selectContigs_usedReads;
	std::set<unsigned> selectContigs_usedPseudoReads;
	std::set<unsigned> selectContigs_contigs2Remove;
	std::set<unsigned> selectContigs_newSupportReads;

	//IterativeAssemblerOptions o;
	std::vector<AssemblyReadInfo> readInfo;
	std::vector<AssembledContig>  contigs;
	//buffs
	std::vector<AssembledContig> tmpContigs;

	/***********************************Walk**************************************/
	bool walk(const std::string &seed, const unsigned wordLength, const str_uint_map_t &wordCount,
			const str_set_uint_map_t &wordReads, const std::set<std::string> &repeatWords,
			std::set<std::string> &unusedWords, AssembledContig &contig);
	//buffs for walk
	std::set<std::string> walk_wordsInContig;
	std::set<unsigned> walk_maxWordReads;
	std::set<unsigned> walk_maxContigWordReads;
	std::set<unsigned> walk_previousWordReads;
	std::set<unsigned> walk_supportReads2Remove;
	std::set<unsigned> walk_rejectReads2Add;
	std::set<unsigned> walk_contigWordReads;
	std::set<unsigned> walk_sharedReads;
	std::set<unsigned> walk_toRemove;
	std::set<unsigned> walk_toAdd;
	std::set<unsigned> walk_sharedReads_alleles;
	std::set<unsigned> walk_toUpdate;
};

// std :: vector <CON_suppR> menta_uesd(std :: vector <string> stringV, int maxAssemblyC);
// void mantaAss_m(std :: map <int,std :: vector <readstr_name>> clupos_str, std :: vector <ConsensusINFO> &consensus_str);
//demo:
/**
 * void main()
 * {
 * 		AssemblyManager am;
 * 		while(1)
 * 		{
 * 			if(END_of_all_ass_event)
 * 				break;
 * 			am.clear();
 * 			for(read in read_list)
 * 				am.addRead(read);
 * 			am.assembley();
 * 			auto &contigs = am.getResults();
 * 			for(auto &contig:contigs)
 * 				std::cout << contig;
 * 		}
 * }
 *
 */
