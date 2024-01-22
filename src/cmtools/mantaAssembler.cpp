//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Xiaoyu Chen
/// \author Ole Schulz-Trieglaff
///

#include "mantaAssembler.hpp"

#define PY_SSIZE_T_CLEAN

#define MAXSEQLEN 4096 

//using namespace std;


// compile with this macro to get verbose output:
//#define DEBUG_ASBL
//#define DEBUG_WALK

/// Adds base @p base to the end (isEnd is true) or start (otherwise) of the contig.
///
/// \return the extended contig.
static std::string addBase(const std::string &contig, const char base,
		const bool isEnd) {
	if (isEnd)
		return contig + base;
	else
		return base + contig;
}

/// Returns a suffix (isEnd is true) or prefix (otherwise) of the input contig with the specified length.
static std::string getEnd(const std::string &contig, const unsigned length,
		const bool isEnd) {
	const unsigned csize(contig.size());
	assert(length <= csize);

	if (isEnd)
		return contig.substr((csize - length), length);
	else
		return contig.substr(0, length);
}

/// Construct a contig from the seed provided
///
/// The contig is extened in both directions until
/// either the there is no sufficient support evidence
/// or the last k-mer is repeatitive (i.e. part of a bubble in the graph)
///
/// \return True if the contig runs into a repeatitive k-mer when extending in either mode
bool AssemblyManager::walk(const std::string &seed,
		const unsigned wordLength, const str_uint_map_t &wordCount,
		const str_set_uint_map_t &wordReads,
		const std::set<std::string> &repeatWords,
		std::set<std::string> &unusedWords, AssembledContig &contig) {

	//load buffs
	std::set<std::string> &wordsInContig = walk_wordsInContig;  wordsInContig.clear();
	std::set<unsigned> &maxWordReads = walk_maxWordReads;
	std::set<unsigned> &maxContigWordReads = walk_maxContigWordReads;
	std::set<unsigned> &previousWordReads = walk_previousWordReads;
	std::set<unsigned> &supportReads2Remove = walk_supportReads2Remove;
	std::set<unsigned> &rejectReads2Add = walk_rejectReads2Add;
	std::set<unsigned> &contigWordReads = walk_contigWordReads;
	std::set<unsigned> &sharedReads = walk_sharedReads;
	std::set<unsigned> &toRemove = walk_toRemove;
	std::set<unsigned> &toAdd = walk_toAdd;
	std::set<unsigned> &sharedReads_alleles = walk_sharedReads_alleles;
	std::set<unsigned> &toUpdate = walk_toUpdate;

	const str_uint_map_t::const_iterator wordCountEnd(wordCount.cend());
	const str_set_uint_map_t::const_iterator wordReadsEnd(wordReads.cend());

	// we start with the seed
	str_set_uint_map_t::const_iterator wordReadsIter(wordReads.find(seed));
	assert(wordReadsIter != wordReadsEnd);
	contig.supportReads = wordReadsIter->second;
	contig.seq = seed;
	unusedWords.erase(seed);

	if (repeatWords.find(seed) != repeatWords.end()) {
		contig.conservativeRange_bgn = (0);
		contig.conservativeRange_end = (wordLength);
		return true;
	}

	// collecting words used to build the contig
	wordsInContig.insert(seed);

	const std::string tmpTrunk = getEnd(seed, wordLength - 1, false);
	// collecting rejecting reads for the seed from the unselected branches
	for (const char symbol : o.alphabet) {
		// the seed itself
		if (symbol == seed[wordLength - 1])
			continue;

		// add rejecting reads from an unselected word/branch
		const std::string newKey(addBase(tmpTrunk, symbol, true));

		wordReadsIter = wordReads.find(newKey);
		if (wordReadsIter == wordReadsEnd)
			continue;
		const std::set<unsigned> &unselectedReads(wordReadsIter->second);

		contig.rejectReads.insert(unselectedReads.begin(),
				unselectedReads.end());
	}

	bool isRepeatFound(false);

	// 0 => walk to the right, 1 => walk to the left
	for (unsigned mode(0); mode < 2; ++mode) {
		const bool isEnd(mode == 0);
		unsigned conservativeEndOffset(0);

		while (true) {
			const std::string previousWord = getEnd(contig.seq, wordLength,
					isEnd);
			const std::string trunk(getEnd(contig.seq, wordLength - 1, isEnd));
			unsigned maxBaseCount(0);
			unsigned maxContigWordReadCount(0);
			char maxBase(o.alphabet[0]);
			std::string maxWord;
			maxWordReads.clear();
			maxContigWordReads.clear();
			previousWordReads.clear();;
			supportReads2Remove.clear();;
			rejectReads2Add.clear();

			for (const char symbol : o.alphabet) {
				const std::string newKey(addBase(trunk, symbol, isEnd));
				const str_uint_map_t::const_iterator wordCountIter(
						wordCount.find(newKey));
				if (wordCountIter == wordCountEnd)
					continue;
				const unsigned currWordCount(wordCountIter->second);

				wordReadsIter = wordReads.find(newKey);
				if (wordReadsIter == wordReadsEnd)
					continue;
				const std::set<unsigned> &currWordReads(wordReadsIter->second);

				// get the shared supporting reads between the contig and the current word
				contigWordReads.clear();
				std::set_intersection(contig.supportReads.begin(),
						contig.supportReads.end(), currWordReads.begin(),
						currWordReads.end(),
						std::inserter(contigWordReads,
								contigWordReads.begin()));

				// get the shared supporting reads across two alleles
				sharedReads.clear();
				std::set_intersection(maxContigWordReads.begin(),
						maxContigWordReads.end(), currWordReads.begin(),
						currWordReads.end(),
						std::inserter(sharedReads, sharedReads.begin()));

				if (contigWordReads.empty())
					continue;

				const unsigned contigWordReadCount(contigWordReads.size());
				if (contigWordReadCount > maxContigWordReadCount) {
					// the old shared reads support an unselected allele if they don't support the new word
					// remove them from the contig's supporting reads
					if (!maxContigWordReads.empty()) {
						toRemove.clear();
						std::set_difference(maxContigWordReads.begin(),
								maxContigWordReads.end(), sharedReads.begin(),
								sharedReads.end(),
								std::inserter(toRemove, toRemove.begin()));
						if (!toRemove.empty())
							supportReads2Remove.insert(toRemove.begin(),
									toRemove.end());
					}

					// the old supporting reads is for an unselected allele if they don't support the new word
					// they become rejecting reads for the currently selected allele
					if (!maxWordReads.empty()) {
						toAdd.clear();
						std::set_difference(maxWordReads.begin(),
								maxWordReads.end(), sharedReads.begin(),
								sharedReads.end(),
								std::inserter(toAdd, toAdd.begin()));
						if (!toAdd.empty())
							rejectReads2Add.insert(toAdd.begin(), toAdd.end());
					}
					// new supporting reads for the currently selected allele
					maxWordReads = currWordReads;

					maxContigWordReadCount = contigWordReadCount;
					maxContigWordReads = contigWordReads;
					maxBaseCount = currWordCount;
					maxBase = symbol;
					maxWord = newKey;
				} else {
					toRemove.clear();
					std::set_difference(contigWordReads.begin(),
							contigWordReads.end(), sharedReads.begin(),
							sharedReads.end(),
							std::inserter(toRemove, toRemove.begin()));
					if (!toRemove.empty())
						supportReads2Remove.insert(toRemove.begin(),
								toRemove.end());

					toAdd.clear();
					std::set_difference(currWordReads.begin(),
							currWordReads.end(), sharedReads.begin(),
							sharedReads.end(),
							std::inserter(toAdd, toAdd.begin()));
					if (!toAdd.empty())
						rejectReads2Add.insert(toAdd.begin(), toAdd.end());
				}
			}

			if (maxBaseCount < o.minCoverage) {

				break;
			}

			// stop walk in the current mode after seeing one repeat word
			if (wordsInContig.find(maxWord) != wordsInContig.end()) {
				isRepeatFound = true;
				break;
			}

			contig.seq = addBase(contig.seq, maxBase, isEnd);

			if ((conservativeEndOffset != 0)
					|| (maxBaseCount < o.minConservativeCoverage))
				conservativeEndOffset += 1;

			// TODO: can add threshold for the count or percentage of shared reads
			{
				// walk backwards for one step at a branching point
				if (maxWordReads != previousWordReads) {
					const char tmpSymbol = (
							isEnd ? previousWord[0] : previousWord[wordLength
											- 1]);
					for (const char symbol : o.alphabet) {
						// the selected branch: skip the backward word itself
						if (symbol == tmpSymbol)
							continue;

						// add rejecting reads from an unselected branch
						const std::string newKey(
								addBase(trunk, symbol, !isEnd));
						// the selected branch: skip the word just extended
						if (newKey == maxWord)
							continue;

						wordReadsIter = wordReads.find(newKey);
						if (wordReadsIter == wordReadsEnd)
							continue;

						const std::set<unsigned> &backWordReads(
								wordReadsIter->second);
						// get the shared supporting reads across two alleles
						sharedReads_alleles.clear();
						std::set_intersection(maxContigWordReads.begin(),
								maxContigWordReads.end(), backWordReads.begin(),
								backWordReads.end(),
								std::inserter(sharedReads_alleles,
										sharedReads_alleles.begin()));

						toUpdate.clear();
						std::set_difference(backWordReads.begin(),
								backWordReads.end(), sharedReads_alleles.begin(),
								sharedReads_alleles.end(),
								std::inserter(toUpdate, toUpdate.begin()));
						if (!toUpdate.empty()) {
							rejectReads2Add.insert(toUpdate.begin(),
									toUpdate.end());
							supportReads2Remove.insert(toUpdate.begin(),
									toUpdate.end());
						}
					}
				}
				previousWordReads = maxWordReads;
				// update rejecting reads
				// add reads that support the unselected allele
				for (const unsigned rd : rejectReads2Add) {
					contig.rejectReads.insert(rd);
				}

				// update supporting reads
				// add reads that support the selected allel
				for (const unsigned rd : maxWordReads) {
					if (contig.rejectReads.find(rd) == contig.rejectReads.end())
						contig.supportReads.insert(rd);
				}
				// remove reads that do NOT support the selected allel anymore
				for (const unsigned rd : supportReads2Remove) {
					contig.supportReads.erase(rd);
				}
			}

			// remove the last word from the unused list, so it cannot be used as the seed in finding the next
			// contig
			unusedWords.erase(maxWord);
			// collect the words used to build the contig
			wordsInContig.insert(maxWord);
		}

		// set conservative coverage range for the contig
		if (mode == 0)
			contig.conservativeRange_end = (conservativeEndOffset);
		else
			contig.conservativeRange_bgn = (conservativeEndOffset);
	}

	contig.conservativeRange_end = ( contig.seq.size() - contig.conservativeRange_end);

	return isRepeatFound;
}

/// Construct k-mer maps
/// k-mer ==> number of reads containing the k-mer
/// k-mer ==> a list of read IDs containg the k-mer
void AssemblyManager::getKmerCounts( const unsigned wordLength, str_uint_map_t &wordCount, str_set_uint_map_t &wordSupportReads) {
	const unsigned readCount(reads.size());
	std::set<std::string> &readWords = getKmerCounts_readWords_BUFF;
	for (unsigned readIndex(0); readIndex < readCount; ++readIndex) {
		// stores the index of a kmer in a read sequence
		const std::string &seq(reads[readIndex]);
		const unsigned readLen(seq.size());
		// this read is unusable for assembly:
		if (readLen < wordLength)
			continue;
		// track all words from the read, including repetitive words
		readWords.clear();
		for (unsigned j(0); j <= (readLen - wordLength); ++j) {
			const std::string word(seq.substr(j, wordLength));
			// filter words with "N" (either directly from input alignment or marked due to low basecall quality:
			if (word.find('N') != std::string::npos) continue;
			readWords.insert(word);
		}

		unsigned wordCountAdd = 1;
		// pseudo reads must have passed coverage check with smaller kmers
		// Assigning minCoverage (instead of 1) to a pseudo read allows the pseudo read to rescue the regions
		// where coverage is as low as minCoverage and reads overlap is small
		if (readInfo[readIndex].isPseudo)
			wordCountAdd = o.minCoverage;
		// total occurrences from this read
		for (const std::string &word : readWords) {
			wordCount[word] += wordCountAdd;
			// record the supporting read
			wordSupportReads[word].insert(readIndex);
		}
	}
}

/// Identify repetitive k-mers
/// i.e. k-mers that form a circular subgraph
///
static unsigned searchRepeats(const std::string &alphabet, const unsigned index, const std::string &word,
		str_pair_uint_map_t &wordIndices, std::vector<std::string> &wordStack,
		std::set<std::string> &repeatWords) {
	// set the depth index for the current word to the smallest unused index
	wordIndices[word] = std::pair<unsigned, unsigned>(index, index);
	unsigned nextIndex = index + 1;
	wordStack.push_back(word);

	const std::string tmp(getEnd(word, word.size() - 1, true));
	for (const char symbol : alphabet) {
		// candidate successor of the current word
		const std::string nextWord(addBase(tmp, symbol, true));

		// homopolymer
		if (word == nextWord) {
			repeatWords.insert(word);
			continue;
		}

		// the successor word does not exist in the reads
		if (wordIndices.count(nextWord) == 0)
			continue;

		const unsigned nextWordIdx = wordIndices[nextWord].first;
		if (nextWordIdx == 0) {
			// the successor word has not been visited
			// recurse on it
			nextIndex = searchRepeats(alphabet, nextIndex, nextWord, wordIndices,
					wordStack, repeatWords);
			// update the current word's lowlink
			const unsigned wordLowLink = wordIndices[word].second;
			const unsigned nextWordLowLink = wordIndices[nextWord].second;
			wordIndices[word].second = std::min(wordLowLink, nextWordLowLink);
		} else {
			const bool isContained(
					std::find(wordStack.begin(), wordStack.end(), nextWord)
							!= wordStack.end());
			if (isContained) {
				// the successor word is in stack and therefore in the current circle of words
				// only update the current word's lowlink
				const unsigned wordLowLink = wordIndices[word].second;
				wordIndices[word].second = std::min(wordLowLink, nextWordIdx);
			}
		}
	}

	// if the current word is a root node,
	if (wordIndices[word].second == index) {
		const std::string &lastWord(wordStack.back());
		// exclude singletons
		bool isSingleton(lastWord == word);
		if (isSingleton) {
			wordStack.pop_back();
		} else {
			// record identified repeat words (i.e. words in the current circle) if the circle is small

			const unsigned lastWordIndex(wordIndices[lastWord].first);
			const bool isSmallCircle((lastWordIndex - index) <= 50);
			while (true) {
				const std::string repeatWd = wordStack.back();
				if (isSmallCircle)
					repeatWords.insert(repeatWd);
				wordStack.pop_back();

				if (repeatWd == word)
					break;
			}
		}
	}

	return nextIndex;
}

void AssemblyManager::getRepeatKmers(const str_uint_map_t &wordCount, std::set<std::string> &repeatWords) {
	str_pair_uint_map_t &wordIndices = getRepeatKmers_wordIndices;
	wordIndices.clear();
	for (const str_uint_map_t::value_type &wdct : wordCount) {
		wordIndices[wdct.first] = std::pair<unsigned, unsigned>(0, 0);
	}
	unsigned index = 1;
	std::vector<std::string> &wordStack = getRepeatKmers_wordStack;
	wordStack.clear();
	for (const str_pair_uint_map_t::value_type &wdidx : wordIndices) {
		const std::string word = wdidx.first;
		const unsigned wordIdx = wdidx.second.first;
		if (wordIdx == 0)
			index = searchRepeats(o.alphabet, index, word, wordIndices, wordStack,
					repeatWords);
	}
}

bool AssemblyManager::buildContigs(const unsigned wordLength) {
	tmpContigs.clear();
	bool isAssemblySuccess(true);
	// counts the number of occurrences for each kmer in all reads
	str_uint_map_t &wordCount = buildContigs_wordCount; wordCount.clear();
	// records the supporting reads for each kmer
	str_set_uint_map_t &wordSupportReads = buildContigs_wordSupportReads; wordSupportReads.clear();
	// identify repeat kmers (i.e. circles from the de bruijn graph)
	std::set<std::string> &repeatWords = buildContigs_repeatWords; repeatWords.clear();
	std::set<std::string> &unusedWords = buildContigs_unusedWords;

	// get counts and supporting reads for each kmer
	getKmerCounts(wordLength, wordCount, wordSupportReads);
	getRepeatKmers(wordCount, repeatWords);

	// track kmers can be used as seeds for searching for the next contig
	unusedWords.clear();
	for (const str_uint_map_t::value_type &wdct : wordCount) {
		// filter out kmers with too few coverage
		if (wdct.second >= o.minCoverage)
			unusedWords.insert(wdct.first);
	}
	// limit the number of contigs generated for the seek of speed
	while ((!unusedWords.empty()) && (tmpContigs.size() < 2 * o.maxAssemblyCount)) {
		std::string maxWord;
		unsigned maxWordCount(0);
		// get the kmers corresponding the highest count
		for (auto &word : unusedWords) {
			assert(wordCount.count(word) > 0);
			const unsigned currWordCount = wordCount.at(word);
			if (currWordCount > maxWordCount) {
				maxWord = word;
				maxWordCount = currWordCount;
			}
		}
		// solve for a best contig in the graph by a heuristic greedy maxflow-ish criteria
		AssembledContig contig;
		bool isRepeatFound = walk(maxWord, wordLength, wordCount,
				wordSupportReads, repeatWords, unusedWords, contig);
		if (isRepeatFound)
			isAssemblySuccess = false;
		tmpContigs.push_back(contig);
	}

	// done with this now
    return isAssemblySuccess;
}

void AssemblyManager::selectContigs(const unsigned normalReadCount) {
	contigs.clear();
	unsigned finalContigCount(0);
	// a set of reads that has been used to construct contigs, including pseudo ones.
	std::set<unsigned> &usedReads = selectContigs_usedReads; usedReads.clear();
	// a set of pseudo reads that has been used to construct contigs
	std::set<unsigned> &usedPseudoReads = selectContigs_usedPseudoReads; usedPseudoReads.clear();
	std::set<unsigned> &contigs2Remove = selectContigs_contigs2Remove;
	std::set<unsigned> &newSupportReads = selectContigs_newSupportReads;

	// contig are selected based on the number of supporting reads that are not pseudo
	while ((!tmpContigs.empty())
			&& (finalContigCount < o.maxAssemblyCount)) {
		// count unused reads that are not pseudo reads
		const unsigned usedNormalReads = usedReads.size() - usedPseudoReads.size();
		const unsigned unusedNormalReads = normalReadCount - usedNormalReads;
		if (unusedNormalReads < o.minUnusedReads)	return;

		unsigned contigIndex(0);
		contigs2Remove.clear();

		AssembledContig selectedContig;
		unsigned selectedContigIndex;
		unsigned maxSupport(0);
		unsigned maxLength(0);
		for (const AssembledContig &contig : tmpContigs) {
			// identify new supporting reads that were not used for the previously identified contigs
			newSupportReads.clear();
			std::set_difference(contig.supportReads.begin(),
					contig.supportReads.end(), usedReads.begin(),
					usedReads.end(),
					std::inserter(newSupportReads, newSupportReads.end()));

			// count the number of new supporting reads that are not pseudo reads
			unsigned newNormalSupport(0);
			for (const unsigned rd : newSupportReads) {
				const AssemblyReadInfo &rinfo(readInfo[rd]);
				if (!rinfo.isPseudo)
					newNormalSupport++;
			}
			if (newNormalSupport < o.minSupportReads) {
				contigs2Remove.insert(contigIndex);
				contigIndex++;
				continue;
			}

			// either more supporting reads that were not used
			// or the same number of supports but longer contig
			const unsigned currNewSupport = newSupportReads.size();
			const unsigned currContigLen = contig.seq.size();
			const bool isBetterContig(
					(currNewSupport > maxSupport)
							|| ((currNewSupport == maxSupport)
									&& (currContigLen > maxLength)));
			if (isBetterContig) {
				selectedContig = contig;
				selectedContigIndex = contigIndex;
				maxSupport = currNewSupport;
				maxLength = currContigLen;
			}
			contigIndex++;
		}

		// no more contigs selected, selection is done.
		if (maxSupport == 0)
			break;
		// select one contig
		contigs.push_back(selectedContig);
		contigs2Remove.insert(selectedContigIndex);
		// remove selected & failed contigs
		std::set<unsigned>::reverse_iterator it;      //逆向迭代器
		for (it = contigs2Remove.rbegin(); it != contigs2Remove.rend(); it++) { //逆向迭代加即为减
			tmpContigs.erase(tmpContigs.begin() + *it);
		}
		// update the info about used reads
		for (const unsigned rd : selectedContig.supportReads) {
			usedReads.insert(rd);
			AssemblyReadInfo &rinfo(readInfo[rd]);
			// read info record the ID of all contigs that the read supports
			if (!rinfo.isUsed)
				rinfo.isUsed = true;
			rinfo.contigIds.push_back(finalContigCount);
			if (rinfo.isPseudo)
				usedPseudoReads.insert(rd);
		}
		finalContigCount++;
	}
}

#define AssemblySuccessfully true
void AssemblyManager::assembley() {
	const unsigned normalReadCount(reads.size());

	assert(o.alphabet.size() > 1);
	readInfo.resize(normalReadCount);//todo::reverse

	for (unsigned wordLength(o.minWordLength); wordLength <= o.maxWordLength;
			wordLength += o.wordStepSize) {

		if (AssemblySuccessfully == buildContigs(wordLength))
			break;

		// remove pseudo reads from the previous iteration
		const unsigned readCount(reads.size());
		for (unsigned readIndex(0); readIndex < readCount; ++readIndex) {
			AssemblyReadInfo &rinfo(readInfo[readIndex]);
			if (rinfo.isPseudo) {
				reads.erase(reads.begin() + readIndex, reads.end());
				readInfo.erase(readInfo.begin() + readIndex,
						readInfo.end());
				break;
			}
		}

		//  Add contigs from the current iteration as pseudo reads
		for (const AssembledContig &contig : tmpContigs) {
			//contig.debug_print(std::cerr);		//debug::
			if (contig.seq.size() > (wordLength + o.wordStepSize)) {
				reads.push_back(contig.seq);
				readInfo.push_back(AssemblyReadInfo(true));
			}
		}
	}

	// greedy selection of contigs to be returned
	selectContigs(normalReadCount);
}
// void mantaAss_m(std :: map <int,std :: vector <readstr_name>> clupos_str, std :: vector <ConsensusINFO> & consensus_str) {
// 	std :: map <int,std :: vector <readstr_name>> :: iterator itpos = clupos_str.begin();
// 	while(itpos != clupos_str.end()){//which pos
// 		//length of the list
// 		int seq_n = itpos->second.size();
// 		AssemblyManager *am_ = new AssemblyManager[1];
// 		AssemblyManager &am = *am_;
// 		am.o.minWordLength = 30;
// 		am.o.maxWordLength = 70;
// 		for(int i = 0; i < seq_n; ++i) am.reads.emplace_back(itpos->second[i].str);
// 		am.assembley();
// 		auto &contigs = am.getResults();
// 		int cons_n = contigs.size();
// 		if(cons_n > 0){
// 			ConsensusINFO css_;
// 			css_.chrID = itpos->second[0].chrID;
// 			css_.all_cov_readnum = seq_n;
// 			css_.ref_pos = itpos->first * 50;
// 			for(auto &contig:contigs) {
// 				CON_suppR s1;
// 				s1.consensus = contig.seq;
// 				s1.suppsnv_readnum = contig.supportReads.size();
// 				css_.con_suppr.push_back(s1);
// 			}
// 			consensus_str.push_back(css_);
// 		}
// 		delete []am_;
// 		itpos++;
// 	}
// }
// std :: vector <CON_suppR> menta_uesd(std :: vector <string> stringV, int maxAssemblyC){
// 	//length of the list
// 	std :: vector <CON_suppR> conresult;
// 	int seq_n = stringV.size();
// 	AssemblyManager *am_ = new AssemblyManager[1];
// 	AssemblyManager &am = *am_;
// 	am.o.minWordLength = 30;
// 	am.o.maxWordLength = 70;
// 	am.o.maxAssemblyCount = maxAssemblyC;
// 	for(int i = 0; i < seq_n; ++i) am.reads.emplace_back(stringV[i]);
// 	am.assembley();
// 	auto &contigs = am.getResults();
// 	int cons_n = contigs.size();
// 	if(cons_n > 0){
// 		for(auto &contig:contigs) {
// 			CON_suppR s1;
// 			s1.consensus = contig.seq;
// 			s1.suppsnv_readnum = contig.supportReads.size();
// 			conresult.push_back(s1);
// 		}
// 	}
// 	delete []am_;
// 	return conresult;
// }

// void mantaAss_m(std :: map <int,poscluINFO> pos_clu, std :: map <int , std :: vector <vrtReault>> ReaultV, std :: vector <ConsensusINFO> & consensus_str){
// 	int searchR = 0;
// 	std :: map <int,poscluINFO> :: iterator itchr = pos_clu.begin();
//     while (itchr != pos_clu.end()){//which chr
// 		std :: map <int , std :: vector <vrtReault>> :: iterator itResult = ReaultV.find(itchr->first);
// 		if(itResult != ReaultV.end()){
// 			std :: map<uint32_t,int> read_line;//start_pos,aff: 0-consensus,1-result
// 			std :: map <int,std :: vector <readstr_name>> :: iterator itpos1 = itchr->second.clupos_str.begin();
// 			while(itpos1 != itchr->second.clupos_str.end()){//which pos
// 				uint32_t startpos = itpos1->first * 50;
// 				read_line[startpos] = 0;
// 				itpos1++;
// 			} 
// 			int idx = 1;
// 			for(auto result1 : itResult->second){
// 				read_line[result1.ref_pos] = idx;
// 				idx  = idx + 1;
// 			}

// 			std :: map <int,std :: vector <readstr_name>> :: iterator itpos2 = itchr->second.clupos_str.begin();
// 			while(itpos2 != itchr->second.clupos_str.end()){//which pos
// 				std :: vector<string> phase_str[4];//0:00 1:10 2:01 3:11
// 				uint32_t stratpos = itpos2->first * 50;
// 				ConsensusINFO css_;
// 				css_.chrID = itchr->first;
// 				css_.ref_pos = stratpos;
// 				std :: map<uint32_t,int> :: iterator itRL = read_line.find(stratpos);
// 				if(itRL != read_line.end()){
// 					std :: vector<readstr_name> consenflag;//consensus reads ,include name && string
// 					int consenflag_score[1000] = {0};
// 					consenflag.insert(consenflag.end(), itpos2->second.begin(), itpos2->second.end());
// 					int seqll = consenflag.size();
// 					int findnum = 0;
// 					for(; findnum < 2 ; findnum++){
// 						itRL++;
// 						uint32_t stoppos = stratpos + 1000;
// 						if(itRL != read_line.end() && itRL->second != 0 && itRL->first < stoppos){
// 							for(int s1 = 0 ; s1 < seqll ; s1++){
// 								auto itF = std :: find(itResult->second[itRL->second].readsname.begin() , itResult->second[itRL->second].readsname.end() , consenflag[s1].name);
// 								if(itF != itResult->second[itRL->second].readsname.end()){
// 									consenflag_score[s1] = consenflag_score[s1] + findnum + 1;
// 								}
// 							}
// 						}
// 						else break;
// 					}
// 					for(int inputidx = 0 ; inputidx < seqll ; inputidx++){
// 						// cout << inputidx << "\t" << consenflag_score[inputidx] << endl;
// 						phase_str[consenflag_score[inputidx]].push_back(consenflag[inputidx].str);
// 					}

// 					if(findnum == 0){
// 						uint32_t refpos = itpos2->first * 50;
// 						std :: vector <string> readstrs;
// 						for (auto str1 : itpos2->second) readstrs.emplace_back(str1.str); 
// 						std :: vector <CON_suppR> conresult = menta_uesd(readstrs, 2);
// 						if(conresult.size() > 0){
// 							css_.con_suppr.insert(css_.con_suppr.begin() , conresult.begin() , conresult.end());
// 							css_.all_cov_readnum = readstrs.size();
// 							consensus_str.push_back(css_);
// 						}
// 					} else if(findnum == 1){
// 						// int seqall_n = 0;
// 						for(int phaseN = 0 ; phaseN < 2 ; phaseN ++){
// 							int seq_n = phase_str[phaseN].size();
// 							if(seq_n > 0){
// 								// seqall_n = seqall_n + seq_n;
// 								std :: vector <CON_suppR> conresult = menta_uesd(phase_str[phaseN], 1);
// 								if(conresult.size() > 0) css_.con_suppr.insert(css_.con_suppr.begin() , conresult.begin() , conresult.end());
// 							}
// 						}
// 						if(css_.con_suppr.size() > 0){
// 							css_.all_cov_readnum = seqll;
// 							consensus_str.push_back(css_);
// 						}
// 					} else{
// 						int type1 = phase_str[0].size() + phase_str[3].size();
// 						int type2 = phase_str[1].size() + phase_str[2].size();
// 						if(type1 > type2){
// 							// int seqall_n = 0;
// 							for(int phaseN = 0 ; phaseN < 4 ; phaseN ++){
// 								if(phaseN == 0 || phaseN == 3){
// 									int seq_n = phase_str[phaseN].size();
// 									if(seq_n > 0){
// 										// seqall_n = seqall_n + seq_n;
// 										std :: vector <CON_suppR> conresult = menta_uesd(phase_str[phaseN], 1);
// 										if(conresult.size() > 0) css_.con_suppr.insert(css_.con_suppr.begin() , conresult.begin() , conresult.end());
// 									}
// 								}
// 							}
// 							if(css_.con_suppr.size() > 0){
// 								css_.all_cov_readnum = seqll;
// 								consensus_str.push_back(css_);
// 							}
// 						} else{
// 							// int seqall_n = 0;
// 							for(int phaseN = 1 ; phaseN < 3 ; phaseN ++){
// 								int seq_n = phase_str[phaseN].size();
// 								if(seq_n > 0){
// 									// seqall_n = seqall_n + seq_n;
// 									std :: vector <CON_suppR> conresult = menta_uesd(phase_str[phaseN], 1);
// 									if(conresult.size() > 0) css_.con_suppr.insert(css_.con_suppr.begin() , conresult.begin() , conresult.end());
// 								}
// 							}
// 							if(css_.con_suppr.size() > 0){
// 								css_.all_cov_readnum = seqll;
// 								consensus_str.push_back(css_);
// 							}
// 						}
// 					}
// 				} else{
// 					uint32_t refpos = itpos2->first * 50;
// 					std :: vector <string> readstrs;
// 					for (auto str1 : itpos2->second) readstrs.emplace_back(str1.str);
// 					std :: vector <CON_suppR> conresult = menta_uesd(readstrs, 2);
// 					if(conresult.size() > 0){
// 						css_.con_suppr.insert(css_.con_suppr.begin() , conresult.begin() , conresult.end());
// 						css_.all_cov_readnum = readstrs.size();
// 						consensus_str.push_back(css_);
// 					}
// 				}
// 				itpos2++;
// 			}
// 		} else{
// 			std :: map <int,std :: vector <readstr_name>> :: iterator itpos = itchr->second.clupos_str.begin();
// 			while(itpos != itchr->second.clupos_str.end()){//which pos
// 				uint32_t refpos = itpos->first * 50;
// 				ConsensusINFO css_;
// 				css_.chrID = itchr->first;
// 				css_.ref_pos = refpos;
// 				std :: vector <string> readstrs;
// 				for (auto str1 : itpos->second) readstrs.emplace_back(str1.str); 
// 				std :: vector <CON_suppR> conresult = menta_uesd(readstrs, 2);
// 				if(conresult.size() > 0){
// 					css_.con_suppr.insert(css_.con_suppr.begin() , conresult.begin() , conresult.end());
// 					css_.all_cov_readnum = readstrs.size();
// 					consensus_str.push_back(css_);
// 				}
// 				itpos++;
// 			}
// 		}
// 		itchr++;
// 	}
// }