/**********************************************************
callDist.h (c) 2021 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 04/24/2024
-------------------------
Provides main functionality
***********************************************************/

#pragma once

#include "Options.h"
#include "ChromData.h"
#include "Distrib.h"
#include "FqReader.h"
#include <unordered_map>

enum optValue {		// options id
	oINPUT,
	oCHROM,
	oDTYPE,
	oDUPL,
	oPR_DIST,
	oPR_STATS,
	oOUTFILE,
	oTIME,
	oSUMM,
	oVERSION,
	oHELP,
};

// Input data type
enum class InpType { FRAG, READ };

// Base length distribution class 
class LenDist
{
	Distrib	_freq;					// length frequency statistics
	RBedReader* _file = nullptr;	// valid in constructor only!

protected:
	// pass through file records
	template<typename T>
	void Pass(T* obj, RBedReader& file) {
		_file = &file;
		file.Pass(*obj);
		_file = nullptr;
	}

	// Gets the file being read
	inline const RBedReader& File() { return *_file; }

	// Adds frag length to the frequency distribution
	inline void AddLen(fraglen len) { _freq.AddVal(len); }

public:
	// Print actual frequency distribution on a new line
	//	@param ctype: combined type of distribution
	//	@param prDistr: if true then print distribution additionally
	void Print(Distrib::eCType ctype, bool prDistr)
	{
		// empty input is checked already in the 'UniBedReader' constructor
		if (_freq.Size())
			_freq.Print(dout, ctype, prDistr);
	}
};

// 'FragDist' represents fragment's length frequency statistics ('fragment distribution')
class FragDist : public LenDist
{
	unordered_map<ULONG, Read> _waits;	// 'waiting list' - pair mate candidate's collection
	vector<UniBedReader::Issue> _issues = { "duplicates" };
	ULONG	_cnt = 0;					// count of items
	chrlen	_pos[2] = {0,0};			// mates start positions ([0] - neg read, [1] - pos read)
	bool	_dupl;						// if TRUE then duplicate frags are allowed
	bool	_uncheckedPE = true;		// if TRUE then reads have not yet been checked for PE

	// Adds frag to the freq distribution
	//	@r1: first read in a pair
	//	@r2: second read in a pair
	inline void AddFrag(const Read& r1, const Read& r2) { AddLen(r1.FragLen(r2)); }

public:
	FragDist(const char* fname, bool prStats) : _dupl(Options::GetBVal(oDUPL))
	{
		RBedReader file(fname, nullptr, BYTE_UNDEF, eOInfo::NONE, false);
		Pass(this, file);

		UniBedReader::PrintItemCount(_cnt, "fragments");
		if (_issues[0].Cnt) {
			if (_dupl)	_issues[0].Action = UniBedReader::eAction::ACCEPT;
			UniBedReader::PrintStats(_cnt, _issues[0].Cnt, _issues, prStats);
		}
	}

	// treats current read
	bool operator()();

	// Closes current chrom, open next one
	inline void operator()(chrid, chrlen, size_t, chrid) {}

	// Closes last chrom
	inline void operator()(chrid, chrlen, size_t, size_t) {}
};

// 'ReadDist' represents Read's length frequency statistics ('Read distribution')
class ReadDist : public LenDist
{
public:
	// Constructor by BAM/BED file
	ReadDist(const char* fname, bool prStats) {
		RBedReader file(fname, nullptr,
			BYTE(Options::GetRDuplLevel(oDUPL)),
			prStats ? eOInfo::STAT : eOInfo::STD,
			false);

		Pass(this, file);
	}

	// treats current read
	inline bool operator()() { AddLen(File().ItemRegion().Length()); return true; }

	// Closes current chrom, open next one
	inline void operator()(chrid, chrlen, size_t, chrid) {}

	// Closes last chrom
	inline void operator()(chrid, chrlen, size_t, size_t) {}
};

// 'FqReadDist' represents 'row' (fastq) Read's length frequency statistics ('Read distribution')
class FqReadDist : public LenDist
{
public:
	// Constructor by FastQ file
	FqReadDist(const char* fname)  {
		ULONG	cnt = 0;		// count of reads

		for (FqReader file(fname); file.GetSequence(); cnt++)
			AddLen(file.ReadLength());
		UniBedReader::PrintItemCount(cnt, FT::ItemTitle(FT::eType::FQ, cnt > 1));
		dout << LF;
	}
};