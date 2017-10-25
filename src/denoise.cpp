#include <iostream>
#include <fstream>
#include <bitset>
#include <string>
#include <vector>
#include <array>
#include <list>
#include <algorithm>
#include <cstring>
#include <string>
#include <numeric>
#include <cstdio>
#include <omp.h>
#include "BooPHF.h"
#include "config.h"


typedef std::bitset<3*readlen> bitset;
typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
typedef boomphf::mphf<  u_int64_t, hasher_t> boophf_t;

class bbhashdict
{
	public:
	boophf_t * bphf;
	uint32_t numkeys;
	uint32_t *startpos;
	uint32_t *read_id;
	bool *empty_bin;
	void findpos(int64_t *dictidx, uint32_t &startposidx);
	void remove(int64_t *dictidx, uint32_t &startposidx, uint32_t current);
	bbhashdict()
	{
		bphf = NULL;
		startpos = NULL;
		read_id = NULL;
	}
	~bbhashdict()
	{
		delete[] startpos;
		delete[] read_id;
		delete bphf;
	}	
};

uint32_t numreads = 0;
std::string param1;
double param2,param3;

std::string infilenumreads;

std::string infile;
std::string infile_quality;
std::string outfile;
std::string outfile_quality;
std::string outdir;


//Some global arrays (some initialized in setglobalarrays())
char revinttochar[8] = {'A','N','G','#','C','#','T','#'};//used in bitsettostring
char inttochar[] = {'A','C','G','T','N'};
char chartorevchar[128];//A-T etc for reverse complement
bitset basemask[readlen][128];//bitset for A,G,C,T,N at each position 
//used in stringtobitset, chartobitset and bitsettostring
bitset positionmask[readlen];//bitset for each position (1 at two bits and 0 elsewhere)
//used in bitsettostring
bitset mask63;//bitset with 63 bits set to 1 (used in bitsettostring for conversion to ullong)
char longtochar[] = {'A','C','G','T','N'};
long chartolong[128];


int *dict_start;
int *dict_end; 
bitset mask[maxshift];
bitset revmask[maxshift];
bitset indexmask[numdict];


double quality_to_p_pbar[42];
double log_1_p[42];
double log_p_3[42];


void denoise(bitset *read, char(*quality)[readlen+1], bbhashdict *dict);

void find_overlapping_reads(bitset &current_bitset, std::string &current_quality, std::vector<std::string> &overlap_reads, std::vector<int> &overlap_shift, 
							std::vector<std::string> &overlap_quality, bitset *read, bbhashdict *dict);


void denoise_read(std::string &current_read, std::string &current_quality, std::string &denoised_read, std::string &denoised_quality ,std::vector<std::string> &overlap_reads, 
		std::vector<int> &overlap_shift, std::vector<std::string> &overlap_quality);

void readDnaFile(bitset *read);

void readQualityFile(char(*quality)[readlen+1]);

void constructdictionary(bitset *read, bbhashdict *dict);

bitset stringtobitset(std::string s);

void generateindexmasks(bitset *indexmask);

void generatemasks(bitset *mask,bitset *revmask);

std::string bitsettostring(bitset b);

std::string reverse_complement(std::string s);

void setglobalarrays();


int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outdir = basedir + "/output/";
	infile = basedir + "/output/input.dna";
	infile_quality = basedir + "/output/input.quality";
	outfile = basedir + "/output/output.dna";
	outfile_quality = basedir + "/output/output.quality";
	infilenumreads = basedir + "/output/numreads.bin";
	
	param1 = std::string(argv[2]);
	param2 = atof(argv[3]);
	param3 = atof(argv[4]);

	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.read((char*)&numreads,sizeof(uint32_t));	
	omp_set_num_threads(num_thr);	
	setglobalarrays();

	std::cout << "Reading file: " << infile << std::endl;

	bitset *read = new bitset [numreads];
	readDnaFile(read);

	char(*quality)[readlen+1] = new char [numreads][readlen+1];
	readQualityFile(quality);

	bbhashdict dict[numdict];
	std::cout << "Constructing dictionaries\n";
	constructdictionary(read,dict);

	std::cout << "Denoising reads\n";
	denoise(read,quality,dict);

	delete[] read;
	delete[] quality;
	std::cout << "Done!\n";
	return 0;
}



void denoise(bitset *read, char(*quality)[readlen+1], bbhashdict *dict)
{
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	uint32_t i, stop;	
	//doing initial setup and first read
	i = uint64_t(tid)*numreads/omp_get_num_threads();//spread out first read equally
	stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
	if(tid == omp_get_num_threads()-1)
		stop = numreads;
	std::ofstream fout(outfile+'.'+std::to_string(tid));
	std::ofstream fout_quality(outfile_quality+'.'+std::to_string(tid));
	while(i < stop)
	{
		bitset current_bitset = read[i];
		std::string current_quality = quality[i];
		std::vector<std::string> overlap_reads;//reads overlapping with current read
		std::vector<int> overlap_shift;//shift of overlapping reads wrt current read - 0 for current read, -ve for reads on left
		std::vector<std::string> overlap_quality;//qualities of overlapping reads (reversed if read was reversed)
		find_overlapping_reads(current_bitset,current_quality,overlap_reads,overlap_shift,overlap_quality,read,dict);
		std::string current_read = bitsettostring(current_bitset);
		std::string denoised_read = current_read, denoised_quality = current_quality;
		denoise_read(current_read, current_quality,denoised_read,denoised_quality,overlap_reads,overlap_shift,overlap_quality);
		fout << denoised_read << "\n";
		fout_quality << denoised_quality << "\n";
		i++;	
	}
	fout.close();
	fout_quality.close();
	}//parallel end

	//combine files produced by the threads
	std::ofstream fout(outfile);
	std::ofstream fout_quality(outfile_quality);

	for(int tid = 0; tid < num_thr; tid++)
	{
		std::ifstream fin(outfile+'.'+std::to_string(tid));
		std::ifstream fin_quality(outfile_quality+'.'+std::to_string(tid));

		fout << fin.rdbuf();
		fout_quality << fin_quality.rdbuf();
		
		fout.clear();
		fout_quality.clear();

		remove((outfile+'.'+std::to_string(tid)).c_str());
		remove((outfile_quality+'.'+std::to_string(tid)).c_str());
	}
	return;
}

void find_overlapping_reads(bitset &current_bitset, std::string &current_quality, std::vector<std::string> &overlap_reads, std::vector<int> &overlap_shift, 
							std::vector<std::string> &overlap_quality, bitset *read, bbhashdict *dict)
{
	return;
}

void denoise_read(std::string &current_read, std::string &current_quality, std::string &denoised_read, std::string &denoised_quality ,std::vector<std::string> &overlap_reads, 
		std::vector<int> &overlap_shift, std::vector<std::string> &overlap_quality)
{
	return;
}


void readDnaFile(bitset *read)
{
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	uint32_t i, stop;	
	//doing initial setup and first read
	i = uint64_t(tid)*numreads/omp_get_num_threads();//spread out first read equally
	stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
	if(tid == omp_get_num_threads()-1)
		stop = numreads;
	std::ifstream f(infile, std::ifstream::in);
	f.seekg(uint64_t(i)*(readlen+1), f.beg);
	std::string s;
	while(i < stop)
	{
		std::getline(f,s);
		read[i] = stringtobitset(s);
		i++;
	}
	f.close();
	}
	return;
}


void readQualityFile(char(*quality)[readlen+1])
{
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	uint32_t i, stop;	
	//doing initial setup and first read
	i = uint64_t(tid)*numreads/omp_get_num_threads();//spread out first read equally
	stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
	if(tid == omp_get_num_threads()-1)
		stop = numreads;
	std::ifstream f(infile_quality, std::ifstream::in);
	f.seekg(uint64_t(i)*(readlen+1), f.beg);
	std::string s;
	while(i < stop)
	{
		f.getline(quality[i],readlen+1);
		i++;
	}
	f.close();
	}
	return;
}

void constructdictionary(bitset *read, bbhashdict *dict)
{
	for(int j = 0; j < numdict; j++)
	{
		uint64_t *ull = new uint64_t[numreads];
		#pragma omp parallel
		{
		bitset b;
		int tid = omp_get_thread_num();
		std::ofstream foutkey(outdir+std::string("keys.bin.")+std::to_string(tid),std::ios::binary);
		uint32_t i, stop;
		i = uint64_t(tid)*numreads/omp_get_num_threads();
		stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
		if(tid == omp_get_num_threads()-1)
			stop = numreads;
		//compute keys and write to file and store in ull
		for(; i < stop; i++)
		{
			b = read[i]&indexmask[j];
			ull[i] = (b>>3*dict_start[j]).to_ullong();
			foutkey.write((char*)&ull[i], sizeof(uint64_t));
		}
		foutkey.close();
		}//parallel end
		
		//deduplicating ull
		std::sort(ull,ull+numreads);
		uint32_t k = 0;
		for (uint32_t i = 1; i < numreads; i++) 
		        if (ull[i] != ull[k])         
				ull[++k] = ull[i];
		dict[j].numkeys = k+1;
		//construct mphf
		auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(ull), static_cast<const u_int64_t*>(ull+dict[j].numkeys));
		double gammaFactor = 5.0;//balance between speed and memory
		dict[j].bphf = new boomphf::mphf<u_int64_t,hasher_t>(dict[j].numkeys,data_iterator,num_thr,gammaFactor,true,false);
	
		delete[] ull;

		//compute hashes for all reads
		#pragma omp parallel
		{
		int tid = omp_get_thread_num();	
		std::ifstream finkey(outdir+std::string("keys.bin.")+std::to_string(tid),std::ios::binary);
		std::ofstream fouthash(outdir+std::string("hash.bin.")+std::to_string(tid)+'.'+std::to_string(j),std::ios::binary);
		uint64_t currentkey,currenthash;
		uint32_t i, stop;
		i = uint64_t(tid)*numreads/omp_get_num_threads();
		stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
		if(tid == omp_get_num_threads()-1)
			stop = numreads;
		for(; i < stop; i++)
		{
			finkey.read((char*)&currentkey, sizeof(uint64_t));
			currenthash = (dict[j].bphf)->lookup(currentkey);
			fouthash.write((char*)&currenthash, sizeof(uint64_t));
		}
		finkey.close();
		remove((outdir+std::string("keys.bin.")+std::to_string(tid)).c_str());
		fouthash.close();
		}//parallel end
		
	}
		
	//for rest of the function, use numdict threads to parallelize
	omp_set_num_threads(std::min(numdict,num_thr));
	#pragma omp parallel
	{
	#pragma omp for	
	for(int j = 0; j < numdict; j++)
	{
		//fill startpos by first storing numbers and then doing cumulative sum
		dict[j].startpos = new uint32_t[dict[j].numkeys+1]();//1 extra to store end pos of last key
		uint64_t currenthash;
		for(int tid = 0; tid < num_thr; tid++)
		{
			std::ifstream finhash(outdir+std::string("hash.bin.")+std::to_string(tid)+'.'+std::to_string(j),std::ios::binary);
			finhash.read((char*)&currenthash,sizeof(uint64_t));
			while(!finhash.eof())
			{
				dict[j].startpos[currenthash+1]++;
				finhash.read((char*)&currenthash,sizeof(uint64_t));
			}
			finhash.close();
		}
	
		dict[j].empty_bin = new bool[dict[j].numkeys]();
		for(uint32_t i = 1; i < dict[j].numkeys; i++)
			dict[j].startpos[i] =  dict[j].startpos[i] +  dict[j].startpos[i-1];

		//insert elements in the dict array
		dict[j].read_id = new uint32_t[numreads];
		uint32_t i = 0;
		for(int tid = 0; tid < num_thr; tid++)
		{
			std::ifstream finhash(outdir+std::string("hash.bin.")+std::to_string(tid)+'.'+std::to_string(j),std::ios::binary);
			finhash.read((char*)&currenthash,sizeof(uint64_t));
			while(!finhash.eof())
			{
				dict[j].read_id[dict[j].startpos[currenthash]++] = i;
				i++;
				finhash.read((char*)&currenthash,sizeof(uint64_t));
			}
			finhash.close();
			remove((outdir+std::string("hash.bin.")+std::to_string(tid)+'.'+std::to_string(j)).c_str());
		}
		
		//correcting startpos array modified during insertion
		for(int64_t i = dict[j].numkeys; i >= 1 ; i--)	
			dict[j].startpos[i] = dict[j].startpos[i-1];
		dict[j].startpos[0] = 0;
	}//for end
	}//parallel end
	omp_set_num_threads(num_thr);
	return;
}

bitset stringtobitset(std::string s)
{
	bitset b;
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return b;
}

void generateindexmasks(bitset *indexmask)
//masks for dictionary positions
{
	for(int j = 0; j < numdict; j++)
		indexmask[j].reset();
	for(int j = 0; j < numdict; j++)
		for(int i = 3*dict_start[j]; i < 3*(dict_end[j]+1); i++)
			indexmask[j][i] = 1;
	return;
}

void generatemasks(bitset *mask,bitset *revmask)
{
	for(int i = 0; i < maxshift; i++)
	{	
		mask[i].reset();
		revmask[i].reset();
		for(int j = 0; j < 3*readlen - 3*i; j++)
			mask[i][j] = 1;
		for(int j = 3*i; j < 3*readlen; j++)
			revmask[i][j] = 1; 	
	}
	return;
}

void bbhashdict::findpos(int64_t *dictidx, uint32_t &startposidx)
{
	dictidx[0] = startpos[startposidx];
	auto endidx = startpos[startposidx+1];
	if(read_id[endidx-1] == numreads)//means exactly one read has been removed
		dictidx[1] = endidx-1;
	else if(read_id[endidx-1] == numreads+1)//means two or more reads have been removed (in this case second last entry stores the number of reads left)
		dictidx[1] = dictidx[0] + read_id[endidx-2];
	else
		dictidx[1] = endidx;//no read deleted
	return;
}

void bbhashdict::remove(int64_t *dictidx, uint32_t &startposidx, uint32_t current)
{
	auto size = dictidx[1] - dictidx[0];
	if(size == 1)//just one read left in bin
	{
		empty_bin[startposidx] = 1;
		return; //need to keep one read to check during matching
	}
	uint32_t pos = std::lower_bound(read_id+dictidx[0],read_id+dictidx[1],current)-(read_id+dictidx[0]);
	
	std::move(read_id+dictidx[0]+pos+1,read_id+dictidx[1],read_id+dictidx[0]+pos);
	auto endidx = startpos[startposidx+1];
	if(dictidx[1] == endidx)//this is first read to be deleted
		read_id[endidx-1] = numreads;
	else if(read_id[endidx-1] == numreads)//exactly one read has been deleted till now
	{
		read_id[endidx-1] = numreads + 1;
		read_id[endidx-2] = size - 1;//number of reads left in bin
	}
	else//more than two reads have been deleted
		read_id[endidx-2]--;

	return;	
}



std::string bitsettostring(bitset b)
{
	char s[readlen+1];
	s[readlen] = '\0';
	unsigned long long ull,rem;
	for(int i = 0; i < 3*readlen/63+1; i++)
	{	
		ull = (b&mask63).to_ullong();
		b>>=63;
		for(int j = 21*i  ; j < 21*i+21 && j < readlen ; j++)
		{
			s[j] = revinttochar[ull%8];	
			ull/=8;
		}
	}
	std::string s1 = s;
	return s1;
}


std::string reverse_complement(std::string s)
{
	std::string s1(s);
	for(int j = 0; j < readlen; j++)
		s1[j] = chartorevchar[s[readlen-j-1]];
	return s1;
}

void setglobalarrays()
{
	chartorevchar['A'] = 'T';
	chartorevchar['C'] = 'G';
	chartorevchar['G'] = 'C';
	chartorevchar['T'] = 'A';
	chartorevchar['N'] = 'N';

	chartolong['A'] = 0;
	chartolong['C'] = 1;
	chartolong['G'] = 2;
	chartolong['T'] = 3;
	chartolong['N'] = 4;
	#if numdict == 1
	{
		dict_start = new int[1];
		dict_end = new int[1];
		dict_start[0] = dict1_start;
		dict_end[0] = dict1_end;
	}
	#elif numdict == 2
	{
		dict_start = new int[2];
		dict_end = new int[2];
		dict_start[0] = dict1_start;
		dict_end[0] = dict1_end;
		dict_start[1] = dict2_start;
		dict_end[1] = dict2_end;
	}
	#elif numdict == 3
	{
		dict_start = new int[3];
		dict_end = new int[3];
		dict_start[0] = dict1_start;
		dict_end[0] = dict1_end;
		dict_start[1] = dict2_start;
		dict_end[1] = dict2_end;
		dict_start[2] = dict3_start;
		dict_end[2] = dict3_end;
	}
	#elif numdict == 4
	{
		dict_start = new int[4];
		dict_end = new int[4];
		dict_start[0] = dict1_start;
		dict_end[0] = dict1_end;
		dict_start[1] = dict2_start;
		dict_end[1] = dict2_end;
		dict_start[2] = dict3_start;
		dict_end[2] = dict3_end;
		dict_start[3] = dict4_start;
		dict_end[3] = dict4_end;
	}
	#endif
	for(int i = 0; i < 63; i++)
		mask63[i] = 1;
	for(int i = 0; i < readlen; i++)
	{
		basemask[i]['A'][3*i] = 0;
		basemask[i]['A'][3*i+1] = 0;
		basemask[i]['A'][3*i+2] = 0;
		basemask[i]['C'][3*i] = 0;
		basemask[i]['C'][3*i+1] = 0;
		basemask[i]['C'][3*i+2] = 1;
		basemask[i]['G'][3*i] = 0;
		basemask[i]['G'][3*i+1] = 1;
		basemask[i]['G'][3*i+2] = 0;
		basemask[i]['T'][3*i] = 0;
		basemask[i]['T'][3*i+1] = 1;
		basemask[i]['T'][3*i+2] = 1;
		basemask[i]['N'][3*i] = 1;
		basemask[i]['N'][3*i+1] = 0;
		basemask[i]['N'][3*i+2] = 0;
		positionmask[i][3*i] = 1;
		positionmask[i][3*i+1] = 1;
		positionmask[i][3*i+2] = 1;
	}		
	generatemasks(mask,revmask);
	generateindexmasks(indexmask);

	for(int i=0; i<42;i++) 
	{	
		double _prob = pow(10.0,-((double)i/10.0));
		quality_to_p_pbar[i] = _prob/(1.0-_prob);
		log_1_p[i] = log10(1-_prob);
		log_p_3[i] = log10(_prob/3);
	}
	return;
}

/*
void reorder(bitset *read, bbhashdict *dict)
{
	omp_lock_t *dict_lock = new omp_lock_t [num_locks];
	omp_lock_t *read_lock = new omp_lock_t [num_locks];
	for(int j = 0; j < num_locks; j++)
	{
		omp_init_lock(&dict_lock[j]);
		omp_init_lock(&read_lock[j]);
	}
	uint8_t _readlen = readlen;//used for writing to binary file
	bitset mask[maxshift];
	bitset revmask[maxshift];
	generatemasks(mask,revmask);
	bitset mask1[numdict];
	generateindexmasks(mask1);
	bool *remainingreads = new bool [numreads];
	std::fill(remainingreads, remainingreads+numreads,1);

	//we go through remainingreads array from behind as that speeds up deletion from bin arrays

	uint32_t firstread = 0, unmatched[num_thr];
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	std::string tid_str = std::to_string(tid);	
	std::ofstream foutRC(outfileRC + '.' + tid_str,std::ofstream::out);
	std::ofstream foutflag(outfileflag + '.' + tid_str,std::ofstream::out);
	std::ofstream foutpos(outfilepos + '.' + tid_str,std::ofstream::out|std::ios::binary);
	std::ofstream foutorder(outfileorder + '.' + tid_str,std::ofstream::out|std::ios::binary);
	std::ofstream foutorder_s(outfileorder + ".singleton." + tid_str,std::ofstream::out|std::ios::binary);
	
	unmatched[tid] = 0;
	bitset ref,revref,b;
	int count[4][readlen];
	int64_t dictidx[2];//to store the start and end index (end not inclusive) in the dict read_id array
	uint32_t startposidx;//index in startpos
	bool flag = 0, done = 0, prev_unmatched = false;
	uint32_t current, prev;
	uint64_t ull;
	//flag to check if match was found or not
	
	int64_t remainingpos = numreads-1;//used for searching next unmatched read when no match is found
	#pragma omp critical
	{//doing initial setup and first read
		current = firstread;	
		firstread += numreads/omp_get_num_threads();//spread out first read equally
		remainingreads[current] = 0;
		unmatched[tid]++;
	}
	#pragma omp barrier
	updaterefcount(read[current],ref,revref,count,true,false,0);
	prev_unmatched = true;	
	prev = current;
	uint32_t numdone= 0;
	while(!done)
	{
//		numdone++;
//		if(numdone%1000000==0)
//			std::cout<<tid<<":"<<numdone<<"\n";
		//delete the read from the corresponding dictionary bins
		for(int l = 0; l < numdict; l++)
		{	b = read[current]&mask1[l];
			ull = (b>>2*dict_start[l]).to_ullong();
			startposidx = dict[l].bphf->lookup(ull);
			//check if any other thread is modifying same dictpos
			omp_set_lock(&dict_lock[startposidx & 0xFFFFFF]);
			dict[l].findpos(dictidx,startposidx);
			dict[l].remove(dictidx,startposidx,current);
			omp_unset_lock(&dict_lock[startposidx & 0xFFFFFF]);
		}
		flag = 0;
		uint32_t k;
		for(int j = 0; j < maxshift; j++)
		{
			//find forward match
			for(int l = 0; l < numdict; l++)
			{
				if(dict_end[l]+j >= readlen)
					continue;
				b = ref&mask1[l];
				ull = (b>>2*dict_start[l]).to_ullong();
				startposidx = dict[l].bphf->lookup(ull);
				if(startposidx >= dict[l].numkeys)//not found
					continue;
				//check if any other thread is modifying same dictpos
				omp_set_lock(&dict_lock[startposidx & 0xFFFFFF]);
				dict[l].findpos(dictidx,startposidx);
				if(dict[l].empty_bin[startposidx])//bin is empty
				{
					omp_unset_lock(&dict_lock[startposidx & 0xFFFFFF]);
					continue;
				}
				uint64_t ull1 = ((read[dict[l].read_id[dictidx[0]]]&mask1[l])>>2*dict_start[l]).to_ullong();
				if(ull == ull1)//checking if ull is actually the key for this bin
				{	
					for (int64_t i = dictidx[1] - 1 ; i >= dictidx[0] && i >= dictidx[1] - maxsearch; i--)
					{
						auto rid = dict[l].read_id[i];
						if((ref^(read[rid]&mask[j])).count()<=thresh)
						{	
							omp_set_lock(&read_lock[rid & 0xFFFFFF]);
							if(remainingreads[rid])
							{
								remainingreads[rid]=0;
								k = rid;
								flag = 1;
							}
							omp_unset_lock(&read_lock[rid & 0xFFFFFF]);
							if(flag == 1)
								break;
						}
					}
				}
				omp_unset_lock(&dict_lock[startposidx & 0xFFFFFF]);
				
				if(flag == 1)
				{
					current = k;
					updaterefcount(read[current],ref,revref,count,false,false,j);
					if(prev_unmatched == true)//prev read not singleton, write it now
					{
						foutRC << 'd';
						foutorder.write((char*)&prev,sizeof(uint32_t));
						foutflag << 0;//for unmatched
						foutpos.write((char*)&_readlen,sizeof(uint8_t));
					}	
					foutRC << 'd';
					foutorder.write((char*)&current,sizeof(uint32_t));
					foutflag << 1;//for matched
					foutpos.write((char*)&j,sizeof(uint8_t));
					
					prev_unmatched = false;
					break;
				}
				
			}
			if(flag==1)
				break;	

			//find reverse match
			for(int l = 0; l < numdict; l++)
			{
				if(dict_start[l] <= j)
					continue;
				b = revref&mask1[l];
				ull = (b>>2*dict_start[l]).to_ullong();
				startposidx = dict[l].bphf->lookup(ull);
				if(startposidx >= dict[l].numkeys)//not found
					continue;
				//check if any other thread is modifying same dictpos
				omp_set_lock(&dict_lock[startposidx & 0xFFFFFF]);
				dict[l].findpos(dictidx,startposidx);
				if(dict[l].empty_bin[startposidx])//bin is empty
				{	
					omp_unset_lock(&dict_lock[startposidx & 0xFFFFFF]);
					continue;
				}
				uint64_t ull1 = ((read[dict[l].read_id[dictidx[0]]]&mask1[l])>>2*dict_start[l]).to_ullong();
				if(ull == ull1)//checking if ull is actually the key for this bin
				{
					for (int64_t i = dictidx[1] - 1 ; i >= dictidx[0] && i >= dictidx[1] - maxsearch; i--)
					{
						auto rid = dict[l].read_id[i];
						if((revref^(read[rid]&revmask[j])).count()<=thresh)
						{	
							omp_set_lock(&read_lock[rid & 0xFFFFFF]);
							if(remainingreads[rid])
							{
								remainingreads[rid]=0;
								k = rid;
								flag = 1;
							}
							omp_unset_lock(&read_lock[rid & 0xFFFFFF]);
							if(flag == 1)
								break;
						}
					}
				}
				omp_unset_lock(&dict_lock[startposidx & 0xFFFFFF]);
				if(flag == 1)
				{
					current = k;
					updaterefcount(read[current],ref,revref,count,false,true,j);
					if(prev_unmatched == true)//prev read not singleton, write it now
					{
						foutRC << 'd';
						foutorder.write((char*)&prev,sizeof(uint32_t));
						foutflag << 0;//for unmatched
						foutpos.write((char*)&_readlen,sizeof(uint8_t));
					}
					foutRC << 'r';
					foutorder.write((char*)&current,sizeof(uint32_t));
					foutflag << 1;//for matched
					foutpos.write((char*)&j,sizeof(uint8_t));
					
					prev_unmatched = false;
					break;
				}
			}	
			if(flag==1)
				break;
		
			revref<<=2;
			ref>>=2;
		}
		if(flag == 0)//no match found
		{
			for(int64_t j = remainingpos; j>=0; j--)
			
				if(remainingreads[j] == 1)
				{
					omp_set_lock(&read_lock[j & 0xFFFFFF]);
					if(remainingreads[j])//checking again inside critical block
					{
						current = j;
						remainingpos = j-1;
						remainingreads[j] = 0;
						flag = 1;
						unmatched[tid]++;
					}
					omp_unset_lock(&read_lock[j & 0xFFFFFF]);
					if(flag == 1)
						break;
				}
					
			if(flag == 0)
			{
				if(prev_unmatched == true)//last read was singleton, write it now
				{
					foutorder_s.write((char*)&prev,sizeof(uint32_t));
				}
				done = 1;//no reads left
			}
			else
			{
				updaterefcount(read[current],ref,revref,count,true,false,0);
				if(prev_unmatched == true)//prev read singleton, write it now
				{
					foutorder_s.write((char*)&prev,sizeof(uint32_t));
				}
				prev_unmatched = true;
				prev = current;
			}
		}
	}//while(!done) end
	
	foutRC.close();
	foutorder.close();
	foutflag.close();
	foutpos.close();
	foutorder_s.close();
//	std::cout << tid << ":Done"<<"\n";
	}//parallel end
		
	delete[] remainingreads;
		
	std::cout << "Reordering done, "<< std::accumulate(unmatched,unmatched+num_thr,0) <<" were unmatched\n";
	return;
}

*/

// void bitsettostring(bitset b,char *s)
// {
// 	unsigned long long ull,rem;
// 	for(int i = 0; i < 2*readlen/64+1; i++)
// 	{	
// 		ull = (b&mask64).to_ullong();
// 		b>>=64;
// 		for(int j = 32*i  ; j < 32*i+32 && j < readlen ; j++)
// 		{
// 			s[j] = revinttochar[ull%4];	
// 			ull/=4;
// 		}
// 	}
// 	return;
// }

// bitset chartobitset(char *s)
// {
// 	bitset b;
// 	for(int i = 0; i < readlen; i++)
// 		b |= basemask[i][s[i]];
// 	return b;
// }

// void reverse_complement(char* s, char* s1)
// {
// 	for(int j = 0; j < readlen; j++)
// 		s1[j] = chartorevchar[s[readlen-j-1]];
// 	return;
// }



/*float alternate_thresh = 0.0;
bitset stringtobitset(std::string s);

std::string bitsettostring(bitset b);

void readsingletons(bitset *read, uint32_t *order_s, std::string *quality_sN);

void constructdictionary(bitset *read, bbhashdict *dict);


std::string reverse_complement(std::string s);

void generateindexmasks(bitset *mask1);

void encode(bitset *read, bbhashdict *dict, uint32_t *order_s, std::string *quality_sN);

void packbits();

std::string generateRef(std::vector<std::array<long,5>> count);

std::vector<std::array<long,5>> buildcontig(std::list<std::string> reads, std::list<long> pos, uint32_t list_size);//using lists to enable faster insertion of singletons
void writecontig(std::vector<std::array<long,5>> count, std::list<long> &pos, std::list<std::string> &reads, std::list<uint32_t> &order, std::list<char> &RC, std::list<std::string> &quality, std::ofstream& f_seq, std::ofstream& f_pos, std::ofstream& f_noise, std::ofstream& f_noisepos, std::ofstream& f_order, std::ofstream &f_RC, std::ofstream &f_order_N_pe, std::ofstream &f_flag_N, std::ofstream &f_quality, std::ofstream &f_quality_N, uint32_t list_size);

void getDataParams();

void setglobalarrays();

void encode(bitset *read, bbhashdict *dict, uint32_t *order_s, std::string *quality_sN)
{
	omp_lock_t *read_lock = new omp_lock_t [numreads_s+numreads_N];
	omp_lock_t *dict_lock = new omp_lock_t [numreads_s+numreads_N];
	for(int j = 0; j < numreads_s+numreads_N; j++)
	{
		omp_init_lock(&read_lock[j]);
		omp_init_lock(&dict_lock[j]);
	}
	bool *remainingreads = new bool [numreads_s+numreads_N];
	std::fill(remainingreads, remainingreads+numreads_s+numreads_N,1);

	bitset mask1[numdict_s];
	generateindexmasks(mask1);
	std::cout<<"Encoding reads\n";
	#pragma omp parallel 
	{
	int tid = omp_get_thread_num();
	std::ifstream f(infile);
	std::ifstream in_flag(infile_flag);
	std::ifstream in_pos(infile_pos,std::ios::binary);
	std::ifstream in_order(infile_order,std::ios::binary);
	std::ifstream in_RC(infile_RC);
	std::ifstream in_quality(infile_quality);
	std::ofstream f_seq(outfile_seq+'.'+std::to_string(tid));
	std::ofstream f_pos(outfile_pos+'.'+std::to_string(tid));
	std::ofstream f_noise(outfile_noise+'.'+std::to_string(tid));
	std::ofstream f_noisepos(outfile_noisepos+'.'+std::to_string(tid));
	std::ofstream f_order(infile_order+'.'+std::to_string(tid),std::ios::binary);
	std::ofstream f_RC(infile_RC+'.'+std::to_string(tid));
	std::ofstream f_order_N_pe(outfile_order_N_pe+'.'+std::to_string(tid),std::ios::binary);
	std::ofstream f_flag_N(outfile_flag_N+'.'+std::to_string(tid),std::ios::binary);
	std::ofstream f_quality(outfile_quality+'.'+std::to_string(tid));
	std::ofstream f_quality_N(outfile_quality_N+'.'+std::to_string(tid));

	uint64_t i, stop;
	int64_t dictidx[2];//to store the start and end index (end not inclusive) in the dict read_id array
	uint32_t startposidx;//index in startpos
	uint64_t ull;
	bool flag = 0;
	//flag to check if match was found or not
	//doing initial setup and first read
	i = uint64_t(tid)*numreads/omp_get_num_threads();//spread out first read equally
	stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
	if(tid == omp_get_num_threads()-1)
		stop = numreads;
	f.seekg(uint64_t(i)*(readlen+1), f.beg);
	in_quality.seekg(uint64_t(i)*(readlen+1), in_quality.beg);
	in_flag.seekg(i, in_flag.beg);
	in_pos.seekg(i*sizeof(uint8_t),in_pos.beg);
	std::vector<std::array<long,5>> count;

// 	char c;
// 	std::vector<std::string> reads;
// 	std::vector<long> pos;
// =======
	in_order.seekg(i*sizeof(uint32_t),in_order.beg);
	in_RC.seekg(i,in_RC.beg);
	std::string current,ref,q;
	bitset forward_bitset,reverse_bitset,b;
	char c,rc;
	std::list<std::string> reads;
	std::list<long> pos;
	std::list<char> RC;
	std::list<uint32_t> order;
	std::list<std::string> quality;
	uint8_t p;
	uint32_t ord, list_size = 0;//list_size variable introduced because list::size() was running very slowly
					// on UIUC machine
	std::list<uint32_t> deleted_rids[numdict_s];
	while(i < stop)
	{
		std::getline(f,current);
		c = in_flag.get();
		rc = in_RC.get();
		std::getline(in_quality,q);
		in_pos.read((char*)&p,sizeof(uint8_t));
		in_order.read((char*)&ord,sizeof(uint32_t));
		if(c=='0'||list_size>10000000)//so that reads vector doesn't get too large
		{
			if(list_size!=0)
			{
				count = buildcontig(reads,pos,list_size);
				//writecontig(count,pos,reads,f_seq,f_pos,f_noise,f_noisepos);
				ref = generateRef(count);
				//try to align the singleton reads to ref
				//first create bitsets from first readlen positions of ref
				forward_bitset = stringtobitset(ref.substr(0,readlen));
				reverse_bitset = stringtobitset(reverse_complement(ref.substr(0,readlen)));
				auto pos_it = pos.begin();
				auto reads_it = reads.begin();
				auto order_it = order.begin();
				auto RC_it = RC.begin();
				auto quality_it = quality.begin();
				long nextpos = 0;
				//storing positions of non-singleton reads, so that singletons 
				//are inserted correctly
				*pos_it = 0; //putting 0 in first pos (originally was readlen and will be converted 
				//to readlen in writecontig
				//convert pos to cumulative list
				long cumsum = 0;
				for(auto it = pos.begin();it!=pos.end();++it)
				{
					*it = cumsum + *it;
					cumsum = *it;
				}	
				for(uint32_t j = 0; j < ref.size()-readlen+1; j++)
				{
					if(j == nextpos)//go through non-singleton reads at this pos
					{
						while(pos_it!= pos.end())
						{
							if(*pos_it == j)
							{
								++pos_it;++reads_it;++order_it;++RC_it;++quality_it;
							}
							else
							{
								nextpos = *pos_it;
								break;
							}
						}
					}	
					//search for singleton reads
					for(int l = 0; l < numdict_s; l++)//forward
					{
						b = forward_bitset&mask1[l];
						ull = (b>>3*dict_start[l]).to_ullong();
					//	std::vector<uint32_t> deleted_rids;
						//if(ull == 0 || ull == ((uint64_t(1)<<(2*(dict_end[l]-dict_start[l]+1)))-1)) 
						//	continue;//added because there were too many of these
						//making the whole thing slow due to locking
						startposidx = dict[l].bphf->lookup(ull);
						if(startposidx >= dict[l].numkeys)//not found
							continue;
						//check if any other thread is modifying same dictpos
						if(!omp_test_lock(&dict_lock[startposidx]))
							continue;
						dict[l].findpos(dictidx,startposidx);
						if(dict[l].empty_bin[startposidx])//bin is empty
						{
							omp_unset_lock(&dict_lock[startposidx]);
							continue;
						}
						uint64_t ull1 = ((read[dict[l].read_id[dictidx[0]]]&mask1[l])>>3*dict_start[l]).to_ullong();
						if(ull == ull1)//checking if ull is actually the key for this bin
						{	
							for (int64_t i = dictidx[1] - 1 ; i >= dictidx[0] && i >= dictidx[1] - maxsearch; i--)
							{
								auto rid = dict[l].read_id[i];
								if((forward_bitset^read[rid]).count()<=thresh_s)
								{	
									omp_set_lock(&read_lock[rid]);
									if(remainingreads[rid])
									{
										remainingreads[rid]=0;
										flag = 1;
									}
									omp_unset_lock(&read_lock[rid]);
								}
								if(flag == 1)//match found
								{
									flag = 0;
									list_size++;
									pos.insert(pos_it,j);
									reads.insert(reads_it,bitsettostring(read[rid]));
									order.insert(order_it,order_s[rid]);
									RC.insert(RC_it,'d');
									quality.insert(quality_it,quality_sN[rid]);
									for(int l1 = 0;l1 < numdict_s; l1++)
										deleted_rids[l1].push_back(rid);
								}
							}
						}
						omp_unset_lock(&dict_lock[startposidx]);
						//delete from dictionaries
						for(int l1= 0; l1 < numdict_s; l1++)
							for(auto it = deleted_rids[l1].begin(); it!=deleted_rids[l1].end();)
							{
								b = read[*it]&mask1[l1];
								ull = (b>>3*dict_start[l1]).to_ullong();
								startposidx = dict[l1].bphf->lookup(ull);
								if(!omp_test_lock(&dict_lock[startposidx]))
								{
									++it;
									continue;
								}
								dict[l1].findpos(dictidx,startposidx);
								dict[l1].remove(dictidx,startposidx,*it);
								it = deleted_rids[l1].erase(it);	
								omp_unset_lock(&dict_lock[startposidx]);
							}
					}
					for(int l = 0; l < numdict_s; l++)//reverse
					{
						b = reverse_bitset&mask1[l];
						ull = (b>>3*dict_start[l]).to_ullong();
						startposidx = dict[l].bphf->lookup(ull);
						if(startposidx >= dict[l].numkeys)//not found
							continue;
						//check if any other thread is modifying same dictpos
						if(!omp_test_lock(&dict_lock[startposidx]))
							continue;
						dict[l].findpos(dictidx,startposidx);
						if(dict[l].empty_bin[startposidx])//bin is empty
						{
							omp_unset_lock(&dict_lock[startposidx]);
							continue;
						}
						uint64_t ull1 = ((read[dict[l].read_id[dictidx[0]]]&mask1[l])>>3*dict_start[l]).to_ullong();
						if(ull == ull1)//checking if ull is actually the key for this bin
						{	
							for (int64_t i = dictidx[1] - 1 ; i >= dictidx[0] && i >= dictidx[1] - maxsearch; i--)
							{
								auto rid = dict[l].read_id[i];
								if((reverse_bitset^read[rid]).count()<=thresh_s)
								{	
									omp_set_lock(&read_lock[rid]);
									if(remainingreads[rid])
									{
										remainingreads[rid]=0;
										flag = 1;
									}
									omp_unset_lock(&read_lock[rid]);
								}
								if(flag == 1)//match found
								{
									flag = 0;
									list_size++;
									pos.insert(pos_it,j);
									reads.insert(reads_it,reverse_complement(bitsettostring(read[rid])));
									order.insert(order_it,order_s[rid]);
									RC.insert(RC_it,'r');
									std::reverse(quality_sN[rid].begin(),quality_sN[rid].end());
									quality.insert(quality_it,quality_sN[rid]);
									for(int l1 = 0;l1 < numdict_s; l1++)
										deleted_rids[l1].push_back(rid);
								}
							}
						}
						omp_unset_lock(&dict_lock[startposidx]);
						//delete from dictionaries
						for(int l1= 0; l1 < numdict_s; l1++)
							for(auto it = deleted_rids[l1].begin(); it!=deleted_rids[l1].end();)
							{
								b = read[*it]&mask1[l1];
								ull = (b>>3*dict_start[l1]).to_ullong();
								startposidx = dict[l1].bphf->lookup(ull);
								if(!omp_test_lock(&dict_lock[startposidx]))
								{
									++it;
									continue;
								}
								dict[l1].findpos(dictidx,startposidx);
								dict[l1].remove(dictidx,startposidx,*it);
								it = deleted_rids[l1].erase(it);	
								omp_unset_lock(&dict_lock[startposidx]);
							}
					}
					if(j != ref.size()-readlen)//not at last position,shift bitsets
					{
						forward_bitset >>= 3;
						forward_bitset |= basemask[readlen-1][ref[j+readlen]];
						reverse_bitset <<= 3;
						reverse_bitset |= basemask[0][chartorevchar[ref[j+readlen]]];
					}	
							
				}
				//convert pos to differences again
				long prevpos = 0;
				for(auto it = pos.begin(); it != pos.end(); ++it)
				{
					*it = *it - prevpos;
					prevpos = prevpos + *it;
				}
				//include singletons in count as well 
				//negligible effect on compression
				//hopefully, helps in denoising
				count = buildcontig(reads,pos,list_size);		
				writecontig(count,pos,reads,order,RC,quality,f_seq,f_pos,f_noise,f_noisepos,f_order,f_RC, f_order_N_pe, f_flag_N, f_quality, f_quality_N, list_size);
			}
			reads = {current};
			pos = {p};
			order = {ord};
			RC = {rc};
			if(rc == 'r')
				std::reverse(q.begin(),q.end());
			quality = {q};
			list_size = 1;
		}
		else
		{
			reads.push_back(current);
			pos.push_back(p);
			order.push_back(ord);
			RC.push_back(rc);
			if(rc == 'r')
				std::reverse(q.begin(),q.end());
			quality.push_back(q);
			list_size++;
		}
		i++;	
					
	}

	count = buildcontig(reads,pos,list_size);
	writecontig(count,pos,reads,order,RC,quality,f_seq,f_pos,f_noise,f_noisepos,f_order,f_RC,f_order_N_pe, f_flag_N, f_quality, f_quality_N, list_size);

	f.close();
	in_flag.close();
	in_pos.close();
	in_order.close();
	in_RC.close();
	in_quality.close();
	f_seq.close();
	f_pos.close();
	f_noise.close();
	f_noisepos.close();
	f_order.close();
	f_order_N_pe.close();
	f_RC.close();
	f_quality.close();
	f_quality_N.close();
	}

	//Combine files produced by the threads
	std::ofstream f_order(infile_order);
	std::ofstream f_order_N_pe(outfile_order_N_pe,std::ios::binary|std::ofstream::app);	
	std::ofstream f_meta(outfile_meta);
	std::ofstream f_quality(outfile_quality);
	std::ofstream f_quality_N(outfile_quality_N);

	for(int tid = 0; tid < num_thr; tid++)
	{
		std::ifstream in_order(infile_order+'.'+std::to_string(tid));
		std::ifstream in_order_N_pe(outfile_order_N_pe+'.'+std::to_string(tid));
		std::ifstream in_quality(outfile_quality+'.'+std::to_string(tid));
		std::ifstream in_quality_N(outfile_quality_N+'.'+std::to_string(tid));

		f_order << in_order.rdbuf();
		f_order_N_pe << in_order_N_pe.rdbuf();
		f_quality << in_quality.rdbuf();
		f_quality_N << in_quality_N.rdbuf();
		
		f_order.clear();//clear error flags if rdbuf empty file	
		f_order_N_pe.clear();//clear error flags if rdbuf empty file	
		f_quality.clear();
		f_quality_N.clear();
		
		remove((infile_order+'.'+std::to_string(tid)).c_str());
		remove((outfile_order_N_pe+'.'+std::to_string(tid)).c_str());
		remove((outfile_quality+'.'+std::to_string(tid)).c_str());
		remove((outfile_quality_N+'.'+std::to_string(tid)).c_str());
	}
	f_meta << readlen << "\n";
	f_meta.close();
	f_order.close();
	f_order_N_pe.close();
	f_quality.close();
	f_quality_N.close();
	//write remaining singleton reads and qualities now
	
	std::ofstream f_singleton(outfile_singleton);
	f_quality.open(outfile_quality,std::ofstream::app);
	f_order.open(infile_order,std::ios::binary|std::ofstream::app);
	f_order_N_pe.open(outfile_order_N_pe,std::ios::binary|std::ofstream::app);
	std::ofstream f_N(infile_N);
  	f_quality_N.open(outfile_quality_N,std::ofstream::app);
	char c = readlen;
	std::string s;
	uint32_t matched_s = numreads_s;
	for (uint32_t i = 0; i < numreads_s; i++)
		if(remainingreads[i] == 1)
		{
			matched_s--;
			f_order.write((char*)&order_s[i],sizeof(uint32_t));
			s = bitsettostring(read[i]);
			f_singleton << s;
			f_quality << quality_sN[i] <<"\n";
		}
	uint32_t matched_N = numreads_N;
	for (uint32_t i = numreads_s; i < numreads_s+numreads_N; i++)
		if(remainingreads[i] == 1)
		{
			matched_N--;
			f_N << bitsettostring(read[i]) << "\n";
			f_order_N_pe.write((char*)&order_s[i],sizeof(uint32_t));
			f_quality_N << quality_sN[i] <<"\n";
		}
	f_order.close();
	f_order_N_pe.close();
	f_N.close();
	f_singleton.close();
	f_quality.close();
	f_quality_N.close();
	delete[] remainingreads;
	packbits();
	std::cout << "Encoding done:\n"; 
	std::cout << matched_s << " singleton reads were aligned\n";
	std::cout << matched_N << " reads with N were aligned\n";
	return;
}

std::vector<std::array<long,5>> buildcontig(std::list<std::string> reads, std::list<long> pos, uint32_t list_size)

{
	auto reads_it = reads.begin();
	std::vector<std::array<long,5>> count(readlen,{0,0,0,0,0});
	for(long i = 0; i < readlen; i++)
		count[i][chartolong[(*reads_it)[i]]] = 1;
	if(list_size == 1)
		return count;

	long prevpos = 0,currentpos;
	auto pos_it = pos.begin();
	++reads_it;
	++pos_it;
	for(; pos_it != pos.end(); ++pos_it,++reads_it)
	{
		count.insert(count.end(),*pos_it,{0,0,0,0,0});
		currentpos = prevpos + *pos_it;
		for(long i = 0; i < readlen; i++)
			count[currentpos+i][chartolong[(*reads_it)[i]]] += 1;
		prevpos = currentpos;
	}
	//set counts of N to 0, because we always
//	for(auto it = count.begin(); it!= count.end(); ++it)
//		(*it)[4] = 0;
	return count;
}

std::string generateRef(std::vector<std::array<long,5>> count)
{
	std::string ref(count.size(),'A');
	for(long i = 0; i < count.size(); i++)
	{
		long max = 0,indmax = 0;
		for(long j = 0; j < 4; j++)
		//not including N because packbits doesn't expect seq to have N
			if(count[i][j]>max)
			{
				max = count[i][j];
				indmax = j;
			}
		ref[i] = longtochar[indmax];
	}
	return ref;
}

void writecontig(std::vector<std::array<long,5>> count, std::list<long> &pos, std::list<std::string> &reads, std::list<uint32_t> &order, std::list<char> &RC, std::list<std::string> &quality, std::ofstream& f_seq, std::ofstream& f_pos, std::ofstream& f_noise, std::ofstream& f_noisepos, std::ofstream& f_order, std::ofstream &f_RC, std::ofstream &f_order_N_pe, std::ofstream &f_flag_N, std::ofstream &f_quality, std::ofstream &f_quality_N, uint32_t list_size)
{
	std::string ref = generateRef(count);
	
	f_seq << ref;
	char c;
	if(list_size == 1)
	{
		f_noise << "\n";
		c = readlen;//(not pos[0] to handle breaks in read sequence due to limit on reads.size() - can't
				//assume  pos[0] = readlen)
		f_pos << c;
		f_order.write((char*)&order.front(),sizeof(uint32_t));
		f_RC << RC.front();
		if(RC.front() == 'r')
		{
			std::string s = quality.front();
			std::reverse(s.begin(),s.end());
			f_quality << s << "\n";	
		}
		else
			f_quality << quality.front() << "\n";
		return;
	}
	long prevj = 0;
	auto pos_it = pos.begin();
	auto reads_it = reads.begin();
	auto order_it = order.begin();
	auto RC_it = RC.begin();
	auto quality_it = quality.begin(); 
	long _read_char_id, _ref_char_id;
	long prevpos,currentpos;
	bool firstread = true;//to specially handle pos for first read in contig
	for(;pos_it!=pos.end(); ++pos_it,++reads_it,++order_it,++RC_it,++quality_it)
	{
		if(firstread == true)
			currentpos = 0;
		else
			currentpos = prevpos + *pos_it;
		long prevj = 0;
		for(long j = 0; j < readlen; j++)
		{
			_read_char_id = chartolong[(*reads_it)[j]];
			_ref_char_id = chartolong[ref[currentpos+j]];
			double alt_thresh_per_pos = alternate_thresh;
			bool allowed_alternate_flag = false;
			//double alt_thresh_per_pos = alternate_thresh*(0.5+ 0.5*(j*1.0/readlen));
			//if(*RC_it == 'r'){
			//    alt_thresh_per_pos = alternate_thresh*(0.5+ 0.5*((readlen-j-1)*1.0/readlen));
			//}
			//double alt_thresh_per_pos = alternate_thresh*(0.4+ 0.6*(j*1.0/readlen));
			// if( _read_char_id == 4)  //Temporary to positions with N
			// 	allowed_alternate_flag = true;
			// else
			for(int i = 0; i < 4; i++)
			{
				if(i != _ref_char_id)
					if(count[currentpos+j][i] >= (int)count[currentpos+j][_ref_char_id]*alt_thresh_per_pos)
						allowed_alternate_flag = true;
			}
			//allowed_alternate_flag = (count[currentpos+j][_read_char_id] >= (int)count[currentpos+j][_ref_char_id]*alt_thresh_per_pos);	
		
			if(((*reads_it)[j] != ref[currentpos+j]) && allowed_alternate_flag)
			{
				f_noise<<enc_noise[ref[currentpos+j]][(*reads_it)[j]];
				c = j-prevj;
				f_noisepos<<c;
				prevj = j;
			}
		}
		f_noise << "\n";
		if(firstread == true)
		{
			c = readlen;
			firstread = false;
		}
		else
			c = *pos_it;
		f_pos << c;
		if(*RC_it == 'r')
			std::reverse((*quality_it).begin(),(*quality_it).end());		
		if((*reads_it).find('N')!=std::string::npos)
		{
			if(preserve_order == "True")
				f_quality_N << *quality_it << "\n";
			else
				f_quality << *quality_it << "\n";	
			f_order_N_pe.write((char*)&(*order_it),sizeof(uint32_t));
			f_flag_N << '1';
		}
		else
		{
			f_quality << *quality_it << "\n";	
			f_order.write((char*)&(*order_it),sizeof(uint32_t));
			f_flag_N << '0';
		}

		
		f_RC << *RC_it;
		prevpos = currentpos;
	}
	return;
}

*/
