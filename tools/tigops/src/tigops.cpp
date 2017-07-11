// mainly used to compute kmer coverage of unitigs
// for more efficient tip removal you could also consider BTRIM from Malfoy's github

#include <tigops.hpp>

#include <unordered_map>
#include <numeric>
#include <sstream> // for setprecision
#include <iomanip> // for setprecision


using namespace std;

// need gcc 4.7 atleast here
//
template<typename T> 
using kmer_hash_t = unordered_map<string, T>;

typedef unordered_map<string, vector<string> > kmer_hash_str_t;

struct Insert_into_hash
{
	public:
		size_t _sizeKmer;
		kmer_hash_t<int> &kmer_hash;
		Insert_into_hash (size_t _sizeKmer, kmer_hash_t<int> &kmer_hash) : _sizeKmer(_sizeKmer), kmer_hash(kmer_hash) {}
		void operator() (const kmer_type &kmer) 
		{
			kmer_type normalized_kmer = std::min(kmer, revcomp(kmer, _sizeKmer));
			string normalized_kmer_string = normalized_kmer.toString(_sizeKmer);

            if (kmer_hash.find(normalized_kmer_string) == kmer_hash.end()) // don't insert/zero-initialize if it was already inserted
            {
			    kmer_hash[normalized_kmer_string] = 0;
            }
		}
};

struct Insert_tig_name_into_hash
{
	public:
		size_t _sizeKmer;
		kmer_hash_str_t &kmer_hash_str;
        string tig_name;
		Insert_tig_name_into_hash (size_t _sizeKmer, kmer_hash_str_t &kmer_hash_str) : _sizeKmer(_sizeKmer), kmer_hash_str(kmer_hash_str) {}
		void operator() (const kmer_type &kmer, int k = 0, bool normalized = false) 
		{
            if (k == 0)
                k = _sizeKmer;
            kmer_type normalized_kmer;
            if (normalized)
                normalized_kmer = kmer;
            else
			    normalized_kmer = std::min(kmer, revcomp(kmer, k));
            string normalized_kmer_str = normalized_kmer.toString(k);
            if (kmer_hash_str.find(normalized_kmer_str) == kmer_hash_str.end())
            {
                vector<string> v;
    			kmer_hash_str[normalized_kmer_str] = v;
            }
    		kmer_hash_str[normalized_kmer_str].push_back(tig_name);
		}
};


struct Increment_existing_hash
{
	public:
		size_t _sizeKmer;
		kmer_hash_t<int> &kmer_hash;
		Increment_existing_hash (size_t _sizeKmer, kmer_hash_t<int> &kmer_hash) : _sizeKmer(_sizeKmer), kmer_hash(kmer_hash) {}
		void operator() (const kmer_type &kmer) 
		{
			kmer_type normalized_kmer = std::min(kmer, revcomp(kmer, _sizeKmer));
			string normalized_kmer_string = normalized_kmer.toString(_sizeKmer);

			unordered_map<string,int>::const_iterator got = kmer_hash.find(normalized_kmer_string);
			if (got != kmer_hash.end())
				kmer_hash[got->first]=got->second+1;
		}
};

template<typename T>
struct Get_from_hash
{
	public:
		size_t _sizeKmer;
		kmer_hash_t<T> &kmer_hash;
		vector<T> values;
		Get_from_hash (size_t _sizeKmer, kmer_hash_t<T> &kmer_hash) : _sizeKmer(_sizeKmer), kmer_hash(kmer_hash) {}

		void start()
		{
			values.clear();
		}

		void operator() (const kmer_type &kmer) 
		{
			kmer_type normalized_kmer = std::min(kmer, revcomp(kmer, _sizeKmer));
			string normalized_kmer_string = normalized_kmer.toString(_sizeKmer);
			
            unordered_map<string,int>::const_iterator got = kmer_hash.find(normalized_kmer_string);
			if (got != kmer_hash.end())
			    values.push_back(got->second);
		}
};



template<typename T>
class Tigop {
public:
  int nb_passes, current_pass;
  kmer_hash_t<T> kmer_hash; // todo optimize to kmer_size kmers
  size_t _sizeKmer;
  Kmer<span>::ModelDirect model;
  Insert_into_hash insert_into_hash;
  Increment_existing_hash increment_existing_hash;
  Get_from_hash<T> get_from_hash;

  // maybe factorize
  kmer_hash_str_t kmer_hash_str; // todo optimize to kmer_size kmers
  Insert_tig_name_into_hash insert_tig_name_into_hash;


  Tigop(size_t _sizeKmer) : current_pass(0), _sizeKmer(_sizeKmer), model(_sizeKmer),
	insert_into_hash(_sizeKmer, kmer_hash), increment_existing_hash(_sizeKmer, kmer_hash), get_from_hash(_sizeKmer, kmer_hash),
    insert_tig_name_into_hash(_sizeKmer, kmer_hash_str)
	{}
  virtual void operator()(Sequence &seq, int seqlen) = 0; //{ printf("empty tigop called\n"); exit(1); }

    void next_pass()
    {
        current_pass++;
    }

    // kmer iteration stuff

    Kmer<span>::ModelDirect::Iterator *itKmer;
    kmer_type get_first_kmer(Sequence seq)
    {
        itKmer = new Kmer<span>::ModelDirect::Iterator (model);
        itKmer->setData (seq.getData());
        itKmer->first();
        return (*itKmer)->value();
    }

    // never call it without get_first_kmer
    kmer_type get_last_kmer()
    {
        for (itKmer->first(); !itKmer->isDone(); itKmer->next())
        {
            // just loop to last
        }
        return (*itKmer)->value();
    }
};


template<typename Functor>
void forAllKmers (Sequence& seq, int seqlen, size_t _sizeKmer, Functor& fct)  
{
		//int nbkmers = seqlen - _sizeKmer + 1;
		kmer_type kmer;

		// iterate over sequence kmers
		Kmer<span>::ModelDirect model (_sizeKmer); //ModelCanonical ModelDirect
		Kmer<span>::ModelDirect::Iterator itKmer (model);
		itKmer.setData (seq.getData());
		for (itKmer.first(); !itKmer.isDone(); itKmer.next())
		{
				kmer = itKmer->value();

				fct(kmer);
		}

}

double getMean( vector<int> &values)
{
    if (values.size() == 0) return 0;
	return accumulate(values.begin(), values.end(), 0.0) / values.size();
}

double getMedian( vector<int> values) // copy the vector
{ // https://stackoverflow.com/questions/2114797/compute-median-of-values-stored-in-vector-c
    if (values.size() == 0) return 0;
    typedef vector<int>::size_type vec_sz;
    vec_sz size = values.size();
    sort(values.begin(), values.end());
    vec_sz mid = size/2;
    return size % 2 == 0 ? (values[mid] + values[mid-1]) / 2 : values[mid];
}

double getSd( vector<int> &v)
{
    if (v.size() == 0) return 0;
    double mean = getMean(v);
    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0); // https://stackoverflow.com/questions/7616511/calculate-mean-and-standard-deviation-from-a-vector-of-samples-in-c-using-boos
    double stdev = std::sqrt(sq_sum / v.size() - mean * mean);
    return stdev;
}



template <typename T>
std::string to_string_p(const T a_value, const int n = 3)
{ // https://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}


string create_cov_id_header(float cov_mean = 0, string sample_name = "", float cov_median = 0, float cov_sd = 0, unsigned int cov_num = 0)
{
    string s;
    if (cov_median == 0)
        s = "_mean_" + to_string_p(cov_mean) + "_ID_" + sample_name; // legacy tigops
    else
        s = "_mean_" + to_string_p(cov_mean) + "_median_" + to_string_p(cov_median) + + "_sd_" + to_string_p(cov_sd) + "_nbkmers_" + to_string_p(cov_num) + "_ID_" + sample_name; // advanced tigops with median/number of kmers
    return s;
}


string reverse_complement(string s)
{
    string t;
    for (unsigned int i = 0; i < s.size(); i ++)
    {
        switch (s[i])
        {
            case 'A':
                t += "T"; break;
            case 'T':
                t += "A"; break;
            case 'G':
                t += "C"; break;
            case 'C':
                t += "G"; break;
            default:
                cout << s[i] << " unexpected in reverse_complement" << endl;
                exit(1);
        }
    }
    return string ( t.rbegin(), t.rend() );
}



class ComputeCoverage: public Tigop<int> {
	public:

		string sample_name;
		IBank* outputBank;

		ComputeCoverage(size_t _sizeKmer, string sample_name, string output_file) : Tigop(_sizeKmer),
		sample_name(sample_name)
	{
		nb_passes = 3;
		outputBank = new BankFasta (output_file);
	}

		virtual void operator()(Sequence &seq, int seqlen) 
		{
			if (current_pass == 0)
				pass0(seq, seqlen);
			else if (current_pass == 1)
				pass1(seq, seqlen);
			else if (current_pass == 2)
				pass2(seq, seqlen, sample_name);
		}

		void pass0(Sequence &seq, unsigned int seqlen)
		{
			if (seqlen < _sizeKmer)
				throw Exception("A tig is smaller than k-mer size");

			forAllKmers(seq, seqlen, _sizeKmer, insert_into_hash); // won't reinsert elements that already exist in the hash

            // handle cases where the k-mer size used to compute coverage is smaller than the k-mer size used to create the graph
			forAllKmers(seq, seqlen, _sizeKmer, increment_existing_hash); // record how many times the kmer was seen in unitigs
		}

        void filter_repeated_kmers()
        {
            // remove kmers that appear in multiple unitigs
            for (auto it = kmer_hash.begin(); it != kmer_hash.end();)
            {
                if (it->second > 1)
                {   
                    //std::cout << "skipping repeated kmer " << it->first << std::endl;
                    it = kmer_hash.erase(it);
                    continue;
                }
                else
                {
                    //std::cout << "kept kmer " << it->first <<  " " << it->second << std::endl;
                    it->second = 0; // reset counts for the next pass
                    it++;
                }
            }
        }

		void pass1(Sequence &seq, unsigned int seqlen)
		{
			forAllKmers(seq, seqlen, _sizeKmer, increment_existing_hash);
		}

		void pass2(Sequence &seq, unsigned int seqlen, string sample_name)
		{
			get_from_hash.start();	
			Sequence outseq (Data::ASCII);
			Data&  data	= outseq.getData();
			//data.resize(seq.getDataSize()); // unsure if needed
			data = seq.getData();

			forAllKmers(seq, seqlen, _sizeKmer, get_from_hash);
			vector<int> &values = get_from_hash.values;
#if 0
			printf("seq: %.*s\n",seq.getDataSize(),seq.getDataBuffer());
			for (vector<int>::iterator it = values.begin(); it != values.end(); it++)
			{
				printf("%d ", *it);
			}
            std::cout << "mean : " << getMean(values) << std::endl;
			printf("\n");
#endif

            if (values.size() == 0)
            {
                std::istringstream iss(seq.getComment());
                string unitig;
                iss >> unitig;
                std::cerr << std::endl << "No unique kmer found for unitig " <<  unitig << " (seq: " << seq.toString() << " length: " << seqlen<< "). Reporting coverage of 0. Use a larger k to avoid this." << std::endl;
            }
            
            double mean = getMean(values);
			double median = getMedian(values);
			double sd = getSd(values);
			stringstream ss;
			ss.precision(1);
            ss << fixed << seq.getComment() << create_cov_id_header(mean, sample_name, median, sd, values.size());
			outseq._comment = ss.str();
			outputBank->insert(outseq);
		}



};

class RemoveTips : public Tigop<int> { 
public:

		unsigned int max_tip_length, coverage_threshold;
		string sample_name;
		IBank* outputBank;
        Kmer<span>::ModelDirect model;

        RemoveTips(size_t _sizeKmer, string output_file) : Tigop(_sizeKmer), model(_sizeKmer)
		{
			nb_passes = 3;
	    		outputBank = new BankFasta (output_file);

            // will be overridden by parameters
            max_tip_length = 100;
			coverage_threshold = 3;
		}

		virtual void operator()(Sequence &seq, int seqlen) 
		{
				if (current_pass == 0)
					pass0(seq, seqlen);
                else if (current_pass == 1 && nb_passes == 3)
                        pass1(seq, seqlen);
                else if ((current_pass == 2 && nb_passes == 3) || (current_pass == 1 && nb_passes==2))
                        pass2(seq, seqlen);
		}


        void pass0(Sequence &seq, unsigned int seqlen)
        {
            if (seqlen < _sizeKmer)
                throw Exception("A tig is smaller than k-mer size");

            // special stuff for first and last kmers
            kmer_type kmer = get_first_kmer(seq);

            insert_into_hash(kmer);

            if (nb_passes == 2)
                increment_existing_hash(kmer); // simulate coverage by incrementing kmers we've seen


            // also if we're going to pass reads: mark incoming and all middle kmers one for coverage computation
            if (seqlen <= max_tip_length && nb_passes == 3)
            {
                // mark possible left neighbors to discover them
                model.iterateIncomingNeighbors(kmer, insert_into_hash);

                forAllKmers(seq, seqlen, _sizeKmer, insert_into_hash);
            }

            // loop until last kmer
            kmer = get_last_kmer();

            insert_into_hash(kmer);

            if (nb_passes == 2)
                increment_existing_hash(kmer); // simulate coverage by incrementing kmers we've seen

            if (seqlen <= max_tip_length && nb_passes == 3)
            {
                // mark possible right neighbors to discover them while passing reads
                model.iterateOutgoingNeighbors(kmer, insert_into_hash);
            }
        }
	
		void pass1(Sequence &seq, int seqlen)
		{
				// populate counts
				forAllKmers(seq, seqlen, _sizeKmer, increment_existing_hash);
		}

        double extract_coverage_information(string header)
        {
            return 0; // TODO
        }

		void pass2(Sequence &seq, unsigned int seqlen)
		{
				Sequence outseq (Data::ASCII);
				Data&  data	= outseq.getData();
				data = seq.getData();
				outseq._comment = seq.getComment();
                string header = seq.getComment();
                
				// is it a tip?

				bool is_short_seq = (seqlen <= max_tip_length);
				if (is_short_seq)
				{

                    double mean;
                    if (nb_passes == 3)
                    {
                        get_from_hash.start();	
                        forAllKmers(seq, seqlen, _sizeKmer, get_from_hash);
                        vector<int> &values = get_from_hash.values;
                        mean = getMean(values);
                    }
                    else
                        mean = extract_coverage_information(header); // will return 0 if there's no coverage information
						
	
                    if (mean <= coverage_threshold) // either no cov information (0) or small enough coverage for it being a tip
                    {
						bool has_left_neighbor = false, has_right_neighbor = false;

                        // special stuff for first and last kmers
                        kmer_type kmer = get_first_kmer(seq);

                        get_from_hash.start();
                        model.iterateIncomingNeighbors(kmer, get_from_hash);
                        vector<int> &values = get_from_hash.values;
                        has_left_neighbor = getMean(values) > 0;

                        // loop until last kmer
                        kmer = get_last_kmer();

                        get_from_hash.start();
                        model.iterateOutgoingNeighbors(kmer, get_from_hash);
                        values = get_from_hash.values;
                        has_right_neighbor = getMean(values) > 0;

                        if ((!has_left_neighbor) || (!has_right_neighbor))
                            return; // no left neighbor, or no right neighbor? that's a tip

                    }

                    // else annotate that non-tip short sequence with its coverage
					stringstream ss;
					ss.precision(1);
                    ss << fixed << seq.getComment() << create_cov_id_header(mean, "0");
					outseq._comment = ss.str();
				}

				// don't save tips; save everything else
				outputBank->insert(outseq);
		}

};

#if 0 // ONGOING
class Compress : public Tigop<int> {
    public:

        IBank* outputBank;

        Compress(size_t _sizeKmer, string output_file) : Tigop(_sizeKmer), model(_sizeKmer)
    {
        nb_passes = 2;
        outputBank = new BankFasta (output_file);

        Kmer<span>::ModelDirect model (_sizeKmer); //ModelCanonical ModelDirect

    }

        virtual void operator()(Sequence &seq, int seqlen)
        {
            if (current_pass == 0)
                pass0(seq, seqlen);
            else if (current_pass == 1)
                pass1(seq, seqlen);
        }

        void pass0(Sequence &seq, int seqlen)
        {
            if (seqlen < _sizeKmer)
                throw Exception("A tig is smaller than k-mer size");


            kmer_type kmer = get_first_kmer(seq);
            insert_into_hash(kmer);

            kmer = get_last_kmer();
            insert_into_hash(kmer);
        }

        void pass1(Sequence &seq, int seqlen)
        {
            Sequence outseq (Data::ASCII);
            Data&  data = outseq.getData();
            data = seq.getData();
            outseq._comment = seq.getComment();

            bool has_single_left_neighbor = false, has_single_right_neighbor = false;

            // special stuff for first and last kmers
            kmer_type kmer = get_first_kmer(seq);

            model.iterateIncomingNeighbors(kmer, increment_existing_hash); // simulate coverage by incrementing kmers we've seen

            get_from_hash.start();

            // special stuff for first and last kmers
            kmer_type kmer = get_first_kmer(seq);

            get_from_hash.start();
            model.iterateIncomingNeighbors(kmer, get_from_hash);
            has_single_left_neighbor = getMean(values) == 1;

            // loop until last kmer
            kmer = get_last_kmer();

            get_from_hash.start();
            model.iterateOutgoingNeighbors(kmer, get_from_hash);
            has_single_right_neighbor = getMean(values) == 1;

            if (has_single_left_neighbor)
            {
                // else annotate that non-tip short sequence with its coverage
                // else annotate that non-tip short sequence with its coverage
                stringstream ss;
                ss.precision(1);
                ss << fixed << seq.getComment() << create_cov_id_header(mean, "0");
                outseq._comment = ss.str();
            }

            // don't save tips; save everything else
            outputBank->insert(outseq);
        }
}
#endif

class Fasta2Fastg : public Tigop<int> { 
public:

		string sample_name;
		IBank* outputBank;

        bool KminusOneMerMethod;
        bool renameHeaders;
        int k;

        // extra info to keep renaming output ids
        unsigned long tig_counter;
        unordered_map<string,int> header_map;
        unordered_map<int,int> length_map;

		Fasta2Fastg(size_t _sizeKmer, string output_file) : Tigop(_sizeKmer)
		{
            KminusOneMerMethod = true;
            k = (KminusOneMerMethod) ? (_sizeKmer - 1) : _sizeKmer;

			nb_passes = 3;
            tig_counter = 1; // start at 1 to avoid a Bandage bug
	    	outputBank = new BankFasta (output_file);

            renameHeaders = true; // will be set in execute()
        }

        string revcomp_header(string header)
        {
            if (header.at(header.length()-1) == '\'')
                return header.substr(0,header.length()-1);
            return header + "'";
        }

		virtual void operator()(Sequence &seq, int seqlen) 
		{
				if (current_pass == 0)
					pass0(seq, seqlen);
				else if (current_pass == 1)
						pass1(seq, seqlen);
		}

		void pass0(Sequence &seq, int seqlen)
		{
//				if (seqlen < _sizeKmer)
//				    throw Exception("A tig is smaller than k-mer size");

                string header = seq.getComment();

                if (renameHeaders)
                {
                    // record header to associate it to a node ID later when saving
                    header_map[header] = tig_counter;
                    length_map[tig_counter] = seqlen;
                }

				// deal with first and last (k-1)-mers of that tig
				Kmer<span>::ModelDirect model (k);
				Kmer<span>::ModelDirect::Iterator itKmer (model);
				itKmer.setData (seq.getData());
				itKmer.first();

				kmer_type kmer = itKmer->value();

                if (KminusOneMerMethod)
                {
                    kmer_type revcomp_kmer = revcomp(kmer, k);
                    string rc_status = (kmer < revcomp_kmer) ? "r" : "f";
    		    	kmer_type normalized_kmer = std::min(kmer, revcomp(kmer, k));
                    insert_tig_name_into_hash.tig_name = rc_status + "L" + header;
                    insert_tig_name_into_hash(normalized_kmer, k, false);
                }
                else
                {
                    insert_tig_name_into_hash.tig_name = header + "'"; // consider that the current tig is revcomp and flag its outgoing neighbors
    				// for all left neighbors, insert the current tig name to 
    				model.iterateIncomingNeighbors(kmer, insert_tig_name_into_hash);
                }

				// loop until last kmer
				for (itKmer.first(); !itKmer.isDone(); itKmer.next())
					kmer = itKmer->value();

                if (KminusOneMerMethod)
                {
                    kmer_type revcomp_kmer = revcomp(kmer, k);
                    string rc_status = (kmer < revcomp_kmer) ? "r" : "f";
    		    	kmer_type normalized_kmer = std::min(kmer, revcomp(kmer, k));
                    insert_tig_name_into_hash.tig_name = rc_status + "R" + header;
                    insert_tig_name_into_hash(normalized_kmer, k, false);
                    //(normalized_kmer.toString(k).compare("AATGCCTGATGCGCTACGCTTATCAGGCCTACA") == 0 || normalized_kmer.toString(k).compare("TGTAGGCCTGATAAGCGTAGCGCATCAGGCATT") == 0)
                }
                else
                {
                    insert_tig_name_into_hash.tig_name = header;
                    // same for right kmers
                    model.iterateOutgoingNeighbors(kmer, insert_tig_name_into_hash);
                }
                
                tig_counter ++;
        }

        string maybeRenameHeader(string header, int seqlen = 0)
        {
            if (!renameHeaders)
                return header;
        
            bool rc = false;   
            if (header.substr(header.size()-1).compare("'") == 0)
            {
                header = header.substr(0,header.size()-1);
                rc = true;
            }

            int header_n = header_map[header];

            string rc_str = (rc? "'":"");
            return "NODE_" + std::to_string( header_n ) + "_length_" +  std::to_string( length_map[header_n] ) + "_cov_0_ID_" + std::to_string( header_map[header] ) + rc_str;
        }

		void pass1(Sequence &seq, int seqlen, bool first_call = true)
		{
				Sequence outseq (Data::ASCII);
				Data&  data	= outseq.getData();
				data = seq.getData();
				outseq._comment = seq.getComment();

                // special stuff for first and last kmers
                Kmer<span>::ModelDirect model (k); //ModelCanonical ModelDirect
                Kmer<span>::ModelDirect::Iterator itKmer (model);
                itKmer.setData (seq.getData());
                itKmer.first();

                kmer_type kmer = itKmer->value();
    			kmer_type normalized_kmer;
                string normalized_kmer_str;

                /* // no need to deal with incoming edges
                normalized_kmer = std::min(kmer, revcomp(kmer, k));
                normalized_kmer_str = normalized_kmer.toString(k);
                vector<string> incoming_tigs = kmer_hash_str[normalized_kmer_str];
                */

                // loop until last kmer
                for (itKmer.first(); !itKmer.isDone(); itKmer.next())
                    kmer = itKmer->value();

    			normalized_kmer = std::min(kmer, revcomp(kmer, k));
                bool rc = (normalized_kmer == kmer);
                normalized_kmer_str = normalized_kmer.toString(k);
                vector<string> outgoing_tigs = kmer_hash_str[normalized_kmer_str];
 
				stringstream ss;
				ss.precision(1);
				ss << fixed << maybeRenameHeader(seq.getComment(), seqlen);
              
                if (KminusOneMerMethod)
                {
                    int nb_neighbors = 0;
                    for (vector<string>::iterator itTigs = outgoing_tigs.begin(); itTigs != outgoing_tigs.end(); itTigs ++)
                    {
                        string node = *itTigs;
                        string nei_header = node.substr(2);
                        bool nei_rc = node.substr(0,1).compare("f");
                        bool nei_direction = node.substr(1,1).compare("L");

                        /* okay several cases:
                         *
                         * current node = xxxx[last kmer (denoted R), forward (denoted f)]
                         * 
                         * candidate neighbors = 
                         *  (1) [first kmer (L), forward (f)]xxxxx
                         *  (2) [first kmer (L), reverse (r)]xxxxx
                         *  (3) xxxxx[last kmer (R), forward (f)]
                         *  (4) xxxxx[last kmer (R), reverse (r)]
                         *
                         *  among these we want only cases (1) and (4)
                         *
                         *
                         *  similarly, if current node = xxxxx[last kmer (R), reverse (r)]
                         *
                         * candidate neighbors are the same case analysis as above
                         *
                         * among these we want only cases (2) and (3) 
                         *
                         * TODO: following megahit's contigs_to_fastg.cpp, it appears that this code could have been greatly simplified if I only indexed the  first kmer and the revcomp of the last kmer.
                         *
                        */

                        if (rc == 0)
                        {
                            if (nei_direction == 0 && nei_rc == 0)
                            {
                                //pass 
                            }
                            else
                            {
                                if (nei_direction != 0 && nei_rc != 0)
                                    nei_header += "'";
                                else
                                    continue;
                            }
                        }
 
                        if (rc == 1)
                        {
                            if (nei_direction == 0 && nei_rc != 0)
                            {
                                // keep
                            }
                            else
                            {
                                if (nei_direction != 0 && nei_rc == 0)
                                    nei_header += "'";
                                else
                                    continue;
                            }
                        }

                       if (nb_neighbors == 0)
                            ss << ":";
                        else 
                            ss << ",";
                        ss << maybeRenameHeader(nei_header) << "";
                        nb_neighbors++;
                    }

                }
                else
                {
                    // this mode is buggy still (see simple_test), this is why i coded k-1-mer method
                    // print outgoing tigs
                    for (vector<string>::iterator itTigs = outgoing_tigs.begin(); itTigs != outgoing_tigs.end(); itTigs ++)
                    {
                        if (itTigs == outgoing_tigs.begin())
                            ss << ":";
                        else
                            ss << ",";
                        ss << revcomp_header(*itTigs) << "";
                    }
                }

                ss  <<";";
			
                outseq._comment = ss.str();

				outputBank->insert(outseq);

                // output the reverse complement node also
                if (first_call)
                {
                    string seq_str = seq.getDataBuffer();
                    seq_str.resize(seq.getDataSize()); // databuffer contains more than just the sequence. hence that resize
                    char * reverse_str = strdup(reverse_complement(seq_str).c_str()); // that strdup is necessary but don't know if/when to free
                    Sequence reverse_seq( reverse_str ); 
                    reverse_seq.setComment(seq.getComment() + "'");
                    pass1(reverse_seq, seqlen, false);
                }
		}

};


template <typename T>
void stream_sequences(string sequences_file, Tigop<T> &tigop)
{
    printf("streaming sequences file %s\n", sequences_file.c_str());

    // We declare a Bank instance defined by a list of filenames
    IBank *b  = Bank::open(sequences_file);

    // We create a sequence iterator for the bank
    Iterator<Sequence>* itSeq = b->iterator();

    //  We create a sequence iterator that notifies listeners every N sequences
    SubjectIterator<Sequence>* iter = new SubjectIterator<Sequence>(itSeq, 1000);

    // We create some listener to be notified every N iterations and attach it to the iterator.
    // Note that we get an estimation of the number of sequences in the bank and give it to the
    // functor in order to compute a percentage. Since this is a crude estimation, it is likely
    // that the final percentage will not be exactly 100%
    iter->addObserver (new ProgressTimer (b->estimateNbItems(), "Iterating sequences"));

    for (iter->first(); !iter->isDone(); iter->next())
    {
        Sequence &seq = iter->item();
		int seqlen = seq.getDataSize();
		if (seqlen == 0) return;

		tigop(seq, seqlen);
    }

    delete iter;
    delete b;
}

/********************************************************************************/

tigops::tigops ()  : Tool ("tigops"), _debug(false)
{
    getParser()->push_front (new OptionOneParam (STR_KMER_SIZE, "kmer size", true));
    getParser()->push_front (new OptionOneParam ("-tigs", "[uni/con]tig file", true));
    getParser()->push_front (new OptionOneParam ("-reads", "reads dataset", false));
    getParser()->push_front (new OptionOneParam ("-h5", "h5 file with counted kmers", false));
    getParser()->push_front (new OptionOneParam ("-name", "give an id (e.g. 1) to the dataset (for compute-coverage)", false));
    getParser()->push_front (new OptionOneParam ("-out", "output tigs file name", true));
    getParser()->push_front (new OptionNoParam ("-rename", "name fastg output ids according to SPAdes/Velvet NODE_.. format", false));
    getParser()->push_front (new OptionNoParam ("-debug", "debugging", false));

    // We add some custom arguments for command line interface
//    getParser()->push_front (new OptionNoParam ("remove-tips", "remove tips", false)); // not ready 
    getParser()->push_front (new OptionNoParam ("coverage", "compute mean kmer coverage per tig", false));
    getParser()->push_front (new OptionNoParam ("fasta2fastg", "create a fastg graph from a fasta file where tigs overlap by exactly (k-1) nucleotides", false));


}

void tigops::execute ()
{
    // We can do here anything we want.

    // We gather some statistics.
    getInfo()->add (1, "input");
    getInfo()->add (1, &LibraryInfo::getInfo());

    /** we get the kmer size chosen by the end user. */
    _sizeKmer = getInput()->getInt (STR_KMER_SIZE);
 
 
    string tigs_file = getInput()->getStr("-tigs");
    string output_file = getInput()->getStr("-out");
    _debug = getParser()->saw("-debug");



    if (getParser()->saw("coverage"))
    {
        string reads_file = getInput()->getStr("-reads");
  	    string sample_name = getInput()->getStr("-name");

	    ComputeCoverage computeCoverage(_sizeKmer, sample_name, output_file);

	    stream_sequences(tigs_file, computeCoverage);
        computeCoverage.filter_repeated_kmers();
	    computeCoverage.next_pass(); 
	    stream_sequences(reads_file, computeCoverage);
	    computeCoverage.next_pass(); 
	    stream_sequences(tigs_file, computeCoverage);
    }

    if (getParser()->saw("fasta2fastg"))
    {

	    Fasta2Fastg fasta2fastg(_sizeKmer, output_file);
        fasta2fastg.renameHeaders = getParser()->saw("-rename");

	    stream_sequences(tigs_file, fasta2fastg);
	    fasta2fastg.next_pass(); 
	    stream_sequences(tigs_file, fasta2fastg);
    }

    if (getParser()->saw("remove-tips"))
    {
        string reads_file = getInput()->getStr("-reads");
	    RemoveTips removeTips(_sizeKmer, output_file);

	    stream_sequences(tigs_file, removeTips);
	    removeTips.next_pass(); 
	    stream_sequences(reads_file, removeTips);
        removeTips.next_pass(); 
	    stream_sequences(tigs_file, removeTips);
    }



}
