#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <iterator>
#include <unordered_map>
#include <utility>
#include "../header/utils.h"

namespace Bio {

    static const std::unordered_map<std::string, char> STANDARD_GENETIC_CODE = {
        {"UUU", 'F'}, {"UUC", 'F'}, {"UUA", 'L'}, {"UUG", 'L'},
        {"UCU", 'S'}, {"UCC", 'S'}, {"UCA", 'S'}, {"UCG", 'S'},
        {"UAU", 'Y'}, {"UAC", 'Y'}, {"UAA", '*'}, {"UAG", '*'}, // Stop codons
        {"UGU", 'C'}, {"UGC", 'C'}, {"UGA", '*'}, {"UGG", 'W'}, // Stop codon

        {"CUU", 'L'}, {"CUC", 'L'}, {"CUA", 'L'}, {"CUG", 'L'},
        {"CCU", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
        {"CAU", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
        {"CGU", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},

        {"AUU", 'I'}, {"AUC", 'I'}, {"AUA", 'I'}, {"AUG", 'M'}, // Start codon (Met)
        {"ACU", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
        {"AAU", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
        {"AGU", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},

        {"GUU", 'V'}, {"GUC", 'V'}, {"GUA", 'V'}, {"GUG", 'V'},
        {"GCU", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
        {"GAU", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
        {"GGU", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
    };

    enum class Format {
                Genbank,
                FASTA
    };

    enum class Alphabet {
        DNA,
        MRNA,
        Protein,
        Unknown
    };

    class Seq {
        private:
            std::string m_Data;    
            Alphabet m_Alphabet;

            void check_valid(const std::string& data, Alphabet alphabet = Alphabet::DNA) {
                /*
                    Check validity of string inputted by user. 
                    If not valid, throw exception
                    If valid, return true
                */

                std::string new_data = utils::strings::toLower(data);
                switch (alphabet)
                {

                case Alphabet::DNA : {
                    const std::string valid_dna_chars = "acgt";
                    for (char c : new_data) {
                        if (valid_dna_chars.find(c) == std::string::npos) {
                            throw std::invalid_argument("Invalid character '" + std::string(1, c) + "' found in DNA sequence.");
                        }
                    }
                    break;
                }

                case Alphabet::MRNA : {
                    const std::string valid_mrna_chars = "acgu";
                    for (char c : new_data) {
                        if (valid_mrna_chars.find(c) == std::string::npos) {
                            throw std::invalid_argument("Invalid character '" + std::string(1, c) + "' found in mRNA sequence.");
                        }
                    }
                    break;
                }

                case Alphabet::Protein : {
                    const std::string validProteinChars = "arndcqegilhkmfpsywv";
                    for (char c : new_data) { 
                        if (validProteinChars.find(c) == std::string::npos) {
                            throw std::invalid_argument("Invalid character '" + std::string(1, c) + "' found in Protein sequence.");
                        }
                    }
                    break;
                }

                case Alphabet::Unknown : {
                    throw std::invalid_argument("Validation for 'Unknown' alphabet is not supported");
                }
                default:
                    throw std::logic_error("Unhandled alphabet type in check_valid");
                    break;
                }
            }

        public:

        /*
            TODO:
                Need to implement:
                    Transcribe
                    Translate
                    Reverse Compliment
        */
            explicit Seq(const std::string& data, Alphabet alphabet = Alphabet::DNA)
                : m_Data(data), m_Alphabet(alphabet) {
                    if (m_Data.empty()) {
                        throw std::invalid_argument("Sequence data cannot be empty!");
                    }
                    check_valid(m_Data, m_Alphabet);
                }

            explicit Seq(const std::vector<unsigned char>& data, Alphabet alphabet = Alphabet::DNA) 
                : m_Data(utils::strings::convertBytesToStr(data)), m_Alphabet(alphabet) {
                    if (data.empty()) {
                        throw std::invalid_argument("Sequence data cannot be empty!");
                    }
                    check_valid(m_Data, m_Alphabet);
                }

            Seq transcribe() const {
                if (m_Alphabet != Alphabet::DNA) {
                    throw std::runtime_error("Transcription is only applicable to DNA sequences!");
                }

                std::string upper_mrna_data = utils::strings::toUpper(m_Data);
                std::string mrna_data = utils::strings::replace(upper_mrna_data, 'T', 'U');

                return Seq(mrna_data, Alphabet::MRNA);
            }

            Seq translate() const {
                if (m_Alphabet != Alphabet::MRNA && m_Alphabet != Alphabet::DNA) {
                    throw std::runtime_error("Translation is only applicable to DNA and mRNA sequences!");
                }

                std::string mrna_data = m_Data;
                std::string upper_mrna_Data = utils::strings::toUpper(mrna_data);

                if (m_Alphabet == Alphabet::DNA) {
                    mrna_data = utils::strings::replace(mrna_data, 'T', 'U');
                }

                std::string protein_data;
                protein_data.reserve(upper_mrna_Data.length() / 3);

                for (size_t i = 0; i + 2 < upper_mrna_Data.length(); i += 3) {
                    std::string codon = upper_mrna_Data.substr(i, 3);

                    auto it = STANDARD_GENETIC_CODE.find(codon);
                    if (it != STANDARD_GENETIC_CODE.end()) {
                        protein_data += it->second;
                    } else {
                        protein_data += 'X';
                        std::cerr << "Warning:\n\tUnknown Codon found: " << codon << std::endl;
                    }
                }

                return Seq(protein_data, Alphabet::Protein);
            }

            Seq reverse_complement() const {
                if (m_Alphabet != Alphabet::DNA) {
                    throw std::runtime_error("Reverse complementing requires a DNA sequence!");
                }
                
                std::string upper_mdata = utils::strings::toUpper(m_Data);
                std::string reversed_data = upper_mdata;
                std::reverse(reversed_data.begin(), reversed_data.end());

                for (char& base : reversed_data) {
                    switch (base) {
                        case 'A': base = 'T'; break;
                        case 'T': base = 'A'; break;
                        case 'C': base = 'G'; break;
                        case 'G': base = 'C'; break;
                        default: break;
                    }
                }

                return Seq(reversed_data, Alphabet::DNA);
            }

            void print() const {
                std::cout << m_Data << "\n";
            }

            size_t len() const {
                return m_Data.length();
            }

            Seq slice(size_t start, size_t end) const {
                if (start > this->m_Data.length()) {
                    throw std::out_of_range("Start index cannot be more than the length of the sequence!");
                }

                if (end > m_Data.length()) {
                    end = m_Data.length();
                }

                if (start > end) {
                    start = end;
                }

                size_t length = end - start;
                std::string slicedVal = m_Data.substr(start, length);
                return Seq(slicedVal, this->m_Alphabet);
            }

            int count(const char& nucleotide) const {
                int counter = 0;
                for (char c : m_Data) {
                    if (c == nucleotide) counter++;
                }

                return counter;
            }

            int find(const std::string& subString) const {
                size_t pos = m_Data.find(subString);
                return (pos != std::string::npos) ? static_cast<int>(pos) : -1;
            }

            int rfind(const std::string& subString) const {
                size_t pos = m_Data.rfind(subString);
                return (pos != std::string::npos) ? static_cast<int>(pos) : -1;
            }

            size_t index(const std::string& subString) const {
                size_t pos = m_Data.find(subString);
                if (pos == std::string::npos) {
                    throw std::out_of_range("Substring '" + subString + "' not found in sequence.");
                }
                return pos;
            }

            size_t rindex(const std::string& subString) const {
                size_t pos = m_Data.rfind(subString);
                if (pos == std::string::npos) {
                    throw std::out_of_range("Substring '" + subString + "' not found in sequence.");
                }
                return pos;
            }

            std::vector<std::pair<size_t, std::string>> search(const std::vector<std::string>& subStrings) const {
                std::vector<std::pair<size_t, std::string>> res;

                if (subStrings.empty()) {
                    return res; // assuming empty results
                }

                size_t current_pos = 0;
                while (current_pos < m_Data.length()) {
                    bool match_found_at_pos = false;
                    size_t best_match_len = 0;
                    std::string found_sub = "";

                    for (const auto& sub : subStrings) {
                        if (sub.empty()) continue;
                        if (current_pos + sub.length() <= m_Data.length()) {
                             if (m_Data.compare(current_pos, sub.length(), sub) == 0) {
                                best_match_len = sub.length();
                                found_sub = sub;
                                match_found_at_pos = true;
                                break;
                            }
                        }
                    }

                    if (match_found_at_pos) {
                        res.push_back({current_pos, found_sub});
                        current_pos += best_match_len;
                    } else {
                        current_pos++;
                    }
                }

                return res;

            }

            // helper functions
            Seq toLower() const {
                return Seq(utils::strings::toLower(this->m_Data), this->m_Alphabet);
            }
            
            Seq toUpper() const {
                return Seq(utils::strings::toUpper(this->m_Data), this->m_Alphabet);
            }

            // string rep => basically `str(Seq);`
            const std::string& str() const { return m_Data; }
            
            // iterator stuff
            std::string::const_iterator begin() const { return m_Data.begin(); }
            std::string::const_iterator end() const { return m_Data.end(); }

            // overloaded operators just to make it very similar enough to the python version (hopefully) 
            bool operator==(const Seq& other) const {
                return this->m_Data == other.m_Data;
            }

            bool operator!=(const Seq& other) const {
                return this->m_Data != other.m_Data;
            }

            Seq operator+(const Seq& other) const {
                if (m_Alphabet != other.m_Alphabet) {
                    throw std::invalid_argument("To concatinate 2 Sequences together, you would need to have them be of the same type of Sequence!");
                }
                return Seq(this->m_Data + other.m_Data, m_Alphabet);
            }

            char operator[](size_t index) const {
                if (index >= this->m_Data.length()) {
                    throw std::out_of_range("Index is out of Sequence bounds!");    
                }
                return this->m_Data[index];
            }
    };

    class SeqRecord {
        private:

            std::unordered_map<std::string, Format> m_FormatMapping = {
                { "genbank", Format::Genbank },
                { "gb", Format::Genbank },
                { "fasta", Format::FASTA }
            };
        
            Seq m_Seq;
            std::string m_Id;
            std::string m_Name;
            std::string m_Description;
            std::string m_Annotation;
            std::string m_Subfeatures;
            Format m_Format;

            Format assignFormat(const std::string& format) {
                std::string lower_format = utils::strings::toLower(format);
                auto it = m_FormatMapping.find(lower_format); 
                if (it == m_FormatMapping.end()) {
                    throw std::invalid_argument("Format is invalid! Please enter either 'genbank' or 'fasta'");
                }
                return it->second;
            }

        public:
            explicit SeqRecord(Seq seqObject, const std::string& id, const std::string& name = "", const std::string& description = "", const std::string& format = "") 
                :   m_Seq(std::move(seqObject)),
                    m_Id(id),
                    m_Name(name),
                    m_Description(description),
                    m_Format(assignFormat(format))
            {
                /* empty :3 */
            }

            // getters, since the vars are private. Assuming other functions will need access to these vars
            const std::string& getId() const { return m_Id; }
            const std::string& getName() const { return m_Name; }
            const std::string& getDescription() const { return m_Description; }
            const Seq& getSeq() const { return m_Seq; }

            void print() const {
                std::cout << "ID: " << m_Id << "\n";
                std::cout << "Name: " << m_Name << "\n";
                std::cout << "Description: " << m_Description << "\n";
                std::cout << "Sequence: ";
                m_Seq.print();
            }

            std::string format(const std::string& format_name) const {
                std::string lower_format = utils::strings::toLower(format_name);
                if (lower_format == "fasta") {
                    std::string formatted_str = ">";
                    formatted_str += m_Id;
                    if (!m_Description.empty()) {
                        formatted_str += " " + m_Description;
                    }
                    formatted_str += "\n";
                    formatted_str += m_Seq.str();
                    return formatted_str;
                } else if (lower_format == "genbank") {
                    // TODO: Implement Genbank formatting (much more complex!)
                    throw std::runtime_error("Genbank formatting not yet implemented.");
                } else {
                    throw std::invalid_argument("Unsupported format for SeqRecord::format: " + format_name);
                }
            }

            void setFormat(std::string seqFormat) {
                m_Format = assignFormat(seqFormat);
            }

            void dir();
            
    };

    class SeqUtils {
        private:
            std::unordered_map<char, float> m_gcValues {
                {'G', 1.0f},
                {'C', 1.0f},
                {'A', 0.0f},
                {'T', 0.0f},
                {'U', 0.0f},
                {'S', 1.0f},
                {'W', 0.0f},
                {'M', 0.5f},
                {'R', 0.5f},
                {'Y', 0.5f},
                {'K', 0.5f},
                {'V', (2.0f / 3.0f)},
                {'B', (2.0f / 3.0f)},
                {'H', (1.0f / 3.0f)},
                {'D', (1.0f / 3.0f)},
                {'B', 0.5f},
                {'N', 0.5f}
            };


        public:

        /*
            We have:
                GC utils:
                    gc_fraction
                    GC123
                    GC
                    GC_skew
                    XGC_skew

                Protein Stuff:
                    seq3 // no clue what this is
                    seq1 // no clue what this is
                Molecular Weight
                Six Frame Translation
        */
    };

     class SeqIO {
        public:
            /*
                TODO:
                    - Write => Sequences(A list or iter of SeqRecord objects OR a single SeqRecord object), handle(file handle / filename as string), Format(lowercase)
                    - Parse => Filename
                    - Read => 

            */
    };

    namespace Align {
        class MultipleSeqAlignment {
            private:
                std::vector<Bio::SeqRecord> m_Records;
            public: 
                MultipleSeqAlignment(const std::vector<Bio::SeqRecord>& records) 
                    : m_Records(records) 
                    {
                        if (!m_Records.empty()) {
                            size_t expected_len = m_Records[0].getSeq().len();

                            for (size_t i = 0; i < m_Records.size(); i++) {
                                if (m_Records[i].getSeq().len() != expected_len) {
                                    throw std::invalid_argument("Sequences in a multiple alignment must have the same length!");
                                }
                            }   
                        }
                    }

                void print() const {
                    for (const auto& record : m_Records) {
                        record.print();
                    }
                }

                const Bio::SeqRecord& operator[](size_t index) const {
                    if (index >= m_Records.size()) {
                        throw std::out_of_range("Alignment index out of bounds");
                    }

                    return m_Records[index];
                }

                std::vector<Bio::SeqRecord>::const_iterator begin() const { return m_Records.begin(); }
                std::vector<Bio::SeqRecord>::const_iterator end() const { return m_Records.end(); }

                size_t size() const { return m_Records.size(); }
                size_t get_alignment_length() const {
                    if (m_Records.empty()) return 0;
                    return m_Records[0].getSeq().len();
                }
        };
    }

    class AlignIO {
        private:

        public:

    };


    namespace pairwise2
    {
        class align {
            // it has a function called alignment function (no clue)
            
            // WHEREIS GLOBAL/LOCAL XX MX AND MS DECLARED???
            // wtf is this bro
            /*
                 def __init__(self, name):
                    """Check to make sure the name of the function is reasonable."""
                    if name.startswith("global"):
                        if len(name) != 8:
                            raise AttributeError("function should be globalXX")
                    elif name.startswith("local"):
                        if len(name) != 7:
                            raise AttributeError("function should be localXX")
                    else:
                        raise AttributeError(name)
                    align_type, match_type, penalty_type = name[:-2], name[-2], name[-1]
                    try:
                        match_args, match_doc = self.match2args[match_type]
                    except KeyError:
                        raise AttributeError(f"unknown match type {match_type!r}")
                    try:
                        penalty_args, penalty_doc = self.penalty2args[penalty_type]
                    except KeyError:
                        raise AttributeError(f"unknown penalty type {penalty_type!r}")

            */
        };
        // some function called format_alignment    
    }
    
   
}
