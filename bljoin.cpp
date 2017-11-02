/*
 * Unix join command, but for biological sequences.
 *
 * This program is used to join together records with matching names
 * in different files (and actually in the same file too). Useful for
 * concatenating gene sequences.
 *
 * Assumes just one of each unique identifier in every file - otherwise
 * there are problems.
 *
 * Also assumes that every file is an alignment, so that it only checks
 * the size of the first record in each file. However, the ends of the
 * sequences are padded to even up the length.
 *
 */

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <seqan/seq_io.h>

#include <tclap/CmdLine.h>

#include <SeqFileInWrapper.h>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::ifstream;
using std::istream;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

using namespace seqan;
using namespace bltools;

vector<string> split(string s, string delim) {
    vector<string> ret;
    string buffer = "";
    bool delim_already_seen = false;
    string current_token = "";
    for(const char& c: s) {
        for(const char& d: delim) {
            if(c == d) {
                if(!delim_already_seen) {
                    ret.push_back(current_token);
                    current_token = "";
                }
                delim_already_seen = true;
            } else {
                delim_already_seen = false;
                current_token += c;
            }
        }
    }
    if(!delim_already_seen) {
        ret.push_back(current_token);
    }
    return ret;
}

int main(int argc, char * argv[]) {
  
  TCLAP::CmdLine cmd("Equivalent of `join' for sequence files", ' ', "0.0");
  TCLAP::SwitchArg no_pad_arg("n", "no-pad",
                              "Do not attempt to pad sequences to the same length",
                              cmd);
  TCLAP::SwitchArg allow_dups_arg("D", "allow-duplicates",
                                  "Allow duplicate names within a file", cmd);
  TCLAP::SwitchArg ignore_case_arg("i", "ignore-case",
                                   "Ignore case when matching", cmd);
  TCLAP::ValueArg<unsigned> field_arg("f", "field",
                                 "Field (1-based) of ID to join on (after splitting); will ignore sequences for which the field is empty; default is to use whole ID",
                                 false, 0, "int", cmd);
  TCLAP::ValueArg<string> pad_char_arg("p", "pad-char",
                                     "Character to use for padding",
                                     false, "-", "char", cmd);
  TCLAP::ValueArg<string> delim_arg("d", "delim",
                                 "Field separator", false, " ", "int", cmd);
  TCLAP::ValueArg<string> separator_arg("s", "separator",
                                        "Separator between joined sequences",
                                        false, "", "string", cmd);
  TCLAP::UnlabeledMultiArg<string> files("FILE(s)", "filenames", false,
                                         "file name(s)", cmd, false);
  cmd.parse(argc, argv);
  bool ignore_case = ignore_case_arg.getValue();
  bool no_pad = no_pad_arg.getValue();
  bool allow_dups = allow_dups_arg.getValue();
  unsigned field = field_arg.getValue();
  string delim = delim_arg.getValue();
  string pad_char = pad_char_arg.getValue();
  string separator = separator_arg.getValue();
  vector<string> infiles = files.getValue();
  if(infiles.size() == 0) infiles.push_back("-");

  CharString id;
  CharString seq_;              // CharString more flexible than Dna5String
  string seq;
  SeqFileInWrapper seq_handle;
  map<string, string> seqs;       // seqs[ID] = joined sequence string
  std::set<string> seqs_in_file;  // for each file, keep track of IDs seen
  unsigned long total_bases = 0;  // total length of joined sequences
  unsigned long seq_size = 0;     // size of the first record in each file
  unsigned long nfiles = 0;       // number of files processed
  vector<unsigned long> seq_lengths;

  for(string& infile: infiles) {
    
    seq_size = 0;
    seqs_in_file.clear();

    try {
        seq_handle.open(infile);
    } catch(Exception const &e) {
      cerr << "Could not open " << infile << endl;
      seq_handle.close();
      return 1;
    }
    
    while(!seq_handle.atEnd()) {

      try {

        readRecord(id, seq_, seq_handle.sqh);
        //nrecs_read++; 

      } catch (Exception const &e) {

        cerr << "Error: " << e.what() << endl;
        seq_handle.close();
        return 1;

      } // End try-catch for record reading.
      
      // Check the size of the first sequence in the file
      if(seq_size == 0) {
        seq_size = length(seq_);
        seq_lengths.push_back(seq_size);
      }

      if(seq_size != length(seq_)) {
        cerr << "Warning " << id << " is not the same size as other seqs"
             << " in the same file " << seq_size << " " <<
            length(seq_) << endl;
      }
      
      seq = string(toCString(seq_));
      
      // Simple method: just use the whole sequence ID to join
      string join_id = string(toCString(id));

      // Fancier: use parts of the ID
      if(ignore_case) {
        // Convert to upper case
        for(unsigned jii = 0; jii < join_id.size(); jii++) {
            join_id.at(jii) = toupper(join_id.at(jii));
        }
      }
      if(field > 0) {
        vector<string> splitted = split(join_id, delim);
        if(splitted.size() < field) {
          // If the field is not found in this record, skip it
          continue;
        }
        join_id = splitted[field-1];
      }


      // Check if this ID has been processed yet
      std::set<string>::iterator its = seqs_in_file.find(join_id);
      if(its != seqs_in_file.end() && !allow_dups) {
          cerr << join_id << " found more than once in " << infile << endl;
          throw("Duplicated ID");
      }
      seqs_in_file.insert(join_id);

      // Test if id is in map, and add if not
      // Also add enough sequence to fill in missed sequences
      map<string, string>::iterator it = seqs.find(join_id);
      if(it != seqs.end()) {
        // Found: do nothing except add padding
        if(nfiles > 0) {
            seqs[join_id] += separator;
        }
      } else {
        if(total_bases > 0 && !no_pad) {
          seqs[join_id].reserve(total_bases + seq_size * 2);
          for(unsigned long si = 0; si < seq_lengths.size()-1; si++) {
            unsigned long sl = seq_lengths[si];
            for(unsigned long j = 0; j < sl; j++) {
              seqs[join_id] += pad_char;
            }
            seqs[join_id] += separator;
          }
        } else {
          seqs[join_id] = "";
        }
      } // End test for existence of ID
      
      // Add the current sequence
      seqs[join_id] += seq;
      
    } // End single file reading loop
    
    total_bases += seq_size;
    
    // Add padding to IDs found in previous files but not this one
    if(!no_pad) {
      unsigned long target_length = total_bases + separator.size() * nfiles;
      for(map<string, string>::iterator it = seqs.begin(); it != seqs.end(); it++) {
        if(it->second.length() < target_length) {
          if(nfiles > 0) {
            it->second += separator;
          }
          for(unsigned long j = it->second.length(); j < target_length; j++) {
            it->second += pad_char;
          }
        }
      }
    }

    if(!seq_handle.close()) {
        cerr << "Problem closing " << infile << endl;
        return 1;
    }

    nfiles++;

  } // End loop over files
  
  // Write the output in fasta format
  for(pair<string, string> item: seqs) {
    cout << ">" << item.first << endl;
    cout << item.second << endl;
  }

  return 0;
}
