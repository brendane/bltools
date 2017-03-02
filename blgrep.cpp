/*
 * Seqan-based program for grepping sequence files. Right now it uses
 * a regex approach to all the matching, but I would like to implement
 * a biological matching approach too.
 *
 * It refuses to read genbank files with a ".gb" extension unless they are
 * piped in. I can't seem to convince SeqFileIn that genbank is the format.
 *
 * TODO:
 *    - Output format based on input format, or at least some choices
 *    - Ability to search within a range of sequence coordinates
 *    - Ability to spit out context around matches or the number of matches
 *    - Ability to split a fasta file into records using a GFF file -
 *      or maybe that should be a different tool?
 *    - Translate search
 *    - Biologically-based search instead of regex
 *
 */

#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include <seqan/find.h>
#include <seqan/seq_io.h>

#include <tclap/CmdLine.h>

using std::cout;
using std::cerr;
using std::endl;
using std::regex;
using std::regex_match;
using std::string;
using std::vector;

using namespace seqan;

int main(int argc, char * argv[]) {
  
  // Regex flags
  std::regex_constants::syntax_option_type regex_flags =
    regex::extended | regex::optimize;
  std::regex_constants::match_flag_type regex_match_flags =
    std::regex_constants::match_any | std::regex_constants::match_not_null;
  
  /*
   * Command line arguments.
   */
  TCLAP::CmdLine cmd("Program to regex search sequence files", ' ', "0.0");
  TCLAP::SwitchArg seq_regex_arg("S", "sequence-regex",
				 "use regex for sequences instead of name; sets -i",
				 cmd);
  TCLAP::SwitchArg invert_regex_arg("v", "invert-match",
				    "Invert matching, like grep -v", cmd);
  TCLAP::SwitchArg ignore_case_arg("i", "ignore-case",
				   "Ignore case in pattern and input", cmd);
  TCLAP::SwitchArg case_sensitive_arg("I", "case-sensitive",
				      "Do not ignore case in pattern and input; only for -S",
				      cmd);
  TCLAP::ValueArg<string> match_type_arg("M", "match-type",
				     "Match type: f=fwd, r=rev, c=compl., R=revcomp; ignored for names",
				     false, "f", "string", cmd);
  TCLAP::UnlabeledValueArg<string> regex_string_arg("PATTERN", "regex pattern",
						    true, "",
						    "regex", cmd, false);
  TCLAP::UnlabeledMultiArg<string> infile_name("FILE(s)", "input file(s) use '-' for stdin or leave blank",
					       false, "file name(s)", cmd, false);
  cmd.parse(argc, argv);
  vector<string> infiles = infile_name.getValue();
  if(infiles.size() == 0) infiles.push_back("-");
  bool seq_regex = seq_regex_arg.getValue();
  bool inverted = invert_regex_arg.getValue();
  string match_type = match_type_arg.getValue();
  if(ignore_case_arg.getValue() || (seq_regex && !case_sensitive_arg.getValue()))
    regex_flags |= regex::icase;
  regex regex_pattern(regex_string_arg.getValue(), regex_flags);
  
  CharString id;
  CharString seq;              // CharString more flexible than Dna5String
  CharString qual;
  SeqFileIn seq_handle;
  
  // Loop over input files
  int nmatched = 0;
  for(string infile: infiles) {

    if(infile == "-") {
      if(!open(seq_handle, std::cin)) {
	cerr << "Could not read file from stdin" << endl;
	return 1;
      }
    } else {
      if(!open(seq_handle, infile.c_str())) {
	cerr << "Could not open " << infile << endl;
	return 1;
      }
    }
  
    bool matched = false;
    while(!atEnd(seq_handle)) {

      try {

	readRecord(id, seq, qual, seq_handle);
	
	// Regex
	if(seq_regex) {

	  if(regex_search(match_type, regex("a"))) {
	    match_type = "frcR";
	  }

	  // Due to some quirks of Seqan, I have to do a number of format
	  // conversions, so this isn't as elegant as I would like.
	  // Also this assumes DNA, not RNA, even though RNA could work
	  // fine. Note that any type of sequence will work with regular
	  // forward matching.
	  matched = false;
	  for(char c: match_type) {
	    switch (c) {
	    case 'f':
	      {
		CharString _seq = seq;
		matched |= regex_search(toCString(_seq), regex_pattern,
					regex_match_flags);
		break;
	      }
	    case 'r':
	      {
		ModifiedString<CharString, ModReverse> rseq(seq);
		CharString _seq(rseq);
		matched |= regex_search(toCString(_seq), regex_pattern,
					regex_match_flags);
		break;
	      }
	    case 'c':
	      {
		Dna5String dseq(seq);
		complement(dseq);
		CharString _seq(dseq);
		matched |= regex_search(toCString(_seq), regex_pattern,
					regex_match_flags);
		break;
	      }
	    case 'R':
	      {
		Dna5String dseq(seq);
		reverseComplement(dseq);
		CharString _seq(dseq);
		matched |= regex_search(toCString(_seq), regex_pattern,
					regex_match_flags);
		break;
	      }
	    }
	  } // End match_type loop

	} else {
	  
	  matched = regex_search(toCString(id), regex_pattern,
				 regex_match_flags);
	} // End regex
      
	// Write out if matched
	if((matched && !inverted) || (!matched && inverted)) {
	  nmatched++;
	  cout << ">" << id << endl << seq << endl;
	}

      } catch (Exception const &e) {

	cout << "Error: " << e.what() << endl;
	return 1;

      } // End try-catch for record reading.

    } // End single file reading loop
    
    close(seq_handle);

  } // End loop over files
  
  if(nmatched) {
    return 0;
  } else {
    return 1;
  }
  
}
