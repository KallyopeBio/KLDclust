/*
 * FileReader.cpp
 *
 *  Created on: Sep 11, 2018
 *      Author: Sarah Christensen
 */

#include "FileReader.h"
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <boost/algorithm/string.hpp>

// Instantiate a new file from file path and specify buffer size
FileReader::FileReader(const std::string file)
{
	// Open stream of file
	m_fileStream.open(file.c_str());

	// If file does not open, print error
	if (!m_fileStream.is_open())
	{
		throw std::runtime_error("ERROR: Could not open the file \"" + file + "\" for reading");
	}

}

//= Returns a tuple (location, location, r score) from input data line
std::tuple<unsigned int, unsigned int, double> FileReader::formatLine(std::string line)
{
	// File must have 3 space separated columns with integer locations and a double for correlation score
	unsigned int loc_a(0);                 // First chromosome location
	unsigned int loc_b(0);				   // Second chromosome location
	double score(0.0);			   // R score for the two locations

	auto it = line.begin();
	bool match = boost::spirit::qi::phrase_parse(it, line.end(),
			// Begin grammar
			boost::spirit::qi::uint_[boost::phoenix::ref(loc_a) = boost::spirit::qi::_1] >>
			boost::spirit::qi::uint_[boost::phoenix::ref(loc_b) = boost::spirit::qi::_1] >>
			boost::spirit::qi::double_[boost::phoenix::ref(score) = boost::spirit::qi::_1]>>
			boost::spirit::qi::eoi,
			// End grammar
			boost::spirit::qi::space);

	if (match){
		return std::make_tuple (loc_a, loc_b, score);
	}
	else {
		throw std::runtime_error("ERROR: Format of line in input file not recognized for reading. Desired space separated format: [unsigned int] [unsigned int] [double]");
	}

}


//= Get next line from file
bool FileReader::getNextRead(unsigned int & loc_a, unsigned int & loc_b, double & score)
{
	// State flag for end of file
	bool eof= false;

	// Initialize values
	loc_a = 0;
	loc_b = 0;
	score = 0;

	// Ensure stream is still open
	if (m_fileStream.is_open())
	{
		std::string dummy;

		// Skip lines where the locations are the same, or the score does not meet the threshold
		while ((loc_a == loc_b)) {
			if (std::getline(m_fileStream, dummy)) {

				std::tuple<unsigned int, unsigned int, double> dummyTuple;
				dummyTuple = formatLine(dummy);

				// Place data into input variables
				if (std::get<0>(dummyTuple) <= std::get<1>(dummyTuple)){
					loc_a = std::get<0>(dummyTuple);
					loc_b = std::get<1>(dummyTuple);
				}
				else {
					loc_a = std::get<1>(dummyTuple);
					loc_b = std::get<0>(dummyTuple);
				}

				score = std::abs(std::get<2>(dummyTuple));
			}
			else {
				m_fileStream.close();
				eof = true;
				break;
			}
		}
	}
	else{
		// We previously read to the end of the file
		eof= true;
	}
	return !eof;
}

// Destroy FileReader making sure to free all resources
FileReader::~FileReader()
{
	if (m_fileStream.is_open())
		{
			m_fileStream.close();
		}
}

// Check if string is file
bool FileReader::is_file(std::string path){
	struct stat s;
	if((stat(path.c_str(),&s) == 0) &&  ( s.st_mode & S_IFREG )){
		// It is a file
		return true;
	}
	else{
		return false;
	}
}

// Check if string is folder
bool FileReader::is_folder(std::string path){
	struct stat s;
	if((stat(path.c_str(),&s) == 0) &&  ( s.st_mode & S_IFDIR )){
		// It is a file
		return true;
	}
	else{
		return false;
	}
}

// Return list of files in folder (or just file if file path)
std::vector<std::string> FileReader::get_files(std::string path){
	// Create vector of input files
	std::vector<std::string> inputFiles;

	// Add input files to this vector
	if (is_file(path)){
		inputFiles.push_back(path);
		return inputFiles;
	}
	else if (is_folder(path)){
		DIR *dirp;
		struct dirent *directory;

		dirp = opendir(path.c_str());
		if (dirp){
			while ((directory = readdir(dirp)) != NULL){
				std::string fileName(path+(directory->d_name));

				// We only add files and ignore sub-directories
				if (is_file(fileName)){
					inputFiles.push_back(fileName);
				}
			}
			closedir(dirp);
		}
		return inputFiles;
	}
	else {
		throw std::runtime_error("An input path is not recognized as a file or a folder");
	}
}
