/*
 * FileReader.h
 *
 *  Created on: Sep 11, 2018
 *      Author: Sarah Christensen
 * Description: Class takes as input a PLINK file and returns a reader to stream over the file line by line
 */

#ifndef FILEREADER_H_
#define FILEREADER_H_

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <assert.h>
#include <climits>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <boost/algorithm/string.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/tuple/tuple.hpp>

class FileReader
{

public:

	// Instantiate a new file instance from file path (with default buffer)
	FileReader(const std::string file);

	// Destructor
	~FileReader();

	// Get next line from file
	// Score threshold is the minimum score needed for line to be read in
	bool getNextRead(unsigned int & loc_a,
			         unsigned int & loc_b,
					 double & score);

	// Check if string is file
	static bool is_file(std::string path);

	// Check if string is folder
	static bool is_folder(std::string path);

	// Return list of files in folder (or just file if file path)
	static std::vector<std::string> get_files(std::string path);

protected:

	// Stream over input file
	std::ifstream m_fileStream;

	//= Returns formated data line from a correlation file with proper types (location, location, r score)
	std::tuple<unsigned int, unsigned int, double> formatLine(std::string line);
};


#endif /* FILEREADER_H_ */
