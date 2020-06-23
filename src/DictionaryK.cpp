/*
 * DictionaryK.cpp
 *
 *  Created on: Sep 13, 2018
 *      Author: Sarah Christensen
 */

#include "DictionaryK.h"

// Instantiate new modularity score input terms from file
DictionaryK::DictionaryK(const std::string file): m_M(0)
{
	// Calculate full K Dictionary from file stream
	FileReader fileStream(file);
	calculateK(fileStream, m_dictK);

	// Get M value on full K dictionary
	calculateM();

	// Normalize K by the square root of 2m to reduce redundant work
	normalize(1.0);
}

// Instantiate new modularity score input terms from file
DictionaryK::DictionaryK(const std::string file, double scaler): m_M(0)
{
	// Calculate full K Dictionary from file stream
	FileReader fileStream(file);
	calculateK(fileStream, m_dictK);

	// Get M value on full K dictionary
	calculateM();

	// Normalize K by the square root of 2m to reduce redundant work
	normalize(scaler);
}

// Instantiate new modularity score input terms from file
DictionaryK::DictionaryK(const std::vector<std::string> fileNames): m_M(0)
{
	// Create vector of dictionaries for each thread
	std::vector<std::unordered_map<unsigned int, double>> dictKVector;
	unsigned int threadCount = omp_get_max_threads();
	dictKVector.resize(threadCount);

	#pragma omp parallel for schedule(static)
	for (unsigned int i = 0; i < fileNames.size(); i++){

		// Get thread number
		unsigned int t = omp_get_thread_num();

		// Get input string
		std::string threadFile(fileNames[i]);
		FileReader fileStream(threadFile);

		// Calculate local dictionary K
		calculateK(fileStream, dictKVector[t]);

	}

	// Reduce dictionaries
	for (unsigned int i=0; i<dictKVector.size(); i++){
		for (std::pair<unsigned int, double> element : dictKVector[i]) {

			// Initialize this crazy type which is what will be returned by dictionary.insert()
			std::pair< std::unordered_map<unsigned int, double>::iterator, bool > f;

			// Try to add first location to the dictionary
			f = m_dictK.insert(std::make_pair(element.first, element.second));

			// Check if the first insertion failed because the key was already present
			if (!f.second){
				// If already present, add in r score to entry
				f.first->second += element.second;
			}
		}
	}

	// Get M value on full K dictionary
	calculateM();

	// Normalize K by the square root of 2m to reduce redundant work
	normalize(1.0);
}

// Instantiate (subset) of dictionary K from file already containing dictionary K values
DictionaryK::DictionaryK(const std::string dictK_file, const std::string m_file, const unsigned int start, const unsigned int stop){

	std::cout << "Reading in value m" << std::endl;
	m_M=0;
	std::ifstream m_in(m_file);
	m_in >> m_M;
	m_in.close();

	std::cout << "Reading in dictionary K" << std::endl;
	unsigned int loc;
	double k_val;
	std::ifstream in(dictK_file);
	while(true){
		in >> loc;
		in >> k_val;
		if( in.eof() ) break;
		// Only add relevant range to dictionary
		if ((loc>= start) && loc <= stop)
			m_dictK[loc]=k_val;
	}
	in.close();
}

// Get dictionary containing degree K_i for each variant i in the input data
void DictionaryK::calculateK(FileReader& fileStream, std::unordered_map<unsigned int, double>& dictK)
{

	// Initialize dummy variables to store lines
	unsigned int loc_a(0);
	unsigned int loc_b(0);
	double r_score(0);

	// Keep track of the number of lines being read
	unsigned int i=0;

	// Loop until the next read is the end of the file
	while(fileStream.getNextRead(loc_a, loc_b, r_score))
	{
		// Verify that r score is in valid range
		assert(r_score>=0.0);
		assert(r_score<=1.0);

		// Initialize this crazy type which is what will be returned by dictionary.insert()
		std::pair< std::unordered_map<unsigned int, double>::iterator, bool > f;

		// Try to add first location to the dictionary
		f = dictK.insert(std::make_pair(loc_a, r_score));

		// Check if the first insertion failed because the key was already present
		if (!f.second){
			// If already present, add in r score to entry
			f.first->second += r_score;
		}

		// Try to add second location to the dictionary
		f = dictK.insert(std::make_pair(loc_b, r_score));

		// Check if the second insertion failed because the key was already present
		if (!f.second){
			// If already present, add in r score to entry
			f.first->second += r_score;
		}

		// Increment i counter
		++i;

	}

	if (i==0)
	{
		throw std::runtime_error("ERROR: No lines read from file when attempting to calculate K");
	}

	return;
}

// Get m corresponding to the total weight of edges in the graph (See m in modularity score equation XX)
void DictionaryK::calculateM()
{
	// Reset m to 0 to be safe
	m_M = 0.0;

	// Iterate over values in K dictionary
	for (auto const& x : m_dictK)
	{
		// Add value from each variant to total
	    m_M += x.second;
	}

	// Divide by 2 since each edge is double counted (appears once at each endpoint)
	m_M = m_M/2.0;

	// Check that the global variable did indeed get calculated and is not equal to the initialization value
	if (m_M==0.0)
	{
		throw std::runtime_error("ERROR: The total weight contained in the input graph is 0");
	}

	return;
}

//= Normalize K dictionary
void DictionaryK::normalize(double scaler)
{
	// Iterate over values in K dictionary
	for (auto & x : m_dictK)
	{
		// Divide value by the square root of 2M
	    x.second = x.second/(sqrt(2.0*m_M*scaler));
	}

	return;
}

//= Allow user to retrieve m
double DictionaryK::getM(){return m_M;}

//= Allow user to retrieve dictionary K containing the degree of each node
std::unordered_map<unsigned int, double> const & DictionaryK::getK(){
	return m_dictK;
}

// Print dictionary K to file
void DictionaryK::printK(const std::string file){

	// Open output file stream
	std::ofstream ofk(file);

	// Loop over the key value pairs and print variant (in BP) space variant degree
	for (std::pair<const unsigned int, double> & kv : m_dictK){
		ofk << kv.first << " " << kv.second << std::endl;
	}
	ofk.close();
}

// Print value m to file
void DictionaryK::printM(const std::string file){
	std::ofstream ofm(file);
	ofm << m_M << std::endl;
	ofm.close();
}



