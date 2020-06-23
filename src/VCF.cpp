/*
 * VCF.cpp
 *
 *  Created on: Nov 6, 2018
 *      Author: sarah
 */
#include "VCF.h"

// Instantiate a new VCF file from file path
VCF::VCF(const std::string file, const unsigned int size, const unsigned int maxLength, const bool phased, const std::string indFile): m_maxLength(maxLength)
{
	// Open stream of file
	std::ifstream fileStream;
	fileStream.open(file.c_str());

	// If file does not open, print error
	if (!fileStream.is_open())
	{
		std::cerr << "ERROR: Could not open the file \"" << file << "\" for reading" << std::endl;
		return;
	}

	// Set column location for where data starts (indexed at 1)
	unsigned int headerCount = 1;
	unsigned int dataRowNum(0);

	// Loop through each line of file
	std::string line;

	while (std::getline(fileStream, line))
	{
		// Check if first character of the line is a header
		if(line.size() && line[0] =='#'){
			// Increase header count
			headerCount++;
			// The only header we are interested in is the one that gives the list of individuals
			if(line.compare(0, 6, "#CHROM") == 0){

				// Find requested individuals in this VCF line
				parseIndividuals(indFile, line);

				// Resize member variables to correct size
				m_variants.resize(size-headerCount);
				m_avgScores.resize(size-headerCount);
				m_counts.resize(size-headerCount);
				m_indScores.resize(size-headerCount);

			}
			// If it is a header line that does not contain the list of individuals, there is no processing
			else continue;
		}

		// Otherwise this is a data line that must be processed
		else {
			parseLine(line, dataRowNum, phased);
			dataRowNum++;
		}
	} // End while loop over file

	// Report on the total number of SNPs
	std::cout << "The total number of SNPs in this file is: " << dataRowNum << std::endl;
	std::cout << "The total number of variant sites in this file is: " << m_variantCounts.size() << std::endl;

	// Remove variants with all missing values, no variation
	fillVariantsRemoved();

	// Close file
	fileStream.close();

}

// Parse individuals included in sample
void VCF::parseIndividuals(const std::string indFile, const std::string line){
	// Put list of relevant individuals into set
	// This potentially should become a different class at some point
	std::cout << "Reading in individual IDs" << std::endl;
	std::unordered_set<std::string> individualsSubsetID;
	std::string ID;
	std::ifstream in(indFile);
	while(true){
		in >> ID;
		if( in.eof() ) break;
		// Add conversion to dictionary
		individualsSubsetID.insert(ID);
	}
	in.close();

	// Match requested IDs to IDs found in VCF file
	std::istringstream iss(line);
	int column = 1;

	// Iterate over column headers
	for(std::string word; iss >> word; ){
		// Check if this individual ID is one of interest
		if (individualsSubsetID.count(word)){
			// If so, add the column number
			m_individualsSubset.insert(column);
		}

		// Increment column
		column++;
	}

	std::cout << "Out of the requested " << individualsSubsetID.size() << " individuals, "<< m_individualsSubset.size() << " individuals found in VCF file" << std::endl;
}

// Parse chromosome location in data line
void VCF::parseChromLoc(const unsigned int loc, const unsigned int index){
	// Add chromosome location to vector where each index is the data row index
	// A location may appear more than once since multiple SNPs can have same BP location
	m_variants[index]=loc;

	// Try to add location to the dictionary of location counts
	{std::pair< std::unordered_map<unsigned int, unsigned int>::iterator, bool > f= m_variantCounts.insert(std::make_pair(loc, 1));

	// Check if the insertion failed because the key was already present
	if (!f.second)
		// If already present, add in r score to entry
		f.first->second += 1;}
}


// Parse a line of data from VCF file assuming that it is not phased (use 0, 1, 2 to represent genotype)
void VCF::parseLine(const std::string line, const unsigned int index, bool phased){

	// Initialize values
	std::istringstream iss(line);
	int column = 1;
	std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int> tmp = std::make_tuple(0, 0, 0, 0, 0);
	m_counts[index]=tmp;

	// Iterate over column headers
	for(std::string data; iss >> data; ){

		// Get CHROM location
		if (column==2){
			parseChromLoc(std::stoul(data), index);
		}

		// We are interested in the individual genotypes in the sample of interest
		if (m_individualsSubset.count(column)){

			if (phased){
				parsePhasedData(data, index, column);
			}
			else{
				parseUnphasedData(data, index, column);
			}

		} // Ends if data row

		// Increment column
		column++;

	} // Ends for loop

	// Print warning message if there are many missing values
	if ((std::get<3>(m_counts[index])/m_individualsSubset.size()) > 0.1){
		std::cerr << "WARNING: The variant at location " << m_variants[index] << " has a fraction of " << (std::get<3>(m_counts[index])/m_individualsSubset.size()) << " missing values" << std::endl;
	}

	// Compute average score for this variant across all present individuals making sure it doesn't round (#1 + 2*#2)/(# Present)
	if (std::get<4>(m_counts[index]) != 0){
		m_avgScores[index] = std::get<1>(m_counts[index])+2*std::get<2>(m_counts[index]);
		m_avgScores[index]/= std::get<4>(m_counts[index]);
	}
}

// Parse a line of data from VCF file assuming that it is not phased (use 0, 1, 2 to represent genotype)
void VCF::parseUnphasedData(const std::string data, const unsigned int index, const unsigned int column){

	// Check that data are available for this individual at this site
	if ((data[0] != '.') && (data[2] != '.')){

		// Get the genotype score for this individual
		// We convert from Char to Int type by subtracting '0'
		int x(data[0]-'0');
		int y(data[2]-'0');
		int score(x+y);

		// Increase count of non-missing observations
		std::get<4>(m_counts[index])++;

		// Make sure each haplotype is binary
		if ((x==0 || x==1) && (y==0 || y==1)){
			// If score is 0 we just increment count to know how many zero observations we have
			if (score == 0){
				std::get<0>(m_counts[index])++;
			}

			// Increment count and add this individual to dictionary of individuals with variant alleles
			// We identify individuals with their column number rather than ID number in data here
			else if (score == 1){
				std::get<1>(m_counts[index])++;
				m_indScores[index].insert({column, 1});
			}

			// Increment count and add this individual to dictionary of individuals with variant alleles
			// We identify individuals with their column number rather than ID number in data here
			else if (score == 2){
				std::get<2>(m_counts[index])++;
				m_indScores[index].insert({column, 2});
			}
		}

		// Otherwise, the variant is multiallelic and will be removed
		else {
			assert(score>2);
			m_variantCounts[m_variants[index]] = score;
		}
	}
	// Otherwise, report that this individual is missing data from this site
	// We use a flag of -1 to indicate the individual is missing
	else {
		std::get<3>(m_counts[index])++;
		m_indScores[index].insert({column, -1});
	}
}

// Parse a line of data from VCF file assuming that it is phased
void VCF::parsePhasedData(const std::string data, const unsigned int index, const unsigned int column){

	// We treat each haplotype as if it is just another individual sampled
	for (unsigned int i=0; i<=2; i=i+2) {
		// Check that data are available for this individual at this site
		if (data[i]!='.'){

			// Increase count of non-missing
			std::get<4>(m_counts[index])++;

			// Increase count of 0 values
			if (data[i] == '0'){
				std::get<0>(m_counts[index])++;
			}

			// Increment count and add this individual to dictionary of individuals with variant alleles
			// We identify individuals with their column number rather than ID number in data here
			// Since data is phase, each individual has a positive and negative ID for each haplotype
			else if (data[i] == '1'){
				std::get<1>(m_counts[index])++;
				if (i==0){m_indScores[index].insert({column, 1});}
				if (i==2){m_indScores[index].insert({(-1)*column, 1});}
			}

			// Otherwise, the score indicates the variant is multi allelic
			else {
				m_variantCounts[m_variants[index]] = 10;
			}
		}

		// Otherwise, report that this individual is missing data from this site
		else{
			std::get<3>(m_counts[index])++;
			if (i==0){m_indScores[index].insert({column, -1});}
			if (i==2){m_indScores[index].insert({(-1)*column, -1});}
		}
	} // Ends for loop
}

// Get the number of individuals in the subset
unsigned int VCF::getObsCount(bool phased){
	// We use the Fisher transformation so that z= [F(r)-F(o)](n-3)^.5
	unsigned int obsCount = m_individualsSubset.size();

	// If phased we get 2 observations from each individual
	if (phased) {
		obsCount *= 2;
	}

	return obsCount;
}

// Add variants to set m_variantsRemoved that should be removed from analysis
void VCF::fillVariantsRemoved(){

	// Initialize values
	unsigned int multi(0);
	unsigned int missing(0);
	unsigned int variation(0);

	// Iterate over all variants in input data
	for (unsigned int i=0; i<m_variants.size(); i++){

		// Remove variant if there is more than one SNP at the site
		if (m_variantCounts[m_variants[i]]!=1){
			m_variantsRemoved.insert(m_variants[i]);
			multi++;
		}

		// Remove variant if all values are missing
		else if (std::get<4>(m_counts[i])==0){
					m_variantsRemoved.insert(m_variants[i]);
					missing++;
		}

		// Remove variant if there is no variation in observed values
		else if ((std::get<0>(m_counts[i]) == std::get<4>(m_counts[i])) ||
					(std::get<1>(m_counts[i]) == std::get<4>(m_counts[i])) ||
						(std::get<2>(m_counts[i]) == std::get<4>(m_counts[i]))){
							m_variantsRemoved.insert(m_variants[i]);
							variation++;
		}

	}
	std::cerr << "WARNING: There are " << m_variantsRemoved.size() << " variants that will be removed for the following reasons" << std::endl;
	std::cerr << "\t More than one SNP or multiallelic at site: " << multi << std::endl;
	std::cerr << "\t All individuals missing for this site: " << missing << std::endl;
	std::cerr << "\t No variation in individuals at this site: " << variation << std::endl;
}

// Calculate theta parameter based on Li and Stephens (2003); EQ 2.8 in Wen and Stephens 2010
double VCF::calculateTheta(){

	// Reset theta back to zero to be safe
	double theta=0;

	// Get the number of individuals in sample times 2
	unsigned int twiceNumInd = 2*m_individualsSubset.size();

	// Initialize sum of inverse
	double sumInv=0;

	// Compute the sum of the finite reciprocal series
	for (unsigned int i=1; i<twiceNumInd; i++){
		sumInv = sumInv + (1/i);
	}

	// It may seem confusing, but now we need to take the inverse of the sum!
	double invSumInv = (1/sumInv);

	// Finally, we can compute theta
	theta = (invSumInv)/(twiceNumInd + invSumInv);

	return theta;
}

// Calculate correlations and write correlations above threshold to file
void VCF::calculateShrinkageCorrelation(std::string outputDir, GeneticMap gmap, unsigned int NE, double minCorr){

	// START: Shrinkage Preprocessing
	// Calculate theta parameter based on Li and Stephens (2003); EQ 2.8 in Wen and Stephens 2010
	double theta(calculateTheta());
	std::cout << "The theta parameter for this population is: " << theta << std::endl;

	// Check that every site in the chromosome is also in the gmap file
	for (unsigned int i=0; i<m_variants.size(); i++){
		try{gmap.getValue(m_variants[i]);}
		catch (...){
			std::cerr << "ERROR: Variant location " << m_variants[i] << " at index " << i << " was not found in genetic map." << std::endl;
			return;
		}
	}
	// END: Shrinkage Preprocessing

	// Create chromosome
	Chromosome chrom(m_variantCounts);

	// Create vector to store variance result for each variant
	std::vector<double> variance;
	variance.resize(m_variants.size());

	// Pre-calculate parameters used across all variants
	double numeratorRR(NE*4.0);
	double denominatorRR(2.0*m_individualsSubset.size());
	double thetaSq((1.0-theta)*(1.0-theta));
	double thetaOver2((theta/2.0)*(1.0-(theta/2.0)));

	// Calculate base case variance for each variant
	unsigned int maxThreads = omp_get_max_threads();
	#pragma omp parallel for schedule(static) num_threads(maxThreads)
	for (unsigned int i=0; i<m_variants.size(); i++){

		// Only include positions that have not been removed
		if (m_variantsRemoved.count(m_variants[i])==0){
			std::tuple<double, unsigned int> tmp(calculateInteraction(i,i));
			unsigned int interactionScore(std::get<0>(tmp));

			// Get the total number of individuals for this variant
			unsigned int n = std::get<4>(m_counts[i]);

			// Now compute variance (no recombination adjustment is applied to diagonal)
			double covariance((interactionScore - (n*m_avgScores[i]*m_avgScores[i]))/(n-1.0));

			// Apply shrinkage adjustment (Eq. 2.6 and 2.7 Wen and Stephens 2010)
			// On diagonal there is an extra term to be added (Second term in Eq. 2.6)
			double shrunkCov(thetaSq*covariance + thetaOver2);

			// Add to variance dictionary
			variance[i]=shrunkCov;

			// We removed this in m_variantsRemoved because it will be in denominator of covariance later
			assert(variance[i]!=0);
		}
	}

	// Create an array of streams
	unsigned int hashNum(maxThreads);
	std::ofstream stream[maxThreads][hashNum];
	for (unsigned int i = 0; i < maxThreads; ++i){
		for (unsigned int j = 0; j < hashNum; ++j){
			std::string outFile(outputDir+std::to_string(i)+"_"+std::to_string(j)+".txt");
			stream[i][j].open(outFile);
			assert(stream[i][j].good());
		}
	}

	std::cout.precision(std::numeric_limits< double >::max_digits10);

	// Now calculate empirical covariance, apply shrinkage, and then convert to correlation
	// We don't need the last variant since there are no variants to right to compute correlation with
	#pragma omp parallel for schedule(static) num_threads(maxThreads)
	for (unsigned int i=0; i<(m_variants.size()-1); i++){
		// Only include positions with a single SNP that has observations and is above minor allele frequency (MAF)
		if (m_variantsRemoved.count(m_variants[i])==0){

			// Get thread number
			unsigned int t = omp_get_thread_num();

			// Find right bound within allowed window length
			unsigned int rightBound = chrom.getRightNeighborBP(m_variants[i], m_maxLength);

			// Initialize variant j to be compared to i
			unsigned int j(i+1);

			// Now consider all variants to the right that are within the given window size
			while ((m_variants[j]<=rightBound) && (j<m_variants.size())){
				// Only include positions with a single SNP and above minor allele frequency (MAF)
				if (m_variantsRemoved.count(m_variants[j])==0) {
					// Calculate recombination rate
					// Equivalent to exp((-df*NE*4)/(2*m)) where m is number of individuals, df is difference in
					double rr=exp(((gmap.getValue(m_variants[i]) - gmap.getValue(m_variants[j]))*numeratorRR)/denominatorRR);

					// Calculate numerator in covariance and get overlap in missing values
					std::tuple<unsigned int, unsigned int> tmp(calculateInteraction(i,j));
					double interactionScore(std::get<0>(tmp));
					unsigned int overlap(std::get<1>(tmp));

					// Compute the total number of individuals present in both samples
					// total - (missing[i] + missing[j] - overlap)  = (present[i]+missing[i]) - (missing[i] + missing[j] - overlap)
					// = present[i]- missing[j] + overlap
					double n = (std::get<4>(m_counts[i])+overlap)-std::get<3>(m_counts[j]);

					// Now compute covariance
					double covariance((interactionScore - (n*m_avgScores[i]*m_avgScores[j]))/(n-1.0));

					// Apply shrinkage to covariance (See Eq. 2.6 and 2.7 in Wen and Stephens 2010)
					double shrunkCov(thetaSq*rr*covariance);

					// Convert to pearson correlation r2
					double r2((shrunkCov*shrunkCov)/(variance[i]*variance[j]));

					// If above threshold, write to file
					assert(r2<=1.0000001);
					if (r2>=minCorr){
						// Hash starting location into file
						unsigned int h(std::hash<unsigned int>{}(m_variants[i])%hashNum);
						stream[t][h] << std::fixed << m_variants[i] << " " << m_variants[j] << " " << r2 << "\n";
					}
				}

				// Increase j index by one
				j++;
			}
		}

	}

	// Close streams
	for (unsigned int i = 0; i < maxThreads; ++i){
		for (unsigned int j=0; j<hashNum; j++){
			assert(stream[i][j].good());
			stream[i][j].close();
		}
	}
}


// Calculate correlations and write correlations above threshold to file
void VCF::calculateMAFCorrelation(std::string outputDir, double minCorr, double maf){

	// Get number of variants removed before MAF threshold
	unsigned int preMAF(m_variantsRemoved.size());

	// Iterate over all variants in input data
	for (unsigned int i=0; i<m_variants.size(); i++){

		// Remove variant if it has low MAF relative to threshold
		if ((m_avgScores[i] <= maf) || (m_avgScores[i] >=(1-maf))){
			m_variantsRemoved.insert(m_variants[i]);
		}
	}

	std::cerr << "WARNING: There are " << m_variantsRemoved.size()-preMAF << " additional variants removed after MAF thresholding" << std::endl;

	// Create chromosome
	Chromosome chrom(m_variantCounts);

	// Create vector to store variance result for each variant
	std::vector<long double> variance;
	variance.resize(m_variants.size());

	// Calculate base case variance for each variant
	unsigned int maxThreads = omp_get_max_threads();
	#pragma omp parallel for schedule(static) num_threads(maxThreads)
	for (unsigned int i=0; i<m_variants.size(); i++){

		// Only include positions that have not been flagged for removal
		if (m_variantsRemoved.count(m_variants[i])==0) {

			// Calculate sum of the product of the individual scores
			std::tuple<double, unsigned int> tmp(calculateInteraction(i,i));
			long double interactionScore(std::get<0>(tmp));

			// Get the total number of individuals for this variant
			long double n = std::get<4>(m_counts[i]);

			// Compute and add variance to dictionary
			variance[i]=(interactionScore - (n*m_avgScores[i]*m_avgScores[i]))/(n-1.0L);

			// We removed this in m_variantsRemoved because it will be in denominator of covariance later
			assert(variance[i]!=0);
		}
	}

	// Create an array of streams
	unsigned int hashNum(maxThreads);
	std::ofstream stream[maxThreads][hashNum];
	for (unsigned int i = 0; i < maxThreads; ++i){
		for (unsigned int j = 0; j < hashNum; ++j){
			std::string outFile(outputDir+std::to_string(i)+"_"+std::to_string(j)+".txt");
			stream[i][j].open(outFile);
			assert(stream[i][j].good());
		}
	}

	std::cout.precision(std::numeric_limits< double >::max_digits10);

	// Now calculate empirical covariance and then convert to correlation
	// We don't need the last variant since there are no variants to right to compute correlation with
	#pragma omp parallel for schedule(static) num_threads(maxThreads)
	for (unsigned int i=0; i<(m_variants.size()-1); i++){

		// Only include positions that have not been flagged for removal
		if (m_variantsRemoved.count(m_variants[i])==0) {

			// Find right bound within allowed window length
			unsigned int rightBound = chrom.getRightNeighborBP(m_variants[i], m_maxLength);

			// Initialize variant j to be compared to i
			unsigned int j(i+1);

			// Now consider all variants to the right that are within the given window size
			while ((m_variants[j]<=rightBound) && (j<m_variants.size())){
				// Only include positions that have not been flagged for removal
				if (m_variantsRemoved.count(m_variants[j])==0) {

					// Calculate numerator in covariance and get overlap in missing values
					std::tuple<unsigned int, unsigned int> tmp(calculateInteraction(i,j));
					long double interactionScore(std::get<0>(tmp));
					unsigned int overlap(std::get<1>(tmp));

					// Compute the total number of individuals present in both samples
					// total - (missing[i] + missing[j] - overlap)  = (present[i]+missing[i]) - (missing[i] + missing[j] - overlap)
					// = present[i]- missing[j] + overlap
					long double n = (std::get<4>(m_counts[i])+overlap)-std::get<3>(m_counts[j]);

					// Now compute covariance
					long double covariance((interactionScore - (n*m_avgScores[i]*m_avgScores[j]))/(n-1.0L));

					// Convert to pearson correlation r2 (in this context equivalent to phi^2 correlation)
					long double r2((covariance*covariance)/(variance[i]*variance[j]));

					assert((r2<=1.01) && (r2>=0));
					if (r2>1.000001){
						std::cerr << "WARNING: Variant " << m_variants[i] << " and " << m_variants[j] << " have a r2 score of " << r2 << " that has been rounded down to 1" << std::endl;
						r2=1;
					}
					if (r2>=minCorr){
						// Get thread number
						unsigned int t = omp_get_thread_num();

						// Hash starting location into file
						unsigned int h(std::hash<unsigned int>{}(m_variants[i])%hashNum);

						// Write to that file
						stream[t][h] << std::fixed << m_variants[i] << " " << m_variants[j] << " " << r2 << "\n";
					}
				}

				// Increase j index by one
				j++;

			}
		}

	}

	// Close streams
	for (unsigned int i = 0; i < maxThreads; ++i){
		for (unsigned int j=0; j<hashNum; j++){
			assert(stream[i][j].good());
			stream[i][j].close();
		}
	}
}

// Calculates numerator for empirical covariance and the number of missing individuals that overlap
std::tuple<unsigned int, unsigned int> VCF::calculateInteraction(unsigned int i, unsigned int j){
	unsigned int interactionScore(0);

	// Initialize variable that keeps track of the overlap in missing individuals
	unsigned int overlap(0);

	for (auto it=m_indScores[i].begin(); it !=m_indScores[i].end(); it++){
		// Check if individual with non-zero score for variant i also has non-zero score for j
		int ind1(it->first);
		int score1(it->second);
		if (m_indScores[j].count(ind1)){
			// Check if this individual was missing from both variants
			if ((score1==-1) && (m_indScores[j][ind1]==-1)){overlap++;}

			// Check if both individuals are not missing for both variants
			if ((score1!=-1) && (m_indScores[j][ind1]!=-1)){
				// If so, this will contribute to the covariance
				interactionScore = interactionScore + (score1*(m_indScores[j][ind1]));
			}
		}
	}
	return std::make_tuple(interactionScore, overlap);
}

