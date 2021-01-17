#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <unordered_map>

static void GetDistributionOfBackers();
static std::unordered_map<double, double> CreateBins(int maxVal, int binLength);
static std::unordered_map<double, double> FindFrequencies(std::unordered_map<int, int> map, std::string fileName);
static std::unordered_map<double,double> Histogram(int maxValue);
static std::vector<int> SortHistKeys(std::unordered_map<double, double> hist);
static std::vector<double> HistValues(std::unordered_map<double, double> hist, std::vector<int> sHistKeys);
static void FindBestDurationOfCampaign();
static void FindSuccessfullProjectTypes();
static void FindSuccessfullDaysMonths();
static void PrintInformation(std::string histogramName);

int Partition(std::vector<int>& ary, int left, int right);
void QuickSort(std::vector<int>& vec, int left, int right);
static double StandardDeviation(std::string fileName);
static void PrintHist();

//double backers = 45957;
//double sDeviation = 688.628479;
//double mean = 69.973192;
// bad global variables! Change later.
std::vector<int> histKeys;
std::vector<double> histValues;
double maxGlobalVal = 0;
std::string histoGlobalName;
int main()
{
	std::vector<double> vec, vec1, pledgeList, backersList;
	std::unordered_map<double, double> histogram;	
	
	std::ifstream myfile("pledged.txt");
	std::ifstream backersfile("backers.txt");
	std::ifstream durationFile("duration.txt");

	std::string line;
	int count = 0;
	int i = 0;
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			if (line == "")
				vec.push_back(i);
			if (line != "")
				pledgeList.push_back(std::stoi(line));
			i++;
			
		}
		myfile.close();
	}

	double sum = 0;
	for (int i = 0; i < pledgeList.size(); i++)
	{
		sum += pledgeList[i];
	}
	double res = sum / pledgeList.size();	
	
	i = 0;
	int maxBinValue = 0;

	if (backersfile.is_open())
	{
		while (getline(backersfile, line))
		{
			if (line == "")
				vec1.push_back(i);
			if (line != "")
			{
				int data = stoi(line);
				if (data > maxBinValue)
					maxBinValue = data;
			}
			i++;
		}
		backersfile.close();
	} 
	double maxBinVal1 = 0;
	if (durationFile.is_open())
	{
		while (getline(durationFile, line))
		{			
			if (line != "")
			{
				double data = stoi(line);
				if (data > maxBinVal1)
					maxBinVal1 = data;
			}
			i++;
		}
		durationFile.close();
	}

	histogram = Histogram(maxBinValue);// backers column
	//histogram = Histogram(maxBinVal1); //duration column max value	
	
	//StandardDeviation("duration.txt");
	PrintHist();
	PrintInformation(histoGlobalName);
	FindBestDurationOfCampaign();
	FindSuccessfullProjectTypes();
	FindSuccessfullDaysMonths();

	return 0;
}
static void PrintHist()
{
	std::ofstream histFile, histFile1;

	histFile.open("histOut1.txt");
	histFile1.open("histOut2.txt");

	if (histFile.is_open() && histFile1.is_open())
	{
		/*for (auto it = histogram.cbegin(); it != histogram.cend(); it++)
		{
			histFile << it->first << "\n";
			histFile1 << it->second << "\n";
		}*/
		if (histKeys.size() == histValues.size())
		{
			for (int i = 0; i < histKeys.size(); i++)
			{
				histFile << histKeys[i] << "\n";
				histFile1 << histValues[i] << "\n";
			}
		}
		else
		{
			histFile << "Error in printing sizes are not the same" << "\n";
			histFile1 << "Error in printing sizes are not the same"  << "\n";
		}
		histFile.close();
		histFile1.close();
	}
}
static std::unordered_map<double, double> CreateBins(int maxVal, int binLength)
{
	//maxVal = 128;
	int binInterval = binLength;
	int steps = (maxVal / binInterval) + 2;
	std::unordered_map<double, double> hist;
	int val = 0;
	
	for (int i = 0; i < steps; i++)
	{				
		hist.emplace(val, 0);
		val += binInterval;	
		
	}
	return hist;
}
static std::unordered_map<double, double> FindFrequencies(std::unordered_map<double, double> hist, std::string fileName)
{
	std::ifstream file(fileName +".txt");
	
	std::string line;	
	histKeys = SortHistKeys(hist);
	//int i = 0;
	if (file.is_open())
	{
		while (getline(file, line))
		{			
			// figure out the frequency for each bin in the histogram
			if (line != "")
			{
				double data = stof(line);
				data = ceil(data);
				
				for (int i = 0; i < histKeys.size(); i++)
				{
					if (data <= histKeys[i])
					{
						hist[histKeys[i]] = hist[histKeys[i]] + 1;					
						break;
					}					
				}
			}
			//i++;	
		}
		file.close();
	}
	histValues = HistValues(hist, histKeys);
	return hist;
}
static std::unordered_map<double,double> Histogram(int maxVal)
{
	std::unordered_map<double,double> hist;
	maxVal = 150;
	hist = CreateBins(maxVal, 10);// set bin size/ bin width.

	histoGlobalName = "backers";
	maxGlobalVal = maxVal;
	hist = FindFrequencies(hist,histoGlobalName);

	return hist;
}

static double StandardDeviation(std::string fileName)
{
	double sDeviation = 0;
	std::ifstream durationFile(fileName + ".txt");

	std::string line;
	int i = 0;
	double mean = 0;
	double sum = 0;
	double variance = 0;
	std::vector<double> durationDataColumnList;
	std::vector<double> outputList;
	if (durationFile.is_open())
	{
		while (getline(durationFile, line))
		{			
			if (line != "")
			{
				double data = stoi(line);
				durationDataColumnList.push_back(data);
				sum += data;
				i++;
			}			
		}
		durationFile.close();
	}
	mean = sum / i;

	for (int i = 0; i < durationDataColumnList.size(); i++)
	{
		double val = pow(abs(durationDataColumnList[i] - mean),2);
		outputList.push_back(val);
	}
	double varianceSum = 0;
	for (int i = 0; i < outputList.size(); i++)
	{
		varianceSum += outputList[i];
	}
	variance = varianceSum / outputList.size();

	sDeviation = sqrt(variance);

	return sDeviation;
}

static std::vector<int> SortHistKeys(std::unordered_map<double, double> hist)
{
	std::vector<int> sortedHistogram;
	

	for (auto it = hist.cbegin(); it != hist.cend(); it++)
	{
		sortedHistogram.push_back(it->first);
	}
	QuickSort(sortedHistogram, 0, sortedHistogram.size() - 1);
	return sortedHistogram;
}

int Partition(std::vector<int>& vec, int left, int right)
{
	int pivot = vec[left];	
	int i = left + 1;

	for (int j = left + 1; j <= right; j++)
	{
		if (vec[j] < pivot)
		{
			int temp = vec[j];
			vec[j] = vec[i];
			vec[i] = temp;
			i = i + 1;
		}
	}	
	int temp = vec[left];
	vec[left] = vec[i - 1];
	vec[i - 1] = temp;
	i = i - 1;
	return i;
}

void QuickSort(std::vector<int>& vec, int left, int right)
{
	if (left >= right)
		return;

	int pivotPostion = Partition(vec, left, right);

	QuickSort(vec, left, pivotPostion - 1);
	QuickSort(vec, pivotPostion + 1, right);

}
static std::vector<double> HistValues(std::unordered_map<double, double> hist, std::vector<int> sHistKeys)
{
	std::vector<double> histValues;

	for (int i = 0; i < sHistKeys.size(); i++)
	{
		histValues.push_back(hist[sHistKeys[i]]);
	}

	return histValues;

}
static double Mean(std::string fileName)
{
	std::ifstream file;
	std::string line;
	std::vector<double> vec;
	file.open(fileName + ".txt");
	if (file.is_open())
	{
		while (getline(file, line))
		{
			
			if (line != "")
			{
				double data = stof(line);
				vec.push_back(data);
				
			}
			
		}
		file.close();
	}
	double sum = 0;
	for (int i = 0; i < vec.size(); i++)
	{
		sum += vec[i];
	}
	double mean = sum / vec.size();
	return mean;
}
static void PrintInformation(std::string histogramName)
{
	std::cout << "The following lines contain information collected from the column: " 
		<< histogramName 
		<< " in the Kickstarter dataset"
		<< std::endl;
	std::cout << "The standard deviation of the " << histogramName << " column is: " << StandardDeviation(histogramName) << std::endl;
	std::cout << "The mean of the " << histogramName << " column is: " << Mean(histogramName) << std::endl;
	std::cout << "The max value of the " << histogramName << " column is: " << maxGlobalVal << std::endl;	
	
	std::cout << "Printing the Bin and Frequency values for the Histogram of the column " << histogramName << std::endl;
	std::cout << "Bins" << "\t" << "Frequencies" << std::endl;
	if (histKeys.size() == histValues.size())
	{
		for (int i = 0; i < histKeys.size(); i++)
		{
			std::cout << histKeys[i] << "\t";
			std::cout << histValues[i] << std::endl;				
		}
	}
	else
	{
		std::cout << "Error in printing sizes are not the same" << "\n";
		std::cout << "Error in printing sizes are not the same" << "\n";
	}
}

static void FindBestDurationOfCampaign()
{
	std::unordered_map<double, double> histogram;
	std::ifstream statusFile("status.txt");
	std::ifstream durationFile("duration.txt");
	if (statusFile.is_open() && durationFile.is_open())
	{
		std::string statusLine, durationLine;

		while (getline(statusFile, statusLine) && getline(durationFile, durationLine))
		{
			if (statusLine != "" && durationLine != "")
			{				
				if (statusLine == "successful")
				{
					double data = stof(durationLine);
					data = ceil(data);
					if (histogram.find(data) != histogram.end())
					{
						histogram[data] += 1;
					}
					else
					{
						histogram.emplace(data, 1);
					}
				}

			}
		}
		statusFile.close();
		durationFile.close();
	}

	

	std::vector<int> days;
	std::vector<double> frequencies;

	days = SortHistKeys(histogram);
	frequencies = HistValues(histogram, days);
	std::ofstream statusOutFile("statusOut.txt");
	std::ofstream durationOutFile("goalOut.txt");

	/*std::cout << "Histogram of successfull goal amount" << std::endl;
	std::cout << "Goal Amount" << "\t" << "Frequency" << std::endl;*/
	std::cout << "Duration Of Campaign" << std::endl;
	std::cout << "Days" << "\t" << "Frequency" << std::endl;
	if (statusOutFile.is_open() && durationOutFile.is_open())
	{
		if (days.size() == frequencies.size())
		{
			for (int i = 0; i < days.size(); i++)
			{
				statusOutFile << days[i] << "\n";
				durationOutFile << frequencies[i] << "\n";	
				std::cout << days[i] << "\t"<<"\t";
				std::cout << frequencies[i] << "\n";
			}
		}
		else
		{
			statusOutFile << "Error in printing sizes are not the same" << "\n";
			durationOutFile << "Error in printing sizes are not the same" << "\n";
		}
		statusOutFile.close();
		durationOutFile.close();	
	}
	int x = 0;	
}
static void FindSuccessfullProjectTypes()
{
	std::unordered_map<std::string, double> histogram;
	std::ifstream statusFile("status.txt");
	std::ifstream dateFile("category.txt");
	if (statusFile.is_open() && dateFile.is_open())
	{
		std::string statusLine, categoryLine;

		while (getline(statusFile, statusLine) && getline(dateFile, categoryLine))
		{
			if (statusLine != "" && categoryLine != "")
			{
				if (statusLine == "successful")
				{
					std::string data = categoryLine;
				
					if (data == "Film &amp; Video")
						data = "Film & Video";
					if (histogram.find(data) != histogram.end())
					{
						histogram[data] += 1;
					}
					else
					{
						histogram.emplace(data, 1);
					}
				}

			}
		}
		statusFile.close();
		dateFile.close();
	}
	
	std::ofstream statusOutFile("statusOut.txt");
	std::ofstream categoryOutFile("categoryOut.txt");

	std::cout << "Histogram of successfull Categories" << std::endl;
	std::cout << "Categories" << "\t" <<"\t" << "Frequency" << std::endl;
	if (statusOutFile.is_open() && categoryOutFile.is_open())
	{
		
		for (auto itr = histogram.cbegin(); itr!= histogram.cend(); itr++)
		{
			statusOutFile << itr->first << "\n";
			categoryOutFile << itr->second << "\n";
			std::cout << itr->first << "\t" << "\t" <<"\t";
			std::cout << itr->second << "\n";
		}		
		
		statusOutFile.close();
		categoryOutFile.close();
	}
	int x = 0;
}
static void FindSuccessfullDaysMonths()
{
	std::unordered_map<std::string, double> histogram, monthHist;
	std::ifstream statusFile("status.txt");
	std::ifstream dateFile("fundedDate.txt");
	if (statusFile.is_open() && dateFile.is_open())
	{
		std::string statusLine, dateLine;

		while (getline(statusFile, statusLine) && getline(dateFile, dateLine))
		{
			if (statusLine != "" && dateLine != "")
			{
				if (statusLine == "successful")
				{
					std::string data = dateLine;
					std::string month = data.substr(8,3);
					data = data.substr(0, 3);
				
					if (histogram.find(data) != histogram.end())
					{
						histogram[data] += 1;
					}
					else
					{
						histogram.emplace(data, 1);
					}
					if (monthHist.find(month) != monthHist.end())
					{
						monthHist[month] += 1;
					}
					else
					{
						monthHist.emplace(month, 1);
					}
				}

			}
		}
		statusFile.close();
		dateFile.close();
	}

	std::ofstream statusOutFile("statusOut.txt");
	std::ofstream dateOutFile("dateOut.txt");
	std::ofstream monthOutFile("monthOut.txt");

	std::cout << "Histogram of Days" << std::endl;
	std::cout << "Days" << "\t" << "\t" << "Frequency" << std::endl;
	if (statusOutFile.is_open() && dateOutFile.is_open())
	{
		for (auto itr = histogram.cbegin(); itr != histogram.cend(); itr++)
		{
			statusOutFile << itr->first << "\n";
			dateOutFile << itr->second << "\n";
			std::cout << itr->first << "\t" << "\t" << "\t";
			std::cout << itr->second << "\n";
		}

		statusOutFile.close();
		dateOutFile.close();
	}
	std::cout << "Histogram of Months" << std::endl;
	std::cout << "Months" << "\t" << "\t" << "Frequency" << std::endl;
	statusOutFile.open("months.txt");
	if (monthOutFile.is_open())
	{
		for (auto itr = monthHist.cbegin(); itr != monthHist.cend(); itr++)
		{
			statusOutFile << itr->first << "\n";
			monthOutFile << itr->second << "\n";
			std::cout << itr->first << "\t" << "\t" << "\t";
			std::cout << itr->second << "\n";
		}		
		monthOutFile.close();
	}
	int x = 0;
}