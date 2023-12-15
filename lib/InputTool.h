#include "ErrorHandler.h"

template <typename... T>
void ReadFile(std::string file_name, T&... args)
{
	CheckInputFile(file_name);
	std::ifstream file(file_name);
	
	((file >> args), ...);
}

double *ReadFileIntoArray(std::string file_name, const int size)
{
	CheckInputFile(file_name);
	double *buff = new double[size];
	std::ifstream file(file_name);
	
	for (int i = 0; i < size; i++) {file >> buff[i];};
	return buff;
}
