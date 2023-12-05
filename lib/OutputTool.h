#pragma once

#include "ErrorHandler.h"

#include <string>
#include <vector>

template <typename... T>
void Print(T... args)
{
	((std::cout << args << " "), ...);
	std::cout << std::endl;
}

void PrintInfo(std::string info_line)
{
	std::cout << OutputColor::green << "INFO: " << OutputColor::reset << info_line << std::endl;
}

int utf8_strlen(const std::string& str)
{
	int len = 0;
	for (int i=0; i < str.length(); i++, len++)
	{
		const unsigned char c = (unsigned char) str[i];
		if      (c>=0   && c<=127) i+=0;
		else if ((c & 0xE0) == 0xC0) i+=1;
		else if ((c & 0xF0) == 0xE0) i+=2;
		else if ((c & 0xF8) == 0xF0) i+=3;
		else return 0;
	}
	return len;
}

void PrintSimpleSeparator(std::string left_edge = "|", 
		std::string body = "-", 
		std::string right_edge = "|",
		const int length = 80)
{
	std::cout << left_edge;
	for (int i = 0; i < length - utf8_strlen(left_edge) - utf8_strlen(right_edge); i++) std::cout << body;
	std::cout << right_edge << std::endl;
	std::cout << OutputColor::reset;
}

void PrintSeparator(std::string text, 
	std::string color = OutputColor::bold_cyan,
	std::string left_edge = "//", 
	std::string body = "-", 
	std::string right_edge = "//",
	const int length = 80)
{
	std::cout << left_edge;
	for (int i = 0; i < length/2 - utf8_strlen(text)/2 - utf8_strlen(left_edge)*2 - 1; i++) std::cout << body;
	std::cout << left_edge << " " << color << text << OutputColor::reset << " " << right_edge;
	for (int i = 0; i < length - length/2 - utf8_strlen(text) + utf8_strlen(text)/2 - utf8_strlen(right_edge)*2 - 1; i++) std::cout << body;
	std::cout << right_edge << std::endl;
}

void PrintEdgedLine(std::string entry1, std::string entry2,
	std::string left_edge = "|",
	std::string right_edge = "|",
	const int length = 80)
{
	std::cout << left_edge << " " << entry1;
	int space_size = length - utf8_strlen(entry1) - utf8_strlen(entry2) - utf8_strlen(left_edge) - utf8_strlen(right_edge) - 2;
	for (int i = 0; i < space_size; i++) std::cout << " ";
	std::cout << entry2 << " " << right_edge << std::endl;
}

void PrintBigSeparator(std::string text,
	std::string color = OutputColor::bold_cyan,
	std::string ul_corner = "╓", 
	std::string ur_corner = "╖",
	std::string horizontal_line = "─",
	std::string vertical_line = "║",
	std::string dl_corner = "╙",
	std::string dr_corner = "╜")
{
	PrintSimpleSeparator(" " + ul_corner, horizontal_line, ur_corner);
	PrintSeparator(text, color, " " + vertical_line, " ", vertical_line);
	PrintSimpleSeparator(" " + dl_corner, horizontal_line, dr_corner);
}
