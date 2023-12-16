#pragma once

namespace Tool
{
	template <typename... Ts> 
	double Maximum(Ts... args)
	{
		constexpr int size = sizeof...(args);
		double entries[size] = {static_cast<double>(args)...};
		double result = entries[0];
		for (double val : entries) if (val > result) result = val;
		return result;
	}

	template <typename... Ts> 
	double Minimum(Ts... args)
	{
		constexpr int size = sizeof...(args);
		double entries[size] = {static_cast<double>(args)...};
		double result = entries[0];
		for (double val : entries) if (val < result) result = val;
		return result;
	}

	template <typename... Ts> 
	double Average(Ts... args)
	{
		constexpr int size = sizeof...(args);
		double entries[size] = {static_cast<double>(args)...};
		double result = 0.;
		for (double val : entries) result += val/static_cast<double>(size);
		return result;
	}
}
