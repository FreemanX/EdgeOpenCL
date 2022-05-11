#include <unordered_map>
#include <iostream>
#include "csv.h"

std::unordered_map<std::string, std::vector<size_t>> b4Map;
std::unordered_map<std::string, std::vector<size_t>> b2Map;
std::unordered_map<std::string, std::vector<size_t>> b1Map;
std::unordered_map<std::string, std::vector<size_t>> l4Map;
std::unordered_map<std::string, std::vector<size_t>> gpuMap;

void printCleanWG(std::string &str_wg) {
  char lastChar = 'a';
  for (char &c : str_wg) {
	if (c == '0' && lastChar == ' ') {
	  std::cout << ']';
	  break;
	}
	std::cout << c;
	lastChar = c;
  }
}

int main() {
  io::CSVReader<6> in("KernelWG_dirty.csv");
  in.read_header(io::ignore_extra_column, "kernel", "b4", "b2", "b1", "l4", "GPU");
  std::string kernel;
  std::string b4;
  std::string b2;
  std::string b1;
  std::string l4;
  std::string gpu;
  std::cout << "kernel,b4,b2,b1,l4,GPU\n";
  while (in.read_row(kernel, b4, b2, b1, l4, gpu)) {
	std::cout << kernel;
	std::cout << ",";
	printCleanWG(b4);
	std::cout << ",";
	printCleanWG(b2);
	std::cout << ",";
	printCleanWG(b1);
	std::cout << ",";
	printCleanWG(l4);
	std::cout << ",";
	printCleanWG(gpu);
	std::cout << "\n";
  }
  return 0;
}