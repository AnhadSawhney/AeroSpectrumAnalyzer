#pragma once
#include <string>
#include <fstream>

class Options {
public:
	int readInt(std::string arg) {
		return std::stoi(query(arg));
	}

	double readDouble(std::string arg) {
		return std::stod(query(arg));
	}

	LPCWSTR readLPCWSTR(std::string arg) {
		std::string s = query(arg);
		return std::wstring(s.begin(), s.end()).c_str();
	}

	Options(std::string f) : file(f, std::ios::in) {
		//file.open(f, std::ios::in);
	};

	void close() {
		file.close();
	}
private:
	std::ifstream file;
	
	std::string query(std::string arg) {
		if (file.is_open()) {
			file.clear();
			file.seekg(0, std::ios::beg);
			std::string line;
			while (getline(file, line) && line.find(arg) != 0) {
				//std::cout << line.c_str();
			}
			return line.substr(line.find('=')+1);
		}
		return "\0";
	}
};