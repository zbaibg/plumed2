#include <iostream>
#include <fstream>
#include <string>

// Simple test program that checks if the PRECAL_DIFF test suite ran successfully
int main() {
    std::ifstream logfile("test_output");
    if (!logfile.is_open()) {
        std::cerr << "ERROR: test_output not found" << std::endl;
        return 1;
    }
    
    std::string line;
    bool found_summary = false;
    
    // Look for the summary section in the log
    while (std::getline(logfile, line)) {
        if (line.find("Summary (first hill sigma):") != std::string::npos) {
            found_summary = true;
            break;
        }
    }
    
    if (!found_summary) {
        std::cerr << "ERROR: Test summary not found in output" << std::endl;
        return 1;
    }
    
    std::cout << "PRECAL_DIFF test completed successfully" << std::endl;
    return 0;
}