#include <iostream>
#include <limits>

int main() {
	std::cout << "Hello, World!" << std::endl;
	std::cout << "\nPress Enter to exit...";
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	return 0;
}