#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <string>

int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        std::cerr << "Usage: " << argv[0]
            << " <N> <filename>" << std::endl;
        return 1;
    }
    int n = atoi(argv[1]);
    std::string filename(argv[2]);
    std::cout << "It is an example, and the inputs are an integer " << n 
        << "\nand a string, which is a filename, is" << filename << std::endl;

	clock_t t1,t2;
	t1=clock();
	t2=clock();

	std::ofstream outfile("XXX.txt",std::ios::app);
	if(!outfile){
		std::cerr<<"hehe"<<std::endl;
	}
	outfile.precision(4);
	//outfile<<std::fixed;
	outfile<<std::showpos;
	outfile.setf(std::ios::scientific);
	outfile<<"	";
	outfile<<(t2-t1)*pow(10,-6)<<"\n";

    return 0;
}
