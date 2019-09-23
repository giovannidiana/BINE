

/*
Last line from the file
*/

#include"../include/ReadLine.h"
#include <fstream>
#include<sstream>

using namespace std;

void rline(string filename,int* ids,int N)
{
    std::ifstream read(filename, std::ios_base::ate);//open file
    std::string line;

    int length = 0;
    int ind;
    int i;

    char c = '\0';

    if( read )
    {
        length = read.tellg();//Get file size

        // loop backward over the file

        for(i = length-2; i > 0; i-- )
        {
            read.seekg(i);
            c = read.get();
            if( c == '\r' || c == '\n' )//new line?
                 break;
        }

        std::getline(read, line);//read last line
        stringstream ss(line);
        i=0;while(ss>>ind){
            if(i>N-1) exit(1);
            ids[i]=ind;
            i++;
        }

    }

    read.close();

}

