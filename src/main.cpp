#include "variational_bayes.cpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unistd.h>

int main(int argc, char *argv[]){
    int i, opt;
    opterr = 0;

    string data_directory = "sim1";
    string sample_name = "sample-1";
    int start_J = 2;
    int iteration = 1;
    vector<string> ref_sig;
    ref_sig.reserve(100);
    string ABS_path = "./";
    string mode = "pe";
    bool with_noise = false; string temp;

    while((opt = getopt(argc, argv, "d:s:j:i:r:p:m:w:")) != -1){
        switch(opt){
            case 'd':
                data_directory = optarg;
                break;
            case 's':
                sample_name = optarg;
                break;
            case 'j':
                start_J = atoi(optarg);
                break;
            case 'i':
                iteration = atoi(optarg);
                break;
            case 'r':
                ref_sig.push_back(optarg);
                break;
            case 'p':
                ABS_path = optarg;
                break;
            case 'm':
                mode = optarg;
                break;
            case 'w':
                temp = optarg;
                if(temp == "with_noise") with_noise = true;
                else with_noise = false;
                break;
            default:
                cout << "Usage: hoge" << endl;
                exit(1);
        }
    }
    run_VB(data_directory, sample_name, start_J, iteration, ref_sig, ABS_path, mode, with_noise);
    cout << "Complete " << sample_name << " J:" << start_J << " Replicated:" << iteration << endl;
}

