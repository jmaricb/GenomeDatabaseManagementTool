//
//  main.cpp
//  MetagenomicClassifier
//
//  Created by Josip Maric on 21/01/2020.
//  Copyright © 2020 Josip Maric. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <map>
#include <set>
#include <vector>
#include <tuple>
#include <dirent.h>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <pthread.h>

using namespace std;

#define NUM_THREADS 8
#define WINDOW 5
#define KMER_LEN 14
int32_t front_zero_mask = 268435455;

char rev_compl(char c) {
    if (c == 'G' || c == 'g') {
        return 'C';
    } else if(c == 'C' || c == 'c') {
        return 'G';
    } else if(c == 'A' || c == 'a') {
        return 'T';
    } else {
        return 'A';
    }
}

void reverse_complement(std::string &read, std::string &r_compl) {
    for(int i = read.size()-1; i >= 0; i--) {
        r_compl.push_back(rev_compl(read[i]));
    }
}

void progress_buffer(int32_t &buffer, char letter) {
    if (letter == 'A' || letter == 'a') {
        buffer |= 0;
    } else if(letter == 'C' || letter == 'c') {
        buffer |= 1;
    } else if(letter== 'G' || letter == 'g') {
        buffer |= 2;
    } else if(letter == 'T' || letter == 't') {
        buffer |= 3;
    }
}

struct thread_data {
    int  thread_id;
    std::vector<std::string> tax_ids;
    std::vector<std::string> dirs;
    uint64_t** index;
};

struct thread_data2 {
    int  thread_id;
    int64_t **index;
};

class Classifier {
public:
    Classifier() {
    }

    void classify_read(std::string &read, uint64_t** index, int number_of_partitions, int* results) {

//        int results[number_of_partitions*64];

        for (int i = 0; i < number_of_partitions*64; i++) {
            results[i] = 0;
        }

        int32_t buffer = 0;
        for (int i = 0; i < KMER_LEN; i++) {
            if (i != 0) {
                buffer = buffer << 2;
            }
            if(read[i] == 'N') {
                continue;
            }
            progress_buffer(buffer, read[i]);

            for(int i = 0; i < number_of_partitions; i++) {
                for(int k = 0; k < 64; k++) {
                    uint64_t mask = 1;
                    mask = mask << (63 - k);
                    uint64_t rez = index[buffer][i] & mask;
                    rez = rez >> (63 - k);
                    results[i * 64 + k] += 1;
                }
            }
        }

        for (int i = 14; i < read.size(); i++) {
            if(read[i] == 'N') {
                continue;
            }
            buffer = buffer << 2;
            progress_buffer(buffer, read[i]);
            buffer &= front_zero_mask;

            for(int i = 0; i < number_of_partitions; i++) {
                for(int k = 0; k < 64; k++) {
                    uint64_t mask = 1;
                    mask = mask << (63 - k);
                    uint64_t rez = index[buffer][i] & mask;
                    rez = rez >> (63 - k);
                    results[i * 64 + k] += rez;
                }
            }
        }
    }
};


void *ExecuteThat(void *threadarg) {
    struct thread_data *my_data;
    my_data = (struct thread_data *) threadarg;
    int thread_id = my_data->thread_id;
    std::vector<std::string> directories = my_data->dirs;

    int block_size = (int) directories.size() / NUM_THREADS;
    int upper_limit = thread_id == NUM_THREADS-1 ? directories.size() : (thread_id + 1) * block_size;

//    std::cout << "Thread: " << thread_id << "  (" << (thread_id * block_size) << ", " <<  upper_limit << ")" << std::endl;

    for(int i = thread_id * block_size; i < upper_limit; i++ ) {
        std::ifstream myfile (directories[i]);
        if (!myfile.is_open()) {
            std::cout << "cannot open " << directories[i] << std::endl;
        }

        int number_of_partition = i / 64;
        int position_in_partition = i % 64;
        uint64_t mask = 1;
        mask = mask << (63 - position_in_partition);

//        std::cout <<"Thread: " << thread_id << " " << directories[i] << " ---  " << (i-thread_id * block_size) << std::endl;

        std::string line;
        int is_start = 0;

        int32_t buffer = 0;
        int counter = 0;

        while ( getline (myfile, line) ) {
            if (line.size() > 0 && line[0] == '>') {
                is_start = 1;
            } else {
                if(is_start) {
                    is_start = 0;
                    buffer = 0;
                    for (int i = 0; i < KMER_LEN; i++) {
                        if (i != 0) {
                            buffer = buffer << 2;
                        }
                        if(line[i] == 'N') {
                            continue;
                        }
                        progress_buffer(buffer, line[i]);
                        my_data->index[buffer][number_of_partition] |= mask;
                    }

                    for (int i = KMER_LEN; i < line.size(); i++) {
                        if(line[i] == 'N') {
                            continue;
                        }
                        buffer = buffer << 2;
                        progress_buffer(buffer, line[i]);
                        buffer &= front_zero_mask;
                        my_data->index[buffer][number_of_partition] |= mask;
                        counter += 1;
                    }
                } else {
                    for (int i = 0; i < line.size(); i++) {
                        if(line[i] == 'N') {
                            continue;
                        }
                        buffer = buffer << 2;
                        progress_buffer(buffer, line[i]);
                        buffer &= front_zero_mask;
                        my_data->index[buffer][number_of_partition] |= mask;
                        counter += 1;
                    }
                }
            }
        }
        myfile.close();
    }

//    std::cout << "Thread with id : " << thread_id << "  ...exiting " << std::endl;
    pthread_exit(NULL);
}

void analyseResults(std::string path1, std::string path3, std::string path2) {
    std::ifstream myfile1 (path1);
    if (!myfile1.is_open()) {
        std::cout << "cannot open " << path1 << std::endl;
    }
    
    std::ifstream myfile3 (path3);
    if (!myfile3.is_open()) {
        std::cout << "cannot open " << path3 << std::endl;
    }


    std::ofstream outfile_test;
    outfile_test.open ("eval.txt");
    
    std::vector<std::vector<int>> results;

    std::string line;
    while ( getline (myfile1, line) ) {
        std::string delimiter = "\t";

        size_t pos = 0;
        std::string token;

        std::vector<int> one_line_rez;

        while ((pos = line.find(delimiter)) != std::string::npos) {
            token = line.substr(0, pos);
            std::string::size_type sz;
            int rez = std::stoi(token, &sz);
            one_line_rez.emplace_back(rez);
            line.erase(0, pos + delimiter.length());
        }

        results.emplace_back(one_line_rez);
    }
    
    std::vector<std::vector<int>> results_comp;

    while ( getline (myfile3, line) ) {
        std::string delimiter = "\t";

        size_t pos = 0;
        std::string token;

        std::vector<int> one_line_rez;

        while ((pos = line.find(delimiter)) != std::string::npos) {
            token = line.substr(0, pos);
            std::string::size_type sz;
            int rez = std::stoi(token, &sz);
            one_line_rez.emplace_back(rez);
            line.erase(0, pos + delimiter.length());
        }

        results_comp.emplace_back(one_line_rez);
    }

    std::cout << "results compl size " << results_comp.size() << std::endl;

    std::vector<int> legend;

    std::ifstream myfile2 (path2);
    if (!myfile2.is_open()) {
        std::cout << "cannot open " << path2 << std::endl;
    }

    while ( getline (myfile2, line) ) {
        std::string delimiter = "\t";

        size_t pos = 0;
        std::string token;
        int index = 0;

        while ((pos = line.find(delimiter)) != std::string::npos) {
            token = line.substr(0, pos);

            if(index == 1) {
                std::string::size_type sz;
                int tax_id = std::stoi(token, &sz);
                legend.emplace_back(tax_id);
                break;
            }

            line.erase(0, pos + delimiter.length());
            index += 1;
        }
    }

    std::cout << "legend size " << legend.size() << std::endl;

    for(int i = 0; i < results.size(); i++) {
        std::vector<int> rez = results[i];
        
        int biggestFive[5];
        for(int i = 0; i < 5; i++) {
            biggestFive[i] = 0;
        }
        int biggestFiveIndex[5];

        for(int j = 0; j < rez.size(); j++) {
            int current = rez[j];
            int curr_index = 0;
            while (curr_index < 5) {
                if(biggestFive[curr_index] < current) {
                    break;
                }
                curr_index += 1;
            }
            if(curr_index < 5) {
                int prev = biggestFive[curr_index];
                int prev_index = biggestFiveIndex[curr_index];
                biggestFive[curr_index] = current;
                biggestFiveIndex[curr_index] = j;
                curr_index += 1;
                while (curr_index < 5) {
                    int tmp = biggestFive[curr_index];
                    int tmp_index = biggestFiveIndex[curr_index];
                    biggestFive[curr_index] = prev;
                    biggestFiveIndex[curr_index] = prev_index;
                    prev = tmp;
                    prev_index = tmp_index;
                    curr_index += 1;
                }
            }
        }
        
        std::vector<int> rez2 = results_comp[i];

        int biggestFive2[5];
        for(int i = 0; i < 5; i++) {
            biggestFive2[i] = 0;
        }
        int biggestFiveIndex2[5];

        for(int j = 0; j < rez2.size(); j++) {
            int current = rez2[j];
            int curr_index = 0;
            while (curr_index < 5) {
                if(biggestFive2[curr_index] < current) {
                    break;
                }
                curr_index += 1;
            }
            if(curr_index < 5) {
                int prev = biggestFive2[curr_index];
                int prev_index = biggestFiveIndex2[curr_index];
                biggestFive2[curr_index] = current;
                biggestFiveIndex2[curr_index] = j;
                curr_index += 1;
                while (curr_index < 5) {
                    int tmp = biggestFive2[curr_index];
                    int tmp_index = biggestFiveIndex2[curr_index];
                    biggestFive2[curr_index] = prev;
                    biggestFiveIndex2[curr_index] = prev_index;
                    prev = tmp;
                    prev_index = tmp_index;
                    curr_index += 1;
                }
            }
        }

        int biggestTen[10];
        int biggestTenIndex[10];

        int firstI = 0;
        int secondI = 0;
        
        while(firstI < 5 && secondI < 5) {
            if(biggestFive[firstI] > biggestFive2[secondI]) {
                biggestTen[firstI + secondI] = biggestFive[firstI];
                biggestTenIndex[firstI + secondI] = biggestFiveIndex[firstI];
                firstI += 1;
            } else {
                biggestTen[firstI + secondI] = biggestFive2[secondI];
                biggestTenIndex[firstI + secondI] = biggestFiveIndex2[secondI];
                secondI += 1;
            }
        }
        
        while (firstI < 5) {
            biggestTen[firstI + secondI] = biggestFive[firstI];
            biggestTenIndex[firstI + secondI] = biggestFiveIndex[firstI];
            firstI += 1;
        }
        
        while (secondI < 5) {
            biggestTen[firstI + secondI] = biggestFive2[secondI];
            biggestTenIndex[firstI + secondI] = biggestFiveIndex2[secondI];
            secondI += 1;
        }
        
        outfile_test << "### " << i << std::endl;
        for(int i = 0; i < 10; i++) {
            outfile_test << biggestTen[i] << "\t" << legend[biggestTenIndex[i]] << std::endl;
        }
    }
    outfile_test.close();
}

int main(int argc, const char * argv[]) {
    
    // READING READS
    
    const char *path2 = argv[2];
    std::ifstream myfile (path2);
    if (!myfile.is_open()) {
        std::cout << "cannot open " << path2 << std::endl;
    }

    std::vector<std::string> reads;

    std::string line;
    while ( getline (myfile, line) ) {
        if (line.rfind(">", 0) == 0) {
            getline (myfile, line);
            reads.emplace_back(line);
        }
    }

    std::cout << "Read " << reads.size() << " reads" << std::endl;
    
    myfile.close();

    // READING DATABASE
    
    const char *path = argv[1];

    std::vector<std::string> directories;
    std::vector<std::string> tax_ids;

    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (path)) != NULL) {
        while ((ent = readdir (dir)) != NULL) {
            std::string test = ent->d_name;

            if(test.compare("..") == 0 || test.compare(".") == 0) {
                continue;
            }

            DIR *dir_tax;
            struct dirent *ent_tax;
            std::string total_path = path;
            total_path = total_path + "/" + ent->d_name;

            if ((dir_tax = opendir (total_path.c_str())) != NULL) {
                while ((ent_tax = readdir (dir_tax)) != NULL) {
                    std::string test_tax = ent_tax->d_name;
                    if(test_tax.compare("..") == 0 || test_tax.compare(".") == 0) {
                        continue;
                    }
                    std::string total_path_tax = total_path + "/" + ent_tax->d_name;

                    directories.emplace_back(total_path_tax);
                    tax_ids.push_back(ent->d_name);
                    break;
                }
                closedir (dir_tax);
            } else {
                perror ("error cannot open");
            }
        }
        closedir (dir);
    } else {
        perror ("error cannot open");
    }
    
    std::cout << "Read " << tax_ids.size() << " tax ids" << std::endl;
    
    // WRITING DATABASE LEGEND
    
    ofstream myfile_legend;
    myfile_legend.open ("legend.txt");

    for(int i = 0; i < directories.size(); i++) {
        myfile_legend << i << "\t" << tax_ids[i] << "\t" << directories[i] << std::endl;
    }

    myfile_legend.close();
    
    std::cout << "Wrote legend file" << std::endl;
    
    // CREATING INDEX
    
    std::cout << "Allocating space for inndex" << std::endl;

    int key_l = KMER_LEN;
    int limit = pow(2, 2*key_l);
    uint64_t** index;
    index = new uint64_t*[limit];
    int number_of_partitions = directories.size() / 64 + 1;
    
//    uint64_t** index2;
//    index2 = new uint64_t*[limit];

//    std::fstream infile_binary;
//    infile_binary.open("index.bin", ios:: binary);
//
//    for(int i = 0; i < limit; i++) {
//        index[i] = new uint64_t[number_of_partitions];
//        infile_binary.read((char *) index[i], sizeof (uint64_t) * number_of_partitions);
//    }
//
//    infile_binary.close();

    for(int i = 0; i < limit; i++) {
        index[i] = new uint64_t[number_of_partitions];
        for(int j = 0; j < number_of_partitions; j++) {
            index[i][j] = 0;
        }
    }

    pthread_t threads[NUM_THREADS];
    struct thread_data td[NUM_THREADS];
    int rc;

    std::cout << "Creating index" << std::endl;

    for( int i = 0; i < NUM_THREADS; i++ ) {
//        cout <<"main() : creating thread, " << i << endl;
        td[i].thread_id = i;
        td[i].dirs = directories;
        td[i].index = index;
        td[i].tax_ids = tax_ids;
        rc = pthread_create(&threads[i], NULL, ExecuteThat, (void *)&td[i]);
        if (rc) {
            cout << "Error:unable to create thread," << rc << endl;
            exit(-1);
        }
    }

    void *status;

    for( int i = 0; i < NUM_THREADS; i++ ) {
        rc = pthread_join(threads[i], &status);
        if (rc) {
            cout << "Error:unable to join," << rc << endl;
            exit(-1);
        }
//        cout << "Main: completed thread id :" << i ;
    }

    std::cout << "Index created" << std::endl;
    
//    std::cout << "Storing index" << std::endl;

//    std::ofstream outfile_binary;
//    outfile_binary.open("index.bin", ios:: binary);
//    for(int i = 0; i < limit; i++) {
//        outfile_binary.write((char *) index[i], sizeof(uint64_t) * number_of_partitions);
//    }
//    outfile_binary.close();

//    std::cout << "Index stored" << std::endl;
        
//    std::cout << "Reading index" << std::endl;
    
//    uint64_t** index2;
//    index2 = new uint64_t*[limit];
//
//    std::fstream infile_binary;
//    infile_binary.open("index.bin", ios:: binary);
//
//    for(int i = 0; i < limit; i++) {
//        index2[i] = new uint64_t[number_of_partitions];
//        infile_binary.read((char *) index2[i], sizeof (uint64_t) * number_of_partitions);
//    }
//
//    infile_binary.close();
    
//    std::cout << "Index read" << std::endl;
    
//    std::ofstream outfile_test;
//    outfile_test.open ("test.txt");
//    for(int j = 0; j < 10; j++) {
//        for (int i = 0; i < number_of_partitions; i++) {
//            outfile_test << index[j][i] << " " << index[j][i] << std::endl;
//        }
//    }
//
//    outfile_test.close();

    std::cout << "Clasiffying" << std::endl;
    
    Classifier classifier = Classifier();
    for (int i = 0; i < reads.size(); i++) {
//        std::cout << "processing read " << i << " with partition size: " << number_of_partitions << std::endl;
        std::string r_compl;
        reverse_complement(reads[i], r_compl);

        int results[number_of_partitions*64];
        int results_comp[number_of_partitions*64];
        
        classifier.classify_read(reads[i], index, number_of_partitions, results);
        classifier.classify_read(r_compl, index, number_of_partitions, results_comp);
        
        std::ofstream myfile_out;
        std::string filename = "rezults_whole.txt";
        myfile_out.open (filename, std::ios_base::app);
        for(int i = 0; i < number_of_partitions * 64; i++) {
            myfile_out << results[i] << "\t";
        }
        myfile_out << std::endl;
        myfile_out.close();
        
        std::ofstream myfile_out2;
        std::string filename2 = "rezults_whole_compl.txt";
        myfile_out2.open (filename2, std::ios_base::app);
        for(int i = 0; i < number_of_partitions * 64; i++) {
            myfile_out2 << results_comp[i] << "\t";
        }
        myfile_out2 << std::endl;
        myfile_out2.close();
    }

    for(int i = 0; i < limit; ++i) {
        delete[] index[i];
    }
    delete[] index;

    std::cout << "Deleted memory" << std::endl;
    
    analyseResults("rezults_whole.txt", "rezults_whole_compl.txt",  "legend.txt");
        
    std::cout << "Finished " << std::endl;

    return 0;
}
 
