#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <mpi.h>

#define BLOCK_LOW(id,p,n) (((id)*(n))/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n) - 1)
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id+1),p,n) - BLOCK_LOW(id,p,n))
#define BLOCK_OWNER(index,p,n) (((p)*((index)+1)-1)/(n))

#define MERGE_TAG 0
#define BUFFER_SIZE 1000
#ifndef PRINT_FLAG
#define PRINT_FLAG 1
#endif

using namespace std;


typedef pair<string, int> PAIR;


struct comparePAIR {
  bool operator()(const PAIR& lhs, const PAIR& rhs) {
    return lhs.second > rhs.second;
  }
};

bool isLetter(char c) {
    if ( (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')) { 
        return true;} else {
        return false;}
}

void wordCountSingle(FILE *file_inptr, int size, map<string, int> &dict) {
    char c;
    bool is_letter = false;
    char word[BUFFER_SIZE];
    int letter_idx = 0, n_chars = 0;

    // wordCountSingle algorithm
    while ((c=fgetc(file_inptr)) != EOF) {
        if (isLetter(c)) {
            is_letter = true;word[letter_idx++] =  tolower(c);    
        } else {
            if (is_letter) {
                is_letter = false;word[letter_idx] = '\0';
                letter_idx = 0;string new_word(word);
                if (dict.find(new_word) == dict.end()) {
                    dict[new_word] = 1;} else {
                    dict[new_word] += 1;
                }}}
        if ((++n_chars) == size)break;
    }
}

