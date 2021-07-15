#include "wordcount.h"

int main(int argc, char *argv[]) {
    int p, id, i, j, file_size, block_low, block_size;
    string filename = "./big_file/big_100.txt";
    FILE *file_inptr;
    map<string, int> dict;
    int n_words, n_chars, word_freq;
    char word[BUFFER_SIZE];
    MPI_Status status;
    double elapsed_time, elapsed_time_total = 0.0;
    int iter;
    int n_iter = 1;

    MPI_Init(&argc, &argv);
    if (argc > 1) {if (id == 0) 
    printf("[info] replacing `filename` with %s\n", argv[1]);
        filename = string(argv[1]);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    for (iter = 0; iter < n_iter; ++ iter) {

        MPI_Barrier(MPI_COMM_WORLD);elapsed_time = -MPI_Wtime();
        file_inptr = fopen(filename.c_str(), "r");fseek(file_inptr, 0, SEEK_END);
        file_size = (int)ftell(file_inptr);
        block_low = BLOCK_LOW(id, p, file_size);
        block_size = BLOCK_SIZE(id, p, file_size);
        
        // find forward until the first non-word character
        if (id) {
            if(isLetter(fgetc(file_inptr))){
                fseek(file_inptr, block_low, SEEK_SET);
                while(isLetter(fgetc(file_inptr))){
                    --block_low;
                    fseek(file_inptr, block_low, SEEK_SET);
                }
            }
            // fseek(file_inptr, block_low-1, SEEK_SET);
            // if (isLetter(fgetc(file_inptr)) && isLetter(fgetc(file_inptr))) {
            //     do {
            //         -- block_low;
            //         fseek(file_inptr, block_low, SEEK_SET);
            //     } while (isLetter(fgetc(file_inptr)));
            // }
        }

        wordCountSingle(file_inptr, BLOCK_LOW(id+1, p, file_size)-block_low-1, dict);

        fclose(file_inptr);
        for (i = 1; i < p; ++ i) {
            if (id == 0) {
            MPI_Recv(&n_words, 1, MPI_INT, i, MERGE_TAG, MPI_COMM_WORLD, &status);
        for (j = 0; j < n_words; ++ j) {
            MPI_Recv(&n_chars, 1, MPI_INT, i, MERGE_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(word, n_chars, MPI_CHAR, i, MERGE_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&word_freq, 1, MPI_INT, i, MERGE_TAG, MPI_COMM_WORLD, &status);
            string new_word(word);
            if (dict.find(new_word) == dict.end()) {dict[new_word] = word_freq;
            } else {dict[new_word] += word_freq;}
        }
        // MPI_Recv(&n_chars, 1, MPI_INT, i, MERGE_TAG, MPI_COMM_WORLD, &status);
        //     MPI_Recv(word, n_chars, MPI_CHAR, i, MERGE_TAG, MPI_COMM_WORLD, &status);
        //     MPI_Recv(&word_freq, 1, MPI_INT, i, MERGE_TAG, MPI_COMM_WORLD, &status);
        //     string new_word(word);
            } else if (id == i) {
        n_words = dict.size();
        MPI_Send(&n_words, 1, MPI_INT, 0, MERGE_TAG, MPI_COMM_WORLD);
        for (auto &x: dict) {
            n_chars = x.first.size() + 1;
            MPI_Send(&n_chars, 1, MPI_INT, 0, MERGE_TAG, MPI_COMM_WORLD);
            MPI_Send(x.first.c_str(), n_chars, MPI_CHAR, 0, MERGE_TAG, MPI_COMM_WORLD);
            word_freq = x.second;
            MPI_Send(&word_freq, 1, MPI_INT, 0, MERGE_TAG, MPI_COMM_WORLD);
        }
            }
        }if (PRINT_FLAG) {
            vector<PAIR> dict_vec(dict.begin(), dict.end());
            sort(dict_vec.begin(), dict_vec.end(), comparePAIR());
            
            if (id == 0) {
                int top_n = 20, i = 0;
                printf("most frequent words: TOP %d\n", top_n);
                for (auto &x: dict_vec) {
                    printf("%-10s%6d\n", x.first.c_str(), x.second);
                    if (++i == top_n)break;
                }
            }
        }MPI_Barrier(MPI_COMM_WORLD);elapsed_time += MPI_Wtime();
        elapsed_time_total += elapsed_time;
    }
    if (id == 0) {printf("Total time: %10.6f\n",
     elapsed_time_total/n_iter);}
    MPI_Finalize();
    return 0;
}
