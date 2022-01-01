/*
 * vector_filterMain.c
 *
 *  Created on: Nov 12, 2013
 *      Author: hxin
 */

//#ifndef BOOST_PP_IS_ITERATING
//#include "print.h"
#include <string>
#include <sys/times.h>
#include <unistd.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SIMD_ED.h"
#include "SHD.h"

#define BATCH_RUN 1000000
#ifndef _MAX_LENGTH_
#define _MAX_LENGTH_ 256
#endif

using namespace std;

//char read[128];
//char ref[128];

char init_all_NULL[128] = "";

//char read_t[128] __aligned;// = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
//char ref_t[128] __aligned;// = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

int main(int argc, char* argv[]) {

    string *read_strs = new string [BATCH_RUN];
    string *ref_strs = new string [BATCH_RUN];
    bool *valid_buff = new bool [BATCH_RUN];

    uint8_t **read0 = new uint8_t* [BATCH_RUN];
    uint8_t **read1 = new uint8_t* [BATCH_RUN];
    uint8_t **ref0 = new uint8_t* [BATCH_RUN];
    uint8_t **ref1 = new uint8_t* [BATCH_RUN];
    size_t *length = new size_t [BATCH_RUN];

    for (int i = 0; i < BATCH_RUN; i++) {
        read0[i] = new uint8_t[_MAX_LENGTH_ / 8];
        read1[i] = new uint8_t[_MAX_LENGTH_ / 8];
        ref0[i] = new uint8_t[_MAX_LENGTH_ / 8];
        ref1[i] = new uint8_t[_MAX_LENGTH_ / 8];
    }
/*
*/
    if (argc < 2) {
        printf("Usage: $>bin error [use SHD, use_levenshtein] \n");
        exit(1);
    }
    int use_SHD = -1; // use default as default.
    if (argc > 2) {
        use_SHD = atoi(argv[2]);
    }
    int use_levenshtein = 1; // use levenshtein as default.
    if (argc > 3) {
        use_levenshtein = atoi(argv[3]);
    }


    int error = atoi(argv[1]);

    size_t lineLength;
    char* tempstr = NULL;

    long long unsigned int passNum = 0;
    long long unsigned int totalNum = 0;

    long long read_size;
    long long read_idx;
    bool stop = false;

    tms start_time;
    tms end_time;
    tms elp_time;

    elp_time.tms_stime = 0;
    elp_time.tms_utime = 0;
    elp_time.tms_cstime = 0;
    elp_time.tms_cutime = 0;

    SIMD_ED ed_obj;
    if (use_levenshtein) {
        bool tmp = (use_SHD == -1? true : use_SHD);
        ed_obj.init_levenshtein(error, ED_GLOBAL, tmp);
    }
    else {
        bool tmp = (use_SHD == -1? false : use_SHD);
        ed_obj.init_affine(error, error * 3, ED_GLOBAL, 2, 3, 1, tmp);
    }

    //ed_obj.init_levenshtein(error, ED_LOCAL, false);
    //ed_obj.init_affine(error, error * 3, ED_GLOBAL, 2, 3, 1, true);

    do {
        //clear past result
//		strncpy(read, init_all_NULL, 128);
//		strncpy(ref, init_all_NULL, 128);

        for (read_size = 0; read_size < BATCH_RUN; read_size++) {

            //get read
            getline(&tempstr, &lineLength, stdin);
            int length_t = strlen(tempstr) - 1;
            length[read_size] = length_t;
            //Get rid of the new line character
            tempstr[length_t] = '\0';

            if (strcmp(tempstr, "end_of_file\0") == 0) {
                stop = true;
                break;
            }
            read_strs[read_size].assign(tempstr);

            //get ref
            getline(&tempstr, &lineLength, stdin);
            length_t = strlen(tempstr) - 1;
            //length[read_size] = length_t;
            //Get rid of the new line character
            tempstr[length_t] = '\0';
            ref_strs[read_size].assign(tempstr);
            valid_buff[read_size] = false;

            if (length[read_size] > 256)
                length[read_size] = 256;

            //cout << "length: " << length[read_size] << endl;

            ed_obj.convert_reads((char*) read_strs[read_size].c_str(), (char*) ref_strs[read_size].c_str(),
                                 length[read_size], read0[read_size], read1[read_size], ref0[read_size], ref1[read_size]);

        }

        times(&start_time);

        for (read_idx = 0; read_idx < read_size; read_idx++) {

            //int length_t = read_strs[read_idx].length();
            int length_t = length[read_idx];
/*
			__m256i read0_YMM = _mm256_loadu_si256((__m256i*) read0[read_idx]);
			__m256i read1_YMM = _mm256_loadu_si256((__m256i*) read1[read_idx]);
			__m256i ref0_YMM = _mm256_loadu_si256((__m256i*) ref0[read_idx]);
			__m256i ref1_YMM = _mm256_loadu_si256((__m256i*) ref1[read_idx]);

#ifdef debug
			cout << "read0: ";
			for (int i = 0; i < _MAX_LENGTH_ / 8; i++) {
				cout << (int) read0[read_idx][i] << " ";
			}
			cout << endl;
			print256_bit(read0_YMM);

			cout << "read1: ";
			for (int i = 0; i < _MAX_LENGTH_ / 8; i++) {
				cout << (int) read1[read_idx][i] << " ";
			}
			cout << endl;
			print256_bit(read1_YMM);

			cout << "ref0: ";
			for (int i = 0; i < _MAX_LENGTH_ / 8; i++) {
				cout << (int) ref0[read_idx][i] << " ";
			}
			cout << endl;
			print256_bit(ref0_YMM);

			cout << "ref1: ";
			for (int i = 0; i < _MAX_LENGTH_ / 8; i++) {
				cout << (int) ref1[read_idx][i] << " ";
			}
			cout << endl;
			print256_bit(ref1_YMM);
#endif
*/


            //valid_buff[read_idx] = true; }
            //ed_obj.load_reads((char*) read_strs[read_idx].c_str(), (char*) ref_strs[read_idx].c_str(), length_t);
            ed_obj.load_reads(read0[read_idx], read1[read_idx], ref0[read_idx], ref1[read_idx], length_t);
            ed_obj.calculate_masks();
            ed_obj.reset();
            ed_obj.run();
            if (ed_obj.check_pass() ) {
                //ed_obj.backtrack();
                //fprintf(stderr, "%.*s\n", 128, ed_obj.get_CIGAR().c_str() );
                valid_buff[read_idx] = true;
            }
/*
*/
/*
			ed_obj.load_reads(read0[read_idx], read1[read_idx], ref0[read_idx], ref1[read_idx], length_t);
			ed_obj.calculate_masks();
			ed_obj.reset();
			ed_obj.run();
			if (!ed_obj.check_pass() ) {
				ed_obj.load_reads((char*) read_strs[read_idx].c_str(), (char*) ref_strs[read_idx].c_str(), length_t);
				ed_obj.calculate_masks();
				ed_obj.reset();
				ed_obj.run();
				if (ed_obj.check_pass() ) {
					fprintf(stdout, "%.*s\n", 128, read_strs[read_idx].c_str() );
					fprintf(stdout, "%.*s\n", 128, ref_strs[read_idx].c_str() );
				}
			}

			else {
				//ed_obj.load_reads(read0[read_idx], read1[read_idx], ref0[read_idx], ref1[read_idx], length_t);
				ed_obj.load_reads((char*) read_strs[read_idx].c_str(), (char*) ref_strs[read_idx].c_str(), length_t);
				ed_obj.calculate_masks();
				ed_obj.reset();
				ed_obj.run();
				if (ed_obj.check_pass() ) {
					fprintf(stdout, "%.*s\n", 128, read_strs[read_idx].c_str() );
					fprintf(stdout, "%.*s\n", 128, ref_strs[read_idx].c_str() );
				}
			}
/*
			else {
				fprintf(stdout, "error!\n");
				fprintf(stdout, "%.*s\n", 128, read_strs[read_idx].c_str() );
				fprintf(stdout, "%.*s\n", 128, ref_strs[read_idx].c_str() );
				exit(1);
			}
			/*
			strncpy(read_t, init_all_NULL, 128);
			strncpy(ref_t, init_all_NULL, 128);

			length = read_strs[read_idx].length();

			if (length > 128)
				length = 128;
			strncpy(read_t, read_strs[read_idx].c_str(), length);

			length = ref_strs[read_idx].length();
			//Get rid of the new line character
			if (length > 128)
				length = 128;
			strncpy(ref_t, ref_strs[read_idx].c_str(), length);

			if (bit_vec_filter_sse1(read_t, ref_t, length, error))
				valid_buff[read_idx] = true;
			*/
        }

        times(&end_time);

        for (read_idx = 0; read_idx < read_size; read_idx++) {

            if (valid_buff[read_idx]) {
                //fprintf(stderr, "%.*s\n", 128, read_strs[read_idx].c_str() );
                //fprintf(stderr, "%.*s\n", 128, ref_strs[read_idx].c_str() );
                passNum++;
            }
            totalNum++;
        }


        elp_time.tms_stime += end_time.tms_stime - start_time.tms_stime;
        elp_time.tms_utime += end_time.tms_utime - start_time.tms_utime;

        if (stop)
            break;

    } while (1);

    fprintf(stderr, "end_of_file\n");
    printf("passNum:\t%lld\n", passNum);
    printf("totalNum:\t%lld\n", totalNum);
    printf("total_time: %f\n", (double) elp_time.tms_utime / sysconf(_SC_CLK_TCK) );

    delete [] read_strs;
    delete [] ref_strs;
    delete [] valid_buff;

    for (int i = 0 ; i < BATCH_RUN; i++) {
        delete [] read0[i];
        delete [] read1[i];
        delete [] ref0[i];
        delete [] ref1[i];
    }
    delete [] read0;
    delete [] read1;
    delete [] ref0;
    delete [] ref1;
    delete [] length;
/*
*/

    return 0;

}

//#endif
