#include <getopt.h>
#include <stdio.h>

#ifndef __OCR__
	#define __OCR__
	#include "ocr.h"
#endif

#ifndef STREAM_TYPE
#	define STREAM_TYPE double
#endif

void printHelp() {
	PRINTF( "Usage:\n"
			"serial_stream [-d db_size] [-e] [-h] [-i num_iter]\n"
			"              [-r] [-s scalar_value] [-v]\n\n"
			"List of Options:\n" 
			"  -d|--db_size     = Size of data blocks a, b, c (1 million by default).\n"
			"  -e|--export      = Exports results to csv file.\n"
			"  -h|--help        = List of options.\n"
			"  -i|--iterations  = Number of iterations (1 by default).\n"
			"* -p|--split       = Splits db_size by given value (split value should be set so db_size % split = 0 and is 1 by default).\n"
			"  -r|--verify      = Verify results.\n"
			"  -s|--scalar      = Set scalar value (3.0 by default).\n"
			"  -v|--verbose     = Verbose output.\n"
			"Options with a * next to them are only available in parallel implementations\n");
	return;
}

int parseOptions(int argc, char ** argv, u64 * db_size, int * efile, u64 * iterations, int * verify, STREAM_TYPE * scalar, int * verbose) {
	char c;
	int help = 0;
	*db_size = 1000000;
	*efile = 0;
	*iterations = 1;
	*verify = 0;
	*scalar = 3.0;
	*verbose = 0;

	while (1) {
		int option_index = 0;
		static struct option long_options[] = {
			{"db_size",     required_argument,  0,              'd'},
			{"export",      no_argument,        0,              'e'},
			{"help",        no_argument,        0,              'h'},
			{"iterations",  required_argument,  0,              'i'},
			{"verify",      no_argument,        0,              'r'},
			{"scalar",      required_argument,  0,              's'},
			{"verbose",     no_argument,        0,              'v'},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "d:ehi:rs:v", long_options, &option_index);
		if (c == -1)
			break;
		switch (c) {
			case 'd':
				sscanf(optarg, "%llu", db_size);
				break;
			case 'e':
				*efile = 1;
				break;
			case 'h':
				printHelp();
				help = 1;
				break;
			case 'i':
				sscanf(optarg, "%llu", iterations);
				break;
			case 'r':
				*verify = 1;
				break;
			case 's':
				sscanf(optarg, "%lf", scalar);
				break;
			case 'v':
				*verbose = 1;
				break;
			case '?':
				PRINTF("Unknown option -%c.\n",optopt);
				break;
		}
	}
	return help;
}

int parseParallelOptions(int argc, char ** argv, u64 * db_size, int * efile, u64 * iterations, u64 * split, u64 * chunk, int * verify, STREAM_TYPE * scalar, int * verbose) {
	char c;
	int help = 0;
	*db_size = 1000000;
	*efile = 0;
	*iterations = 1;
	*split = 1;
	*verify = 0;
	*scalar = 3.0;
	*verbose = 0;

	while (1) {
		int option_index = 0;
		static struct option long_options[] = {
			{"db_size",     required_argument,  0,              'd'},
			{"export",      no_argument,        0,              'e'},
			{"help",        no_argument,        0,              'h'},
			{"iterations",  required_argument,  0,              'i'},
			{"split",       required_argument,  0,              'p'},
			{"verify",      no_argument,        0,              'r'},
			{"scalar",      required_argument,  0,              's'},
			{"verbose",     no_argument,        0,              'v'},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "d:ehi:n:rs:v", long_options, &option_index);
		if (c == -1)
			break;
		switch (c) {
			case 'd':
				sscanf(optarg, "%llu", db_size);
				break;
			case 'e':
				*efile = 1;
				break;
			case 'h':
				printHelp();
				help = 1;
				break;
			case 'i':
				sscanf(optarg, "%llu", iterations);
				break;
			case 'p':
				sscanf(optarg, "%llu", split);
				break;
			case 'r':
				*verify = 1;
				break;
			case 's':
				sscanf(optarg, "%lf", scalar);
				break;
			case 'v':
				*verbose = 1;
				break;
			case '?':
				PRINTF("Unknown option -%c.\n",optopt);
				break;
		}
	}
	*chunk = *db_size / *split;
	return help;
}