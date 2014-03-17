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
			"serial_stream [-d db_size] [-e file_name] [-h] [-i num_iter]\n"
			"              [-r] [-s scalar_value] [-v]\n\n"
			"List of Options:\n" 
			"-d|--db_size     = Size of data blocks a, b, c (1 million by default).\n"
			"-e|--export      = Exports results to csv file.\n"
			"-h|--help        = List of options.\n"
			"-i|--iterations  = Number of iterations (1 by default).\n"
			"-r|--verify      = Verify results.\n"
			"-s|--scalar      = Set scalar value (3.0 by default).\n"
			"-v|--verbose     = Verbose output.\n");
	return;
}

int parseOptions(int argc, char ** argv, u64 * db_size, char * efile, u64 * iterations, STREAM_TYPE * scalar, int * verify, int * verbose) {
	char c;
	int help = 0;
	efile = NULL;
	*iterations = 1;
	*db_size = 1000000;
	*scalar = 3.0;
	*verify = 0;
	*verbose = 0;

	while (1) {
		int option_index = 0;
		static struct option long_options[] = {
			{"db_size",     required_argument,  0,              'd'},
			{"export_to",   required_argument,  0,              'e'},
			{"help",        no_argument,        0,              'h'},
			{"iterations",  required_argument,  0,              'i'},
			{"scalar",      required_argument,  0,              's'},
			{"verify",      no_argument,        0,              'r'},
			{"verbose",     no_argument,        0,              'v'},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "e:hi:n:rs:v", long_options, &option_index);
		if (c == -1)
			break;
		switch (c) {
			case 'd':
				sscanf(optarg, "%llu", db_size);
				break;
			case 'e':
				sscanf(optarg, "%s", efile);
				break;
			case 'h':
				printHelp();
				help = 1;
				break;
			case 'i':
				sscanf(optarg, "%d", iterations);
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