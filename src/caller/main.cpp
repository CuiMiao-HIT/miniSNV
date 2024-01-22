// verson : 20230725

#include "../cmtools/desc.hpp"
#include "callerstep1A3.hpp"

using namespace std;

int usage(void) {
	fprintf(stderr, "\n");
    fprintf(stderr, "Program:   %s\n", PACKAGE_NAME); 
    fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION); 
    fprintf(stderr, "Contact:   %s\n", CONTACT); 

    fprintf(stderr, "Usage:    ./%s <Command> [Options]\n", PACKAGE_NAME); 
	fprintf(stderr, "Command: ");
	// fprintf(stderr, "\tindex	    index reference sequence\n");
	fprintf(stderr, "\t\tcaller	    call snp\n");
	fprintf(stderr, "\n");

    return 2;
}

int main(int argc, char *argv[]){
    int r = 0;
	if (argc < 2) return usage();
	// if (strcmp(argv[1], "index") == 0)	r = run_index(argc, argv);
	else if (strcmp(argv[1], "caller") == 0){
		if (strcmp(argv[2], "callerStep1") == 0){
			r = run_callerStep1(argc, argv);
			return r;
		} else if (strcmp(argv[2], "callerStep3") == 0){
			r = run_callerStep3(argc, argv);
			return r;
		}
	}
	else if (strcmp(argv[1], "--help") == 0) return usage();
	else {
		fprintf(stderr, "[Waring!!!] wrong command: '%s'\n", argv[1]);
		return 2;
	}
    return 0;
}
