#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io.h"

#define DIM	3


// Write in simple text format for Python:
void write_txt_(
	char *Filename, 
	INT *N,
	REAL *X,
	int namelen
) {
	int k, n = *N;
	char filename[256];
	strncpy(filename, Filename, namelen);
	filename[namelen]='\0';
printf("Writing %d records to %s\n",n,filename);fflush(stdout);
	FILE *out = fopen(filename, "w");
	k = 0;
	for (int i=0; i<n; i++) {
		fprintf(out,"%g",X[k++]);
		for (int j=1; j<DIM; j++) {
			fprintf(out,"\t%g",X[k++]);
		}
		fprintf(out,"\n");
	}
	fclose(out);
}

// Write in simple text format for Qhull:
void write_dat_(
	char *Filename, 
	INT *N,
	REAL *X,
	int namelen
) {
	INT k, n = *N;
	char filename[256];
	strncpy(filename, Filename, namelen);
	filename[namelen]='\0';
	printf("Writing %d records to %s\n",n,filename);fflush(stdout);
	FILE *out = fopen(filename, "w");
   fprintf(out,"%d\n%d\n",DIM,n);fflush(out);
	k = 0;
	for (INT i=0; i<n; i++) {
		fprintf(out,"%g",X[k++]);
		for (int j=1; j<DIM; j++) {
			fprintf(out,"\t%g",X[k++]);
		}
		fprintf(out,"\n");
	}
	fclose(out);
}

