#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "common.h"
#include "cmdline.h"
#include "processbam.h"
#include "processfq.h"
#include "config.h"
#include "external.h"
#include "tardis.h"
#include "bamonly.h"
#include "variants.h"
#include "splitread.h"
#include "vh_main.h"
#include "sonic/sonic.h"
#include "processrefgen.h"
#include "free.h"

FILE *logFile = NULL;

#define MAXCHAR 10000
#define EXONS 750993
#define TOKENS 11
#define CHR1 66348
#define CHR2 54739
#define CHR3 46233
#define CHR4 26279
#define CHR5 31829
#define CHR6 34319
#define CHR7 32873
#define CHR8 30008
#define CHR9 27685
#define CHR10 33021
#define CHR11 36311
#define CHR12 37747
#define CHR13 14013
#define CHR14 20497
#define CHR15 25045
#define CHR16 26591
#define CHR17 35271
#define CHR18 11648
#define CHR19 31536
#define CHR20 16018
#define CHR21 9026
#define CHR22 14360
#define CHRX 23334
#define CHRY 4854
#define CHRM 2
#define CHROMS 26

int compare(const void *s1, const void *s2)
{
  exon_info *e1 = *(exon_info **)s1;
  exon_info *e2 = *(exon_info **)s2;

  // printf("%d\n", e1->start - e2->start);
  return e1->start - e2->start;
}

int main( int argc, char** argv)
{
	bam_info** in_bams;
	parameters* params;
	configuration* cfg;
	int return_value;
	char* username;
	int i;
	time_t rawtime;
	struct tm * timeinfo;
	char *log_file_path;

	/* For loading exons from refFlat.txt */
	FILE *fp;
	char str[MAXCHAR];
	char *saveptr1, *saveptr2;
	char* filename;
	char * start;
	char * end;
	exon_info **in_exons;
	char ** tokens;
	int this_exon_code;
	/**/

	time ( &rawtime);
	timeinfo = localtime( &rawtime);

	print_quote();

	/* Set program parameters */
	init_params( &params);
	/* Parse command line arguments */	
	return_value = parse_command_line( argc, argv, params);


	if( return_value == 0){
		/*        fprintf(stderr, "parse_command_line returned 0.\n"); */
		exit(EXIT_SUCCESS);
	}
	else if( return_value != 1){
		/*      fprintf(stderr, "parse_command_line returned %d.\n", return_value); */
		exit( return_value);
	}

	/* Keeping simple logs in tardis.log file */
	log_file_path = (char *) getMem(sizeof(char) * (12+strlen(params->outprefix)+strlen(params->outdir)));
	sprintf(log_file_path, "%s%s-%s", params->outdir, params->outprefix, "tardis.log");
	logFile = safe_fopen (log_file_path, "w");
	fprintf( logFile, "#CreationDate=%d.%d.%d\n\n", timeinfo->tm_year + 1900, timeinfo->tm_mon + 1, timeinfo->tm_mday);
	free(log_file_path);

	fprintf( logFile, "Command line:\n \t%s ", argv[0]);
	for ( i = 1; i < argc; i++)
		fprintf( logFile, "%s ", argv[i]);
	fprintf( logFile, "\n\n");	

	/* Load configuration file (created if it does not exist) */
	cfg = ( configuration*) getMem( sizeof( configuration));
	load_config( cfg, params);

	/* make_sonic is standalone. Execute and return.  */
	if ( params->make_sonic)
		return sonic_build(params->ref_genome, params->gaps, params->reps, params->dups, params->sonic_info, params->sonic_file);

	/* Load SONIC */
	params->this_sonic = sonic_load(params->sonic_file);

	if (params->last_chr < params->first_chr)
		params->last_chr = params->this_sonic->number_of_chromosomes - 1;

	if (! compare_sonic_ref(params))
	{
		fprintf(stderr, "Reference FASTA and SONIC file do not match. Check if you are using the same version.\n");
		exit(EXIT_WRONG_SONIC);	    
	}


	if ( TARDIS_DEBUG)
		print_params( params);

	print_quote();

	/* Read BAM files and calculate the median/avg/std of fragment sizes per library */
	in_bams = ( bam_info**) getMem( (params->num_bams) * sizeof( bam_info*));
	for( i = 0; i < params->num_bams; i++)
	{
		in_bams[i] = ( bam_info*) getMem( sizeof( bam_info));
		in_bams[i]->sample_name = NULL;
		load_bam( params, cfg, in_bams[i], params->bam_file_list[i], params->alt_mapping, i, params->ref_genome);
	}

	/* Load genes */
	/* For loading exons from refFlat.txt */
	fprintf(stderr, "Loading genes ...");
	filename = "refFlat.txt";
	fp = fopen(filename, "r");
	if (fp == NULL){
		fprintf(stderr, "Could not open file %s",filename);
	}
	else {
		chrom_count = (int *) getMem( 26 * sizeof(int));
		chrom_count[0] = 61406;
		chrom_count[1] = CHR1;
		chrom_count[2] = CHR2;
		chrom_count[3] = CHR3;
		chrom_count[4] = CHR4;
		chrom_count[5] = CHR5;
		chrom_count[6] = CHR6;
		chrom_count[7] = CHR7;
		chrom_count[8] = CHR8;
		chrom_count[9] = CHR9;
		chrom_count[10] = CHR10;
		chrom_count[11] = CHR11;
		chrom_count[12] = CHR12;
		chrom_count[13] = CHR13;
		chrom_count[14] = CHR14;
		chrom_count[15] = CHR15;
		chrom_count[16] = CHR16;
		chrom_count[17] = CHR17;
		chrom_count[18] = CHR18;
		chrom_count[19] = CHR19;
		chrom_count[20] = CHR20;
		chrom_count[21] = CHR21;
		chrom_count[22] = CHR22;
		chrom_count[23] = CHRX;
		chrom_count[24] = CHRY;
		chrom_count[25] = CHRM;


		// in_exons = ( exon_info**) getMem( (EXONS) * sizeof( exon_info*));
		int * chrom_index = (int *) getMem( 26 * sizeof(int));
		for (int l = 0; l < 26; l++)
			chrom_index[l] = 0;

		exon_info ***in_exons = ( exon_info***) getMem( (CHROMS) * sizeof( exon_info**));
		for (int l = 0; l < CHROMS; l++) {
			in_exons[l] = ( exon_info**) getMem( (chrom_count[l]) * sizeof( exon_info*));
		}

		tokens = (char **) getMem((TOKENS) * sizeof(char *));
		this_exon_code = 0;

		while (fgets(str, MAXCHAR, fp) != NULL) {
			int k = 0;
			char * token = strtok(str, "	");

	   		// loop through the string to extract all other tokens
			while( token != NULL ) {
				tokens[k] = NULL;
				set_str( &(tokens[k]), token);

				if(k==10) {

					start = strtok_r(tokens[9], ",", &saveptr1);
					end = strtok_r(tokens[10], ",", &saveptr2);
	   				
					for(int j = 0; j < atoi(tokens[8]); j++)
					{
						int index;
					
						if(strcmp(tokens[2] + 3, "X") == 0)
							index = 23;
						else if(strcmp(tokens[2] + 3, "Y") == 0)
							index = 24;
						else if(strcmp(tokens[2] + 3, "M") == 0)
							index = 25;
						else if(strlen(tokens[2]) < 6)
							index = atoi(tokens[2] + 3);
						else
							index = 0;

						in_exons[index][chrom_index[index]] = ( exon_info*) getMem( sizeof( exon_info));

						in_exons[index][chrom_index[index]]->gene_id = NULL;
						set_str( &(in_exons[index][chrom_index[index]]->gene_id), tokens[0]);

						in_exons[index][chrom_index[index]]->transcript_id = NULL;
						set_str( &(in_exons[index][chrom_index[index]]->transcript_id), tokens[1]);

						in_exons[index][chrom_index[index]]->chr = NULL;
						set_str( &(in_exons[index][chrom_index[index]]->chr), tokens[2]);
						
						in_exons[index][chrom_index[index]]->start = atoi(start);
						
						in_exons[index][chrom_index[index]]->end = atoi(end);
						
						in_exons[index][chrom_index[index]]->exon_code = this_exon_code;
						
						if (strcmp(tokens[3], "+") == 0)
							in_exons[index][chrom_index[index]]->strand = 0;
						else
							in_exons[index][chrom_index[index]]->strand = 1;

						char exon_index[5];
						sprintf(exon_index, "%d", j);

						char * str3 = (char *) malloc(4 + strlen(tokens[0]) + strlen(tokens[1]) + strlen(exon_index) + strlen(tokens[2]));
						strcpy(str3, tokens[0]);
						strcat(str3, "_");
						strcat(str3, tokens[1]);
						strcat(str3, "_");
						strcat(str3, exon_index);
						strcat(str3, "_");
						strcat(str3, tokens[2]);

						in_exons[index][chrom_index[index]]->exon_id = NULL;
						set_str( &(in_exons[index][chrom_index[index]]->exon_id), str3);

						free(str3);

						if(j!= atoi(tokens[8])) {
							start = strtok_r(NULL, ",", &saveptr1);
							end = strtok_r(NULL, ",", &saveptr2);
						}
						this_exon_code++;
						chrom_index[index] = chrom_index[index] + 1;
					}
					break;
				}
				k++;
				token = strtok(NULL, "	");
			}
		}
		fclose(fp);
	}
	//qsort(in_exons, EXONS, sizeof(exon_info *), compare);
	for (int l = 1; l < CHROMS; l++) {
		qsort(in_exons[l], chrom_count[l], sizeof(exon_info *), compare);
	}

	// for (int i = 0; i < EXONS; i++)
	// 	in_exons[i]->exon_code = i;
	fprintf(stderr, "Loaded genes");

	/* Passing the flags to VH */
	ten_x_flag = params->ten_x;
	output_hs_flag = params->output_hs;

	print_quote();
	fprintf( stderr, "\n\tRun. Run, you clever boy... And remember.\n");

	if ( !params->no_soft_clip)
		init_hash_count( params);


	if ( !params->histogram_only)
	{
		/* Sensitive Mode */
		if( running_mode == SENSITIVE)
		{
			fprintf( stderr, "\nTARDIS (v%s) is running in Sensitive Mode - using mrFAST...\n\n", TARDIS_VERSION);
			fprintf( logFile, "(Running in sensitive mode - using mrFAST)\n\n");

			/* If you already have the correct divet file */
			if( params->skip_mrfast == 0)
			{
				/* Extract FASTQs of discordants, OEAs, and orphans */
				for( i = 0; i < params->num_bams; i++)
					create_fastq( in_bams[i], params->bam_file_list[i], params);

				fprintf( stderr, "All FASTQ files ready for remapping.\n");

				/* Remap with mrFAST */
				return_value = remap_mrfast( params, in_bams, cfg);
				if( return_value != RETURN_SUCCESS)
					return EXIT_EXTERNAL_PROG_ERROR;
			}
			/* Read depth calculation */
			return_value = run_rd( in_bams, params);
			if( return_value != RETURN_SUCCESS)
				return EXIT_EXTERNAL_PROG_ERROR;

			/* Read pair calculation */
			return_value = run_vh( params, in_bams);
			if( return_value != RETURN_SUCCESS)
				return EXIT_EXTERNAL_PROG_ERROR;

			clean_up_temp_files(params);
			free_sensitive( in_bams, params);
		}
		/* Quick Mode */
		else
		{
			fprintf( stderr, "\nTARDIS (v%s) is running in Quick Mode\n\n", TARDIS_VERSION);
			fprintf( logFile, "(Running in quick mode)\n\n");
			return_value = bamonly_run( params, in_bams, in_exons);
			if ( return_value != RETURN_SUCCESS)
				return EXIT_EXTERNAL_PROG_ERROR;

			clean_up_temp_files(params);
			free_quick( in_bams, params);
		}
	}

	username = ( char*) getMem( MAX_SEQ * sizeof( char));
	getlogin_r( username, (MAX_SEQ - 1));
	fprintf( stderr, "\n%s, before I go, I just want to tell you: you were fantastic. Absolutely fantastic. And you know what? So was I.\n", username);

	fclose( logFile);

	free( username);
	return EXIT_SUCCESS;
}
