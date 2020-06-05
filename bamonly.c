/*
 * bo_bamhandler.c
 *
 *  Created on: Aug 23, 2016
 *      Author: tardis
 */

#include "bamonly.h"
#include "common.h"
#include <math.h>
#include "variants.h"
#include "vh/vh_common.h"
#include "vh/vh_maximalCluster.h"
#include "vh/vh_setcover.h"
#include "vh_createMaxClusterMEI.h"
#include "vh_createMaxClusterNUMT.h"
#include "splitread.h"
#include "free.h"
#include "mappings.h"

#define EXONS 750993
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
#define CHROX 23334
#define CHROY 4854
#define CHRM 2
#define CHROMS 26

long del_cnt_bam = 0;
long ins_cnt_bam = 0;
long inv_cnt_bam = 0;
long mei_cnt_bam = 0;
long numt_cnt_bam = 0;
long tandup_cnt_bam = 0;
long sr_cnt_bam = 0;
long alt_cnt_bam = 0;

long cnt_concordant = 0, cnt_discordant = 0, cnt_unmapped = 0, cnt_total_reads;

long total_read_count = 0;


char **allReadNameList;
int * chrom_count;

int is_retro_cnv_deprecated(const char * chromosome_name, int start, int end, exon_info** in_exons) {
 	char str[3] = "chr";
	for (int i = 0; i < EXONS; i++)
	{
		int length = 3 + strlen(chromosome_name);
		char * chr = malloc(length * sizeof(char));
		for(int j = 0; j < 3; j++)
			chr[j] = str[j];
		for(int j = 3; j < length; j++)
			chr[j] = chromosome_name[j - 3];
		// printf("\n%s\n", chr);
		if(strcmp(in_exons[i]->chr, chr) == 0 ) 

		{
			if(start >= in_exons[i]->start && end <= in_exons[i]->end) {
				free(chr);
				return i;
			}
		}
		free(chr);
	}
	return -1;
}

int this_interval_intersects(int pos_start, int pos_end, int start, int end){
  /* all in */
  if (pos_start >= start && pos_end < end)
    return 1;

  /* all cover */
  else if (pos_start <= start && pos_end > end)
    return 1;

  /* left */
  else if (pos_start <= start && pos_end >= start)
    return 1;

  /* right */
  else if (pos_start <= end && pos_end > end)
    return 1;

  /* no hit.  */

  return 0;
}

int is_retro_cnv_deprecated2(const char * chromosome_name, int pos_start, int pos_end, exon_info** in_exons){

  int start;
  int end;
  int med;
  int interval_count;

  char str[3] = "chr";
  int length = 3 + strlen(chromosome_name);
  char * chr = malloc(length * sizeof(char));
  for(int j = 0; j < 3; j++)
	  chr[j] = str[j];
  for(int j = 3; j < length; j++)
	  chr[j] = chromosome_name[j - 3];
  // printf("\n%s\n", chr);
  
  start = 0;
  end = EXONS - 1;

  med = (start + end) / 2;

    
  while (1){

    if (start > end)
      return -1;

    if (strcmp(in_exons[med]->chr, chr) == 0 && this_interval_intersects(pos_start, pos_end, in_exons[med]->start, in_exons[med]->end))
      return in_exons[med]->exon_code;

    /* no hit. search is exhausted */
    if (start == med || end == med){
      if (strcmp(in_exons[start]->chr, chr) == 0 && this_interval_intersects(pos_start, pos_end, in_exons[start]->start, in_exons[start]->end))
		return in_exons[start]->exon_code;
      else if (strcmp(in_exons[end]->chr, chr) == 0 && this_interval_intersects(pos_start, pos_end, in_exons[end]->start, in_exons[end]->end))
		return in_exons[end]->exon_code;      
      return -1;
    }

    /* no hit, search left half */
    else if (pos_start < in_exons[med]->start){
      end = med;
      med = (start + end) / 2;
    }

    /* no hit, search right half */
    else {
      start = med;
      med = (start + end) / 2;      
    }
      
  }

  return -1;
}

exon_info * is_retro_cnv(const char * chromosome_name, int pos_start, int pos_end, exon_info*** in_exons){

	exon_info * potential_retro_cnv = NULL;
	int index = -1;
	if(strcmp(chromosome_name, "X") == 0)
		index = 23;
	else if(strcmp(chromosome_name, "Y") == 0)
		index = 24;
	else if(strcmp(chromosome_name, "MT") == 0)
		index = 25;
	else if((strlen(chromosome_name) < 3))
		index = atoi(chromosome_name);

	if(index == -1)
		return potential_retro_cnv;
	int start;
	int end;
	int med;
	int interval_count;

	// char str[3] = "chr";
	// int length = 3 + strlen(chromosome_name);
	// char * chr = malloc(length * sizeof(char));
	// for(int j = 0; j < 3; j++)
	// 	chr[j] = str[j];
	// for(int j = 3; j < length; j++)
	// 	chr[j] = chromosome_name[j - 3];
  // printf("\n%s\n", chr);

	start = 0;
	end = chrom_count[index] - 1;

	med = (start + end) / 2;


	while (1){

		if (start > end)
			return potential_retro_cnv;


		if (this_interval_intersects(pos_start, pos_end, in_exons[index][med]->start, in_exons[index][med]->end)) {
			potential_retro_cnv = ( exon_info*) getMem( sizeof( exon_info));
			potential_retro_cnv->gene_id = NULL;
			set_str( &(potential_retro_cnv->gene_id), in_exons[index][med]->gene_id);
			potential_retro_cnv->transcript_id = NULL;
			set_str( &(potential_retro_cnv->transcript_id), in_exons[index][med]->transcript_id);
			potential_retro_cnv->chr = NULL;
			set_str( &(potential_retro_cnv->chr), in_exons[index][med]->chr);
			potential_retro_cnv->start = in_exons[index][med]->start;
			potential_retro_cnv->end = in_exons[index][med]->end;
			potential_retro_cnv->exon_code = in_exons[index][med]->exon_code;	
			potential_retro_cnv->strand = in_exons[index][med]->strand;
			potential_retro_cnv->exon_id = NULL;
			set_str( &(potential_retro_cnv->exon_id), in_exons[index][med]->exon_id);
			return potential_retro_cnv;
		}

    /* no hit. search is exhausted */
		if (start == med || end == med){
			if (this_interval_intersects(pos_start, pos_end, in_exons[index][start]->start, in_exons[index][start]->end)) {
				potential_retro_cnv = ( exon_info*) getMem( sizeof( exon_info));
				potential_retro_cnv->gene_id = NULL;
				set_str( &(potential_retro_cnv->gene_id), in_exons[index][start]->gene_id);
				potential_retro_cnv->transcript_id = NULL;
				set_str( &(potential_retro_cnv->transcript_id), in_exons[index][start]->transcript_id);
				potential_retro_cnv->chr = NULL;
				set_str( &(potential_retro_cnv->chr), in_exons[index][start]->chr);
				potential_retro_cnv->start = in_exons[index][start]->start;
				potential_retro_cnv->end = in_exons[index][start]->end;
				potential_retro_cnv->exon_code = in_exons[index][start]->exon_code;	
				potential_retro_cnv->strand = in_exons[index][start]->strand;
				potential_retro_cnv->exon_id = NULL;
				set_str( &(potential_retro_cnv->exon_id), in_exons[index][start]->exon_id);
				return potential_retro_cnv;
			}
			else if (this_interval_intersects(pos_start, pos_end, in_exons[index][end]->start, in_exons[index][end]->end)) {
				potential_retro_cnv = ( exon_info*) getMem( sizeof( exon_info));
				potential_retro_cnv->gene_id = NULL;
				set_str( &(potential_retro_cnv->gene_id), in_exons[index][end]->gene_id);
				potential_retro_cnv->transcript_id = NULL;
				set_str( &(potential_retro_cnv->transcript_id), in_exons[index][end]->transcript_id);
				potential_retro_cnv->chr = NULL;
				set_str( &(potential_retro_cnv->chr), in_exons[index][end]->chr);
				potential_retro_cnv->start = in_exons[index][end]->start;
				potential_retro_cnv->end = in_exons[index][end]->end;
				potential_retro_cnv->exon_code = in_exons[index][end]->exon_code;	
				potential_retro_cnv->strand = in_exons[index][end]->strand;
				potential_retro_cnv->exon_id = NULL;
				set_str( &(potential_retro_cnv->exon_id), in_exons[index][end]->exon_id);
				return potential_retro_cnv;
			}
			return potential_retro_cnv;
		}

    /* no hit, search left half */
		else if (pos_start < in_exons[index][med]->start){
			end = med;
			med = (start + end) / 2;
		}

    /* no hit, search right half */
		else {
			start = med;
			med = (start + end) / 2;      
		}

	}

	return potential_retro_cnv;
}

void findUniqueReads( bam_info** in_bam, parameters *params, char *outputread)
{
	int i, j, k;
	int totalCountRead = 0;
	long read_name_count;

	softClip *softClipPtr;
	discordantMappingMEI *discordantReadPtrMEI;
	discordantMappingNUMT *discordantReadPtrNUMT;
	discordantMapping *discordantReadPtr;
	alternativeMapping *discordantReadPtrAlt;

	FILE *fileOutputReadName;

	if( debug_mode)
		fileOutputReadName = safe_fopen ( outputread, "w");

	allReadNameList = ( char **) getMem ( ( total_read_count + 1) * sizeof ( char *));
	if( allReadNameList == NULL)
		fprintf( stderr,"Memory problem\n");

	read_name_count = 0;

	/* Put the names of all the reads in allReadNameList and write them to .name file in sorted order */
	for( i = 0; i < params->num_bams; i++)
	{
		for( j = 0; j < in_bam[i]->num_libraries; j++)
		{
			for( k = 0; k < NHASH; k++)
			{
				discordantReadPtr = in_bam[i]->libraries[j]->mappings_discordant[k];
				while( discordantReadPtr != NULL)
				{
					allReadNameList[read_name_count] = NULL;
					set_str( &(allReadNameList[read_name_count]), discordantReadPtr->readName);
					read_name_count++;
					discordantReadPtr = discordantReadPtr->next;
				}
			}

			if( params->alt_mapping != 0)
			{
				for( k = 0; k < NHASH; k++)
				{
					discordantReadPtrAlt = in_bam[i]->libraries[j]->mappings_alternative[k];
					while( discordantReadPtrAlt != NULL)
					{
						allReadNameList[read_name_count] = NULL;
						set_str( &(allReadNameList[read_name_count]), discordantReadPtrAlt->readName);
						read_name_count++;
						discordantReadPtrAlt = discordantReadPtrAlt->next;
					}
				}
			}

			discordantReadPtrMEI = in_bam[i]->libraries[j]->listMEI_Mapping;
			while( discordantReadPtrMEI != NULL)
			{
				allReadNameList[read_name_count] = NULL;
				set_str( &(allReadNameList[read_name_count]), discordantReadPtrMEI->readName);
				read_name_count++;
				discordantReadPtrMEI = discordantReadPtrMEI->next;
			}

			discordantReadPtrNUMT = in_bam[i]->libraries[j]->listNUMT_Mapping;
			while( discordantReadPtrNUMT != NULL)
			{
				allReadNameList[read_name_count] = NULL;
				set_str( &(allReadNameList[read_name_count]), discordantReadPtrNUMT->readName);
				read_name_count++;
				discordantReadPtrNUMT = discordantReadPtrNUMT->next;
			}

			if( !params->no_soft_clip)
			{
				softClipPtr = in_bam[i]->libraries[j]->listSoftClip;
				while(softClipPtr!=NULL)
				{
					allReadNameList[read_name_count] = NULL;
					set_str( &(allReadNameList[read_name_count]), softClipPtr->readName);
					read_name_count++;
					softClipPtr = softClipPtr->next;
				}
			}
		}
	}
	qsort( allReadNameList, read_name_count, sizeof( char *), vh_cmprReadNameStr);

	for( i = 0; i < read_name_count; i++)
	{
		if ( i == 0 || strcmp( allReadNameList[i], allReadNameList[i-1]) != 0)
			totalCountRead++;
	}
	if( debug_mode)
		fprintf( fileOutputReadName, "%i\n", totalCountRead);

	/* Write the read names to read_names structure for set_cover */
	read_names = (readEl *) getMem( ( read_name_count + 1) * sizeof( readEl));
	j = 0;

	for( i = 0; i < read_name_count; i++)
	{
		if( i == 0 || strcmp( allReadNameList[i], allReadNameList[i-1]) != 0)
		{
			if( debug_mode)
				fprintf( fileOutputReadName, "%s\n", allReadNameList[i]);

			read_names[j].readName = NULL;
			set_str( &(read_names[j].readName), allReadNameList[i]);
			read_names[j].readCovered = 0;
			read_names[j].readId = j;
			read_names[j].libId = -1;
			read_names[j].indId = -1;
			read_names[j].next = NULL;
			j++;
		}
	}
	read_names_count = j;

	if( debug_mode)
		fclose( fileOutputReadName);

	for( i = 0; i < total_read_count; i++)
	{
		if( allReadNameList[i] != NULL)
			free( allReadNameList[i]);
	}
	if( allReadNameList != NULL)
		free( allReadNameList);

	total_read_count = 0;
}


void discordant_mapping( library_properties *library, parameters* params, bam_alignment_region* bam_align, int svType, int chrID)
{
	int len, h;

	discordantMapping *newEl;
	newEl = ( discordantMapping *) getMem( sizeof( discordantMapping));

	/* Need to be put into into divet row */
	if ( params->ten_x)
		newEl->ten_x_barcode = bam_align->ten_x_barcode;

	/* This is the first read */
	if( bam_align->pos_left < bam_align->pos_right)
	{
		newEl->readName = NULL;
		set_str( &(newEl->readName), bam_align->read_name);

		/* Get the name of the chromosome */
		len = strlen( params->this_sonic->chromosome_names[chrID]);
		newEl->chromosome_name = ( char *) getMem( sizeof ( char) * len + 1);
		strcpy( newEl->chromosome_name, params->this_sonic->chromosome_names[chrID]);

		newEl->pos1 = bam_align->pos_left;
		newEl->pos2 = bam_align->pos_right;
		newEl->mQual1 = bam_align->qual;
		newEl->mQual2 = 0;

		newEl->cigar = ( uint32_t*) getMem( 20 * sizeof( uint32_t));
		memcpy( newEl->cigar, bam_align->cigar, sizeof( bam_align->cigar));

		newEl->n_cigar = bam_align->n_cigar;

		if( ( svType == RPDEL || svType == RPINS ) &&
				( bam_align->flag & BAM_FREVERSE) == 0 && ( bam_align->flag & BAM_FMREVERSE) != 0)
		{
			newEl->orient1 = FORWARD;
			newEl->orient2 = REVERSE;
			if( svType == RPDEL)
			{
				newEl->svType = DELETION;
				del_cnt_bam++;
			}
			if( svType == RPINS)
			{
				newEl->svType = INSERTION;
				ins_cnt_bam++;
			}
		}
		else if( ( svType == RPTDUP) && ( bam_align->flag & BAM_FREVERSE) != 0 && ( bam_align->flag & BAM_FMREVERSE) == 0)
		{
			newEl->orient1 = REVERSE;
			newEl->orient2 = FORWARD;
			newEl->svType = TANDEMDUP;
			tandup_cnt_bam++;
		}
		else if( ( svType == RPTDUP) && ( bam_align->flag & BAM_FREVERSE) == 0 && ( bam_align->flag & BAM_FMREVERSE) != 0)
		{
			newEl->orient1 = FORWARD;
			newEl->orient2 = REVERSE;
			newEl->svType = TANDEMDUP;
			tandup_cnt_bam++;
			fprintf(stderr, "ERORRRRRR\n");
		}
		else if( svType == RPINV)
		{
			if( ( bam_align->flag & BAM_FREVERSE) != 0 && ( bam_align->flag & BAM_FMREVERSE) != 0)
			{
				newEl->orient1 = REVERSE;
				newEl->orient2 = REVERSE;
			}
			else if( ( bam_align->flag & BAM_FREVERSE) == 0 && ( bam_align->flag & BAM_FMREVERSE) == 0)
			{
				newEl->orient1 = FORWARD;
				newEl->orient2 = FORWARD;
			}
			inv_cnt_bam++;
			newEl->svType = INVERSION;
		}
		else
		{
			fprintf(stderr, "Error in add_discordant_FR() - flag:%016x pos1=%d, pos2=%d SV=%d\n", bam_align->flag, bam_align->pos_left, bam_align->pos_right, svType);
			return;
		}

		newEl->editDistance = bam_align->edit_distance;
		newEl->editDistance_left = bam_align->edit_distance;

		newEl->pos1_End = newEl->pos1 + library->read_length;

		if( bam_cigar_opchr( bam_align->cigar[bam_align->n_cigar - 1]) == 'S')
			newEl->pos1_End = newEl->pos1_End - bam_cigar_oplen( bam_align->cigar[bam_align->n_cigar - 1]);

		if( bam_cigar_opchr( bam_align->cigar[0]) == 'S')
			newEl->pos1 = newEl->pos1 + bam_cigar_oplen( bam_align->cigar[0]);

		newEl->pos1_End--;
		//newEl->pos1++;
		newEl->pos2_End = newEl->pos2 + library->read_length;

		newEl->flag = bam_align->flag;

		h = vh_getHash( newEl->readName);

		newEl->next = library->mappings_discordant[h];
		library->mappings_discordant[h] = newEl;

	}
	/* This is the second read */
	else if( bam_align->pos_left > bam_align->pos_right)
	{
		/* Check primary mappings */
		h = vh_getHash( bam_align->read_name);
		newEl = library->mappings_discordant[h];

		while( newEl != NULL)
		{
			if( strncmp( newEl->readName, bam_align->read_name, strlen( bam_align->read_name)) == 0)
			{
				newEl->editDistance = newEl->editDistance_left + bam_align->edit_distance;

				newEl->pos2_End = newEl->pos2 + library->read_length;
				if( bam_cigar_opchr( bam_align->cigar[bam_align->n_cigar - 1]) == 'S')
					newEl->pos2_End = newEl->pos2_End - bam_cigar_oplen( bam_align->cigar[bam_align->n_cigar - 1]);

				if( bam_cigar_opchr( bam_align->cigar[0]) == 'S')
					newEl->pos2 = newEl->pos2 + bam_cigar_oplen( bam_align->cigar[0]);

				newEl->pos2_End--;
				//newEl->pos2++;
				newEl->mQual2 = bam_align->qual;
			}
			newEl = newEl->next;
		}
	}
}

void discordant_mapping_MEI( library_properties *library, parameters* params,
		bam_alignment_region* bam_align, char* mei_subclass, char* mei_class, int MEI_Type, int chrID)
{
	int len;
	int countOp;

	discordantMappingMEI *newEl;
	newEl = ( discordantMappingMEI *) getMem( sizeof( discordantMappingMEI));

	/* Need to be put into into divet row */
	if ( params->ten_x)
		newEl->ten_x_barcode = bam_align->ten_x_barcode;

	newEl->readName = NULL;
	set_str( &(newEl->readName), bam_align->read_name);

	/* Get the name of the chromosome */
	len = strlen( params->this_sonic->chromosome_names[chrID]);
	newEl->chromosome_name = NULL;
	set_str( &(newEl->chromosome_name), params->this_sonic->chromosome_names[chrID]);

	/* Get the name of the mei subclass */
	len = strlen( mei_subclass);
	newEl->MEI_subclass = NULL;

	set_str( &(newEl->MEI_subclass), mei_subclass);

	/* Get the name of the mei class */
	len = strlen( mei_class);
	newEl->MEI_class = NULL;

	set_str( &(newEl->MEI_class), mei_class);

	newEl->pos = bam_align->pos_left;
	newEl->qual = bam_align->qual;

	if( ( bam_align->flag & BAM_FREVERSE) != 0)
		newEl->orient = REVERSE;
	else
		newEl->orient = FORWARD;

	mei_cnt_bam++;
	newEl->MEI_Type = MEI_Type;

	newEl->pos_End = newEl->pos;
	for( countOp = 0; countOp < bam_align->n_cigar; countOp++)
	{
		if( bam_cigar_opchr( bam_align->cigar[countOp]) == 'S')
			newEl->pos_End = newEl->pos_End + bam_cigar_oplen( bam_align->cigar[countOp]);
	}

	newEl->next = library->listMEI_Mapping;
	library->listMEI_Mapping = newEl;
}


void discordant_mapping_NUMT( library_properties *library, parameters* params,
		bam_alignment_region* bam_align, int numt_type, int chrID)
{
	int len;
	int countOp;

	discordantMappingNUMT *newEl;
	newEl = ( discordantMappingNUMT *) getMem( sizeof( discordantMappingNUMT));

	/* Need to be put into into divet row */
	if ( params->ten_x)
		newEl->ten_x_barcode = bam_align->ten_x_barcode;

	newEl->readName = NULL;
	set_str( &(newEl->readName), bam_align->read_name);

	/* Get the name of the chromosome */
	len = strlen( params->this_sonic->chromosome_names[chrID]);
	newEl->chromosome_name = NULL;
	set_str( &(newEl->chromosome_name), params->this_sonic->chromosome_names[chrID]);

	newEl->pos = bam_align->pos_left;
	newEl->qual = bam_align->qual;

	numt_cnt_bam++;
	newEl->NUMT_Type = numt_type;

	if( ( bam_align->flag & BAM_FREVERSE) != 0)
		newEl->orientation = REVERSE;
	else
		newEl->orientation = FORWARD;

	newEl->pos_End = newEl->pos + library->read_length;
	/*newEl->pos_End = newEl->pos;
	for( countOp = 0; countOp < bam_align->n_cigar; countOp++)
	{
		if( bam_cigar_opchr( bam_align->cigar[countOp]) == 'S')
			newEl->pos_End = newEl->pos_End + bam_cigar_oplen( bam_align->cigar[countOp]);
	}
	 */

	if( bam_cigar_opchr( bam_align->cigar[bam_align->n_cigar - 1]) == 'S')
		newEl->pos_End = newEl->pos_End - bam_cigar_oplen( bam_align->cigar[bam_align->n_cigar - 1]);

	if( bam_cigar_opchr( bam_align->cigar[0]) == 'S')
		newEl->pos = newEl->pos + bam_cigar_oplen( bam_align->cigar[0]);

	newEl->next = library->listNUMT_Mapping;
	library->listNUMT_Mapping = newEl;
}


int find_numt_bam( bam_alignment_region* bam_align, char *chromosome_name_left, char *chromosome_name_right)
{
	/* If the read maps to a chromosome and the right pair maps to MT */
	if( ( strstr( chromosome_name_left, "MT") == NULL || strstr( chromosome_name_left, "chrM") == NULL) &&
			( strstr( chromosome_name_right, "MT") != NULL || strstr( chromosome_name_right, "chrM") != NULL))
	{
		if( ( bam_align->flag & BAM_FMREVERSE) != 0)
			return NUMTR;
		else
			return NUMTF;
	}
	else
		return NOTNUMT;
}


int find_mei_bam( parameters *params, exon_info*** in_exons, char *chromosome_name, char** mei_subclass, char** mei_class, int start, int end, int flag)
{
	int ind, len;
	int return_type = NOTMEI;

	sonic_repeat *repeat_item;
	//int exon;
	exon_info * exon;

	/* Check if the right end is inside the annotated transposon */
	// fprintf( stderr, "\n Checking chromosome %s with start %d and end %d and length %d", chromosome_name, start, end, end - start);
	//repeat_item = sonic_is_mobile_element( params->this_sonic, chromosome_name, start, end, params->mei );
	exon = is_retro_cnv(chromosome_name, start, end, in_exons);
	// if( repeat_item == NULL) {
	// 	fprintf( stderr, "\n NOTMEI\n");
	// 	return NOTMEI;
	// }
	if( exon == NULL) {
		// fprintf( stderr, "\n NOTMEI\n");
		return NOTMEI;
	}
	// else {
	// 	fprintf( stderr, "\n strand: %d, repeat_type: %s, repeat_class: %s, repeat_start: %d, repeat_end: %d, mei_code: %d ", repeat_item->strand, repeat_item->repeat_type,
	// 		repeat_item->repeat_class, repeat_item->repeat_start, repeat_item->repeat_end,
	// 		repeat_item->mei_code);
	// }

	(*mei_subclass) = NULL;
	//set_str( mei_subclass, repeat_item->repeat_type);
	set_str( mei_subclass, exon->exon_id);

	(*mei_class) = NULL;
	//set_str( mei_class, repeat_item->repeat_class);
	set_str( mei_class, exon->exon_id);

	/* NOTE: SONIC keeps repeat_class as SINE/Alu, LINE/L1, etc. */
	// if( ( ( flag & BAM_FMREVERSE) == 0 && repeat_item->strand == SONIC_STRAND_REV)
	// 		|| ( ( flag & BAM_FMREVERSE) != 0 && repeat_item->strand == SONIC_STRAND_FWD))
		//return_type = repeat_item->mei_code * 2;
	if( ( ( flag & BAM_FMREVERSE) == 0 && exon->strand == SONIC_STRAND_REV)
			|| ( ( flag & BAM_FMREVERSE) != 0 && exon->strand == SONIC_STRAND_FWD))
		return_type = exon->exon_code;
	else
		//return_type = (repeat_item->mei_code * 2) + 1;
		return_type = exon->exon_code + EXONS;
	// fprintf( stderr, "exon_type: %d\n", return_type);
	return return_type;
}

int read_mapping( library_properties *library, exon_info*** in_exons, parameters* params, bam1_t* bam_alignment, int32_t *bamToRefIndex, bam_alignment_region* bam_align)
{
	int svType, meiType = NOTMEI, numtType = NOTNUMT, left_end_id, right_end_id, i, insLen;
	char* mei_subclass, *mei_class;

	/* Make sure the chromosome name is within the correct range */
	if( bam_align->chrID_left >= 0 && bam_align->chrID_left < params->this_sonic->number_of_chromosomes &&
			bam_align->chrID_right >= 0 && bam_align->chrID_right < params->this_sonic->number_of_chromosomes)
	{
		left_end_id = bamToRefIndex[bam_align->chrID_left];
		right_end_id = bamToRefIndex[bam_align->chrID_right];

		if( right_end_id >= params->this_sonic->number_of_chromosomes || right_end_id < 0 ||
				left_end_id >= params->this_sonic->number_of_chromosomes || left_end_id < 0)
			return -1;
	}
	else
		return -1;

	/* Since MT is circular, we need to eliminate the read-pairs at both ends of the chromosome */
	insLen = abs( bam_align->pos_left - bam_align->pos_right);

	/* Find the SVs */
	if( sonic_is_satellite( params->this_sonic, params->this_sonic->chromosome_names[left_end_id], bam_align->pos_left, bam_align->pos_left + library->read_length) == 0
			&& sonic_is_satellite( params->this_sonic, params->this_sonic->chromosome_names[right_end_id], bam_align->pos_right, bam_align->pos_right + library->read_length) == 0
			&& insLen < (params->this_sonic->chromosome_lengths[left_end_id]) - (2 * library->conc_max) && bam_align->qual > params->mq_threshold
			&& is_proper( bam_align->flag) && !( ( left_end_id == right_end_id) && ( insLen < MIN_INS_LEN)))
		/* Remove the paired-ends that both end overlap each other */
	{
		svType = is_concordant_bamonly( bam_align->pos_left, bam_align->pos_right, bam_align->flag,
				bam_align->isize, library->conc_min, library->conc_max);

		/* For statistics */
		cnt_total_reads++;

		if( svType == RPCONC)
			cnt_concordant++;
		else if( svType == RPUNMAPPED)
			cnt_unmapped++;
		else
			cnt_discordant++;

		if( svType != RPCONC && svType != RPUNMAPPED)
		{
			meiType = find_mei_bam( params, in_exons, params->this_sonic->chromosome_names[right_end_id], &mei_subclass, &mei_class, bam_align->pos_right,
					bam_align->pos_right + library->read_length, bam_align->flag);

			numtType = find_numt_bam( bam_align, params->this_sonic->chromosome_names[left_end_id], params->this_sonic->chromosome_names[right_end_id]);

			/* NUMT */
			if( numtType != NOTNUMT)
				discordant_mapping_NUMT( library, params, bam_align, numtType, left_end_id);
			/* MEI */
			else if( meiType != NOTMEI && ( left_end_id != right_end_id ||
					abs( bam_align->pos_left - bam_align->pos_right) > MIN_MEI_DISTANCE))
				discordant_mapping_MEI( library, params, bam_align, mei_subclass, mei_class, meiType, left_end_id);
			/* Deletion, Insertion or Tandem Duplication */
			else if( svType != RPINV && left_end_id == right_end_id)
				discordant_mapping( library, params, bam_align, svType, left_end_id);
			/* Inversion */
			else if( left_end_id == right_end_id && abs( bam_align->pos_left - bam_align->pos_right) < MAX_INV_LEN)
				discordant_mapping( library, params, bam_align, svType, left_end_id);
		}
		/* Soft Clipping */
		if( !params->no_soft_clip && bam_align->n_cigar < MAX_NUM_CIGAR && params->alt_mapping == 0)
		{
			/* We need to have soft clip at the beginning, end, or both with length >MIN_SOFTCLIP_LEN */
			if( ( bam_cigar_opchr( bam_align->cigar[0]) == 'S' && bam_cigar_oplen( bam_align->cigar[0]) > MIN_SOFTCLIP_LEN) ||
					( bam_cigar_opchr( bam_align->cigar[bam_align->n_cigar - 1]) == 'S' &&
							bam_cigar_oplen( bam_align->cigar[bam_align->n_cigar - 1]) > MIN_SOFTCLIP_LEN))
				addSoftClip( library, bam_align, bam_alignment, params->this_sonic->chromosome_names[left_end_id]);
		}

		return RETURN_SUCCESS;
	}
	return -1;
}

int read_bam( bam_info* in_bam, exon_info*** in_exons, parameters* params)
{
	/* Variables */
	int i, chr_index_bam, return_type, ed, len, lib_index;
	char *xa_str = NULL, *library_name = NULL;
	int32_t bamToRefIndex[in_bam->bam_header->n_targets];
	char fai_file[MAX_SEQ];

	sprintf( fai_file,"%s.fai", params->ref_genome);
	hts_set_fai_filename( in_bam->bam_file, fai_file);
	//hts_set_fai_filename(in_bam->bam_file, "/home/tardis/tardis/human_g1k_v37.fasta.fai");

	bam1_core_t bam_alignment_core;
	bam1_t* bam_alignment = bam_init1();

	/* The array is used to map the chromosome indices in bam file to the ones in reference genome */
	for( i = 0; i < in_bam->bam_header->n_targets; i++)
		bamToRefIndex[i] = sonic_refind_chromosome_index( params->this_sonic, in_bam->bam_header->target_name[i]);


	while( sam_itr_next( in_bam->bam_file, in_bam->iter, bam_alignment) > 0)
	{
		bam_alignment_core = bam_alignment->core;

		if( in_bam->num_libraries > 1)
		{
			/* Get library index */
			set_str( &library_name, bam_aux_get( bam_alignment, "RG"));
			lib_index = find_library_index( in_bam, library_name + 1);
		}
		else
			lib_index = 0;

		return_type = primary_mapping( in_bam, in_exons, params, lib_index, bam_alignment, bamToRefIndex);
		if( return_type == -1)
			continue;

		/* Increase the read depth and read count for RD filtering */
		in_bam->read_depth_per_chr[bam_alignment_core.pos]++;
		in_bam->read_count++;

		/* Get alternative mapping locations - XA field */
		if( params->alt_mapping != 0)
		{
			xa_str = NULL;
			set_str( &xa_str, bam_aux_get( bam_alignment, "XA"));

			if( xa_str != NULL)
				find_alt_mappings( in_bam, params, lib_index, bam_alignment, ( xa_str + 1));
		}
	}
	bam_destroy1( bam_alignment);
	fprintf( stderr, "\n[%s] %li DEL, %li INV, %li INS, %li TANDUP, %li MEI, %li NUMT, %li Split Read and %li XA regions found in BAM/CRAM.\n",
			in_bam->sample_name, del_cnt_bam, inv_cnt_bam, ins_cnt_bam, tandup_cnt_bam, mei_cnt_bam, numt_cnt_bam, sr_cnt_bam, alt_cnt_bam);
	fprintf( logFile, "\n[%s] %li DEL, %li INV, %li INS, %li TANDUP, %li MEI, %li Split Read and %li XA regions found in BAM/CRAM.\n",
			in_bam->sample_name, del_cnt_bam, inv_cnt_bam, ins_cnt_bam, tandup_cnt_bam, mei_cnt_bam, sr_cnt_bam, alt_cnt_bam);

	if (del_cnt_bam + inv_cnt_bam + ins_cnt_bam + tandup_cnt_bam + mei_cnt_bam + numt_cnt_bam + sr_cnt_bam + alt_cnt_bam == 0)
		return 1;
	return 0;
}


void bamonly_vh_clustering( bam_info** in_bams, exon_info*** in_exons, parameters *params)
{
	int i, bam_index, chr_index, chr_index_bam, return_value, not_in_bam = 0, invdup_location, interdup_location;
	int total_sv = 0, total_sv_lowqual = 0, divet_row_count;
	char outputfile[MAX_SEQ];
	char outputread[MAX_SEQ];
	char svfile[MAX_SEQ];
	FILE *fpVcf = NULL;
	int skip_chromosome = 0;
	int any_in_bam = 0;

	sprintf( outputread, "%s%s_name.log", params->outdir, params->outprefix);
	sprintf( outputfile, "%s%s_clusters.log", params->outdir, params->outprefix);

	/* Print all structural variations in .vcf format */
	sprintf( svfile, "%s%s.vcf", params->outdir, params->outprefix);
	fprintf( stderr, "Output file: %s\n", svfile);
	fpVcf = safe_fopen( svfile,"w");

	print_vcf_header( fpVcf, in_bams, params);

	/* Open the .clusters file if running in debug mode */
	if( debug_mode)
		fileOutput = safe_fopen( outputfile, "w");

	for( chr_index = params->first_chr; chr_index <= params->last_chr; chr_index++)
	{
		if( strstr( params->this_sonic->chromosome_names[chr_index], "GL000220") != NULL)
			continue;


		for( bam_index = 0; bam_index < params->num_bams; bam_index++)
		{
			/* HTS implementation */
			in_bams[bam_index]->bam_file = safe_hts_open( params->bam_file_list[bam_index], "r");

			/* Read in BAM header information */
			in_bams[bam_index]->bam_header = sam_hdr_read(in_bams[bam_index]->bam_file);

			/* Load the bam index file */
			in_bams[bam_index]->bam_file_index = sam_index_load( in_bams[bam_index]->bam_file, params->bam_file_list[bam_index]);
			if( in_bams[bam_index]->bam_file_index == NULL)
			{
				fprintf( stderr, "\nError: BAM/CRAM Index cannot be loaded (sam_index_load): %s\n", params->bam_file_list[bam_index]);
				exit( EXIT_BAM_INDEX);
			}

			chr_index_bam = find_chr_index_bam( params->this_sonic->chromosome_names[chr_index], in_bams[bam_index]->bam_header);
			not_in_bam = 0;
			if( chr_index_bam == -1)
			{
				fprintf( stderr, "\nCannot find chromosome name %s in BAM/CRAM %s", params->this_sonic->chromosome_names[chr_index], in_bams[bam_index]->sample_name);
				not_in_bam = 1;
				continue;
			}
			else
				any_in_bam = 1;

			in_bams[bam_index]->iter = sam_itr_queryi( in_bams[bam_index]->bam_file_index, chr_index_bam, 0, params->this_sonic->chromosome_lengths[chr_index]);
			if( in_bams[bam_index]->iter == NULL)
			{
				fprintf( stderr, "Error: Iterator cannot be loaded (bam_itr_queryi)\n");
				exit( EXIT_ITERATOR_LOAD);
			}

			fprintf( stderr, "\n                                                        ");
			fflush( stderr);
			fprintf( stderr, "\nReading BAM [%s] - Chromosome: %s", in_bams[bam_index]->sample_name, in_bams[bam_index]->bam_header->target_name[chr_index_bam]);
			fflush( stderr);

			/*
			if ( chr_index != params->first_chr)
			free_rd_per_chr(in_bams[bam_index], params, chr_index); */
			/* Initialize the read depth and read count */
			init_rd_per_chr( in_bams[bam_index], params, chr_index);

			/* Read bam file for this chromosome */
			skip_chromosome = read_bam( in_bams[bam_index], in_exons, params);


			if ( skip_chromosome)
			{
				/* Close the BAM file */
				return_value = hts_close( in_bams[bam_index]->bam_file);
				if( return_value != 0)
				{
					fprintf( stderr, "Error closing BAM file\n");
					exit( EXIT_BAM_CLOSE);
				}
				/* Free the bam related files */
				sam_itr_destroy( in_bams[bam_index]->iter);
				bam_hdr_destroy( in_bams[bam_index]->bam_header);
				hts_idx_destroy( in_bams[bam_index]->bam_file_index);
				continue;
			}

			//fprintf( stderr, "%d altLeftPrimRight, %d primRightAltLeft, %d altLeftAltRight, %d primLeftAltRight\n", altLeftPrimRight, primRightAltLeft, altLeftAltRight, primLeftAltRight);

			if( !params->no_soft_clip)
			{
				/* Count the number of softclip reads which are clustering for each read */
				fprintf( stderr, "\nReading reference genome");
				readReferenceSeq( params, chr_index);

				fprintf( stderr, "\nRunning split read mapping..");
				countNumSoftClipInCluster( params, in_bams[bam_index], chr_index);
				fprintf( stderr, "..");
				mapSoftClipToRef( in_bams[bam_index], params, chr_index);
				fprintf( stderr, "..");
			}

			/* Mean value (mu) calculation */
			calc_mean_per_chr( params, in_bams[bam_index], chr_index);

			/* Close the BAM file */
			return_value = hts_close( in_bams[bam_index]->bam_file);
			if( return_value != 0)
			{
				fprintf( stderr, "Error closing BAM file\n");
				exit( EXIT_BAM_CLOSE);
			}
			/* Free the bam related files */
			sam_itr_destroy( in_bams[bam_index]->iter);
			bam_hdr_destroy( in_bams[bam_index]->bam_header);
			hts_idx_destroy( in_bams[bam_index]->bam_file_index);

			/* Increase the total read count for this chromosome and make the number of SVs 0 for this bam */
			total_read_count += del_cnt_bam + inv_cnt_bam + ins_cnt_bam + tandup_cnt_bam + mei_cnt_bam + numt_cnt_bam + sr_cnt_bam + alt_cnt_bam;
			del_cnt_bam = 0, ins_cnt_bam = 0, inv_cnt_bam = 0, mei_cnt_bam = 0, numt_cnt_bam = 0, tandup_cnt_bam = 0, sr_cnt_bam = 0, alt_cnt_bam = 0;
		}

		if( not_in_bam == 1 || total_read_count == 0)
			continue;

		divet_row_count = load_Divet_bam( in_bams, params, chr_index);


		/* Make all clusters NULL */
		for( i = 0; i < MaxClusterCount; i++)
			clusters_all[i] = NULL;

		// /* Deletion */
		// fprintf( stderr, "\nPreparing Deletion clusters");
		// vh_initializeReadMapping_Deletion( params->this_sonic, chr_index);
		// fprintf( stderr, "..");
		// fflush( stderr);
		// vh_createDeletionClusters( params->this_sonic->chromosome_lengths[chr_index]);
		// fprintf( stderr, "..");
		// fflush( stderr);
		// vh_finalizeReadMapping( params->this_sonic->chromosome_names[chr_index], params->this_sonic->chromosome_lengths[chr_index]);
		// fprintf( stderr, "..");
		// fflush( stderr);

		// /* Inversion */
		// fprintf( stderr, "\nPreparing Inversion clusters");
		// vh_initializeReadMapping_Inversion( params->this_sonic, chr_index);
		// fprintf( stderr, "..");
		// fflush( stderr);
		// vh_createInversionClusters( params->this_sonic->chromosome_lengths[chr_index]);
		// fprintf( stderr, "..");
		// fflush( stderr);
		// vh_finalizeReadMapping( params->this_sonic->chromosome_names[chr_index], params->this_sonic->chromosome_lengths[chr_index]);
		// fprintf( stderr, "..");
		// fflush( stderr);

		// /* Insertion */
		// fprintf( stderr, "\nPreparing Insertion clusters");
		// vh_initializeReadMapping_Insertion( params->this_sonic, chr_index);
		// fprintf( stderr, "..");
		// fflush( stderr);
		// vh_createInsertionClusters( params->this_sonic->chromosome_lengths[chr_index]);
		// fprintf( stderr, "..");
		// fflush( stderr);
		// vh_finalizeReadMapping( params->this_sonic->chromosome_names[chr_index], params->this_sonic->chromosome_lengths[chr_index]);
		// fprintf( stderr, "..");
		// fflush( stderr);

		// /* Tandem Duplication */
		// fprintf( stderr, "\nPreparing Tandem Duplication clusters");
		// vh_initializeReadMapping_TDup( params->this_sonic, chr_index);
		// fprintf( stderr, "..");
		// fflush( stderr);
		// vh_createTDupClusters( params->this_sonic->chromosome_lengths[chr_index]);
		// fprintf( stderr, "..");
		// fflush( stderr);
		// vh_finalizeReadMapping( params->this_sonic->chromosome_names[chr_index], params->this_sonic->chromosome_lengths[chr_index]);
		// fprintf( stderr, "..");
		// fflush( stderr);

		/* Mei */
		
			fprintf( stderr, "\nPreparing MEI clusterss");
			int a = initializeReadMapping_MEI( in_bams, params, chr_index);
			if(a == 0)
				fprintf( stderr, "\nHELLO QUICK\n");
			else
				fprintf( stderr, "\nHELLO SENSITIVE\n");
			fprintf( stderr, "..");
			//fflush( stderr);
			MEICluster_Region( params, chr_index);
			fprintf( stderr, "..");
			//fflush( stderr);
			vh_finalizeReadMapping_Mei( params->this_sonic->chromosome_lengths[chr_index]);
			fprintf( stderr, "..");
			//fflush( stderr);

		// /* NumT */
		// if( ( strcmp( params->this_sonic->chromosome_names[chr_index], "MT") != 0)
		// 		&& ( strcmp( params->this_sonic->chromosome_names[chr_index], "chrM") != 0))
		// {
		// 	fprintf( stderr, "\nPreparing NUMT clusters");
		// 	initializeReadMapping_NUMT( in_bams, params, chr_index);
		// 	fprintf( stderr, "..");
		// 	fflush( stderr);
		// 	NUMTCluster_Region( params, chr_index);
		// 	fprintf( stderr, "..");
		// 	fflush( stderr);
		// 	vh_finalizeReadMapping_NUMT( params->this_sonic->chromosome_lengths[chr_index]);
		// 	fprintf( stderr, "..");
		// 	fflush( stderr);
		// }

		// if( params->no_interdup == 0)
		// {
		// 	/* Interspersed Direct Duplication */
		// 	fprintf( stderr, "\nPreparing Interspersed Duplication clusters");
		// 	for( interdup_location = 0; interdup_location <= RIGHTSIDE; interdup_location++)
		// 	{
		// 		vh_initializeReadMapping_InterDup( params->this_sonic, chr_index, interdup_location);
		// 		fprintf( stderr, ".");
		// 		fflush( stderr);
		// 		vh_createInterDupClusters( params->this_sonic->chromosome_lengths[chr_index], interdup_location);
		// 		fprintf( stderr, ".");
		// 		fflush( stderr);
		// 		vh_finalizeReadMapping_InterDup( params->this_sonic->chromosome_names[chr_index], params->this_sonic->chromosome_lengths[chr_index]);
		// 		fprintf( stderr, ".");
		// 		fflush( stderr);
		// 	}

		// 	/* Interspersed Inverted Duplication */
		// 	fprintf( stderr, "\nPreparing Interspersed Duplication (Inverted) clusters");
		// 	for( invdup_location = 0; invdup_location <= RIGHTSIDE; invdup_location++)
		// 	{
		// 		vh_initializeReadMapping_InvDup( params->this_sonic, chr_index, invdup_location);
		// 		fprintf( stderr, ".");
		// 		fflush( stderr);
		// 		vh_createInvDupClusters( params->this_sonic->chromosome_lengths[chr_index], invdup_location);
		// 		fprintf( stderr, ".");
		// 		fflush( stderr);
		// 		vh_finalizeReadMapping_InvDup( params->this_sonic->chromosome_names[chr_index], params->this_sonic->chromosome_lengths[chr_index]);
		// 		fprintf( stderr, ".");
		// 		fflush( stderr);
		// 	}
		// }

		fprintf( stderr, "\n");
		findUniqueReads( in_bams, params, outputread);
		

		/* Free the mappings and libraries */
		free_mappings( in_bams, params);
		free_libraries();
		
		/* Apply Set-Cover */
		vh_setcover( in_bams, params, fpVcf);
		
		total_sv += sv_count;
		total_sv_lowqual += sv_lowqual_count;

		if( !params->no_soft_clip)
		{
			/* Free the hash */
			// free_HashIndex();
			free_hash_table( params);
		}

		free_the_rest( in_bams, params);
	}


	fprintf( stderr, "\n");
	fclose( fpVcf);

	if ( !any_in_bam)
	{
		fprintf( stderr, "No chromosome is found in the input BAM/CRAM. Your SONIC and/or reference FASTA files do not match the reference used to generate the alignment file(s).\n");
		fprintf( stderr, "Check the reference genome versions and make sure you use the correct SONIC and reference FASTA files.\n");
		exit(EXIT_WRONG_SONIC);
	}


	if( debug_mode)
		fclose( fileOutput);

	fflush( stderr);
	fprintf( stderr, "\n");

	fprintf( stderr, "TARDIS is complete. Found %d SVs", total_sv);

	fprintf( logFile, "\n\nTotal read count = %li; %li Concordant, %li Unmapped, %li Discordant\n",
			cnt_total_reads, cnt_concordant, cnt_unmapped, cnt_discordant);
	fprintf( stderr, "\n\nTotal read count = %li; %li Concordant, %li Unmapped, %li Discordant\n",
			cnt_total_reads, cnt_concordant, cnt_unmapped, cnt_discordant);

	print_sv_stats();
}

int bamonly_run( parameters *params, bam_info ** in_bams, exon_info *** in_exons)
{
	int rd_del_filtered, bam_index;
	int sv_total, i, len;

	/* Initialize and read bam file */
	fprintf( stderr, "Processing BAM file for read pair and read depth filtering\n"
			"(Mapping Quality Threshold: %d; RP Support Threshold: %d)\n\n"
			, params->mq_threshold, params->rp_threshold);

	fprintf( logFile,"\n--> Processing BAM file for read pair and read depth filtering\n"
			"(Mapping Quality Threshold: %d; RP Support Threshold: %d)\n\n"
			, params->mq_threshold, params->rp_threshold);

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
	chrom_count[23] = CHROX;
	chrom_count[24] = CHROY;
	chrom_count[25] = CHRM;

	bamonly_vh_clustering( in_bams, in_exons, params);

	return RETURN_SUCCESS;
}
