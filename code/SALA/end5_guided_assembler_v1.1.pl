#!/usr/bin/env perl

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
no warnings 'utf8';
use warnings;
use strict;
use File::Path;
use File::Copy;
use File::Basename;
use File::Spec::Functions qw(rel2abs abs2rel);
use Time::HiRes qw( time );
use threads ('stack_size' => 64*4096);
use threads::shared;
use Getopt::Long 'HelpMessage';
use List::Util qw (sum shuffle min max);
use Cwd 'abs_path';
use AutoLoader qw/AUTOLOAD/;
#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
=head1 SYNOPSIS

 Description:
   This tool assigns ONT CAGE reads (as query) to a set of reference transcripts (e.g. GENCODE) using a 5'end centric approach. It will take a set of 
   user-defined of confident 5'end clusters (or de novo defined by clustering) and 3'end clusters (or de novo defined by clustering) and assign 
   the query reads to the reference transcripts in a 5'end centric approach with the following step:

   1. A query read is classified as complete if both of its 5' and 3' end overlap a confident 5' and 3' end cluster, otherwise as incomplete;
   2. An incomplete read without a confident 5' cluster will be flagged and excluded from downstream procedures;
   3. A complete query read will be assigned to a reference transcript if it shares the same 1) 5'end cluster, 2) 3'end cluster 
      and 3) internal splicing structure (i.e. same splicing juctions or both unspliced);
   4. An incomplete query read will be assigned to a reference transcript if it shares the same 5'end cluster and a partial  
      internal splicing structure (i.e. contains part if the reference transcript plicing juctions or unspliced 
      but overlap with reference transcript 1st exon);
   5. All unassigned query reads will be flagged as novel;
   6. Novel complete query reads with the same 1) 5'end cluster, 2) 3'end cluster and 3) internal splicing structure (i.e. same splicing 
      juctions or both unspliced) will be collapsed as a novel transcript model will new ID assigned;
   7. Novel incomplete query reads will be assigned to the novel transcript models if it shares the same 5'end cluster and a partial  
      internal splicing structure (i.e. contains part of the reference transcript splicing juctions or unspliced but overlap with 
      novel transcript models 1st exon);
   8. All remaining unassigned novel incomplete query reads will be grouped by their 1) 5'end clusters and 2) internal splicing structure 
      (i.e. same splicing juctions or unspliced) and each group will be collapsed as a novel transcript model with a new ID.
   9. The 5'end all transcript models will be adjusted to the summit of the 5'end clusters (de novo or user defined).
  10. The 3'end complete transcript models will be adjusted to the summit of the 3'end clusters (de novo or user defined).
  11. The 3'end incomplete transcript models will be adjusted to the furtherest 3'end of its query reads.

 Usage:
   end5_guided_assembler_v1.1.pl [options] --qry_bed_bgz --ref_bed_bgz --out_dir
   
   --qry_bed_bgz                <required> [path]    bed 12 of the ONT CAGE reads, 4th column must be read ID and in bgz format, 
                                                     for multiple query bed, user can supply a list of path in plain text format, one line one path
   --ref_bed_bgz                <required> [path]    bed 12 of the reference transcript models, 4th column must be transcript ID and in bgz format
   --out_dir                    <required> [path]    output directory
   --chrom_size_path            <required> [path]    a txt file contains the chromsome size in format of chrom\tsize
   --chrom_fasta_path           <required> [path]    genome fasta file
   --conf_end5_bed_bgz          <required> [path]    a bed bgz 12 file contains the 5'end clusters, summit must be provide in the thick end column
   --conf_end3_bed_bgz          <required> [path]    a bed bgz 12 file contains the 3'end clusters, summit must be provide in the thick end column
   --signal_end5_bed_bgz        <required> [path]    a single nucleotide piled up end5 signal bed (ctss bed file) used to define conf_end5_bed_bgz
   --signal_end3_bed_bgz        <required> [path]    a single nucleotide piled up end3 signal bed (ctes bed file) used to define conf_end3_bed_bgz
   --out_prefix                 (optional) [string]  output files prefix, if not defined, qry_bed_bgz filename will be used
   --novel_model_prefix         (optional) [string]  prefix of the novel transcript models [default=ONTC]
   --min_qry_score              (optional) [integer] the minimum score in the query bed file (assumes MAPQ) to be taken for assembly [default=10]
   --conf_end3_merge_flank      (optional) [integer] the flanking distance (on each side) of the 3'end clusters used to merge as a end3 region.
                                                     Use '-1' to turn off. [default=50]
   --conf_end5_merge_flank      (optional) [integer] the flanking distance (on each side) of the 5'end clusters used to merge as a end5 region.
                                                     Use '-1' to turn off. [default=50]
   --conf_end3_add_ref          (optional) [yes/no]  to add reference 3'end into the user defined confident 3'end clusters or not. if yes, the ref 3'end 
                                                     will bed extended by conf_end3_merge_flank nt and merged with confident 3'end clusters
   --min_exon_length            (optional) [integer] minimum length of an exon in a transcript to be considered as valid. If a transcript contains
                                                     an exon shorter than min_exon_length, the transcript will be discarded [default=1]
   --min_transcript_length      (optional) [integer] minimum length of a transcript (including intron) to be considered as valid. If a transcript 
                                                     is shorter than min_transcript_length, the transcript will be discarded [default=50]
   --filter_conf_end5           (optional) [yes/no]  to filter out query reads that is out the original ranges in the conf_end5_bed_bgz
                                                     will bed extended by conf_end3_merge_flank nt and merged with confident 3'end clusters
   --trnscpt_set_end_priority   (optional) [string]  Priority of methods to determine the ends of transcript set? 
                                                     1) based on "summit" : the signal summit in confident end3/end5 clusters, in signal_end*_bed_bgz
                                                     2) based on "commonest" : the observed position that is the most frequent in transcripts of the set
                                                     3) based on "longest": the observed position that is the furtherest in transcripts of the set
                                                     for (1) and (2), there is chances of causing conflicts in the transcript set ranges (e.g. 3'end is
                                                     more downstream than the 5'end in the transcript set). (3) is guranteed to be conflict free.
                                                     use a colon (:) delimited string to indiciate priority e.g. "summit:commonest:longest"
                                                     [default=summit:commonest:longest]
   --enforce_qry_original_end   (optional) [yes/no]  Overrides 
   --print_trnscrptID           (optional) [yes/no]  Print out the transcript ID or not 
   --doubtful_end_merge_dist    (optional) [integer] Distance to merge incomplete ends as groups [default=100]
   --doubtful_end_avoid_summit  (optional) [yes/no]  Overrides --trnscpt_set_end_priority from using "summit" 
   --retain_no_qry_ref_bound_set(optional) [yes/no]  report the bound set or not if the bound set is not detected from the query reads
   --min_summit_dist_split      (optional) [integer] When splitting an end cluster into two, the minimum distance between two summits [default=50]
   --min_size_split             (optional) [integer] When splitting an end cluster into two, the minimum size of the cluster [default=100]
   --min_frac_split             (optional) [integer] When splitting an end cluster into two, the minimum fraction of signal from the two summits [default=0.2]
   --max_thread                 (optional) [integer] number of threads to be used [default=5]
   --bedtools_bin               (optional) [path]    path to the binary of bedtools, if not provided, "bedtools" will be called
   --tabix_bin                  (optional) [path]    path to the binary of tabix, if not provided, "tabix" will be called
   --bgzip_bin                  (optional) [path]    path to the binary of bgzip, if not provided, "bgzip" will be called

 Dependencies:
   perl
   bedtools
   tabix
   bgzip
   

=head1 VERSION

1.1   -debut
=cut
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	notCalledBefore
#
#	notCalledBefore
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
my $scriptDirPath = dirname(rel2abs($0));
my $scriptAbsPath = abs_path($0);
my ($curntTimeStamp) = &timeStamp();#->3732
my $ARGVStr = join "\n", (&currentTime(), $scriptAbsPath, @ARGV);#->871
my $globalReadmeHsh_ref = {};
our $tmplog_fh;
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#
#<section ID="startingTasks" num="0">
my ($qry_bed_bgz, $ref_bed_bgz, $chrom_size_path, $chrom_fasta_path, $novel_model_prefix, $conf_end3_merge_flank, $conf_end5_merge_flank, $doubtful_end_merge_dist, $doubtful_end_avoid_summit, $use_ref_only_end_pos, $conf_end5_bed_bgz, $conf_end3_bed_bgz, $conf_junction_bed, $signal_end5_bed_bgz, $signal_end3_bed_bgz, $conf_end3_add_ref, $conf_end5_add_ref, $min_size_split, $min_frac_split, $min_summit_dist_split, $min_qry_score, $retain_no_qry_ref_bound_set, $trnscpt_set_end_priority_hsh_ref, $min_transcript_length, $min_exon_length, $print_trnscrptID, $bedtools_bin, $tabix_bin, $bgzip_bin, $max_thread, $min_output_qry_count, $out_prefix, $out_dir) = &readParameters();#->3438
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineHardCodedParam
#
#<section ID="defineHardCodedParam" num="1">
my $paramTag = "$out_prefix";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_defineout_dirPath
#
#<section ID="defineout_dirPath" num="2">
my @mkDirAry;
my $result_dir = "$out_dir/$paramTag"; push @mkDirAry, $result_dir;
my $result_tmp_dir = "$result_dir/tmp/"; push @mkDirAry, $result_tmp_dir;
my $result_bed_dir = "$result_dir/bed/"; push @mkDirAry, $result_bed_dir;
my $result_log_dir = "$result_dir/log/"; push @mkDirAry, $result_log_dir;
my $result_script_dir = "$result_dir/script/"; push @mkDirAry, $result_script_dir;
foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}

open $tmplog_fh, ">", "$result_dir/00_screen_log.$curntTimeStamp.log.txt";
&logCalledCMDAndScript($ARGVStr, $result_script_dir, $scriptAbsPath);#->1907
&printStartOrFinishMessage("startMessage");#->3115
&reportAndLogStatus("Parameters:", 10, "\n");#->3614
&reportAndLogStatus("novel_model_prefix = $novel_model_prefix", 10, "\n");#->3614
&reportAndLogStatus("conf_end3_merge_flank = $conf_end3_merge_flank", 10, "\n");#->3614
&reportAndLogStatus("conf_end5_merge_flank = $conf_end5_merge_flank", 10, "\n");#->3614
&reportAndLogStatus("conf_end3_add_ref = $conf_end3_add_ref", 10, "\n");#->3614
&reportAndLogStatus("conf_end3_add_ref = $conf_end5_add_ref", 10, "\n");#->3614
&reportAndLogStatus("print_trnscrptID = $print_trnscrptID", 10, "\n");#->3614
&reportAndLogStatus("min_size_split = $min_size_split", 10, "\n");#->3614
&reportAndLogStatus("min_frac_split = $min_frac_split", 10, "\n");#->3614
&reportAndLogStatus("min_summit_dist_split = $min_summit_dist_split", 10, "\n");#->3614
&reportAndLogStatus("use_ref_only_end_pos = $use_ref_only_end_pos", 10, "\n");#->3614
&reportAndLogStatus("retain_no_qry_ref_bound_set = $retain_no_qry_ref_bound_set", 10, "\n");#->3614
foreach my $end_type (sort {$trnscpt_set_end_priority_hsh_ref->{$a} <=> $trnscpt_set_end_priority_hsh_ref->{$b}} keys %{$trnscpt_set_end_priority_hsh_ref}) {
	my $priority = $trnscpt_set_end_priority_hsh_ref->{$end_type};
	&reportAndLogStatus("trnscpt_set_end_priority $priority = $end_type", 10, "\n");#->3614
};

&reportAndLogStatus("doubtful_end_merge_dist = $doubtful_end_merge_dist", 10, "\n");#->3614
&reportAndLogStatus("doubtful_end_avoid_summit = $doubtful_end_avoid_summit", 10, "\n");#->3614
&reportAndLogStatus("min_transcript_length = $min_transcript_length", 10, "\n");#->3614
&reportAndLogStatus("max_thread = $max_thread", 10, "\n");#->3614
&reportAndLogStatus("min_qry_score = $min_qry_score", 10, "\n");#->3614
&reportAndLogStatus("min_output_qry_count = $min_output_qry_count", 10, "\n");#->3614
&reportAndLogStatus("Start processing......", 10, "\n");#->3614

my $remove_ref_transcriptID_version = 'no'; #---[2023/07/21 13:21] will remove all the reference transcript ID version (i.e. ".\d+") if yes
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_process
#
#<section ID="process" num="3">
my ($qry_bed_ary_ref) = &defineQryBedAry($qry_bed_bgz);#->1281
my ($chrom_info_hsh_ref) = &defineChromInfo($chrom_size_path, $result_tmp_dir, $qry_bed_ary_ref, $ref_bed_bgz);#->889
&processPerChromosome($chrom_size_path, $novel_model_prefix, $conf_end3_merge_flank, $conf_end5_merge_flank, $doubtful_end_merge_dist, $doubtful_end_avoid_summit, $conf_end5_bed_bgz, $conf_end3_bed_bgz, $signal_end5_bed_bgz, $signal_end3_bed_bgz, $conf_junction_bed, $conf_end3_add_ref, $conf_end5_add_ref, $bedtools_bin, $tabix_bin, $bgzip_bin, $max_thread, $result_tmp_dir, $chrom_info_hsh_ref, $trnscpt_set_end_priority_hsh_ref, $min_transcript_length, $min_exon_length, $print_trnscrptID, $min_size_split, $min_frac_split, $min_summit_dist_split, $use_ref_only_end_pos, $retain_no_qry_ref_bound_set, $chrom_fasta_path);#->3216
&poolChromLog($chrom_info_hsh_ref, $result_log_dir, $paramTag);#->2286
my ($output_path_hsh_ref) = &definePoolOutputPath($result_bed_dir, $result_log_dir, $paramTag);#->1180
&poolChromBedAndTable($chrom_info_hsh_ref, $output_path_hsh_ref, $result_tmp_dir, $tabix_bin, $max_thread);#->2199
&filterModelAndSetWithQry($bgzip_bin, $tabix_bin, $paramTag, $result_bed_dir, $output_path_hsh_ref, $min_output_qry_count);#->1714
#system "rm -rf $result_tmp_dir";
system "chmod -R 755 $result_dir";

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_finishingTasks
#
#<section ID="finishingTasks" num="4">
&printOutputFileListAndReadme($ARGVStr, $paramTag, $out_dir);#->3000
&printStartOrFinishMessage("finishMessage");#->3115
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

#====================================================================================================================================================#
#List of subroutines by category
#
#	general [n=5]:
#		currentTime, logCalledCMDAndScript, printStartOrFinishMessage
#		readParameters, timeStamp
#
#	log [n=1]:
#		reportAndLogStatus
#
#	output [n=1]:
#		printOutputFileListAndReadme
#
#	time [n=1]:
#		timeStamp
#
#	unassigned [n=36]:
#		assignConfidentEnds, assignDoubtfulEnds, assignEnd3CompleteSetToEnd3CompleteModels
#		assignEnd3IncompleteSetToEnd3CompleteModels, assignEnd3IncompleteSetToEnd3IncompleteModels, assignOrphanSetToExistingModels
#		checkRefTrncptIDNumbering, convertBoundToBed, defineChromInfo
#		defineEndRegion, defineEndSummit, definePoolOutputPath
#		defineQryBedAry, defineRefEndRegBed, defineTrnscptModels
#		defineTrnscptSets, extractTranscriptBounds, filterModelAndSetWithQry
#		getEndNum, getJunctSeq, mergeInConfEndRegion
#		mergeRefRegionAndInConfEndRegion, poolChromBedAndTable, poolChromLog
#		printInfoBedBothEnd, printInfoBedJunction, printInfoModel
#		printInfoModelSetPair, printInfoSet, printInfoTranscript
#		printLogTxt, printModelBed, printTranscriptBed
#		processPerChromosome, readConfJunction, retainNoQryRefBoundSet
#
#====================================================================================================================================================#

sub assignConfidentEnds {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $bedtools_bin, $chrom, $conf_end5_bed_bgz, $confident_end3_bed, $confident_end5_bed, $end3_info_hsh_ref, $end5_info_hsh_ref, $input_bed_hsh_ref, $report_tag, $tabix_bin, $threadID, $trnscpt_info_hsh_ref
#	output: 
#	toCall: &assignConfidentEnds($tabix_bin, $end5_info_hsh_ref, $end3_info_hsh_ref, $bedtools_bin, $trnscpt_info_hsh_ref, $input_bed_hsh_ref, $confident_end5_bed, $confident_end3_bed, $threadID, $chrom, $conf_end5_bed_bgz, $report_tag);
#	calledInLine: 3325
#....................................................................................................................................................#
	my ($tabix_bin, $end5_info_hsh_ref, $end3_info_hsh_ref, $bedtools_bin, $trnscpt_info_hsh_ref, $input_bed_hsh_ref, $confident_end5_bed, $confident_end3_bed, $threadID, $chrom, $conf_end5_bed_bgz, $report_tag) = @_;
	
	my $info_end_hsh_ref = {
		'end5' => {
			'end_bed_key' => 'end5_bed',
			'end_hsh_ref' => $end5_info_hsh_ref,
			'region_bed' => $confident_end5_bed,
			'prefix' => 'F',
		},
		'end3' => {
			'end_bed_key' => 'end3_bed',
			'end_hsh_ref' => $end3_info_hsh_ref,
			'region_bed' => $confident_end3_bed,
			'prefix' => 'T',
		},
	};
	
	foreach my $end (keys %{$info_end_hsh_ref}) {
		my $region_bed = $info_end_hsh_ref->{$end}{'region_bed'};
		my $prefix = $info_end_hsh_ref->{$end}{'prefix'};
		my $end_bed_key = $info_end_hsh_ref->{$end}{'end_bed_key'};
		my $end_hsh_ref = $info_end_hsh_ref->{$end}{'end_hsh_ref'};
		foreach my $input_ID (keys %{$input_bed_hsh_ref}) {
			my $ref_qry = $input_bed_hsh_ref->{$input_ID}{'ref_qry'};
			my $end_bed = $input_bed_hsh_ref->{$input_ID}{$end_bed_key};
			my $bedtools_cmd = "$bedtools_bin intersect -sorted -s -wo -a $end_bed -b $region_bed | cut -f 4,10 |";
			my $num_proc = 0;
			open BEDTOOLS, "$bedtools_cmd";
			while (<BEDTOOLS>) {
				#iPSC_rep1_run1_11847564	F019368
				#iPSC_rep1_run1_16692563	F019368
				#iPSC_rep1_run1_6633243	F019369
				#iPSC_rep1_run1_6165830	F019370
				#iPSC_rep1_run1_3774910	F019371
				chomp;
				my ($trnscpt_ID, $endID) = split /\t/;
				$num_proc++;
				&reportAndLogStatus("$report_tag Assigning $end >> $num_proc transcripts assigned in $input_ID", 10, "\n") if $num_proc%100000 == 0;#->3614
				$trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{$prefix} = $endID;
				$end_hsh_ref->{$endID}{'CT'}{$ref_qry}++;
			}
			close BEDTOOLS;
			&reportAndLogStatus("$report_tag Finished assigning $end to $input_ID ## $num_proc transcripts assigned", 10, "\n");#->3614
		}
	}

	return ();
}
sub assignDoubtfulEnds {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $bedtools_bin, $bgzip_bin, $chrom, $chrom_num, $doubtful_end3_bed, $doubtful_end5_bed, $doubtful_end_merge_dist, $end3_info_hsh_ref, $end5_info_hsh_ref, $log_ary_ref, $report_tag, $trnscpt_info_hsh_ref
#	output: 
#	toCall: &assignDoubtfulEnds($bedtools_bin, $bgzip_bin, $trnscpt_info_hsh_ref, $doubtful_end_merge_dist, $end5_info_hsh_ref, $end3_info_hsh_ref, $report_tag, $doubtful_end5_bed, $doubtful_end3_bed, $chrom, $chrom_num, $log_ary_ref);
#	calledInLine: 3328
#....................................................................................................................................................#
	my ($bedtools_bin, $bgzip_bin, $trnscpt_info_hsh_ref, $doubtful_end_merge_dist, $end5_info_hsh_ref, $end3_info_hsh_ref, $report_tag, $doubtful_end5_bed, $doubtful_end3_bed, $chrom, $chrom_num, $log_ary_ref) = @_;
	
	my $info_end_hsh_ref = {
		'end5' => {
			'end_hsh_ref' => $end5_info_hsh_ref,
			'end_bed' => $doubtful_end5_bed,
			'prefix' => 'F',
		},
		'end3' => {
			'end_hsh_ref' => $end3_info_hsh_ref,
			'end_bed' => $doubtful_end3_bed,
			'prefix' => 'T',
		},
	};
	
	my $count_hsh_ref;
	foreach my $end (qw/end3 end5/) {
		$count_hsh_ref->{$end} = 0;
	}
	
	my $delimiter = '_:_';
	my $num_proc = 0;
	foreach my $end (sort keys %{$info_end_hsh_ref}) {
		my $end_bed = $info_end_hsh_ref->{$end}{'end_bed'};
		my $prefix = $info_end_hsh_ref->{$end}{'prefix'};
		my $end_hsh_ref = $info_end_hsh_ref->{$end}{'end_hsh_ref'};

		my $num_doubful = 0;
		foreach my $ref_qry (sort keys %{$trnscpt_info_hsh_ref}) {
			foreach my $trnscpt_ID (sort keys %{$trnscpt_info_hsh_ref->{$ref_qry}}) {
				next if defined $trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{$prefix};
				$num_doubful++;
			}
		}

		if ($num_doubful == 0) {
			&reportAndLogStatus("No doubtful end$end. Skipping doubtful end assignment.", 10, "\n");
			next;
		}
		open ENDBED, "| sort -k2,2n -k6,6 | $bedtools_bin merge -s -d $doubtful_end_merge_dist -c 4,5,6 -o distinct,sum,distinct -i stdin | $bgzip_bin -c >$end_bed";
		
		foreach my $ref_qry (sort keys %{$trnscpt_info_hsh_ref}) {
			foreach my $trnscpt_ID (sort keys %{$trnscpt_info_hsh_ref->{$ref_qry}}) {
				#my ($exon_num, $strand, $chromStart, $chromEnd, $exon1_start, $exon1_end) = @{$trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'I'}};
				next if defined $trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{$prefix};
				$num_proc++;
				&reportAndLogStatus("$report_tag Collecting doubtful $end >> $num_proc transcripts processed", 10, "\n") if $num_proc%10000 == 0;#->3614
				my $strand = $trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'I'}[0];
				my $chromStart = $trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'I'}[1];
				my $chromEnd = $trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'I'}[2];
				my $end_pos;
				if ($end eq 'end5') {
					if ($strand eq '+') {
						$end_pos = $chromStart + 1;
					} elsif ($strand eq '-') {
						$end_pos = $chromEnd;
					}
				} elsif ($end eq 'end3') {
					if ($strand eq '+') {
						$end_pos = $chromEnd;
					} elsif ($strand eq '-') {
						$end_pos = $chromStart + 1;
					}
				} else {
					die;
				}
				$end_pos = 1 if $end_pos == 0;
				print ENDBED join "", (join "\t", ($chrom, $end_pos-1, $end_pos, $trnscpt_ID.$delimiter.$end_pos.$delimiter.$ref_qry, 1, $strand)), "\n";
			}
		}
		close ENDBED;
		&reportAndLogStatus("$report_tag Finished collecting doubtful $end >> $num_proc transcripts processed", 10, "\n");#->3614

		my $doubtful_num = 0;
		open ENDBED, "gzip -dc $end_bed |";
		while (<ENDBED>) {
			#chr1	31096	31109	ENST00000469289.1:31109,ENST00000473358.1:31097	2	+
			chomp;
			$doubtful_num++;
			&reportAndLogStatus("$report_tag Assigning doubtful $end >> $doubtful_num end group processed", 10, "\n") if $doubtful_num%1000 == 0;#->3614
			my $endID = 'X'.$prefix.$chrom_num.$doubtful_num;
			my (undef, $chromStart, $chromEnd, $ID_str, $score, $strand) = split /\t/;
			my $pos_hsh_ref = {};
			$end_hsh_ref->{$endID}{'score'} = $score;
			$end_hsh_ref->{$endID}{'ST'} = $strand;
			$end_hsh_ref->{$endID}{'CS'} = $chromStart;
			$end_hsh_ref->{$endID}{'CE'} = $chromEnd;
			$end_hsh_ref->{$endID}{'OS'} = $chromStart; #---[2023/07/24 19:28] observed start
			$end_hsh_ref->{$endID}{'OE'} = $chromEnd; #---[2023/07/24 19:28] observed end
			$end_hsh_ref->{$endID}{'ref_only'} = 'N';

			foreach my $trnscpt_ID_pos_ref_qry (split /,/, $ID_str) {
				my ($trnscpt_ID, $pos, $ref_qry) = split /$delimiter/, $trnscpt_ID_pos_ref_qry;
				$end_hsh_ref->{$endID}{'CT'}{$ref_qry}++;
				$pos_hsh_ref->{$pos}++;
				$trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{$prefix} = $endID;
			}
			my $max_count = 0;
			my $pos_count_hsh_ref = {};
			foreach my $pos (keys %{$pos_hsh_ref}) {
				my $count = $pos_hsh_ref->{$pos};
				$max_count = $count if $count > $max_count;
				push @{$pos_count_hsh_ref->{$count}}, $pos;
			}
			
			my @summit_ary = sort {$a <=> $b} @{$pos_count_hsh_ref->{$max_count}};
			my $summit;
			if (($strand eq '+' and $end eq 'end5') or ($strand eq '-' and $end eq 'end3')) {
				$summit = $summit_ary[0];
			} elsif (($strand eq '+' and $end eq 'end3') or ($strand eq '-' and $end eq 'end5')) {
				$summit = $summit_ary[-1];
			} else {
				die;
			}
			$end_hsh_ref->{$endID}{'summit'} = $summit;
			$count_hsh_ref->{$end}++;
		}
		close ENDBED;
		&reportAndLogStatus("$report_tag Finished assigning doubtful $end >> $doubtful_num end group processed", 10, "\n");#->3614
	}

	foreach my $end (qw/end5 end3/) {
		my $count = $count_hsh_ref->{$end};
		my $log_str = "doubtful $end region num";
		push @{$log_ary_ref}, [$log_str, $count];
	}
	
	return ();
}
sub assignEnd3CompleteSetToEnd3CompleteModels {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: defineTrnscptModels|1382
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $end5_ID, $end5_info_hsh_ref, $novel_model_chrom_prefix, $novel_model_num, $set_info_hsh_ref, $within_end5_model_hsh_ref
#	output: $novel_model_num
#	toCall: my ($novel_model_num) = &assignEnd3CompleteSetToEnd3CompleteModels($within_end5_model_hsh_ref, $set_info_hsh_ref, $end5_info_hsh_ref, $novel_model_chrom_prefix, $novel_model_num, $end5_ID);
#	calledInLine: 1405
#....................................................................................................................................................#
	my ($within_end5_model_hsh_ref, $set_info_hsh_ref, $end5_info_hsh_ref, $novel_model_chrom_prefix, $novel_model_num, $end5_ID) = @_;
	
	#---[2023/07/21 14:41] assign complete models
	foreach my $set_ID (keys %{$end5_info_hsh_ref->{$end5_ID}{'set'}}) {
		next if $set_info_hsh_ref->{$set_ID}{'BS'} =~ m/XT/; #---[2023/07/21 14:59] skip if incomplete 3'end
		my $model_ID;
		my $novelty;
		if (exists $set_info_hsh_ref->{$set_ID}{'T'}{'ref'}) {
			my ($ref_trnscpt_ID) = sort keys %{$set_info_hsh_ref->{$set_ID}{'T'}{'ref'}};
			$model_ID = $ref_trnscpt_ID;
			$novelty = 'ref';
		} else {
			#ENST00000685199.1
			$novel_model_num = sprintf "%.9d", $novel_model_num;
			$model_ID = $novel_model_chrom_prefix.$novel_model_num.".1";
			$novelty = 'new';
			$novel_model_num++;
		}
		
		if ($end5_ID =~ m/XF/) {
			$within_end5_model_hsh_ref->{$model_ID}{'IS'}{$set_ID}++; #---[2023/07/21 14:46] incomplete set
		}
		
		$within_end5_model_hsh_ref->{$model_ID}{'CM'} = $set_info_hsh_ref->{$set_ID}{'CM'}; #---[2023/07/21 14:46] complete
		$within_end5_model_hsh_ref->{$model_ID}{'set'} = $set_ID; #---[2023/07/21 14:46] principle cound
		$within_end5_model_hsh_ref->{$model_ID}{'N'} = $novelty;

		$set_info_hsh_ref->{$set_ID}{'AM'} = 'Y'; #---[2023/07/24 14:30] set as model
		$set_info_hsh_ref->{$set_ID}{'model'}{$model_ID}++; #---[2023/07/21 14:55] model ID
	}
	
	return ($novel_model_num);
}
sub assignEnd3IncompleteSetToEnd3CompleteModels {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: defineTrnscptModels|1382
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $end5_ID, $end5_info_hsh_ref, $set_info_hsh_ref, $strand, $within_end5_model_hsh_ref
#	output: 
#	toCall: &assignEnd3IncompleteSetToEnd3CompleteModels($strand, $set_info_hsh_ref, $within_end5_model_hsh_ref, $end5_info_hsh_ref, $end5_ID);
#	calledInLine: 1408
#....................................................................................................................................................#
	my ($strand, $set_info_hsh_ref, $within_end5_model_hsh_ref, $end5_info_hsh_ref, $end5_ID) = @_;

	my $debug_set_ID = 'S084101';
	foreach my $qry_set_ID (keys %{$end5_info_hsh_ref->{$end5_ID}{'set'}}) {
		next if exists $set_info_hsh_ref->{$qry_set_ID}{'model'}; #---[2023/07/21 16:54] only for anything no assigned;
		
		#---[2023/07/21 14:49] only non-model incomplete set left
		my $trimmed_qry_incomplete_bound_str = $set_info_hsh_ref->{$qry_set_ID}{'BS'};
		$trimmed_qry_incomplete_bound_str =~ s/\_XT.+$//; #---[2023/07/21 14:59] trimm the incomplete 3'end
		#&reportAndLogStatus("$trimmed_qry_incomplete_bound_str", 10, "\n") if $qry_set_ID eq $debug_set_ID;#->3614
		
		foreach my $model_ID (keys %{$within_end5_model_hsh_ref}) {
			my $model_set_ID = $within_end5_model_hsh_ref->{$model_ID}{'set'};
			my $model_bound_str = $set_info_hsh_ref->{$model_set_ID}{'BS'};
		
			my $assign = 'N';
			if ($trimmed_qry_incomplete_bound_str =~ m/J/) { #---[2023/07/21 15:35] spliced set
				$assign = 'Y' if ($model_bound_str =~ m/^$trimmed_qry_incomplete_bound_str/);#---[2023/07/21 15:00] the model_bound_str contains the whole trimmed_incomplete_bound_str
			} else {#---[2023/07/21 15:39] if qry_set is not-spliced : automatically implied the model_set is spliced since within an end5 there is only one set is not splice
				my ($qry_exon1_start, $qry_exon1_end) = @{$set_info_hsh_ref->{$qry_set_ID}{'E1'}};
				my ($model_exon1_start, $model_exon1_end) = @{$set_info_hsh_ref->{$model_set_ID}{'E1'}};
				$assign = 'Y' if $model_exon1_start <= $qry_exon1_start and $model_exon1_end >= $qry_exon1_end;
			}

			if ($assign eq 'Y') {
				$set_info_hsh_ref->{$qry_set_ID}{'AM'} = 'N'; #---[2023/07/22 15:11] As Model no
				$set_info_hsh_ref->{$qry_set_ID}{'model'}{$model_ID}++;
				$within_end5_model_hsh_ref->{$model_ID}{'IS'}{$qry_set_ID}++; #---[2023/07/21 14:46] incomplete set
			}
		}
	}	

	return ();
}
sub assignEnd3IncompleteSetToEnd3IncompleteModels {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: defineTrnscptModels|1382
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $end5_ID, $end5_info_hsh_ref, $novel_model_chrom_prefix, $novel_model_num, $set_info_hsh_ref, $strand, $within_end5_model_hsh_ref
#	output: $novel_model_num
#	toCall: my ($novel_model_num) = &assignEnd3IncompleteSetToEnd3IncompleteModels($within_end5_model_hsh_ref, $set_info_hsh_ref, $end5_info_hsh_ref, $novel_model_chrom_prefix, $novel_model_num, $strand, $end5_ID);
#	calledInLine: 1411
#....................................................................................................................................................#
	my ($within_end5_model_hsh_ref, $set_info_hsh_ref, $end5_info_hsh_ref, $novel_model_chrom_prefix, $novel_model_num, $strand, $end5_ID) = @_;

	my $unassigned_set_hsh_ref = {};
	foreach my $set_ID (keys %{$end5_info_hsh_ref->{$end5_ID}{'set'}}) {
		$unassigned_set_hsh_ref->{$set_ID}++ if not exists $set_info_hsh_ref->{$set_ID}{'model'};
	}

	#my $debug_set_ID = 'S072847';
	#my $debug_set_ID = 'S075873';
	my $debug_set_ID = 'S075213';
	foreach my $ref_set_ID (keys %{$unassigned_set_hsh_ref}) {
		
		#---[2023/07/21 17:43] test if the ref_set_ID is the longest and make new model if yes
		my $longest = 'Y';
		my $qry_assign_to_ref_hsh_ref = {};
		my $trimmed_ref_incomplete_bound_str = $set_info_hsh_ref->{$ref_set_ID}{'BS'};
		$trimmed_ref_incomplete_bound_str =~ s/\_XT.+$//; #---[2023/07/21 14:59] trimm the incomplete 3'end

		my $ref_size = $set_info_hsh_ref->{$ref_set_ID}{'CE'} - $set_info_hsh_ref->{$ref_set_ID}{'CS'};
		foreach my $qry_set_ID (keys %{$unassigned_set_hsh_ref}) {
			next if $qry_set_ID eq $ref_set_ID;
			my $trimmed_qry_incomplete_bound_str = $set_info_hsh_ref->{$qry_set_ID}{'BS'};
			$trimmed_qry_incomplete_bound_str =~ s/\_XT.+$//; #---[2023/07/21 14:59] trimm the incomplete 3'end
			my $qry_size = $set_info_hsh_ref->{$qry_set_ID}{'CE'} - $set_info_hsh_ref->{$qry_set_ID}{'CS'};

			#&reportAndLogStatus("##################################", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614
			#&reportAndLogStatus("ref_set_ID=$ref_set_ID qry_set_ID=$qry_set_ID", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614
			#&reportAndLogStatus("$trimmed_ref_incomplete_bound_str $trimmed_qry_incomplete_bound_str", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614

			my $ref_assign_to_qry = 'N';
			my $qry_assign_to_ref = 'N';
			if ($trimmed_ref_incomplete_bound_str =~ m/J/ and $trimmed_qry_incomplete_bound_str =~ m/J/) { #---[2023/07/21 15:35] spliced set
			#if ($trimmed_ref_incomplete_bound_str =~ m/J/) { #---[2023/07/21 15:35] spliced set
				if ($trimmed_ref_incomplete_bound_str ne $trimmed_qry_incomplete_bound_str) {
					$ref_assign_to_qry = 'Y' if ($trimmed_qry_incomplete_bound_str =~ m/^$trimmed_ref_incomplete_bound_str/);#---[2023/07/21 15:00] the qry contains the whole ref str
					$qry_assign_to_ref = 'Y' if ($trimmed_ref_incomplete_bound_str =~ m/^$trimmed_qry_incomplete_bound_str/);#---[2023/07/21 15:00] the ref contains the whole qry str
				} else {#---[2023/07/24 18:11] same splicing structure but different 3'end
					#die "DEBUG: $ref_set_ID vs $qry_set_ID ref_size == qry_size\n" if $ref_size == $qry_size;
					if ($ref_size >= $qry_size) {
						$qry_assign_to_ref = 'Y';
					} else {
						$ref_assign_to_qry = 'Y';
					}
				}

			} else {
			
				my $ref_within_qry = 'N';
				my $qry_within_ref = 'N';
				my ($qry_exon1_start, $qry_exon1_end) = @{$set_info_hsh_ref->{$qry_set_ID}{'E1'}};
				my ($ref_exon1_start, $ref_exon1_end) = @{$set_info_hsh_ref->{$ref_set_ID}{'E1'}};
				
				#---[2023/08/18 1:57] adjust the 5'end to the same position (for commonest 5end might be different )
				if ($strand eq '+') {
					$qry_exon1_start = $ref_exon1_start;
				} elsif ($strand eq '-') {
					$qry_exon1_end = $qry_exon1_end;
				}
				
				$qry_within_ref = 'Y' if $ref_exon1_start <= $qry_exon1_start and $ref_exon1_end >= $qry_exon1_end;
				$ref_within_qry = 'Y' if $qry_exon1_start <= $ref_exon1_start and $qry_exon1_end >= $ref_exon1_end;
				#&reportAndLogStatus("ref $ref_exon1_start $ref_exon1_end", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614
				#&reportAndLogStatus("qry $qry_exon1_start $qry_exon1_end", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614
				#&reportAndLogStatus("ref_within_qry->$ref_within_qry qry_within_ref->$qry_within_ref", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614

				if ($trimmed_ref_incomplete_bound_str !~ m/J/ and $trimmed_qry_incomplete_bound_str =~ m/J/) {#---[2023/07/21 15:39] if ref_set is not-spliced set: automatically implied the qry_set is spliced since within an end5 there is only one set is not spliced
						$ref_assign_to_qry = 'Y' if $ref_within_qry eq 'Y';
						#$ref_assign_to_qry = 'Y' if $ref_within_qry eq 'N'; #---[2023/07/25 20:03] to intentionally create bug
				} elsif ($trimmed_ref_incomplete_bound_str =~ m/J/ and $trimmed_qry_incomplete_bound_str !~ m/J/) {#---[2023/07/21 15:39] if ref_set is not-spliced set: automatically implied the qry_set is spliced since within an end5 there is only one set is not spliced
						$qry_assign_to_ref = 'Y' if $qry_within_ref eq 'Y';
				} elsif ($trimmed_ref_incomplete_bound_str !~ m/J/ and $trimmed_qry_incomplete_bound_str !~ m/J/) {#---[2023/07/21 15:39] if ref_set is not-spliced set: automatically implied the qry_set is spliced since within an end5 there is only one set is not spliced
						$ref_assign_to_qry = 'Y' if $ref_within_qry eq 'Y';
						$qry_assign_to_ref = 'Y' if $qry_within_ref eq 'Y';
				} else {
					die "DEBUG: $trimmed_ref_incomplete_bound_str $trimmed_qry_incomplete_bound_str\n";
				}
			}

			#&reportAndLogStatus("ref_assign_to_qry->$ref_assign_to_qry qry_assign_to_ref->$qry_assign_to_ref", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614

			#---[2023/07/21 17:43] if it can be assigned 
			$longest = 'N' if ($ref_assign_to_qry eq 'Y');
			$qry_assign_to_ref_hsh_ref->{$qry_set_ID}++ if $qry_assign_to_ref eq 'Y';
		}
		
		#&reportAndLogStatus("longest=$longest", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614
		
		if ($longest eq 'Y') {
			my $model_ID;
			my $novelty;
			if (exists $set_info_hsh_ref->{$ref_set_ID}{'T'}{'ref'}) {
				($model_ID) = sort keys %{$set_info_hsh_ref->{$ref_set_ID}{'T'}{'ref'}};
				$novelty = 'ref'; #---[2023/07/21 18:03] should not contain ref if added ref 3'end as confident end3
			} else {
				#ENST00000685199.1
				$novel_model_num = sprintf "%.9d", $novel_model_num;
				$model_ID = $novel_model_chrom_prefix.$novel_model_num.".1";
				$novelty = 'new';
				$novel_model_num++;
			}
			$within_end5_model_hsh_ref->{$model_ID}{'CM'} =  $set_info_hsh_ref->{$ref_set_ID}{'CM'}; #---[2023/07/21 14:46] complete
			$within_end5_model_hsh_ref->{$model_ID}{'set'} = $ref_set_ID; #---[2023/07/21 14:46] principle cound
			$within_end5_model_hsh_ref->{$model_ID}{'N'} = $novelty;

			$set_info_hsh_ref->{$ref_set_ID}{'model'}{$model_ID}++; #---[2023/07/21 14:55] model ID
			$set_info_hsh_ref->{$ref_set_ID}{'AM'} = 'Y';

			foreach my $qry_set_ID (keys %{$qry_assign_to_ref_hsh_ref}) {
				#&reportAndLogStatus("qry_set_ID $qry_set_ID assigned", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614
				$set_info_hsh_ref->{$qry_set_ID}{'AM'} = 'N';
				$set_info_hsh_ref->{$qry_set_ID}{'model'}{$model_ID}++;
				$within_end5_model_hsh_ref->{$model_ID}{'IS'}{$qry_set_ID}++; #---[2023/07/21 14:46] incomplete set
			}
		}
	}

	return ($novel_model_num);
}
sub assignOrphanSetToExistingModels {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: defineTrnscptModels|1382
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $end5_ID, $end5_info_hsh_ref, $novel_model_chrom_prefix, $novel_model_num, $set_info_hsh_ref, $strand, $within_end5_model_hsh_ref
#	output: 
#	toCall: &assignOrphanSetToExistingModels($within_end5_model_hsh_ref, $set_info_hsh_ref, $end5_info_hsh_ref, $novel_model_chrom_prefix, $novel_model_num, $strand, $end5_ID);
#	calledInLine: 1414
#....................................................................................................................................................#
	my ($within_end5_model_hsh_ref, $set_info_hsh_ref, $end5_info_hsh_ref, $novel_model_chrom_prefix, $novel_model_num, $strand, $end5_ID) = @_;

	#---for the sets can be assigned to another set but the set itself is not the longest

	my $unassigned_set_hsh_ref = {};
	my $assigned_set_hsh_ref = {};
	foreach my $set_ID (keys %{$end5_info_hsh_ref->{$end5_ID}{'set'}}) {
		if (not exists $set_info_hsh_ref->{$set_ID}{'model'}) {
			$unassigned_set_hsh_ref->{$set_ID}++;
		} else {
			$assigned_set_hsh_ref->{$set_ID}++;
		}
	}

	#my $debug_set_ID = 'S072847';
	#my $debug_set_ID = 'S075873';
	my $debug_set_ID = 'S075213';

	foreach my $ref_set_ID (keys %{$unassigned_set_hsh_ref}) {
		
		my $trimmed_ref_incomplete_bound_str = $set_info_hsh_ref->{$ref_set_ID}{'BS'};
		$trimmed_ref_incomplete_bound_str =~ s/\_XT.+$//; #---[2023/07/21 14:59] trimm the incomplete 3'end

		my $ref_size = $set_info_hsh_ref->{$ref_set_ID}{'CE'} - $set_info_hsh_ref->{$ref_set_ID}{'CS'};
		foreach my $qry_set_ID (sort keys %{$assigned_set_hsh_ref}) {
			my $trimmed_qry_incomplete_bound_str = $set_info_hsh_ref->{$qry_set_ID}{'BS'};
			$trimmed_qry_incomplete_bound_str =~ s/\_XT.+$//; #---[2023/07/21 14:59] trimm the incomplete 3'end

			my $ref_assign_to_qry = 'N';

			if ($trimmed_ref_incomplete_bound_str =~ m/J/ and $trimmed_qry_incomplete_bound_str =~ m/J/) { #---[2023/07/21 15:35] spliced set
					$ref_assign_to_qry = 'Y' if ($trimmed_qry_incomplete_bound_str =~ m/^$trimmed_ref_incomplete_bound_str/);#---[2023/07/21 15:00] the qry contains the whole ref str

			} else {
			
				my $ref_within_qry = 'N';
				my $qry_within_ref = 'N';
				my ($qry_exon1_start, $qry_exon1_end) = @{$set_info_hsh_ref->{$qry_set_ID}{'E1'}};
				my ($ref_exon1_start, $ref_exon1_end) = @{$set_info_hsh_ref->{$ref_set_ID}{'E1'}};
				$qry_within_ref = 'Y' if $ref_exon1_start <= $qry_exon1_start and $ref_exon1_end >= $qry_exon1_end;
				$ref_within_qry = 'Y' if $qry_exon1_start <= $ref_exon1_start and $qry_exon1_end >= $ref_exon1_end;

				$ref_assign_to_qry = 'Y' if $ref_within_qry eq 'Y';

				#&reportAndLogStatus("##################################", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614
				#&reportAndLogStatus("ref_set_ID=$ref_set_ID qry_set_ID=$qry_set_ID", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614
				#&reportAndLogStatus("$trimmed_ref_incomplete_bound_str $trimmed_qry_incomplete_bound_str", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614
				#&reportAndLogStatus("ref $ref_exon1_start $ref_exon1_end", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614
				#&reportAndLogStatus("qry $qry_exon1_start $qry_exon1_end", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614
				#&reportAndLogStatus("ref_within_qry->$ref_within_qry", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614
				#&reportAndLogStatus("ref_assign_to_qry->$ref_assign_to_qry", 10, "\n") if $ref_set_ID eq $debug_set_ID;#->3614
				
				
			}


			if ($ref_assign_to_qry eq 'Y') {
				foreach my $model_ID (keys %{$set_info_hsh_ref->{$qry_set_ID}{'model'}}) {
					#---[2023/07/25 19:23] for all the models qry_set_ID assigned too
					$within_end5_model_hsh_ref->{$model_ID}{'IS'}{$ref_set_ID}++; #---[2023/07/21 14:46] incomplete set
					$set_info_hsh_ref->{$ref_set_ID}{'AM'} = 'N';
					$set_info_hsh_ref->{$ref_set_ID}{'model'}{$model_ID}++;
				}
			}
		}
	}
	
	foreach my $set_ID (keys %{$end5_info_hsh_ref->{$end5_ID}{'set'}}) {
		die "DEBUG: set_ID $set_ID has no models assigned\n" if (not exists $set_info_hsh_ref->{$set_ID}{'model'});
	}

	return ();
}
sub checkRefTrncptIDNumbering {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $chrom_num, $novel_model_chrom_prefix, $trnscpt_info_hsh_ref
#	output: $inital_novel_model_num
#	toCall: my ($inital_novel_model_num) = &checkRefTrncptIDNumbering($chrom_num, $novel_model_chrom_prefix, $trnscpt_info_hsh_ref);
#	calledInLine: 3334
#....................................................................................................................................................#
	my ($chrom_num, $novel_model_chrom_prefix, $trnscpt_info_hsh_ref) = @_;
	
	my $inital_novel_model_num = 0;
	foreach my $ref_trnscpt_ID (keys %{$trnscpt_info_hsh_ref->{'ref'}}) {
		if ($ref_trnscpt_ID =~ m/^$novel_model_chrom_prefix(\d+)/) {
			my $ref_novel_model_num = int($1);
			$inital_novel_model_num = $ref_novel_model_num + 1 if $ref_novel_model_num > $inital_novel_model_num;
		}
	}
	
	return ($inital_novel_model_num);
}
sub convertBoundToBed {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: defineTrnscptSets|1431
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $RR, $bound_str_ary_ref, $chrom, $chromStartEnd_hsh_ref, $doubtful_end_avoid_summit, $junct_info_hsh_ref, $min_transcript_length, $qry_count, $set_ID, $test_endtype_ary_ref
#	output: $CE_end_type, $CS_end_type, $bed_line_ary_ref, $exon1_end, $exon1_start
#	toCall: my ($bed_line_ary_ref, $CS_end_type, $CE_end_type, $exon1_start, $exon1_end) = &convertBoundToBed($junct_info_hsh_ref, $bound_str_ary_ref, $chromStartEnd_hsh_ref, $set_ID, $min_transcript_length, $test_endtype_ary_ref, $qry_count, $chrom, $doubtful_end_avoid_summit, $RR);
#	calledInLine: 1547
#....................................................................................................................................................#
	my ($junct_info_hsh_ref, $bound_str_ary_ref, $chromStartEnd_hsh_ref, $set_ID, $min_transcript_length, $test_endtype_ary_ref, $qry_count, $chrom, $doubtful_end_avoid_summit, $RR) = @_;
	#e.g. 0-end5,1-J,2-J,3-end3,4-strand 1 to 2 are J

	my @intron_bound_ary = ();
	my $strand = $bound_str_ary_ref->[-1];
	foreach my $index (1..($#{$bound_str_ary_ref}-2)) {
		my $junct_ID = $bound_str_ary_ref->[$index];
		my ($intron_start, $intron_end) = split /:/, $junct_info_hsh_ref->{$junct_ID}{'JS'};
		push @intron_bound_ary, ($intron_start, $intron_end);
	}
	@intron_bound_ary = sort {$a <=> $b} @intron_bound_ary;

	my $CS_end_type;
	my $CE_end_type;
	if ($RR eq 'Y') {
		$CS_end_type = 'original';
		$CE_end_type = 'original';
	} else {
		#---[2023/07/24 23:29] check for conflicts by end_type sorted in priority
		foreach my $test_CS_end_type (@{$test_endtype_ary_ref}) {
			if (($doubtful_end_avoid_summit eq 'yes' and $test_CS_end_type eq 'summit') and 
				(($bound_str_ary_ref->[0] =~ m/XF/ and $strand eq '+') or ($bound_str_ary_ref->[-2] =~ m/XT/ and $strand eq '-'))) { #---[2023/08/15 14:17] prevent doubtful ends to assign to submit
				next;
			}
			foreach my $test_CE_end_type (@{$test_endtype_ary_ref}) {
				if (($doubtful_end_avoid_summit eq 'yes' and $test_CE_end_type eq 'summit') and 
					(($bound_str_ary_ref->[0] =~ m/XF/ and $strand eq '-') or ($bound_str_ary_ref->[-2] =~ m/XT/ and $strand eq '+'))) { #---[2023/08/15 14:17] prevent doubtful ends to assign to submit
					next;
				}
				my $test_CS = $chromStartEnd_hsh_ref->{'CS'}{$test_CS_end_type};
				my $test_CE = $chromStartEnd_hsh_ref->{'CE'}{$test_CE_end_type};
				my $transcript_length = $test_CE - $test_CS;
				next if $transcript_length < $min_transcript_length;
				if (@intron_bound_ary == 0) {
					$CS_end_type = $test_CS_end_type;
					$CE_end_type = $test_CE_end_type;
				} else {
					my @tmp_ary = ($test_CS, @intron_bound_ary, $test_CE);
					$CS_end_type = $test_CS_end_type if $tmp_ary[1] > $tmp_ary[0];
					$CE_end_type = $test_CE_end_type if $tmp_ary[-1] > $tmp_ary[-2];
				}
				last if (defined $CS_end_type and defined $CE_end_type);
			}
			last if (defined $CS_end_type and defined $CE_end_type);
		}
	die "DEBUG: failed to choose endtype: set_ID $set_ID " if not defined $CS_end_type or not defined $CE_end_type;
	}
	
	my @rngAry = ($chromStartEnd_hsh_ref->{'CS'}{$CS_end_type}, @intron_bound_ary, $chromStartEnd_hsh_ref->{'CE'}{$CE_end_type});
	my $exon1_start;
	my $exon1_end;
	if ($strand eq '+') {
		($exon1_start, $exon1_end) = ($rngAry[0], $rngAry[1]);
	} else {
		($exon1_start, $exon1_end) = ($rngAry[-2], $rngAry[-1]);
	}
	my $chromStart = $rngAry[0];
	my $chromEnd = $rngAry[-1];
	my $score = $qry_count;
	my $thickStart = $chromStart;
	my $thickEnd = $chromEnd;
	my $blockCount = 0;
	my @blockSizesAry = ();
	my @blockStartsAry = ();
	for (my $i=0; $i < $#rngAry; $i += 2) {
		$blockCount++;
		my ($start, $end) = ($rngAry[$i], $rngAry[$i+1]);
		my $blockSize = $end - $start;
		my $blockStart = $start - $chromStart;
		push @blockSizesAry, $blockSize;
		push @blockStartsAry, $blockStart;
	}
	my $blockSizes = join ",", @blockSizesAry;
	my $blockStarts = join ",", @blockStartsAry;
	
	my $itemRgb = '228,26,28';
	$itemRgb = '55,126,184' if ($strand eq '-');
	
	my $bed_line_ary_ref = [$chrom, $chromStart, $chromEnd, $set_ID, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts];
	
	return ($bed_line_ary_ref, $CS_end_type, $CE_end_type, $exon1_start, $exon1_end);
}
sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: printStartOrFinishMessage|3115, reportAndLogStatus|3614
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 2_defineout_dirPath|153, 4_finishingTasks|214
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 126, 3131, 3135, 3140, 3144, 3630, 3631
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub defineChromInfo {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: >none
#	primaryAppearInSection: 3_process|197
#	secondaryAppearInSection: >none
#	input: $chrom_size_path, $qry_bed_ary_ref, $ref_bed_bgz, $result_tmp_dir
#	output: $chrom_info_hsh_ref
#	toCall: my ($chrom_info_hsh_ref) = &defineChromInfo($chrom_size_path, $result_tmp_dir, $qry_bed_ary_ref, $ref_bed_bgz);
#	calledInLine: 201
#....................................................................................................................................................#
	my ($chrom_size_path, $result_tmp_dir, $qry_bed_ary_ref, $ref_bed_bgz) = @_;
	
	my $chrom_info_hsh_ref = {};
	&reportAndLogStatus("Reading chromosome size", 10, "\n");#->3614
	
	open (CHROMSIZE, "<", $chrom_size_path);
	while (<CHROMSIZE>) {
		chomp;
		my ($chrom, $size) = split /\t/;
		#next if $chrom ne 'chr1' and $chrom ne 'chr3' and $chrom ne 'chr5';
		#next if $chrom ne 'chr1';
		my $chrom_num = $chrom;
		$chrom_num =~ s/chr//;
		$chrom_num = 23 if $chrom_num =~ m/^X/;
		$chrom_num = 24 if $chrom_num eq 'Y';
		$chrom_num = 25 if $chrom_num eq 'M';
		$chrom_num = sprintf "%.2d", $chrom_num;
		my $chrom_dir = "$result_tmp_dir/$chrom";
		system ("mkdir -pm 755 $chrom_dir");

		$chrom_info_hsh_ref->{$chrom}{'size'} = $size;
		$chrom_info_hsh_ref->{$chrom}{'chrom_num'} = $chrom_num;
		$chrom_info_hsh_ref->{$chrom}{'tmp_txt'} = "$chrom_dir/$chrom.tmp.txt";
		$chrom_info_hsh_ref->{$chrom}{'log_txt'} = "$chrom_dir/$chrom.log.txt";

		$chrom_info_hsh_ref->{$chrom}{'ref_in_conf_merge_end3_reg_bed'} = "$chrom_dir/$chrom.ref_in_conf_merge.end3.reg.bed";
		$chrom_info_hsh_ref->{$chrom}{'ref_in_conf_merge_end5_reg_bed'} = "$chrom_dir/$chrom.ref_in_conf_merge.end5.reg.bed";
		$chrom_info_hsh_ref->{$chrom}{'split_merge_in_conf_end3_bed'} = "$chrom_dir/$chrom.split_merge_in_conf.end3.bed";
		$chrom_info_hsh_ref->{$chrom}{'split_merge_in_conf_end5_bed'} = "$chrom_dir/$chrom.split_merge_in_conf.end5.bed";
		$chrom_info_hsh_ref->{$chrom}{'confident_end5_bed'} = "$chrom_dir/$chrom.confident.end5.bed.bgz";
		$chrom_info_hsh_ref->{$chrom}{'confident_end3_bed'} = "$chrom_dir/$chrom.confident.end3.bed.bgz";
		$chrom_info_hsh_ref->{$chrom}{'doubtful_end5_bed'} = "$chrom_dir/$chrom.doubtful.end5.bed.bgz";
		$chrom_info_hsh_ref->{$chrom}{'doubtful_end3_bed'} = "$chrom_dir/$chrom.doubtful.end3.bed.bgz";
		$chrom_info_hsh_ref->{$chrom}{'ref_end3_reg_bed'} = "$chrom_dir/$chrom.ref.end3.reg.bed";
		$chrom_info_hsh_ref->{$chrom}{'ref_end5_reg_bed'} = "$chrom_dir/$chrom.ref.end5.reg.bed";
		$chrom_info_hsh_ref->{$chrom}{'ref_end3_pileup_bed'} = "$chrom_dir/$chrom.ref.end3.pileup.bed";
		$chrom_info_hsh_ref->{$chrom}{'ref_end5_pileup_bed'} = "$chrom_dir/$chrom.ref.end5.pileup.bed";
		$chrom_info_hsh_ref->{$chrom}{'all_end5_bed'} = "$chrom_dir/$chrom.all.end5.bed.bgz";
		$chrom_info_hsh_ref->{$chrom}{'all_end3_bed'} = "$chrom_dir/$chrom.all.end3.bed.bgz";

		$chrom_info_hsh_ref->{$chrom}{'set_bed'} = "$chrom_dir/$chrom.set.bed.bgz";
		$chrom_info_hsh_ref->{$chrom}{'model_bed'} = "$chrom_dir/$chrom.model.bed.bgz";
		$chrom_info_hsh_ref->{$chrom}{'ref_trnscpt_bed'} = "$chrom_dir/$chrom.trnscpt.ref.bed.bgz";
		$chrom_info_hsh_ref->{$chrom}{'qry_trnscpt_bed'} = "$chrom_dir/$chrom.trnscpt.qry.bed.bgz";
		$chrom_info_hsh_ref->{$chrom}{'junct_bed'} = "$chrom_dir/$chrom.junct.bed.bgz";
		$chrom_info_hsh_ref->{$chrom}{'discard_trnscpt_bed'} = "$chrom_dir/$chrom.discard.trnscpt.bed.bgz";

		$chrom_info_hsh_ref->{$chrom}{'splicing_site_bed'} = "$chrom_dir/$chrom.splicing_site.bed";
		$chrom_info_hsh_ref->{$chrom}{'splicing_site_fasta'} = "$chrom_dir/$chrom.splicing_site.fasta";

		$chrom_info_hsh_ref->{$chrom}{'trnscpt_info_tsv'} = "$chrom_dir/$chrom.transcript.info.tsv.gz";
		$chrom_info_hsh_ref->{$chrom}{'model_info_tsv'} = "$chrom_dir/$chrom.model.info.tsv.gz";
		$chrom_info_hsh_ref->{$chrom}{'set_info_tsv'} = "$chrom_dir/$chrom.set.info.tsv.gz";
		$chrom_info_hsh_ref->{$chrom}{'end5_info_tsv'} = "$chrom_dir/$chrom.end5.info.tsv.gz";
		$chrom_info_hsh_ref->{$chrom}{'end3_info_tsv'} = "$chrom_dir/$chrom.end3.info.tsv.gz";
		$chrom_info_hsh_ref->{$chrom}{'junct_info_tsv'} = "$chrom_dir/$chrom.junct.info.tsv.gz";
		$chrom_info_hsh_ref->{$chrom}{'model_set_pair_info_tsv'} = "$chrom_dir/$chrom.model_set_pair.info.tsv.gz";

		my $input_bed_hsh_ref = {
			'ref' => {
				'ref_qry' => 'ref',
				'suffix' => 'dummy', #---[2023/07/22 13:30] dummy and will not be used
				'trnscpt_bed' => $ref_bed_bgz,
				'anno_trnscpt_bed' => "$chrom_dir/$chrom.ref.anno.trnscpt.bed.bgz",
				'end5_bed' => "$chrom_dir/$chrom.ref.end5.bed.bgz",
				'end3_bed' => "$chrom_dir/$chrom.ref.end3.bed.bgz",
			},
		};
		foreach my $index (0..$#{$qry_bed_ary_ref}) {
			my $input_ID = "qry_$index";
			#&reportAndLogStatus("storting $input_ID $chrom path", 10, "\n");#->3614
			my $suffix = $index;
			my $ref_qry = 'qry';
			my $trnscpt_bed = $qry_bed_ary_ref->[$index];
			my $end5_bed = "$chrom_dir/$chrom.$input_ID.end5.bed.gz";
			my $end3_bed = "$chrom_dir/$chrom.$input_ID.end3.bed.gz";
			my $anno_trnscpt_bed = "$chrom_dir/$chrom.$input_ID.anno.trnscpt.bed.bgz";
			$input_bed_hsh_ref->{$input_ID}{'ref_qry'} = $ref_qry;
			$input_bed_hsh_ref->{$input_ID}{'suffix'} = $suffix;
			$input_bed_hsh_ref->{$input_ID}{'trnscpt_bed'} = $trnscpt_bed;
			$input_bed_hsh_ref->{$input_ID}{'end5_bed'} = $end5_bed;
			$input_bed_hsh_ref->{$input_ID}{'end3_bed'} = $end3_bed;
			$input_bed_hsh_ref->{$input_ID}{'anno_trnscpt_bed'} = $anno_trnscpt_bed;
		}

		$chrom_info_hsh_ref->{$chrom}{'input_bed_hsh_ref'} = $input_bed_hsh_ref;
		
		&reportAndLogStatus("chrom $chrom [#$chrom_num] = $size nt", 10, "\n");#->3614
	}
	close CHROMSIZE;

	return ($chrom_info_hsh_ref);
}
sub defineEndRegion {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: defineEndSummit|1073, defineRefEndRegBed|1320, mergeInConfEndRegion|1932, mergeRefRegionAndInConfEndRegion|2082, reportAndLogStatus|3614
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $bedtools_bin, $bgzip_bin, $chrom, $chrom_info_hsh_ref, $chrom_num, $chrom_size_path, $conf_end3_add_ref, $conf_end3_bed_bgz, $conf_end3_merge_flank, $conf_end5_add_ref, $conf_end5_bed_bgz, $conf_end5_merge_flank, $confident_end3_bed, $confident_end5_bed, $input_bed_hsh_ref, $log_ary_ref, $min_frac_split, $min_size_split, $min_summit_dist_split, $report_tag, $signal_end3_bed_bgz, $signal_end5_bed_bgz, $size, $tabix_bin, $threadID, $use_ref_only_end_pos
#	output: $end3_info_hsh_ref, $end5_info_hsh_ref
#	toCall: my ($end5_info_hsh_ref, $end3_info_hsh_ref) = &defineEndRegion($chrom_info_hsh_ref, $input_bed_hsh_ref, $conf_end5_add_ref, $conf_end3_add_ref, $bedtools_bin, $tabix_bin, $bgzip_bin, $conf_end5_merge_flank, $conf_end3_merge_flank, $conf_end5_bed_bgz, $conf_end3_bed_bgz, $confident_end5_bed, $confident_end3_bed, $signal_end5_bed_bgz, $signal_end3_bed_bgz, $chrom_size_path, $chrom, $chrom_num, $threadID, $log_ary_ref, $size, $min_size_split, $min_frac_split, $min_summit_dist_split, $use_ref_only_end_pos, $report_tag);
#	calledInLine: 3312
#....................................................................................................................................................#
	my ($chrom_info_hsh_ref, $input_bed_hsh_ref, $conf_end5_add_ref, $conf_end3_add_ref, $bedtools_bin, $tabix_bin, $bgzip_bin, $conf_end5_merge_flank, $conf_end3_merge_flank, $conf_end5_bed_bgz, $conf_end3_bed_bgz, $confident_end5_bed, $confident_end3_bed, $signal_end5_bed_bgz, $signal_end3_bed_bgz, $chrom_size_path, $chrom, $chrom_num, $threadID, $log_ary_ref, $size, $min_size_split, $min_frac_split, $min_summit_dist_split, $use_ref_only_end_pos, $report_tag) = @_;
	
	my $ref_bed_bgz = $input_bed_hsh_ref->{'ref'}{'trnscpt_bed'};
	my $end5_info_hsh_ref = {};
	my $end3_info_hsh_ref = {};
	my $ref_in_conf_merge_end3_reg_bed = $chrom_info_hsh_ref->{$chrom}{'ref_in_conf_merge_end3_reg_bed'};
	my $ref_in_conf_merge_end5_reg_bed = $chrom_info_hsh_ref->{$chrom}{'ref_in_conf_merge_end5_reg_bed'};
	my $split_merge_in_conf_end3_bed = $chrom_info_hsh_ref->{$chrom}{'split_merge_in_conf_end3_bed'};
	my $split_merge_in_conf_end5_bed = $chrom_info_hsh_ref->{$chrom}{'split_merge_in_conf_end5_bed'};
	my $ref_end3_pileup_bed = $chrom_info_hsh_ref->{$chrom}{'ref_end3_pileup_bed'};
	my $ref_end5_pileup_bed = $chrom_info_hsh_ref->{$chrom}{'ref_end5_pileup_bed'};
	my $ref_end3_reg_bed = $chrom_info_hsh_ref->{$chrom}{'ref_end3_reg_bed'};
	my $ref_end5_reg_bed = $chrom_info_hsh_ref->{$chrom}{'ref_end5_reg_bed'};

	my $end_hsh_ref = {
		'end5' => {
			'conf_bed_bgz' => $conf_end5_bed_bgz,
			'conf_merge_flank' => $conf_end5_merge_flank,
			'split_merge_in_conf_bed' => $split_merge_in_conf_end5_bed,
			'ref_reg_bed' => $ref_end5_reg_bed,
			'ref_pileup_bed' => $ref_end5_pileup_bed,
			'conf_add_ref' => $conf_end5_add_ref,
			'ref_in_conf_merge_reg_bed' => $ref_in_conf_merge_end5_reg_bed,
			'end_prefix' => 'F',
			'confident_end_bed' => $confident_end5_bed,
			'signal_bed_bgz' => $signal_end5_bed_bgz,
			'end_info_hsh_ref' => $end5_info_hsh_ref,
		},
		'end3' => {
			'conf_bed_bgz' => $conf_end3_bed_bgz,
			'conf_merge_flank' => $conf_end3_merge_flank,
			'split_merge_in_conf_bed' => $split_merge_in_conf_end3_bed,
			'ref_reg_bed' => $ref_end3_reg_bed,
			'ref_pileup_bed' => $ref_end3_pileup_bed,
			'conf_add_ref' => $conf_end3_add_ref,
			'ref_in_conf_merge_reg_bed' => $ref_in_conf_merge_end3_reg_bed,
			'end_prefix' => 'T',
			'confident_end_bed' => $confident_end3_bed,
			'signal_bed_bgz' => $signal_end3_bed_bgz,
			'end_info_hsh_ref' => $end3_info_hsh_ref,
		},
	};
	
	foreach my $end (sort keys %{$end_hsh_ref}) {
		&reportAndLogStatus("$report_tag Defining $end regions", 10, "\n");#->3614
		my $conf_bed_bgz = $end_hsh_ref->{$end}{'conf_bed_bgz'};
		my $conf_merge_flank = $end_hsh_ref->{$end}{'conf_merge_flank'};
		my $split_merge_in_conf_bed = $end_hsh_ref->{$end}{'split_merge_in_conf_bed'};
		my $ref_reg_bed = $end_hsh_ref->{$end}{'ref_reg_bed'};
		my $ref_pileup_bed = $end_hsh_ref->{$end}{'ref_pileup_bed'};
		my $conf_add_ref = $end_hsh_ref->{$end}{'conf_add_ref'};
		my $ref_in_conf_merge_reg_bed = $end_hsh_ref->{$end}{'ref_in_conf_merge_reg_bed'};
		my $end_prefix = $end_hsh_ref->{$end}{'end_prefix'};
		my $confident_end_bed = $end_hsh_ref->{$end}{'confident_end_bed'};
		my $signal_bed_bgz = $end_hsh_ref->{$end}{'signal_bed_bgz'};
		my $end_info_hsh_ref = $end_hsh_ref->{$end}{'end_info_hsh_ref'};
		&mergeInConfEndRegion($tabix_bin, $conf_bed_bgz, $bedtools_bin, $conf_merge_flank, $size, $min_size_split, $min_frac_split, $min_summit_dist_split, $chrom, $split_merge_in_conf_bed);#->1932
		&defineRefEndRegBed($end, $conf_add_ref, $ref_reg_bed, $bedtools_bin, $size, $chrom, $ref_bed_bgz, $tabix_bin, $conf_merge_flank, $ref_pileup_bed);#->1320
		&mergeRefRegionAndInConfEndRegion($split_merge_in_conf_bed, $ref_reg_bed, $bedtools_bin, $ref_in_conf_merge_reg_bed, $end_prefix, $chrom_num);#->2082
		&defineEndSummit($bedtools_bin, $bgzip_bin, $ref_in_conf_merge_reg_bed, $signal_bed_bgz, $ref_pileup_bed, $conf_bed_bgz, $confident_end_bed, $chrom, $end_info_hsh_ref, $use_ref_only_end_pos);#->1073
		my $count = keys %{$end_info_hsh_ref};
		my $log_str = "confident $end region num";
		push @{$log_ary_ref}, [$log_str, $count];
	}

	return ($end5_info_hsh_ref, $end3_info_hsh_ref);
}
sub defineEndSummit {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: defineEndRegion|994
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $bedtools_bin, $bgzip_bin, $chrom, $conf_bed_bgz, $confident_end_bed, $end_info_hsh_ref, $ref_in_conf_merge_reg_bed, $ref_pileup_bed, $signal_bed_bgz, $use_ref_only_end_pos
#	output: 
#	toCall: &defineEndSummit($bedtools_bin, $bgzip_bin, $ref_in_conf_merge_reg_bed, $signal_bed_bgz, $ref_pileup_bed, $conf_bed_bgz, $confident_end_bed, $chrom, $end_info_hsh_ref, $use_ref_only_end_pos);
#	calledInLine: 1064
#....................................................................................................................................................#
	my ($bedtools_bin, $bgzip_bin, $ref_in_conf_merge_reg_bed, $signal_bed_bgz, $ref_pileup_bed, $conf_bed_bgz, $confident_end_bed, $chrom, $end_info_hsh_ref, $use_ref_only_end_pos) = @_;
	
	open MERGEREG, "<", $ref_in_conf_merge_reg_bed;
	while (<MERGEREG>) {
		chomp;
		my (undef, $chromStart, $chromEnd, $end_ID, $ref_only_digit, $strand) = split /\t/; 
		$end_info_hsh_ref->{$end_ID}{'CS'} = $chromStart;
		$end_info_hsh_ref->{$end_ID}{'CE'} = $chromEnd;
		$end_info_hsh_ref->{$end_ID}{'ST'} = $strand;
		$end_info_hsh_ref->{$end_ID}{'score'} = 0;
		$end_info_hsh_ref->{$end_ID}{'in_conf'} = 'N';
		my $ref_only = 'N';
		$ref_only = 'Y' if $ref_only_digit == 1;
		$end_info_hsh_ref->{$end_ID}{'ref_only'} = $ref_only;
	}
	close MERGEREG;
	
	open BEDTOOLS, "$bedtools_bin intersect -wo -s -sorted -a $ref_in_conf_merge_reg_bed -b $signal_bed_bgz | cut -f 4 |";
	while (<BEDTOOLS>) {
		chomp;
		my ($end_ID) = split /\t/;
		$end_info_hsh_ref->{$end_ID}{'in_conf'} = 'Y';
	}
	close BEDTOOLS;

	open BEDTOOLS, "$bedtools_bin intersect -wo -s -sorted -a $ref_in_conf_merge_reg_bed -b $signal_bed_bgz | cut -f 4,9,11 |";
	while (<BEDTOOLS>) {
		chomp;
		my ($end_ID, $pos, $signal_count) = split /\t/;
		my $ref_only = $end_info_hsh_ref->{$end_ID}{'ref_only'};
		$end_info_hsh_ref->{$end_ID}{'signal_count'}{$signal_count}{$pos}++;
		$end_info_hsh_ref->{$end_ID}{'score'} += $signal_count;
		push @{$end_info_hsh_ref->{$end_ID}{'pos'}}, $pos if $ref_only eq 'N' and $use_ref_only_end_pos eq 'yes'; #---[2023/08/21 21:20] use ref only pos if the end if ref only
	}
	close BEDTOOLS;

	open BEDTOOLS, "$bedtools_bin intersect -wo -s -sorted -a $ref_in_conf_merge_reg_bed -b $ref_pileup_bed | cut -f 4,9,11 |";
	while (<BEDTOOLS>) {
		chomp;
		my ($end_ID, $pos, $ref_count) = split /\t/;
		$end_info_hsh_ref->{$end_ID}{'ref_count'}{$ref_count}{$pos}++;
		push @{$end_info_hsh_ref->{$end_ID}{'pos'}}, $pos;
	}
	close BEDTOOLS;

	open ENDBED, "| sort -k1,1 -k2,2n -k6,6 | $bgzip_bin -c >$confident_end_bed";
	foreach my $end_ID (keys %{$end_info_hsh_ref}) {
		my $score = $end_info_hsh_ref->{$end_ID}{'score'};
		my $strand = $end_info_hsh_ref->{$end_ID}{'ST'};
		my $chromStart = $end_info_hsh_ref->{$end_ID}{'CS'};
		my $chromEnd = $end_info_hsh_ref->{$end_ID}{'CE'};
		my @signal_pos_ary = sort @{$end_info_hsh_ref->{$end_ID}{'pos'}};
		my $signalStart = $signal_pos_ary[0];
		my $signalEnd = $signal_pos_ary[-1];
		$end_info_hsh_ref->{$end_ID}{'OS'} = $signalStart;
		$end_info_hsh_ref->{$end_ID}{'OE'} = $signalEnd;

		my $itemRgb = '228,26,28';
		$itemRgb = '55,126,184' if ($strand eq '-');
		
		if ($end_info_hsh_ref->{$end_ID}{'in_conf'} eq 'Y') {
			die "signal_count does not exists in end_ID $end_ID\n" if not exists $end_info_hsh_ref->{$end_ID}{'signal_count'}; 
			foreach my $signal_count (sort {$b <=> $a} keys %{$end_info_hsh_ref->{$end_ID}{'signal_count'}}) {
				my @pos_ary = sort {$a <=> $b} keys %{$end_info_hsh_ref->{$end_ID}{'signal_count'}{$signal_count}};
				my $summit = $pos_ary[0];
				$summit = $pos_ary[-1] if ($strand eq '-');
				$end_info_hsh_ref->{$end_ID}{'summit'} = $summit;
				last;
			}
		} else {
			die "ref_count does not exists in end_ID $end_ID\n" if not exists $end_info_hsh_ref->{$end_ID}{'ref_count'}; 
			foreach my $ref_count (sort {$b <=> $a} keys %{$end_info_hsh_ref->{$end_ID}{'ref_count'}}) {
				my @pos_ary = sort {$a <=> $b} keys %{$end_info_hsh_ref->{$end_ID}{'ref_count'}{$ref_count}};
				my $summit = $pos_ary[0];
				$summit = $pos_ary[-1] if ($strand eq '-');
				$end_info_hsh_ref->{$end_ID}{'summit'} = $summit;
				last;
			}
		}
		
		my $thickEnd = $end_info_hsh_ref->{$end_ID}{'summit'};
		my $thickStart = $thickEnd - 1;
		my $blockCount = 1;
		my $blockSizes = $chromEnd - $chromStart;
		my $blockStarts = 0;
		print ENDBED join "", (join "\t", ($chrom, $chromStart, $chromEnd, $end_ID, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts)), "\n";

		delete $end_info_hsh_ref->{$end_ID}{'signal_count'};
		delete $end_info_hsh_ref->{$end_ID}{'ref_count'};
		delete $end_info_hsh_ref->{$end_ID}{'pos'};
	}
	close ENDBED;

	return ();
}
sub definePoolOutputPath {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 3_process|197
#	secondaryAppearInSection: >none
#	input: $paramTag, $result_bed_dir, $result_log_dir
#	output: $output_path_hsh_ref
#	toCall: my ($output_path_hsh_ref) = &definePoolOutputPath($result_bed_dir, $result_log_dir, $paramTag);
#	calledInLine: 204
#....................................................................................................................................................#
	my ($result_bed_dir, $result_log_dir, $paramTag) = @_;
	
	my $output_path_hsh_ref = {
		'junct_bed' => {
			'size_group' => 4,
			'type' => 'bed',
			'path' => "$result_bed_dir/$paramTag.junct.bed.bgz",
		},
		'set_bed' => {
			'size_group' => 5,
			'type' => 'bed',
			'path' => "$result_bed_dir/$paramTag.set.bed.bgz",
		},
		'model_bed' => {
			'size_group' => 5,
			'type' => 'bed',
			'path' => "$result_bed_dir/$paramTag.model.bed.bgz",
		},
		'all_end5_bed' => {
			'size_group' => 0,
			'type' => 'bed',
			'path' => "$result_bed_dir/$paramTag.end5.bed.bgz",
		},
		'all_end3_bed' => {
			'size_group' => 0,
			'type' => 'bed',
			'path' => "$result_bed_dir/$paramTag.end3.bed.bgz",
		},
		'ref_trnscpt_bed' => {
			'size_group' => 0,
			'type' => 'bed',
			'path' => "$result_bed_dir/$paramTag.trnscpt.ref.bed.bgz",
		},
		'qry_trnscpt_bed' => {
			'size_group' => 9,
			'type' => 'bed',
			'path' => "$result_bed_dir/$paramTag.trnscpt.qry.bed.bgz",
		},
		'discard_trnscpt_bed' => {
			'size_group' => 0,
			'type' => 'bed',
			'path' => "$result_bed_dir/$paramTag.trnscpt.discard.bed.bgz",
		},
		'trnscpt_info_tsv' => {
			'size_group' => 9,
			'type' => 'info',
			'path' => "$result_log_dir/$paramTag.trnscpt.info.tsv.gz",
			'header' => ['trnscpt_ID', 'ref_qry', 'completeness', 'set_ID', 'end5_ID', 'end3_ID', 'model_ID_str', 'chrom', 'trnscpt_chromStart', 'trnscpt_chromEnd', 'set_chromStart', 'set_chromEnd', 'strand', 'bound_str'],
		},
		'model_info_tsv' => {
			'size_group' => 4,
			'type' => 'info',
			'path' => "$result_log_dir/$paramTag.model.info.tsv.gz",
			'header' => ['model_ID', 'loc', 'completeness', 'full_set_ID', 'end5_conf', 'end3_conf', 'junct_conf', 'end5_endtype', 'end3_endtype', 'novelty', 'strand', 'partial_set_count', 'partial_set_ID_str', 'full_ref_count', 'full_qry_count', 'partial_ref_count', 'partial_qry_count', 'full_set_bound_str'],
		},
		'set_info_tsv' => {
			'size_group' => 7,
			'type' => 'info',
			'path' => "$result_log_dir/$paramTag.set.info.tsv.gz",
			'header' => ['set_ID', 'loc', 'completeness', 'represent_model', 'end5_conf', 'end3_conf', 'junct_conf', 'end5_endtype', 'end3_endtype', 'ref_count', 'qry_count', 'model_count', 'model_ID_str', 'strand', 'bound_str', 'ref_trnscrptID', 'qry_trnscrptID'],
		},
		'end5_info_tsv' => {
			'size_group' => 1,
			'type' => 'info',
			'path' => "$result_log_dir/$paramTag.end5.info.tsv.gz",
			'header' => ['endID', 'confidence', 'ref_only', 'ref_count', 'qry_count', 'set_count', 'model_count', 'chrom', 'summit', 'strand', 'chromStart', 'chromEnd', 'obsStart', 'obsEnd', 'set_str', 'model_str'],
		},
		'end3_info_tsv' => {
			'size_group' => 1,
			'type' => 'info',
			'path' => "$result_log_dir/$paramTag.end3.info.tsv.gz",
			'header' => ['endID', 'confidence', 'ref_only', 'ref_count', 'qry_count', 'set_count', 'model_count', 'chrom', 'summit', 'strand', 'chromStart', 'chromEnd', 'obsStart', 'obsEnd', 'set_str', 'model_str'],
		},
		'junct_info_tsv' => {
			'size_group' => 5,
			'type' => 'info',
			'path' => "$result_log_dir/$paramTag.junct.info.tsv.gz",
			'header' => ['junct_ID', 'canonical', 'splicing_site', 'chrom', 'intron_start', 'intron_end', 'strand', 'ref_count', 'qry_count'],
		},
		'model_set_pair_info_tsv' => {
			'size_group' => 3,
			'type' => 'info',
			'path' => "$result_log_dir/$paramTag.model_set_pair.info.tsv.gz",
			'header' => ['model_ID', 'set_ID', 'loc', 'end5_endtype', 'end3_endtype', 'full_partial', 'junct_num', 'length', 'ref_count', 'qry_count', 'bound_str'],
		},
	};

	return ($output_path_hsh_ref);
}
sub defineQryBedAry {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: >none
#	primaryAppearInSection: 3_process|197
#	secondaryAppearInSection: >none
#	input: $qry_bed_bgz
#	output: $qry_bed_ary_ref
#	toCall: my ($qry_bed_ary_ref) = &defineQryBedAry($qry_bed_bgz);
#	calledInLine: 200
#....................................................................................................................................................#
	my ($qry_bed_bgz) = @_;
	
	my $qry_bed_ary_ref = [];
	if ($qry_bed_bgz =~ m/\.bed\.bgz$/) {
		push @{$qry_bed_ary_ref}, $qry_bed_bgz;
	} elsif ($qry_bed_bgz =~ m/\.txt$/) {
		open QRY, "<", $qry_bed_bgz;
		while (<QRY>) {
			chomp;
			next if $_ =~ m/^#/;
			my ($indiv_qry_bed_bgz) = split /\t/;
			if ($indiv_qry_bed_bgz =~ m/\.bed\.bgz$/) {
				push @{$qry_bed_ary_ref}, $indiv_qry_bed_bgz;
			} else {
				die "Quitting: indiv_qry_bed_bgz $indiv_qry_bed_bgz must be in *.bed.bgz format\n";
			}
		}
		close QRY;
	} else {
		die "Quitting: qry_bed_bgz $qry_bed_bgz must be in *.bed.bgz format or *.txt format\n";
	}
	
	my $num_qry = @{$qry_bed_ary_ref};
	&reportAndLogStatus("$num_qry query bed store", 10, "\n");#->3614

	return ($qry_bed_ary_ref);
}
sub defineRefEndRegBed {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: defineEndRegion|994
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $bedtools_bin, $chrom, $conf_add_ref, $conf_merge_flank, $end, $ref_bed_bgz, $ref_pileup_bed, $ref_reg_bed, $size, $tabix_bin
#	output: 
#	toCall: &defineRefEndRegBed($end, $conf_add_ref, $ref_reg_bed, $bedtools_bin, $size, $chrom, $ref_bed_bgz, $tabix_bin, $conf_merge_flank, $ref_pileup_bed);
#	calledInLine: 1062
#....................................................................................................................................................#
	my ($end, $conf_add_ref, $ref_reg_bed, $bedtools_bin, $size, $chrom, $ref_bed_bgz, $tabix_bin, $conf_merge_flank, $ref_pileup_bed) = @_;
	
	system ("touch $ref_reg_bed");
	my $tabix_cmd  = "$tabix_bin $ref_bed_bgz $chrom | cut -f 2,3,4,6 |";
	my $region_bedtools_merge_cmd = "| sort -k2,2n -k6,6 | $bedtools_bin merge -d -1 -s -c 4,5,6 -o sum,sum,distinct | awk 'BEGIN {FS=OFS=\"\\t\"} { \$4=\"R\"NR; print }' >$ref_reg_bed";
	my $pileup_bedtools_merge_cmd = "| sort -k2,2n -k6,6 | $bedtools_bin merge -d -1 -s -c 4,5,6 -o sum,sum,distinct >$ref_pileup_bed";
	open REGIONBEDTOOLS, $region_bedtools_merge_cmd;
	open PILEUPBEDTOOLS, $pileup_bedtools_merge_cmd;
	open TABIX, $tabix_cmd;
	while (<TABIX>) {
		chomp;
		my ($chromStart, $chromEnd, $trnscpt_ID, $strand) = split /\t/;
		
		my $end5_pos_end;
		my $end3_pos_end;
		if ($strand eq '+') {
			$end5_pos_end = $chromStart+1;
			$end3_pos_end = $chromEnd;
		} elsif ($strand eq '-') {
			$end5_pos_end = $chromEnd;
			$end3_pos_end = $chromStart+1;
		} else {
			die;
		}
		
		my $pos_end;
		if ($end eq 'end5') {
			$pos_end = $end5_pos_end;
		} elsif ($end eq 'end3') {
			$pos_end = $end3_pos_end;
		}
		my $pos_start = $pos_end - 1;
		my $regStart = $pos_start - $conf_merge_flank;
		my $regEnd = $pos_end + $conf_merge_flank;
		
		$regEnd = $size if $regEnd > $size;
		$regStart = 1 if $regStart < 1;
		
		print PILEUPBEDTOOLS join "", (join "\t", ($chrom, $pos_start, $pos_end, 1, 1, $strand)), "\n";
		
		if ($conf_add_ref eq 'yes') {
			print REGIONBEDTOOLS join "", (join "\t", ($chrom, $regStart, $regEnd, 1, 1, $strand)), "\n";
		}
	}
	close TABIX;
	close PILEUPBEDTOOLS;
	close REGIONBEDTOOLS;

	return ();
}
sub defineTrnscptModels {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: assignEnd3CompleteSetToEnd3CompleteModels|447, assignEnd3IncompleteSetToEnd3CompleteModels|492, assignEnd3IncompleteSetToEnd3IncompleteModels|538, assignOrphanSetToExistingModels|668, reportAndLogStatus|3614
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $end5_info_hsh_ref, $inital_novel_model_num, $novel_model_chrom_prefix, $report_tag, $set_info_hsh_ref
#	output: $model_info_hsh_ref
#	toCall: my ($model_info_hsh_ref) = &defineTrnscptModels($end5_info_hsh_ref, $set_info_hsh_ref, $novel_model_chrom_prefix, $inital_novel_model_num, $report_tag);
#	calledInLine: 3342
#....................................................................................................................................................#
	my ($end5_info_hsh_ref, $set_info_hsh_ref, $novel_model_chrom_prefix, $inital_novel_model_num, $report_tag) = @_;
	
	my $model_info_hsh_ref = {};
	my $novel_model_num = $inital_novel_model_num;
	my $num_proc = 0;
	foreach my $end5_ID (sort keys %{$end5_info_hsh_ref}) {
		my $within_end5_model_hsh_ref = {};
		my $strand = $end5_info_hsh_ref->{$end5_ID}{'ST'};
		$num_proc++;
		&reportAndLogStatus("$report_tag Defining models >> $num_proc 5'end regions processed", 10, "\n") if $num_proc%1000 == 0;#->3614

		#---[2023/07/21 14:41] assign complete models
		($novel_model_num) = &assignEnd3CompleteSetToEnd3CompleteModels($within_end5_model_hsh_ref, $set_info_hsh_ref, $end5_info_hsh_ref, $novel_model_chrom_prefix, $novel_model_num, $end5_ID);#->447

		#---[2023/07/21 14:41] assign incomplete transcripts to complete models
		&assignEnd3IncompleteSetToEnd3CompleteModels($strand, $set_info_hsh_ref, $within_end5_model_hsh_ref, $end5_info_hsh_ref, $end5_ID);#->492
		
		#---[2023/07/21 14:41] to identify the incomplete models
		($novel_model_num) = &assignEnd3IncompleteSetToEnd3IncompleteModels($within_end5_model_hsh_ref, $set_info_hsh_ref, $end5_info_hsh_ref, $novel_model_chrom_prefix, $novel_model_num, $strand, $end5_ID);#->538

		#---[2023/07/21 14:41] to assign the sets to aleady defined models, only for those what are the longest and not compatible to the longest set defined in the previous step but compatible to some other sets in the model
		&assignOrphanSetToExistingModels($within_end5_model_hsh_ref, $set_info_hsh_ref, $end5_info_hsh_ref, $novel_model_chrom_prefix, $novel_model_num, $strand, $end5_ID);#->668

		foreach my $model_ID (keys %{$within_end5_model_hsh_ref}) {
			$model_info_hsh_ref->{$model_ID}{'CM'} = $within_end5_model_hsh_ref->{$model_ID}{'CM'}; #---[2023/07/21 14:46] complete
			$model_info_hsh_ref->{$model_ID}{'set'} = $within_end5_model_hsh_ref->{$model_ID}{'set'}; #---[2023/07/21 14:46] principle cound
			$model_info_hsh_ref->{$model_ID}{'N'} = $within_end5_model_hsh_ref->{$model_ID}{'N'}; #---[2023/07/21 14:46] novelty
			foreach my $set_ID (keys %{$within_end5_model_hsh_ref->{$model_ID}{'IS'}}) {
				$model_info_hsh_ref->{$model_ID}{'IS'}{$set_ID} = $within_end5_model_hsh_ref->{$model_ID}{'IS'}{$set_ID}; #---[2023/07/21 14:46] incomplete set
			}
		}
	}

	&reportAndLogStatus("$report_tag Defining models >> $num_proc 5'end regions processed", 10, "\n");#->3614

	return ($model_info_hsh_ref);
}
sub defineTrnscptSets {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: convertBoundToBed|777, getEndNum|1783, reportAndLogStatus|3614, retainNoQryRefBoundSet|3636
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $bgzip_bin, $chrom, $chrom_num, $doubtful_end_avoid_summit, $end3_info_hsh_ref, $end5_info_hsh_ref, $junct_info_hsh_ref, $min_transcript_length, $report_tag, $retain_no_qry_ref_bound_set, $set_bed, $trnscpt_info_hsh_ref, $trnscpt_set_end_priority_hsh_ref
#	output: $set_info_hsh_ref
#	toCall: my ($set_info_hsh_ref) = &defineTrnscptSets($trnscpt_info_hsh_ref, $junct_info_hsh_ref, $chrom_num, $end5_info_hsh_ref, $end3_info_hsh_ref, $report_tag, $trnscpt_set_end_priority_hsh_ref, $set_bed, $min_transcript_length, $chrom, $bgzip_bin, $doubtful_end_avoid_summit, $retain_no_qry_ref_bound_set);
#	calledInLine: 3331
#....................................................................................................................................................#
	my ($trnscpt_info_hsh_ref, $junct_info_hsh_ref, $chrom_num, $end5_info_hsh_ref, $end3_info_hsh_ref, $report_tag, $trnscpt_set_end_priority_hsh_ref, $set_bed, $min_transcript_length, $chrom, $bgzip_bin, $doubtful_end_avoid_summit, $retain_no_qry_ref_bound_set) = @_;
	
	my $set_info_hsh_ref = {};
	my $tmp_hsh_ref = {};
	my $num_proc = 0;
	my $test_endtype_ary_ref = [sort {$trnscpt_set_end_priority_hsh_ref->{$a} <=> $trnscpt_set_end_priority_hsh_ref->{$b}} keys %{$trnscpt_set_end_priority_hsh_ref}];
	my $endtype_prefix_hsh_ref = {
		'commonest' => 'C',
		'longest' => 'L',
		'summit' => 'S',
		'original' => 'O',
	};
	
	foreach my $ref_qry (sort keys %{$trnscpt_info_hsh_ref}) {
		foreach my $trnscpt_ID (sort keys %{$trnscpt_info_hsh_ref->{$ref_qry}}) {
			$num_proc++;
			&reportAndLogStatus("$report_tag Defining sets >> $num_proc transcripts processed", 10, "\n") if $num_proc%10000 == 0;#->3614

			my ($strand, $chromStart, $chromEnd) = @{$trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'I'}};
			my $end5_ID = $trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'F'};
			my $end3_ID = $trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'T'};
			my @junct_ary = @{$trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'J'}};
			my @bound_str_ary = ();
			if ($strand eq '+') {
				@bound_str_ary = ($end5_ID, @junct_ary, $end3_ID, $strand);
			} elsif ($strand eq '-') {
				@bound_str_ary = ($strand, $end3_ID, @junct_ary, $end5_ID);
				@bound_str_ary = reverse @bound_str_ary;
			} else {
				die;
			}
			
			my $bound_str = join "_", @bound_str_ary;
			if (not exists $tmp_hsh_ref->{$bound_str}) {
				$tmp_hsh_ref->{$bound_str}++;;
				my $set_num = keys %{$tmp_hsh_ref};
				my $set_ID = "S".$chrom_num.$set_num;
				$tmp_hsh_ref->{$bound_str} = $set_ID;
				$set_info_hsh_ref->{$set_ID}{'BS'} = $bound_str; #---[2023/07/21 12:34] set string
				$set_info_hsh_ref->{$set_ID}{'ST'} = $strand; #---[2023/07/21 12:34] set string
				$set_info_hsh_ref->{$set_ID}{'AM'} = undef;
				$set_info_hsh_ref->{$set_ID}{'RR'} = 'N'; #---[2023/08/21 23:43] retain reference bounds
				if ($end3_ID =~ m/^X/ or $end5_ID =~ m/^X/) {
					$set_info_hsh_ref->{$set_ID}{'CM'} = 'N';
				} else {
					$set_info_hsh_ref->{$set_ID}{'CM'} = 'Y';
				}
			}
			
			my $set_ID = $tmp_hsh_ref->{$bound_str};
			$set_info_hsh_ref->{$set_ID}{'CS_hsh'}{$chromStart}++; #---[2023/07/21 12:34] complete
			$set_info_hsh_ref->{$set_ID}{'CE_hsh'}{$chromEnd}++; #---[2023/07/21 12:34] complete
			$trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'S'} = $set_ID;
			$end5_info_hsh_ref->{$end5_ID}{'set'}{$set_ID}++;
			$end3_info_hsh_ref->{$end3_ID}{'set'}{$set_ID}++;
			
			$set_info_hsh_ref->{$set_ID}{'T'}{$ref_qry}{$trnscpt_ID}++; #---[2023/07/21 12:57] transcript 
			undef $trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'J'};
		}
	}
	&reportAndLogStatus("$report_tag Finished defining sets ## $num_proc transcripts processed", 10, "\n");#->3614
	
	if ($retain_no_qry_ref_bound_set eq 'yes') {
		&reportAndLogStatus("$report_tag Start retaining ref set bounds without qry support", 10, "\n");#->3614
		my $set_num = keys %{$tmp_hsh_ref};
		my ($end3_num, $end5_num) = &getEndNum($end3_info_hsh_ref, $end5_info_hsh_ref, $chrom_num);#->1783
		&retainNoQryRefBoundSet($report_tag, $set_info_hsh_ref, $trnscpt_info_hsh_ref, $end5_info_hsh_ref, $end3_info_hsh_ref, $end3_num, $end5_num, $set_num, $chrom_num);#->3636
	}
	
	#---[2023/07/21 15:49] summarizing farthest first exon and tail position
	$num_proc = 0;
	open SETBED, "| sort -k2,2n -k6,6 | $bgzip_bin -c >$set_bed";
	foreach my $set_ID (keys %{$set_info_hsh_ref}) {
		$num_proc++;
		&reportAndLogStatus("$report_tag Printing sets >> $num_proc sets printed", 10, "\n") if $num_proc%10000 == 0;#->3614
		my $strand = $set_info_hsh_ref->{$set_ID}{'ST'};
		my $RR = $set_info_hsh_ref->{$set_ID}{'RR'};
		
		my $bound_str_ary_ref = [];
		@{$bound_str_ary_ref} = split /_/, $set_info_hsh_ref->{$set_ID}{'BS'}; #---[2023/07/21 12:34] set string
		my $end5_ID = $bound_str_ary_ref->[0];
		my $end3_ID = $bound_str_ary_ref->[-2];
		my $chromStartEnd_hsh_ref = {};
		
		if ($RR eq 'Y') {
			$chromStartEnd_hsh_ref->{'CS'}{'original'} = $set_info_hsh_ref->{$set_ID}{'OR'}[0];
			$chromStartEnd_hsh_ref->{'CE'}{'original'} = $set_info_hsh_ref->{$set_ID}{'OR'}[1];
		}
		
		if ($strand eq '+') {
			$chromStartEnd_hsh_ref->{'CS'}{'summit'} = $end5_info_hsh_ref->{$end5_ID}{'summit'}-1;
			$chromStartEnd_hsh_ref->{'CE'}{'summit'} = $end3_info_hsh_ref->{$end3_ID}{'summit'};
		} else {
			$chromStartEnd_hsh_ref->{'CS'}{'summit'} = $end3_info_hsh_ref->{$end3_ID}{'summit'}-1;
			$chromStartEnd_hsh_ref->{'CE'}{'summit'} = $end5_info_hsh_ref->{$end5_ID}{'summit'};
		}
		
		($chromStartEnd_hsh_ref->{'CS'}{'longest'}) = sort {$a <=> $b} keys %{$set_info_hsh_ref->{$set_ID}{'CS_hsh'}};
		($chromStartEnd_hsh_ref->{'CE'}{'longest'}) = sort {$b <=> $a} keys %{$set_info_hsh_ref->{$set_ID}{'CE_hsh'}};

		($chromStartEnd_hsh_ref->{'CS'}{'commonest'}) = sort {$set_info_hsh_ref->{$set_ID}{'CS_hsh'}{$b} <=> $set_info_hsh_ref->{$set_ID}{'CS_hsh'}{$a}} keys %{$set_info_hsh_ref->{$set_ID}{'CS_hsh'}};
		($chromStartEnd_hsh_ref->{'CE'}{'commonest'}) = sort {$set_info_hsh_ref->{$set_ID}{'CE_hsh'}{$b} <=> $set_info_hsh_ref->{$set_ID}{'CE_hsh'}{$a}} keys %{$set_info_hsh_ref->{$set_ID}{'CE_hsh'}};
		
		my $qry_count = 0;
		$qry_count = keys %{$set_info_hsh_ref->{$set_ID}{'T'}{'qry'}} if exists $set_info_hsh_ref->{$set_ID}{'T'}{'qry'};
		my ($bed_line_ary_ref, $CS_end_type, $CE_end_type, $exon1_start, $exon1_end) = &convertBoundToBed($junct_info_hsh_ref, $bound_str_ary_ref, $chromStartEnd_hsh_ref, $set_ID, $min_transcript_length, $test_endtype_ary_ref, $qry_count, $chrom, $doubtful_end_avoid_summit, $RR);#->777
		print SETBED join "", (join "\t", (@{$bed_line_ary_ref})), "\n";

		$set_info_hsh_ref->{$set_ID}{'E1'} = [$exon1_start, $exon1_end];
		$set_info_hsh_ref->{$set_ID}{'CS'} = $chromStartEnd_hsh_ref->{'CS'}{$CS_end_type};
		$set_info_hsh_ref->{$set_ID}{'CE'} = $chromStartEnd_hsh_ref->{'CE'}{$CE_end_type};
		
		if ($strand eq '+') {
			$set_info_hsh_ref->{$set_ID}{'T5'} = $endtype_prefix_hsh_ref->{$CS_end_type};
			$set_info_hsh_ref->{$set_ID}{'T3'} = $endtype_prefix_hsh_ref->{$CE_end_type};
		} else {
			$set_info_hsh_ref->{$set_ID}{'T5'} = $endtype_prefix_hsh_ref->{$CE_end_type};
			$set_info_hsh_ref->{$set_ID}{'T3'} = $endtype_prefix_hsh_ref->{$CS_end_type};
		}

		delete $set_info_hsh_ref->{$set_ID}{'CS_hsh'};
		delete $set_info_hsh_ref->{$set_ID}{'CE_hsh'};
	}
	close SETBED;
	&reportAndLogStatus("$report_tag Finised printing sets >> $num_proc sets printed", 10, "\n");#->3614
	
	return ($set_info_hsh_ref);
}
sub extractTranscriptBounds {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $bgzip_bin, $chrom, $chrom_num, $conf_junct_hsh_ref, $discard_trnscpt_bed, $input_bed_hsh_ref, $junct_seq_hsh_ref, $log_ary_ref, $min_exon_length, $min_transcript_length, $report_tag, $tabix_bin, $threadID
#	output: $junct_info_hsh_ref, $trnscpt_info_hsh_ref
#	toCall: my ($junct_info_hsh_ref, $trnscpt_info_hsh_ref) = &extractTranscriptBounds($tabix_bin, $chrom_num, $threadID, $chrom, $input_bed_hsh_ref, $report_tag, $min_transcript_length, $min_exon_length, $log_ary_ref, $discard_trnscpt_bed, $bgzip_bin, $conf_junct_hsh_ref, $junct_seq_hsh_ref);
#	calledInLine: 3321
#....................................................................................................................................................#
	my ($tabix_bin, $chrom_num, $threadID, $chrom, $input_bed_hsh_ref, $report_tag, $min_transcript_length, $min_exon_length, $log_ary_ref, $discard_trnscpt_bed, $bgzip_bin, $conf_junct_hsh_ref, $junct_seq_hsh_ref) = @_;

	my $junct_info_hsh_ref = {};
	my $trnscpt_info_hsh_ref = {};
	my $tmp_hsh_ref = {};

	my $count_hsh_ref = {};
	foreach my $ref_qry (qw/ref qry/) {
		foreach my $remove_retain (qw/removed retained/) {
			$count_hsh_ref->{$ref_qry}{$remove_retain} = 0;
		}
	}
	
	open DISCARDTRNSCPT, "| sort -k2,2n -k6,6 | $bgzip_bin -c >$discard_trnscpt_bed";
	foreach my $input_ID (keys %{$input_bed_hsh_ref}) {
		my $removed_num = 0;
		&reportAndLogStatus("$report_tag Extracting transcript bounds in $input_ID", 10, "\n");#->3614
		my $ref_qry = $input_bed_hsh_ref->{$input_ID}{'ref_qry'};
		my $suffix = $input_bed_hsh_ref->{$input_ID}{'suffix'};
		my $trnscpt_bed = $input_bed_hsh_ref->{$input_ID}{'trnscpt_bed'};
		my $end5_bed = $input_bed_hsh_ref->{$input_ID}{'end5_bed'};
		my $end3_bed = $input_bed_hsh_ref->{$input_ID}{'end3_bed'};
		my $num_proc = 0;
		my $tabix_cmd  = "$tabix_bin $trnscpt_bed $chrom |";
		open END5, "| sort -k2,2n -k6,6 | $bgzip_bin -c >$end5_bed";
		open END3, "| sort -k2,2n -k6,6 | $bgzip_bin -c >$end3_bed";
		open TABIX, "$tabix_cmd";
		while (<TABIX>) {
			chomp;
			my (undef, $chromStart, $chromEnd, $trnscpt_ID, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\t/;
			#$trnscpt_ID = $trnscpt_ID."-".$suffix if $ref_qry eq 'qry';
			$num_proc++;
			my $valid = 'yes';
			foreach my $exon_size (split /,/, $blockSizes) {
				$valid = 'no' if $exon_size < $min_exon_length;
			}
			
			$valid = 'no' if ($chromEnd-$chromStart < $min_transcript_length);
			
			if ($valid eq 'no') {
				$removed_num++;
				$count_hsh_ref->{$ref_qry}{'removed'}++;
				print DISCARDTRNSCPT join "", (join "\t", ($chrom, $chromStart, $chromEnd, $trnscpt_ID, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts)), "\n";
				next;
			}
			$count_hsh_ref->{$ref_qry}{'retained'}++;
			&reportAndLogStatus("$report_tag Extracting bounds >> $num_proc transcript bounds extracted in $input_ID", 10, "\n") if $num_proc%10000 == 0;#->3614
			my @rngAry = ();
			my @blockSizes_ary = split /,/, $blockSizes;
			my @blockStarts_ary = split /,/, $blockStarts;
			foreach my $i (0..$#blockStarts_ary) {
				my $offset = $blockStarts_ary[$i];
				my $exon_start = $chromStart+$offset;
				my $exon_end = $exon_start + $blockSizes_ary[$i];
				push @rngAry, ($exon_start, $exon_end);
			}
		
			$trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'J'} = [];
		
			for (my $i=1; $i < $#rngAry; $i += 2) {
				my ($intron_start, $intron_end) = ($rngAry[$i], $rngAry[$i+1]);
				my $junct_str = join ":", ($intron_start, $intron_end, $strand);
				if (not exists $tmp_hsh_ref->{$junct_str}) {
					$tmp_hsh_ref->{$junct_str}++;
					my $acceptor = $junct_seq_hsh_ref->{$junct_str}{'acceptor'};
					my $donor = $junct_seq_hsh_ref->{$junct_str}{'donor'};
					my $splicing_site = $donor."-".$acceptor;
					my $canonical = 'N';
					#GT/AG, GC/AG and AT/AC as canonical
					$canonical = 'Y' if $splicing_site eq 'GT-AG' or $splicing_site eq 'GC-AG' or $splicing_site eq 'AT-AC';
					my $junct_num = keys %{$tmp_hsh_ref};
					my $junct_ID = "J".$chrom_num.$junct_num;
					$junct_ID = 'X'.$junct_ID if not exists $conf_junct_hsh_ref->{$junct_str} and $canonical eq 'N';
					$junct_info_hsh_ref->{$junct_ID}{'JS'} = $junct_str; #---[2023/07/21 12:34] junction string
					$junct_info_hsh_ref->{$junct_ID}{'SS'} = $splicing_site; #---[2023/07/21 12:34] junction string
					$junct_info_hsh_ref->{$junct_ID}{'CN'} = $canonical; #---[2023/07/21 12:34] junction string
					$tmp_hsh_ref->{$junct_str} = $junct_ID;
				}
				my $junct_ID = $tmp_hsh_ref->{$junct_str};
				$junct_info_hsh_ref->{$junct_ID}{'CT'}{$ref_qry}++; #---[2023/07/21 12:34] count
				push @{$trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'J'}}, $junct_ID;
			}

			my $end5;
			my $end3;

			if ($strand eq '+') {
				$end5 = $chromStart + 1;
				$end3 = $chromEnd;
			
			} elsif ($strand eq '-') {
				$end5 = $chromEnd;
				$end3 = $chromStart + 1;

			} else {
				die "transcript $trnscpt_ID strand is $strand. Strand must be + or -\n";
			}
			print END5 join "", (join "\t", ($chrom, $end5-1, $end5, $trnscpt_ID, 1, $strand)), "\n";
			print END3 join "", (join "\t", ($chrom, $end3-1, $end3, $trnscpt_ID, 1, $strand)), "\n";

			$trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'T'} = undef; #---[2023/07/24 2:46] three prime as X, will be rewritten in later intersection
			$trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'F'} = undef; #---[2023/07/24 2:46] five prime as X, will be rewritten in later intersection
			$trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'I'} = [$strand, $chromStart, $chromEnd]; #---[2023/07/21 12:47] information
		
		}
		close END5;
		close END3;
		&reportAndLogStatus("$report_tag Finished extracting bounds ## $num_proc transcript bounds extracted in $input_ID", 10, "\n");#->3614
		&reportAndLogStatus("$report_tag Finished extracting bounds ## $removed_num transcript were removed by length $min_transcript_length in $input_ID", 10, "\n");#->3614
		
	}
	close DISCARDTRNSCPT;

	foreach my $ref_qry (qw/ref qry/) {
		foreach my $remove_retain (qw/removed retained/) {
			my $count = $count_hsh_ref->{$ref_qry}{$remove_retain};
			my $log_str = "$ref_qry transcript $remove_retain";
			push @{$log_ary_ref}, [$log_str, $count];
		}
	}

	foreach my $junct_ID (keys %{$junct_info_hsh_ref}) {
		foreach my $ref_qry (qw/ref qry/) {
			if (not exists $junct_info_hsh_ref->{$junct_ID}{'CT'}{$ref_qry}) {
				$junct_info_hsh_ref->{$junct_ID}{'CT'}{$ref_qry} = 0;
			}
		}
	}

	return ($junct_info_hsh_ref, $trnscpt_info_hsh_ref);
}
sub filterModelAndSetWithQry {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: >none
#	primaryAppearInSection: 3_process|197
#	secondaryAppearInSection: >none
#	input: $bgzip_bin, $min_output_qry_count, $output_path_hsh_ref, $paramTag, $result_bed_dir, $tabix_bin
#	output: 
#	toCall: &filterModelAndSetWithQry($bgzip_bin, $tabix_bin, $paramTag, $result_bed_dir, $output_path_hsh_ref, $min_output_qry_count);
#	calledInLine: 206
#....................................................................................................................................................#
	my ($bgzip_bin, $tabix_bin, $paramTag, $result_bed_dir, $output_path_hsh_ref, $min_output_qry_count) = @_;
	
	my $all_set_bed_path = $output_path_hsh_ref->{'set_bed'}{'path'};
	my $all_set_info_path = $output_path_hsh_ref->{'set_info_tsv'}{'path'};
	my $qry_set_bed_path = "$result_bed_dir/$paramTag.set.with_qry.bed.bgz";
	
	&reportAndLogStatus("Filtering set with qry_count >= $min_output_qry_count", 10, "\n");#->3614
	my $qry_set_hsh_ref = {};
	open SETINFO, "gzip -dc $all_set_info_path | cut -f 1,11 | ";
	<SETINFO>; #'set_ID', 'loc', 'completeness', 'represent_model', 'end5_conf', 'end3_conf', 'junct_conf', 'end5_endtype', 'end3_endtype', 'ref_count', 'qry_count', 'model_count', 'model_ID_str', 'strand', 'bound_str', 'ref_trnscrptID', 'qry_trnscrptID'
	while (<SETINFO>) {
		chomp;
		my ($set_ID, $qry_count) = split /\t/;
		$qry_set_hsh_ref->{$set_ID} = $qry_count if $qry_count >= $min_output_qry_count;
	}
	close SETINFO;
	
	open QRYSETBED, "| $bgzip_bin -c >$qry_set_bed_path";
	open ALLSETBED, "$bgzip_bin -dc $all_set_bed_path |";
	while (<ALLSETBED>) {
		chomp;
		my @bed_ary = split /\t/;
		print QRYSETBED join "", (join "\t", (@bed_ary)), "\n" if exists $qry_set_hsh_ref->{$bed_ary[3]};
	}
	close ALLSETBED;
	close QRYSETBED;
	system "$tabix_bin -p bed $qry_set_bed_path";

	my $all_model_bed_path = $output_path_hsh_ref->{'model_bed'}{'path'};
	my $all_model_info_path = $output_path_hsh_ref->{'model_info_tsv'}{'path'};
	my $qry_model_bed_path = "$result_bed_dir/$paramTag.model.with_qry.bed.bgz";
	
	&reportAndLogStatus("Filtering model with qry_count >= $min_output_qry_count", 10, "\n");#->3614
	my $qry_model_hsh_ref = {};
	open MODELINFO, "gzip -dc $all_model_info_path | cut -f 1,15,17 | ";
	<MODELINFO>; #'model_ID', 'loc', 'completeness', 'full_set_ID', 'end5_conf', 'end3_conf', 'junct_conf', 'end5_endtype', 'end3_endtype', 'novelty', 'strand', 'partial_set_count', 'partial_set_ID_str', 'full_ref_count', 'full_qry_count', 'partial_ref_count', 'partial_qry_count', 'full_set_bound_str'
	while (<MODELINFO>) {
		chomp;
		my ($model_ID, $full_qry_count, $partial_qry_count) = split /\t/;
		my $qry_count = $full_qry_count + $partial_qry_count;
		$qry_model_hsh_ref->{$model_ID} = $qry_count if $qry_count >= $min_output_qry_count;
	}
	close MODELINFO;
	
	open QRYMODELBED, "| $bgzip_bin -c >$qry_model_bed_path";
	open ALLMODELBED, "$bgzip_bin -dc $all_model_bed_path |";
	while (<ALLMODELBED>) {
		chomp;
		my @bed_ary = split /\t/;
		print QRYMODELBED join "", (join "\t", (@bed_ary)), "\n" if exists $qry_model_hsh_ref->{$bed_ary[3]};
	}
	close ALLMODELBED;
	close QRYMODELBED;
	system "$tabix_bin -p bed $qry_model_bed_path";

	return ();
}
sub getEndNum {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: defineTrnscptSets|1431
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $chrom_num, $end3_info_hsh_ref, $end5_info_hsh_ref
#	output: $end3_num, $end5_num
#	toCall: my ($end3_num, $end5_num) = &getEndNum($end3_info_hsh_ref, $end5_info_hsh_ref, $chrom_num);
#	calledInLine: 1507
#....................................................................................................................................................#
	my ($end3_info_hsh_ref, $end5_info_hsh_ref, $chrom_num) = @_;
	
	my $end3_num = 0;
	my $end5_num = 0;
	foreach my $end3_ID (keys %{$end3_info_hsh_ref}) {
		if ($end3_ID =~ m/T$chrom_num(\d+)/) {
			$end3_num = $1 if $1 >$end3_num;
		} else {
			die "end3_ID wrong format\n";
		}
	}

	foreach my $end5_ID (keys %{$end5_info_hsh_ref}) {
		if ($end5_ID =~ m/F$chrom_num(\d+)/) {
			$end5_num = $1 if $1 >$end5_num;
		} else {
			die "end5_ID wrong format\n";
		}
	}

	return ($end3_num, $end5_num);
}
sub getJunctSeq {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $bedtools_bin, $chrom, $chrom_fasta_path, $chrom_num, $input_bed_hsh_ref, $report_tag, $splicing_site_bed, $splicing_site_fasta, $tabix_bin, $threadID
#	output: $junct_seq_hsh_ref
#	toCall: my ($junct_seq_hsh_ref) = &getJunctSeq($tabix_bin, $chrom_num, $threadID, $chrom, $input_bed_hsh_ref, $splicing_site_bed, $splicing_site_fasta, $chrom_fasta_path, $bedtools_bin, $report_tag);
#	calledInLine: 3318
#....................................................................................................................................................#
	my ($tabix_bin, $chrom_num, $threadID, $chrom, $input_bed_hsh_ref, $splicing_site_bed, $splicing_site_fasta, $chrom_fasta_path, $bedtools_bin, $report_tag) = @_;

	my $tmp_hsh_ref = {};
	my $junct_seq_hsh_ref = {};
	
	foreach my $input_ID (keys %{$input_bed_hsh_ref}) {
		&reportAndLogStatus("$report_tag Extracting junction seq in $input_ID", 10, "\n");#->3614
		my $ref_qry = $input_bed_hsh_ref->{$input_ID}{'ref_qry'};
		my $suffix = $input_bed_hsh_ref->{$input_ID}{'suffix'};
		my $trnscpt_bed = $input_bed_hsh_ref->{$input_ID}{'trnscpt_bed'};
		my $num_proc = 0;
		my $tabix_cmd  = "$tabix_bin $trnscpt_bed $chrom |";
		open TABIX, "$tabix_cmd";
		while (<TABIX>) {
			chomp;
			my (undef, $chromStart, $chromEnd, $trnscpt_ID, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\t/;
			#$trnscpt_ID = $trnscpt_ID."-".$suffix if $ref_qry eq 'qry';
			$num_proc++;
			&reportAndLogStatus("$report_tag Extracting junction seq >> $num_proc transcript junction extracted in $input_ID", 10, "\n") if $num_proc%100000 == 0;#->3614
			my @rngAry = ();
			my @blockSizes_ary = split /,/, $blockSizes;
			my @blockStarts_ary = split /,/, $blockStarts;
			foreach my $i (0..$#blockStarts_ary) {
				my $offset = $blockStarts_ary[$i];
				my $exon_start = $chromStart+$offset;
				my $exon_end = $exon_start + $blockSizes_ary[$i];
				push @rngAry, ($exon_start, $exon_end);
			}
		
			for (my $i=1; $i < $#rngAry; $i += 2) {
				my ($intron_start, $intron_end) = ($rngAry[$i], $rngAry[$i+1]);
				my $junct_str = join ":", ($intron_start, $intron_end, $strand);
				if (not exists $tmp_hsh_ref->{$junct_str}) {
					$tmp_hsh_ref->{$junct_str}++;
				}
			}
		}
		close TABIX;
		&reportAndLogStatus("$report_tag Finished junction seq ## $num_proc transcript junction extracted in $input_ID", 10, "\n");#->3614
	}

	&reportAndLogStatus("$report_tag Getting junction seq fasta", 10, "\n");#->3614
	open SITEBED, "| sort -k2,2n -k6,6 >$splicing_site_bed";
	foreach my $junct_str (sort keys %{$tmp_hsh_ref}) {
		my ($junct_start, $junct_end, $strand) = split /:/, $junct_str;
		my $left_start = $junct_start;
		my $left_end = $left_start + 2;
		my $right_end = $junct_end;
		my $right_start = $right_end - 2;
		
		my ($donor_start, $donor_end, $acceptor_start, $acceptor_end) = ($left_start, $left_end, $right_start, $right_end);
		($acceptor_start, $acceptor_end, $donor_start, $donor_end) = ($donor_start, $donor_end, $acceptor_start, $acceptor_end) if $strand eq '-';
		print SITEBED join "", (join "\t", ($chrom, $donor_start, $donor_end, $junct_str.".donor", '1', $strand)), "\n";
		print SITEBED join "", (join "\t", ($chrom, $acceptor_start, $acceptor_end, $junct_str.".acceptor", '1', $strand)), "\n";
	}
	close SITEBED;
	
	system "$bedtools_bin getfasta -s -name -fi $chrom_fasta_path -bed $splicing_site_bed -fo $splicing_site_fasta";
	
	open SSFASTA, "<", $splicing_site_fasta;
	while (<SSFASTA>) {
		#>chr12_2795244_2797137_+.donnor::chr12:2795244-2795246(+)
		#GT
		#>chr12_2795244_2797137_+.acceptor::chr12:2797135-2797137(+)
		#AG
		chomp;
		if ($_ =~ m/^>(.+)\.(\w+)::/) {
			my $junct_str = $1;
			my $don_acc = $2;
			chomp (my $seq = <SSFASTA>);
			$seq = uc($seq);
			$junct_seq_hsh_ref->{$junct_str}{$don_acc} = $seq;
		}
	}
	close SSFASTA;

	return ($junct_seq_hsh_ref);
}
sub logCalledCMDAndScript {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 2_defineout_dirPath|153
#	secondaryAppearInSection: >none
#	input: $ARGVStr, $result_script_dir, $scriptAbsPath
#	output: 
#	toCall: &logCalledCMDAndScript($ARGVStr, $result_script_dir, $scriptAbsPath);
#	calledInLine: 165
#....................................................................................................................................................#
	my ($ARGVStr, $result_script_dir, $scriptAbsPath) = @_;


	my $cpScriptPath = "$result_script_dir/script.ran.pl";
	my $calledCMDPath = "$result_script_dir/called.cmd.txt";
	system "cp -f $scriptAbsPath $cpScriptPath";
	system "chmod 0444 $cpScriptPath"; #---[07/03/2014 18:02] make it read-only to make sure there'll be accodental change of parameters
	open CALLEDCMD, ">", $calledCMDPath;
	print CALLEDCMD join "", ($ARGVStr), "\n";
	close CALLEDCMD;
	
	return ();
}
sub mergeInConfEndRegion {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: defineEndRegion|994
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $bedtools_bin, $chrom, $conf_bed_bgz, $conf_merge_flank, $min_frac_split, $min_size_split, $min_summit_dist_split, $size, $split_merge_in_conf_bed, $tabix_bin
#	output: 
#	toCall: &mergeInConfEndRegion($tabix_bin, $conf_bed_bgz, $bedtools_bin, $conf_merge_flank, $size, $min_size_split, $min_frac_split, $min_summit_dist_split, $chrom, $split_merge_in_conf_bed);
#	calledInLine: 1061
#....................................................................................................................................................#
	my ($tabix_bin, $conf_bed_bgz, $bedtools_bin, $conf_merge_flank, $size, $min_size_split, $min_frac_split, $min_summit_dist_split, $chrom, $split_merge_in_conf_bed) = @_;
	
	my $merge_info_hsh_ref = {};
	my $group_info_hsh_ref = {};
	my $summit_info_hsh_ref = {};
	my $cluster_dist = $conf_merge_flank*2;

	my $bedtools_cmd  = "$tabix_bin $conf_bed_bgz $chrom | $bedtools_bin cluster -d $cluster_dist -s -i stdin | cut -f 2,3,4,5,6,8,13 |";
	open BEDTOOLS, "$bedtools_cmd";
	while (<BEDTOOLS>) {
		#chr1	248838711	248838831	chr1_248838761_248838781_-	9	-	248838767	248838768	55,126,184	1	20	0	7236
		#chr1	248838734	248839030	chr1_248838784_248838980_-	24	-	248838784	248838785	55,126,184	1	196	0	7236
		#chr1	248839009	248839206	chr1_248839059_248839156_+	4	+	248839059	248839060	228,26,28	1	97	0	7236
		#chr1	248858091	248858287	chr1_248858141_248858237_-	9	-	248858226	248858227	55,126,184	1	96	0	7237
		#chr1	248858222	248858358	chr1_248858272_248858308_+	6	+	248858307	248858308	228,26,28	1	36	0	7237
		#chr1	248858265	248858403	chr1_248858315_248858353_-	6	-	248858315	248858316	55,126,184	1	38	0	7237
		chomp;
		my ($chromStart, $chromEnd, $reg_ID, $score, $strand, $summit, $group_ID) = split /\t/;
		$group_info_hsh_ref->{$group_ID}{'reg_ID'}{$reg_ID}{'info'} = [$summit, $score, $chromStart, $chromEnd];
		push @{$group_info_hsh_ref->{$group_ID}{'bounds'}}, ($chromStart, $chromEnd);
		push @{$group_info_hsh_ref->{$group_ID}{'summit'}}, $summit;
		$summit_info_hsh_ref->{$strand}{$summit} = $reg_ID;
		$group_info_hsh_ref->{$group_ID}{'ST'} = $strand;
		$group_info_hsh_ref->{$group_ID}{'score'} += $score;
		#&reportAndLogStatus("$chromStart $chromEnd $reg_ID $score $strand $summit $group_ID", 10, "\n");#->3614
	}
	close BEDTOOLS;

	my $merge_ID = 0;
	foreach my $group_ID (keys %{$group_info_hsh_ref}) {
		my $total_score = $group_info_hsh_ref->{$group_ID}{'score'};
		my $strand = $group_info_hsh_ref->{$group_ID}{'ST'};
		#&reportAndLogStatus("$group_ID $total_score", 10, "\n");#->3614

		@{$group_info_hsh_ref->{$group_ID}{'bounds'}} = sort {$a <=> $b} @{$group_info_hsh_ref->{$group_ID}{'bounds'}};
		@{$group_info_hsh_ref->{$group_ID}{'summit'}} = sort {$a <=> $b} @{$group_info_hsh_ref->{$group_ID}{'summit'}};
		my $groupStart = $group_info_hsh_ref->{$group_ID}{'bounds'}[0];
		my $groupEnd = $group_info_hsh_ref->{$group_ID}{'bounds'}[-1];
		my $summitStart = $group_info_hsh_ref->{$group_ID}{'summit'}[0];
		my $summitEnd = $group_info_hsh_ref->{$group_ID}{'summit'}[-1];
		my $groupSize = $summitEnd - $summitStart;
		my @bound_ary = ($groupStart, $groupEnd);
		if ($groupSize > $min_size_split) {
			foreach my $reg_ID (keys %{$group_info_hsh_ref->{$group_ID}{'reg_ID'}}) {
				my $score = ${$group_info_hsh_ref->{$group_ID}{'reg_ID'}{$reg_ID}{'info'}}[1];
				my $frac = $score/$total_score;
				if ($frac >= $min_frac_split) {
					my $summit = ${$group_info_hsh_ref->{$group_ID}{'reg_ID'}{$reg_ID}{'info'}}[0];
					push @bound_ary, $summit;
				}
			}
		}

		#0,1
		#0,1,2
		#0,1,2,3
		#0,1,2,3,4
		#0,1,2,3,4,5
		
		@bound_ary = sort {$a <=> $b} @bound_ary;

		my $max_index = $#bound_ary;
		if ($max_index == 1) {#---[2023/08/15 21:51] no split or no reg_ID >$min_frac_split
			$merge_ID++;
			my $chromStart = $groupStart - $conf_merge_flank;
			my $chromEnd = $groupEnd + $conf_merge_flank;
			$chromStart = 1 if $chromStart < 1;
			$chromEnd = $size if $chromEnd > $size;
			$merge_info_hsh_ref->{$merge_ID}{'chromStart'} = $chromStart;
			$merge_info_hsh_ref->{$merge_ID}{'chromEnd'} = $chromEnd;
			$merge_info_hsh_ref->{$merge_ID}{'ST'} = $strand;
			$merge_info_hsh_ref->{$merge_ID}{'score'} = $total_score;
		} else {
			if ($max_index >= 3) {#---[2023/08/15 21:51] at least 2 reg_ID >$min_frac_split
				my @valid_summit_ary = ();
				foreach my $i (1..($max_index-2)) {
					my $curnt_index = $i;
					my $next_index = $i + 1;
					my $curnt_summit = $bound_ary[$curnt_index];
					my $next_summit = $bound_ary[$next_index];
					push @valid_summit_ary, $curnt_summit if $next_summit - $curnt_summit > $min_summit_dist_split;
				}
				@bound_ary = ($bound_ary[0], @valid_summit_ary, $bound_ary[-2], $bound_ary[-1]);
			}
		
			$max_index = $#bound_ary;
			foreach my $i (1..($max_index-1)) {#---[2023/08/15 21:51] has split or at least 1 reg_ID >$min_frac_split
				my $curnt_index = $i;
				my $prior_index = $i - 1;
				my $next_index = $i + 1;
				my $chromStart;
				my $chromEnd;
				my $curnt_summit = $bound_ary[$curnt_index];
				my $curnt_reg_ID = $summit_info_hsh_ref->{$strand}{$curnt_summit};
				my (undef, $curnt_score, $curnt_chromStart, $curnt_chromEnd) = @{$group_info_hsh_ref->{$group_ID}{'reg_ID'}{$curnt_reg_ID}{'info'}};

				if ($prior_index == 0) {
					$chromStart = $groupStart - $conf_merge_flank;
					$chromStart = 1 if $chromStart < 1;
				} else {
					my $prior_summit = $bound_ary[$prior_index];
					my $prior_reg_ID = $summit_info_hsh_ref->{$strand}{$prior_summit};
					my (undef, $prior_score, $prior_chromStart, $prior_chromEnd) = @{$group_info_hsh_ref->{$group_ID}{'reg_ID'}{$prior_reg_ID}{'info'}};
					$chromStart = sprintf "%.0f", $curnt_chromStart - ($curnt_chromStart-$prior_chromEnd)/2;
				}

				if ($next_index == $max_index) {
					$chromEnd = $groupEnd + $conf_merge_flank;
					$chromEnd = $size if $chromEnd > $size;
				} else {
					my $next_summit = $bound_ary[$next_index];
					my $next_reg_ID = $summit_info_hsh_ref->{$strand}{$next_summit};
					my (undef, $next_score, $next_chromStart, $next_chromEnd) = @{$group_info_hsh_ref->{$group_ID}{'reg_ID'}{$next_reg_ID}{'info'}};
					$chromEnd = sprintf "%.0f", $curnt_chromEnd + ($next_chromStart-$curnt_chromEnd)/2;
				}
				$merge_ID++;
				$merge_info_hsh_ref->{$merge_ID}{'ST'} = $strand;
				$merge_info_hsh_ref->{$merge_ID}{'score'} = $curnt_score;
				$merge_info_hsh_ref->{$merge_ID}{'chromStart'} = $chromStart;
				$merge_info_hsh_ref->{$merge_ID}{'chromEnd'} = $chromEnd;
			}
		}
	}

	open MERGESPLITBED, "| sort -k2,2n | awk 'BEGIN {FS=OFS=\"\\t\"} { \$4=\"C\"NR; print }' >$split_merge_in_conf_bed";
	foreach my $merge_ID (sort keys %{$merge_info_hsh_ref}) {
		my $chromStart = $merge_info_hsh_ref->{$merge_ID}{'chromStart'};
		my $chromEnd = $merge_info_hsh_ref->{$merge_ID}{'chromEnd'};
		my $strand = $merge_info_hsh_ref->{$merge_ID}{'ST'};
		my $score = $merge_info_hsh_ref->{$merge_ID}{'score'};
		
		print MERGESPLITBED join "", (join "\t", ($chrom, $chromStart, $chromEnd, $merge_ID, $score, $strand)), "\n";
	}
	close MERGESPLITBED;

	return ();

}
sub mergeRefRegionAndInConfEndRegion {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: defineEndRegion|994
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $bedtools_bin, $chrom_num, $end_prefix, $ref_in_conf_merge_reg_bed, $ref_reg_bed, $split_merge_in_conf_bed
#	output: 
#	toCall: &mergeRefRegionAndInConfEndRegion($split_merge_in_conf_bed, $ref_reg_bed, $bedtools_bin, $ref_in_conf_merge_reg_bed, $end_prefix, $chrom_num);
#	calledInLine: 1063
#....................................................................................................................................................#
	my ($split_merge_in_conf_bed, $ref_reg_bed, $bedtools_bin, $ref_in_conf_merge_reg_bed, $end_prefix, $chrom_num) = @_;
	
	my $group_info_hsh_ref = {};
	open BEDTOOLS, "cat $split_merge_in_conf_bed $ref_reg_bed | sort -k2,2n -k6,6 | $bedtools_bin cluster -s -d -1 -i stdin |";
	while (<BEDTOOLS>) {
		#chr1	9994	10226	C1	6	+	1
		#chr1	11817	11918	R1	1	+	2
		#chr1	11958	12059	R2	1	+	3
		#chr1	12095	12311	C2	5	+	4
		#chr1	14458	14671	C4	26	+	5
		#chr1	20203	20381	C9	11	+	6
		#chr1	20381	20595	C10	5	+	7
		#chr1	28657	28888	C11	8	+	8
		my ($chrom, $chromStart, $chromEnd, $reg_ID, $score, $strand, $group_ID) = split /\t/;
		$group_info_hsh_ref->{$group_ID}{'reg_ID'}{$reg_ID}{'info'} = [$chrom, $chromStart, $chromEnd, $reg_ID, $score, $strand];
		$group_info_hsh_ref->{$group_ID}{'strand'} = $strand;
		$group_info_hsh_ref->{$group_ID}{'chrom'} = $chrom;
		push @{$group_info_hsh_ref->{$group_ID}{'all_bounds'}}, ($chromStart, $chromEnd);
		
		if ($reg_ID =~ m/C/) {
			push @{$group_info_hsh_ref->{$group_ID}{'C_bounds'}}, ($chromStart, $chromEnd);
			$group_info_hsh_ref->{$group_ID}{'C_reg_ID'}{$reg_ID} = $chromStart;
		}
		
	}
	close BEDTOOLS;
	
	my $split_info_hsh_ref = {};
	my $split_num = 0;
	foreach my $group_ID (sort {$a <=> $b} keys %{$group_info_hsh_ref}) {
		my $strand = $group_info_hsh_ref->{$group_ID}{'strand'};
		my $chrom = $group_info_hsh_ref->{$group_ID}{'chrom'};
		my @reg_ID_ary = keys %{$group_info_hsh_ref->{$group_ID}{'reg_ID'}};
		if (@reg_ID_ary == 1) {
			my $reg_ID = $reg_ID_ary[0];
			my (undef, $chromStart, $chromEnd, undef, $score, undef) = @{$group_info_hsh_ref->{$group_ID}{'reg_ID'}{$reg_ID_ary[0]}{'info'}};
			$split_num++;
			my $split_ID = $split_num;
			my $ref_only = 'N';
			$ref_only = 'Y' if ($reg_ID =~ m/R/);
			$split_info_hsh_ref->{$split_ID}{'chrom'} = $chrom;
			$split_info_hsh_ref->{$split_ID}{'chromStart'} = $chromStart;
			$split_info_hsh_ref->{$split_ID}{'chromEnd'} = $chromEnd;
			$split_info_hsh_ref->{$split_ID}{'ST'} = $strand;
			$split_info_hsh_ref->{$split_ID}{'ref_only'} = $ref_only;
		} else {
			die "C_reg_ID does not exists\n" if not exists $group_info_hsh_ref->{$group_ID}{'C_reg_ID'}; #---[2023/08/16 21:51] C_reg_ID must exists since R and C themselves are initially non-overlapping
			@{$group_info_hsh_ref->{$group_ID}{'all_bounds'}} = sort {$a <=> $b} @{$group_info_hsh_ref->{$group_ID}{'all_bounds'}};
			my @extended_bound_ary = sort {$a <=> $b} @{$group_info_hsh_ref->{$group_ID}{'C_bounds'}};
			$extended_bound_ary[0] = $group_info_hsh_ref->{$group_ID}{'all_bounds'}[0]; #---[2023/08/16 22:02] take the leftmost bounds
			$extended_bound_ary[-1] = $group_info_hsh_ref->{$group_ID}{'all_bounds'}[-1]; #---[2023/08/16 22:02] take the rightmost bounds
			
			#0,1
			#0,1,2,3
			#0,1,2,3,4,5
			my $max_index = $#extended_bound_ary;
			for (my $i=0; $i<$max_index; $i=$i+2) {
				my $curnt_start = $extended_bound_ary[$i];
				my $curnt_end = $extended_bound_ary[$i+1];
				my $last_end_index = $i - 1;
				my $next_start_index = $i + 2;
				my $chromEnd;
				my $chromStart;
				
				if ($last_end_index > 0) {
					my $last_end = $extended_bound_ary[$last_end_index];
					$chromStart = sprintf "%.0f", $curnt_start - ($curnt_start-$last_end)/2;
				} else {
					$chromStart = $curnt_start;
				}
				
				if ($next_start_index < $max_index) {
					my $next_start = $extended_bound_ary[$next_start_index];
					$chromEnd = sprintf "%.0f", $curnt_end + ($next_start-$curnt_end)/2;
				} else {
					$chromEnd = $curnt_end;
				}

				$split_num++;
				my $split_ID = $split_num;
				$split_info_hsh_ref->{$split_ID}{'chrom'} = $chrom;
				$split_info_hsh_ref->{$split_ID}{'chromStart'} = $chromStart;
				$split_info_hsh_ref->{$split_ID}{'chromEnd'} = $chromEnd;
				$split_info_hsh_ref->{$split_ID}{'ST'} = $strand;
				$split_info_hsh_ref->{$split_ID}{'ref_only'} = 'N';
			}
		}
	}
	
	open REFINCONFMERGE, "| sort -k2,2n | awk 'BEGIN {FS=OFS=\"\\t\"} { \$4=\"$end_prefix$chrom_num\"NR; print }' >$ref_in_conf_merge_reg_bed";
	foreach my $split_ID (sort keys %{$split_info_hsh_ref}) {
		my $chrom = $split_info_hsh_ref->{$split_ID}{'chrom'};
		my $chromStart = $split_info_hsh_ref->{$split_ID}{'chromStart'};
		my $chromEnd = $split_info_hsh_ref->{$split_ID}{'chromEnd'};
		my $strand = $split_info_hsh_ref->{$split_ID}{'ST'};
		my $ref_only = $split_info_hsh_ref->{$split_ID}{'ref_only'};
		my $ref_only_digit = 0;
		$ref_only_digit = 1 if $ref_only eq 'Y';
		print REFINCONFMERGE join "", (join "\t", ($chrom, $chromStart, $chromEnd, $split_ID, $ref_only_digit, $strand)), "\n";
	}
	close REFINCONFMERGE;
	
	
	return ();
}
sub poolChromBedAndTable {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: >none
#	primaryAppearInSection: 3_process|197
#	secondaryAppearInSection: >none
#	input: $chrom_info_hsh_ref, $max_thread, $output_path_hsh_ref, $result_tmp_dir, $tabix_bin
#	output: 
#	toCall: &poolChromBedAndTable($chrom_info_hsh_ref, $output_path_hsh_ref, $result_tmp_dir, $tabix_bin, $max_thread);
#	calledInLine: 205
#....................................................................................................................................................#
	my ($chrom_info_hsh_ref, $output_path_hsh_ref, $result_tmp_dir, $tabix_bin, $max_thread) = @_;

	my %fileForThrHsh = ();
	my $threadID = 1;
	
	foreach my $out_file (sort {$output_path_hsh_ref->{$b}{'size_group'} <=> $output_path_hsh_ref->{$a}{'size_group'}} keys %{$output_path_hsh_ref}) {
		$threadID = 1 if $threadID > $max_thread;
		$threadID = sprintf "%.2d", $threadID;
		push @{$fileForThrHsh{$threadID}} , $out_file;
		$threadID++;
	}
	
	my %threadHsh =();
	foreach my $threadID (sort {$a <=> $b} keys %fileForThrHsh) {
		my $fileForThrAry_ref = $fileForThrHsh{$threadID};
		($threadHsh{$threadID}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781, the 
	
			sub {
				my ($fileForThrAry_ref) = @_;
				
				foreach my $out_file (@{$fileForThrAry_ref}) {
					#---[2/19/16 13:13] do something on the item here
					my $type = $output_path_hsh_ref->{$out_file}{'type'};
					my $path = $output_path_hsh_ref->{$out_file}{'path'};

					my $report_tag = join " ", (sprintf("%-15s","Thread-$threadID $out_file"), ":");
					&reportAndLogStatus("$report_tag Pooling $type files", 10, "\n");#->3614

					if ($type eq 'info') {
						my $tmp_header_path = "$result_tmp_dir/tmp.$out_file.header.gz";
						my @header_ary = @{$output_path_hsh_ref->{$out_file}{'header'}};
						open HEADER, "| gzip -c >$tmp_header_path";
						print HEADER join "", (join "\t", (@header_ary)), "\n";
						close HEADER;
						my @file_ary = ($tmp_header_path);
						foreach my $chrom (sort keys %{$chrom_info_hsh_ref}) {
							my $chrom_file_path = $chrom_info_hsh_ref->{$chrom}{$out_file};
							push @file_ary, $chrom_file_path;
						}
						my $file_str = join " ", @file_ary;
						system "cat $file_str >$path";
					} elsif ($type eq 'bed') {
						my @file_ary = ();
						foreach my $chrom (sort keys %{$chrom_info_hsh_ref}) {
							my $chrom_file_path = $chrom_info_hsh_ref->{$chrom}{$out_file};
							push @file_ary, $chrom_file_path;
						}
						my $file_str = join " ", @file_ary;
						system "cat $file_str >$path";
						system "$tabix_bin -p bed $path";
					}
				}
				
				return ();
			}
			,($fileForThrAry_ref)
		);
	}
	
	my $data_hsh_ref = {};
	while (keys %threadHsh) {
		my $num_thread = keys %threadHsh;
		&reportAndLogStatus("--------------------------- $num_thread thread running ---------------------------", 10, "\n");#->3614
		foreach my $threadID (keys %threadHsh) {
			if (not $threadHsh{$threadID}->is_running()) {
				$threadHsh{$threadID}->join();
				delete $threadHsh{$threadID};
				&reportAndLogStatus(">>>>>>>>>>>>>>>>>>>>>>>> thread $threadID finished <<<<<<<<<<<<<<<<<<<<<<<<<<<", 10, "\n");#->3614
			}
		}
		sleep 5;
	}
	
	return ();
}
sub poolChromLog {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: >none
#	primaryAppearInSection: 3_process|197
#	secondaryAppearInSection: >none
#	input: $chrom_info_hsh_ref, $paramTag, $result_log_dir
#	output: 
#	toCall: &poolChromLog($chrom_info_hsh_ref, $result_log_dir, $paramTag);
#	calledInLine: 203
#....................................................................................................................................................#
	my ($chrom_info_hsh_ref, $result_log_dir, $paramTag) = @_;

	&reportAndLogStatus("Pooling log txt from all chromosome", 10, "\n");	#->3614
	my $pool_log_path = "$result_log_dir/$paramTag.stats.log.txt";
	my $count_hsh_ref = {};
	foreach my $chrom (sort keys %{$chrom_info_hsh_ref}) {
		my $log_txt = $chrom_info_hsh_ref->{$chrom}{'log_txt'};
		die "log_txt of chrom does not exists. Quitting.\n" if not -s $log_txt;
		#0	ref transcript removed	3
		#1	ref transcript retained	21574
		#2	qry transcript removed	0
		#3	qry transcript retained	227936
		#4	confident end5 region num	9377
		#5	confident end3 region num	18888
		#6	doubtful $end region num	0
		#7	doubtful $end region num	0
		#8	ref transcript complete	9355
		#9	ref transcript incomplete	12219
		#10	ref transcript with confident end5	9355
		#11	ref transcript with doubtful end5	12219
		#12	ref transcript with confident end3	21574
		open CHROMLOG, "<", $log_txt;
		while (<CHROMLOG>) {
			chomp;
			my ($index, $log_str, $count) = split /\t/;
			$count_hsh_ref->{'all'}{$index}{$log_str} += $count;
			$count_hsh_ref->{$chrom}{$index}{$log_str} = $count;
		}
		close CHROMLOG;
	}
	
	open LOG, ">", $pool_log_path;
	print LOG join "", (join "\t", ('chrom', 'stats', 'count')), "\n";
	foreach my $chrom (sort keys %{$count_hsh_ref}) {
		foreach my $index (sort {$a <=> $b} keys %{$count_hsh_ref->{$chrom}}) {
			foreach my $log_str (keys %{$count_hsh_ref->{$chrom}{$index}}) {
				my $count = $count_hsh_ref->{$chrom}{$index}{$log_str};
				print LOG join "", (join "\t", ($chrom, $log_str, $count)), "\n";
			}
		}
	}
	close LOG;

	return ();
}
sub printInfoBedBothEnd {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $all_end3_bed, $all_end5_bed, $bgzip_bin, $chrom, $end3_info_hsh_ref, $end3_info_tsv, $end5_info_hsh_ref, $end5_info_tsv, $report_tag, $set_info_hsh_ref
#	output: 
#	toCall: &printInfoBedBothEnd($set_info_hsh_ref, $end3_info_tsv, $all_end3_bed, $end3_info_hsh_ref, $end5_info_tsv, $all_end5_bed, $end5_info_hsh_ref, $chrom, $report_tag, $bgzip_bin);
#	calledInLine: 3363
#....................................................................................................................................................#
	my ($set_info_hsh_ref, $end3_info_tsv, $all_end3_bed, $end3_info_hsh_ref, $end5_info_tsv, $all_end5_bed, $end5_info_hsh_ref, $chrom, $report_tag, $bgzip_bin) = @_;
	
	my $info_end_hsh_ref = {
		'end5' => {
			'end_hsh_ref' => $end5_info_hsh_ref,
			'end_bed' => $all_end5_bed,
			'info_tsv' => $end5_info_tsv,
		},
		'end3' => {
			'end_hsh_ref' => $end3_info_hsh_ref,
			'end_bed' => $all_end3_bed,
			'info_tsv' => $end3_info_tsv,
		},
	};

	foreach my $end (sort keys %{$info_end_hsh_ref}) {
		my $end_hsh_ref = $info_end_hsh_ref->{$end}{'end_hsh_ref'};
		my $end_bed = $info_end_hsh_ref->{$end}{'end_bed'};
		my $info_tsv = $info_end_hsh_ref->{$end}{'info_tsv'};
		open ENDINFO, "| sort -k10,10n | gzip -c >$info_tsv";
		open ENDBED, "| sort -k2,2n -k6,6 | $bgzip_bin -c >$end_bed";
		my $num_proc = 0;
		foreach my $endID (sort keys %{$end_hsh_ref}) {
			$num_proc++;
			&reportAndLogStatus("$report_tag Printing $end txt and bed >> $num_proc $end region processed", 10, "\n") if $num_proc%10000==0;#->3614
			my $ref_count = 0;
			my $qry_count = 0;
			my $set_count = 0;
			my $set_str = "__na";
			my $model_count = 0;
			my $model_str = "__na";
			my $ref_only = $end_hsh_ref->{$endID}{'ref_only'};
			
			if (exists $end_hsh_ref->{$endID}{'CT'}{'ref'}) {
				$ref_count = $end_hsh_ref->{$endID}{'CT'}{'ref'};
			}
			if (exists $end_hsh_ref->{$endID}{'CT'}{'qry'}) {
				$qry_count = $end_hsh_ref->{$endID}{'CT'}{'qry'};
			}
			
			if (exists $end_hsh_ref->{$endID}{'set'}) {
				$set_count = keys %{$end_hsh_ref->{$endID}{'set'}};
				if ($set_count > 0) {
					my $model_hsh_ref = {};
					$set_str = join ";", sort keys %{$end_hsh_ref->{$endID}{'set'}};
					foreach my $set_ID (keys %{$end_hsh_ref->{$endID}{'set'}}) {
						foreach my $model_ID (keys %{$set_info_hsh_ref->{$set_ID}{'model'}}) {
							$model_hsh_ref->{$model_ID}++;
						}
					}
					$model_str = join ";", sort keys %{$model_hsh_ref};
					$model_count = keys %{$model_hsh_ref};
				}
			}
			
			my $summit = $end_hsh_ref->{$endID}{'summit'};
			my $score = $end_hsh_ref->{$endID}{'score'};
			my $strand = $end_hsh_ref->{$endID}{'ST'};
			my $chromStart = $end_hsh_ref->{$endID}{'CS'};
			my $chromEnd = $end_hsh_ref->{$endID}{'CE'};
			my $obsStart = $end_hsh_ref->{$endID}{'OS'};
			my $obsEnd = $end_hsh_ref->{$endID}{'OE'};
			my $thickEnd = $summit;
			my $thickStart = $thickEnd - 1;
			my $blockCount = 1;
			my $blockSizes = $chromEnd-$chromStart;
			my $blockStarts = 0;
			#my $itemRgb = '77,175,74';
			my $confidence = 'confident';
			if ($endID =~ m/X/) {
				#$itemRgb = '255,127,0';
				$confidence = 'doubtful';
			}
			
			my $itemRgb = '228,26,28';
			$itemRgb = '55,126,184' if ($strand eq '-');
			
			print ENDINFO join "", (join "\t", ($endID, $confidence, $ref_only, $ref_count, $qry_count, $set_count, $model_count, $chrom, $summit, $strand, $chromStart, $chromEnd, $obsStart, $obsEnd, $set_str, $model_str)), "\n";
			print ENDBED join "", (join "\t", ($chrom, $chromStart, $chromEnd, $endID, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts)), "\n";
		}
		&reportAndLogStatus("$report_tag Finished printing $end txt and bed ## $num_proc $end region processed", 10, "\n");#->3614
	}
	close ENDINFO;
	close ENDBED;

	return ();
}
sub printInfoBedJunction {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $bgzip_bin, $chrom, $junct_bed, $junct_info_hsh_ref, $junct_info_tsv
#	output: 
#	toCall: &printInfoBedJunction($bgzip_bin, $junct_info_hsh_ref, $junct_info_tsv, $junct_bed, $chrom);
#	calledInLine: 3338
#....................................................................................................................................................#
	my ($bgzip_bin, $junct_info_hsh_ref, $junct_info_tsv, $junct_bed, $chrom) = @_;
	
	open JUNCTINFO, "| sort -k3,3n | gzip -c >$junct_info_tsv";
	open JUNCTBED, "| sort -k2,2n -k6,6 | $bgzip_bin -c >$junct_bed";
	#print JUNCTINFO join "", (join "\t", ('junct_ID', 'chrom', 'intron_start', 'intron_end', 'strand', 'ref_count', 'qry_count')), "\n";
	foreach my $junct_ID (sort keys %{$junct_info_hsh_ref}) {
		my $ref_count = 0;
		my $qry_count = 0;
		my $splicing_site = $junct_info_hsh_ref->{$junct_ID}{'SS'}; #---[2023/07/21 12:34] junction string
		my $canonical = $junct_info_hsh_ref->{$junct_ID}{'CN'}; #---[2023/07/21 12:34] junction string

		my ($intron_start, $intron_end, $strand) = split /:/, $junct_info_hsh_ref->{$junct_ID}{'JS'}; #---[2023/07/21 12:34] junction string
		
		if (exists $junct_info_hsh_ref->{$junct_ID}{'CT'}{'ref'}) {
			$ref_count = $junct_info_hsh_ref->{$junct_ID}{'CT'}{'ref'};
		}
		if (exists $junct_info_hsh_ref->{$junct_ID}{'CT'}{'qry'}) {
			$qry_count = $junct_info_hsh_ref->{$junct_ID}{'CT'}{'qry'};
		}
		
		print JUNCTBED join "", (join "\t", ($chrom, $intron_start, $intron_end, $junct_ID."|".$splicing_site, $qry_count, $strand)), "\n";
		print JUNCTINFO join "", (join "\t", ($junct_ID, $canonical, $splicing_site, $chrom, $intron_start, $intron_end, $strand, $ref_count, $qry_count)), "\n";
	}
	close JUNCTINFO;
	close JUNCTBED;

	return ();
}
sub printInfoModel {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $chrom, $log_ary_ref, $model_info_hsh_ref, $model_info_tsv, $set_info_hsh_ref
#	output: 
#	toCall: &printInfoModel($model_info_hsh_ref, $model_info_tsv, $set_info_hsh_ref, $chrom, $log_ary_ref);
#	calledInLine: 3351
#....................................................................................................................................................#
	my ($model_info_hsh_ref, $model_info_tsv, $set_info_hsh_ref, $chrom, $log_ary_ref) = @_;

	my $count_hsh_ref = {};
	
	foreach my $all_ref_qry (qw/all ref qry/) {
		foreach my $completeness (qw/complete incomplete/) {
			$count_hsh_ref->{$all_ref_qry}{'model'}{$completeness} = 0;
		}
		foreach my $end (qw/end5 end3/) {
			foreach my $confidence (qw/confident doubtful/) {
				$count_hsh_ref->{$all_ref_qry}{'model_end'}{$end}{$confidence} = 0;
			}
			foreach my $endtype (qw/C L S/) {
				$count_hsh_ref->{$all_ref_qry}{'model_endtype'}{$end}{$endtype} = 0;
			}
		}
	}
	
	open INFO, "| gzip -c >$model_info_tsv";
	#print INFO join "", (join "\t", ('model_ID', 'loc', 'completeness', 'full_set_ID', 'end5_conf', 'end3_conf', 'junct_conf', 'end5_endtype', 'end3_endtype', 'novelty', 'strand', 'partial_set_count', 'partial_set_ID_str', 'full_ref_count', 'full_qry_count', 'partial_ref_count', 'partial_qry_count', 'full_set_bound_str')), "\n";
	foreach my $model_ID (sort keys %{$model_info_hsh_ref}) {
		my $completeness = $model_info_hsh_ref->{$model_ID}{'CM'}; #---[2023/07/21 12:34] completeness: not complete
		my $full_set_ID = $model_info_hsh_ref->{$model_ID}{'set'};#---[2023/07/21 12:34] the setID
		my $full_set_bound_str = $set_info_hsh_ref->{$full_set_ID}{'BS'};#---[2023/07/21 12:34] the bound str
		my $end3_endtype = $set_info_hsh_ref->{$full_set_ID}{'T3'};
		my $end5_endtype = $set_info_hsh_ref->{$full_set_ID}{'T5'};
		my $novelty = $model_info_hsh_ref->{$model_ID}{'N'};
		my $strand = $set_info_hsh_ref->{$full_set_ID}{'ST'};
		my $chromStart = $set_info_hsh_ref->{$full_set_ID}{'CS'};
		my $chromEnd = $set_info_hsh_ref->{$full_set_ID}{'CE'};
		my $loc = $chrom.":".$chromStart."-".$chromEnd;

		my $end5_conf = 'Y';
		my $end3_conf = 'Y';
		my $junct_conf = 'Y';

		$end5_conf = 'N' if $full_set_bound_str =~ m/XF/;
		$end3_conf = 'N' if $full_set_bound_str =~ m/XT/;
		$junct_conf = 'N' if $full_set_bound_str =~ m/XJ/;

		my $full_ref_count = 0;
		my $full_qry_count = 0;
		my $partial_ref_count = 0;
		my $partial_qry_count = 0;
		my $partial_set_count = 0;
		my $partial_set_ID_str = '__na';
		
		if (exists $set_info_hsh_ref->{$full_set_ID}{'T'}{'ref'}) {
			$full_ref_count = keys %{$set_info_hsh_ref->{$full_set_ID}{'T'}{'ref'}};
		}
		if (exists $set_info_hsh_ref->{$full_set_ID}{'T'}{'qry'}) {
			$full_qry_count = keys %{$set_info_hsh_ref->{$full_set_ID}{'T'}{'qry'}};
		}

		if (exists $model_info_hsh_ref->{$model_ID}{'IS'}) {
			$partial_set_count = keys %{$model_info_hsh_ref->{$model_ID}{'IS'}}; #---[2023/07/21 14:55] partial set count
			$partial_set_ID_str = join ";", sort keys %{$model_info_hsh_ref->{$model_ID}{'IS'}}; #---[2023/07/21 14:55] set ID
			foreach my $partial_set_ID (keys %{$model_info_hsh_ref->{$model_ID}{'IS'}}) {
				if (exists $set_info_hsh_ref->{$partial_set_ID}{'T'}{'ref'}) {
					$partial_ref_count = keys %{$set_info_hsh_ref->{$partial_set_ID}{'T'}{'ref'}};
				}
				if (exists $set_info_hsh_ref->{$partial_set_ID}{'T'}{'qry'}) {
					$partial_qry_count = keys %{$set_info_hsh_ref->{$partial_set_ID}{'T'}{'qry'}};
				}
			}
		}
		print INFO join "", (join "\t", ($model_ID, $loc, $completeness, $full_set_ID, $end5_conf, $end3_conf, $junct_conf, $end5_endtype, $end3_endtype, $novelty, $strand, $partial_set_count, $partial_set_ID_str, $full_ref_count, $full_qry_count, $partial_ref_count, $partial_qry_count, $full_set_bound_str)), "\n";
		my $ref_count = $full_ref_count + $partial_ref_count;
		my $qry_count = $full_qry_count + $partial_qry_count;

		if ($completeness eq 'Y') {
			$count_hsh_ref->{'all'}{'model'}{'complete'}++;
			$count_hsh_ref->{'qry'}{'model'}{'complete'}++ if $qry_count > 0;
			$count_hsh_ref->{'ref'}{'model'}{'complete'}++ if $ref_count > 0;
		} else {
			$count_hsh_ref->{'all'}{'model'}{'incomplete'}++;
			$count_hsh_ref->{'qry'}{'model'}{'incomplete'}++ if $qry_count > 0;
			$count_hsh_ref->{'ref'}{'model'}{'incomplete'}++ if $ref_count > 0;
		}

		$count_hsh_ref->{'all'}{'model_endtype'}{'end3'}{$end3_endtype}++;
		$count_hsh_ref->{'qry'}{'model_endtype'}{'end3'}{$end3_endtype}++ if $qry_count > 0;
		$count_hsh_ref->{'ref'}{'model_endtype'}{'end3'}{$end3_endtype}++ if $ref_count > 0;

		$count_hsh_ref->{'all'}{'model_endtype'}{'end5'}{$end5_endtype}++;
		$count_hsh_ref->{'qry'}{'model_endtype'}{'end5'}{$end5_endtype}++ if $qry_count > 0;
		$count_hsh_ref->{'ref'}{'model_endtype'}{'end5'}{$end5_endtype}++ if $ref_count > 0;
		
		
		if ($full_set_bound_str =~ /XT/) {
			$count_hsh_ref->{'all'}{'model_end'}{'end3'}{'doubtful'}++;
			$count_hsh_ref->{'qry'}{'model_end'}{'end3'}{'doubtful'}++ if $qry_count > 0;
			$count_hsh_ref->{'ref'}{'model_end'}{'end3'}{'doubtful'}++ if $ref_count > 0;
		} else {
			$count_hsh_ref->{'all'}{'model_end'}{'end3'}{'confident'}++;
			$count_hsh_ref->{'qry'}{'model_end'}{'end3'}{'confident'}++ if $qry_count > 0;
			$count_hsh_ref->{'ref'}{'model_end'}{'end3'}{'confident'}++ if $ref_count > 0;
		}
		if ($full_set_bound_str =~ /XF/) {
			$count_hsh_ref->{'all'}{'model_end'}{'end5'}{'doubtful'}++;
			$count_hsh_ref->{'qry'}{'model_end'}{'end5'}{'doubtful'}++ if $qry_count > 0;
			$count_hsh_ref->{'ref'}{'model_end'}{'end5'}{'doubtful'}++ if $ref_count > 0;
		} else {
			$count_hsh_ref->{'all'}{'model_end'}{'end5'}{'confident'}++;
			$count_hsh_ref->{'qry'}{'model_end'}{'end5'}{'confident'}++ if $qry_count > 0;
			$count_hsh_ref->{'ref'}{'model_end'}{'end5'}{'confident'}++ if $ref_count > 0;
		}

	}
	close INFO;

	my $endtype_prefix_hsh_ref = {
		'C' => 'commonest',
		'L' => 'longest',
		'S' => 'summit',
	};

	foreach my $all_ref_qry (qw/all ref qry/) {
		foreach my $completeness (qw/complete incomplete/) {
			my $count = $count_hsh_ref->{$all_ref_qry}{'model'}{$completeness};
			my $log_str = "$all_ref_qry model $completeness";
			push @{$log_ary_ref}, [$log_str, $count];
		}
		foreach my $end (qw/end5 end3/) {
			foreach my $confidence (qw/confident doubtful/) {
				my $count = $count_hsh_ref->{$all_ref_qry}{'model_end'}{$end}{$confidence};
				my $log_str = "$all_ref_qry model with $confidence $end";
				push @{$log_ary_ref}, [$log_str, $count];
			}
			foreach my $endtype (qw/C L S/) {
				my $full_endtype = $endtype_prefix_hsh_ref->{$endtype};
				my $count = $count_hsh_ref->{$all_ref_qry}{'model_endtype'}{$end}{$endtype};
				my $log_str = "$all_ref_qry model $end adjusted by $full_endtype";
				push @{$log_ary_ref}, [$log_str, $count];
			}
		}
	}
	
	return ();
}
sub printInfoModelSetPair {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $chrom, $model_info_hsh_ref, $model_set_pair_info_tsv, $set_info_hsh_ref
#	output: 
#	toCall: &printInfoModelSetPair($model_info_hsh_ref, $model_set_pair_info_tsv, $set_info_hsh_ref, $chrom);
#	calledInLine: 3354
#....................................................................................................................................................#
	my ($model_info_hsh_ref, $model_set_pair_info_tsv, $set_info_hsh_ref, $chrom) = @_;
	
	open INFO, "| gzip -c >$model_set_pair_info_tsv";
	#print INFO join "", (join "\t", ('model_ID', 'set_ID', 'loc', 'end5_endtype', 'end3_endtype', 'full_partial', 'junct_num', 'length', 'ref_count', 'qry_count', 'bound_str')), "\n";
	foreach my $model_ID (sort keys %{$model_info_hsh_ref}) {
		my $set_ID = $model_info_hsh_ref->{$model_ID}{'set'};#---[2023/07/21 12:34] the setID
		my $end3_endtype = $set_info_hsh_ref->{$set_ID}{'T3'};
		my $end5_endtype = $set_info_hsh_ref->{$set_ID}{'T5'};
		my $setID_hsh_ref = {
			$set_ID => 'full',
		};
		if (exists $model_info_hsh_ref->{$model_ID}{'IS'}) {
			foreach my $set_ID (keys %{$model_info_hsh_ref->{$model_ID}{'IS'}}) {
				$setID_hsh_ref->{$set_ID} = 'partial';
			}
		}

		foreach my $set_ID (keys %{$setID_hsh_ref}) {
			my $bound_str = $set_info_hsh_ref->{$set_ID}{'BS'};#---[2023/07/21 12:34] the bound str
			my $chromStart = $set_info_hsh_ref->{$set_ID}{'CS'};
			my $chromEnd = $set_info_hsh_ref->{$set_ID}{'CE'};
			my $loc = $chrom.":".$chromStart."-".$chromEnd;
			my $junct_num = $bound_str =~ tr/J//;
			my $length = $chromEnd-$chromStart;
			my $full_partial = $setID_hsh_ref->{$set_ID};
			my $ref_count = 0;
			my $qry_count = 0;
		
			if (exists $set_info_hsh_ref->{$set_ID}{'T'}{'ref'}) {
				$ref_count = keys %{$set_info_hsh_ref->{$set_ID}{'T'}{'ref'}};
			}
			if (exists $set_info_hsh_ref->{$set_ID}{'T'}{'qry'}) {
				$qry_count = keys %{$set_info_hsh_ref->{$set_ID}{'T'}{'qry'}};
			}
			
			print INFO join "", (join "\t", ($model_ID, $set_ID, $loc, $end5_endtype, $end3_endtype, $full_partial, $junct_num, $length, $ref_count, $qry_count, $bound_str)), "\n";
		}
	}
	close INFO;

	return ();
}
sub printInfoSet {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $chrom, $log_ary_ref, $print_trnscrptID, $set_info_hsh_ref, $set_info_tsv
#	output: 
#	toCall: &printInfoSet($set_info_hsh_ref, $set_info_tsv, $chrom, $print_trnscrptID, $log_ary_ref);
#	calledInLine: 3348
#....................................................................................................................................................#
	my ($set_info_hsh_ref, $set_info_tsv, $chrom, $print_trnscrptID, $log_ary_ref) = @_;
	
	my $count_hsh_ref = {};
	
	foreach my $all_ref_qry (qw/all ref qry/) {
		foreach my $completeness (qw/complete incomplete/) {
			$count_hsh_ref->{$all_ref_qry}{'set'}{$completeness} = 0;
		}
		foreach my $end (qw/end5 end3/) {
			foreach my $confidence (qw/confident doubtful/) {
				$count_hsh_ref->{$all_ref_qry}{'set_end'}{$end}{$confidence} = 0;
			}
			foreach my $endtype (qw/C L S/) {
				$count_hsh_ref->{$all_ref_qry}{'set_endtype'}{$end}{$endtype} = 0;
			}
		}
	}
	
	open INFO, "| gzip -c >$set_info_tsv";
	#print INFO join "", (join "\t", ('set_ID', 'loc', 'completeness', 'represent_model', 'end5_conf', 'end3_conf', 'junct_conf', 'end5_endtype', 'end3_endtype', 'ref_count', 'qry_count', 'model_count', 'model_ID_str', 'strand', 'bound_str', 'ref_trnscrptID', 'qry_trnscrptID')), "\n";
	foreach my $set_ID (sort keys %{$set_info_hsh_ref}) {

		my $completeness = $set_info_hsh_ref->{$set_ID}{'CM'}; #---[2023/07/21 12:34] completeness: not complete
		my $represent_model = $set_info_hsh_ref->{$set_ID}{'AM'};#---[2023/07/21 12:34] to be tested as a model or not downstream
		my $bound_str = $set_info_hsh_ref->{$set_ID}{'BS'};
		my $strand = $set_info_hsh_ref->{$set_ID}{'ST'};
		my $chromStart = $set_info_hsh_ref->{$set_ID}{'CS'};
		my $chromEnd = $set_info_hsh_ref->{$set_ID}{'CE'};
		my $end3_endtype = $set_info_hsh_ref->{$set_ID}{'T3'};
		my $end5_endtype = $set_info_hsh_ref->{$set_ID}{'T5'};
		my $end5_conf = 'Y';
		my $end3_conf = 'Y';
		my $junct_conf = 'Y';

		$end5_conf = 'N' if $bound_str =~ m/XF/;
		$end3_conf = 'N' if $bound_str =~ m/XT/;
		$junct_conf = 'N' if $bound_str =~ m/XJ/;

		my $loc = $chrom.":".$chromStart."-".$chromEnd;
		my $ref_count = 0;
		my $qry_count = 0;
		my $model_ID_str = join ";", sort keys %{$set_info_hsh_ref->{$set_ID}{'model'}}; #---[2023/07/21 14:55] model ID
		my $model_count = keys %{$set_info_hsh_ref->{$set_ID}{'model'}}; #---[2023/07/21 14:55] model ID;
		my $ref_trnscrptID = '__na';
		my $qry_trnscrptID = '__na';
		
		if (exists $set_info_hsh_ref->{$set_ID}{'T'}{'ref'}) {
			$ref_count = keys %{$set_info_hsh_ref->{$set_ID}{'T'}{'ref'}};
			if ($print_trnscrptID eq 'yes') {
				$ref_trnscrptID = join ";", keys %{$set_info_hsh_ref->{$set_ID}{'T'}{'ref'}};
			};
 
		}
		if (exists $set_info_hsh_ref->{$set_ID}{'T'}{'qry'}) {
			$qry_count = keys %{$set_info_hsh_ref->{$set_ID}{'T'}{'qry'}};
			if ($print_trnscrptID eq 'yes') {
				$qry_trnscrptID = join ";", keys %{$set_info_hsh_ref->{$set_ID}{'T'}{'qry'}};
			};
		}
		
		print INFO join "", (join "\t", ($set_ID, $loc, $completeness, $represent_model, $end5_conf, $end3_conf, $junct_conf, $end5_endtype, $end3_endtype, $ref_count, $qry_count, $model_count, $model_ID_str, $strand, $bound_str, $ref_trnscrptID, $qry_trnscrptID)), "\n";
		
		if ($completeness eq 'Y') {
			$count_hsh_ref->{'all'}{'set'}{'complete'}++;
			$count_hsh_ref->{'qry'}{'set'}{'complete'}++ if $qry_count > 0;
			$count_hsh_ref->{'ref'}{'set'}{'complete'}++ if $ref_count > 0;
		} else {
			$count_hsh_ref->{'all'}{'set'}{'incomplete'}++;
			$count_hsh_ref->{'qry'}{'set'}{'incomplete'}++ if $qry_count > 0;
			$count_hsh_ref->{'ref'}{'set'}{'incomplete'}++ if $ref_count > 0;
		}

		$count_hsh_ref->{'all'}{'set_endtype'}{'end3'}{$end3_endtype}++;
		$count_hsh_ref->{'qry'}{'set_endtype'}{'end3'}{$end3_endtype}++ if $qry_count > 0;
		$count_hsh_ref->{'ref'}{'set_endtype'}{'end3'}{$end3_endtype}++ if $ref_count > 0;

		$count_hsh_ref->{'all'}{'set_endtype'}{'end5'}{$end5_endtype}++;
		$count_hsh_ref->{'qry'}{'set_endtype'}{'end5'}{$end5_endtype}++ if $qry_count > 0;
		$count_hsh_ref->{'ref'}{'set_endtype'}{'end5'}{$end5_endtype}++ if $ref_count > 0;
		
		
		if ($bound_str =~ /XT/) {
			$count_hsh_ref->{'all'}{'set_end'}{'end3'}{'doubtful'}++;
			$count_hsh_ref->{'qry'}{'set_end'}{'end3'}{'doubtful'}++ if $qry_count > 0;
			$count_hsh_ref->{'ref'}{'set_end'}{'end3'}{'doubtful'}++ if $ref_count > 0;
		} else {
			$count_hsh_ref->{'all'}{'set_end'}{'end3'}{'confident'}++;
			$count_hsh_ref->{'qry'}{'set_end'}{'end3'}{'confident'}++ if $qry_count > 0;
			$count_hsh_ref->{'ref'}{'set_end'}{'end3'}{'confident'}++ if $ref_count > 0;
		}
		if ($bound_str =~ /XF/) {
			$count_hsh_ref->{'all'}{'set_end'}{'end5'}{'doubtful'}++;
			$count_hsh_ref->{'qry'}{'set_end'}{'end5'}{'doubtful'}++ if $qry_count > 0;
			$count_hsh_ref->{'ref'}{'set_end'}{'end5'}{'doubtful'}++ if $ref_count > 0;
		} else {
			$count_hsh_ref->{'all'}{'set_end'}{'end5'}{'confident'}++;
			$count_hsh_ref->{'qry'}{'set_end'}{'end5'}{'confident'}++ if $qry_count > 0;
			$count_hsh_ref->{'ref'}{'set_end'}{'end5'}{'confident'}++ if $ref_count > 0;
		}

	}
	close INFO;

	my $endtype_prefix_hsh_ref = {
		'C' => 'commonest',
		'L' => 'longest',
		'S' => 'summit',
	};

	foreach my $all_ref_qry (qw/all ref qry/) {
		foreach my $completeness (qw/complete incomplete/) {
			my $count = $count_hsh_ref->{$all_ref_qry}{'set'}{$completeness};
			my $log_str = "$all_ref_qry set $completeness";
			push @{$log_ary_ref}, [$log_str, $count];
		}
		foreach my $end (qw/end5 end3/) {
			foreach my $confidence (qw/confident doubtful/) {
				my $count = $count_hsh_ref->{$all_ref_qry}{'set_end'}{$end}{$confidence};
				my $log_str = "$all_ref_qry set with $confidence $end";
				push @{$log_ary_ref}, [$log_str, $count];
			}
			foreach my $endtype (qw/C L S/) {
				my $full_endtype = $endtype_prefix_hsh_ref->{$endtype};
				my $count = $count_hsh_ref->{$all_ref_qry}{'set_endtype'}{$end}{$endtype};
				my $log_str = "$all_ref_qry set $end adjusted by $full_endtype";
				push @{$log_ary_ref}, [$log_str, $count];
			}
		}
	}

	return ();
}
sub printInfoTranscript {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $chrom, $log_ary_ref, $set_info_hsh_ref, $trnscpt_info_hsh_ref, $trnscpt_info_tsv
#	output: 
#	toCall: &printInfoTranscript($set_info_hsh_ref, $trnscpt_info_hsh_ref, $trnscpt_info_tsv, $chrom, $log_ary_ref);
#	calledInLine: 3345
#....................................................................................................................................................#
	my ($set_info_hsh_ref, $trnscpt_info_hsh_ref, $trnscpt_info_tsv, $chrom, $log_ary_ref) = @_;
	my $count_hsh_ref = {};
	
	foreach my $ref_qry (qw/ref qry/) {
		foreach my $completeness (qw/complete incomplete/) {
			$count_hsh_ref->{$ref_qry}{'trnscpt'}{$completeness} = 0;
		}
		foreach my $end (qw/end5 end3/) {
			foreach my $confidence (qw/confident doubtful/) {
				$count_hsh_ref->{$ref_qry}{'trnscpt_end'}{$end}{$confidence} = 0;
			}
		}
	}
	open TRNSCPTINFO, "| gzip -c >$trnscpt_info_tsv";
	#print TRNSCPTINFO join "", (join "\t", ('trnscpt_ID', 'ref_qry', 'completeness', 'set_ID', 'end5_ID', 'end3_ID', 'model_ID_str', 'chrom', 'trnscpt_chromStart', 'trnscpt_chromEnd', 'set_chromStart', 'set_chromEnd', 'strand', 'bound_str')), "\n";
	foreach my $set_ID (sort keys %{$set_info_hsh_ref}) {
		my $completeness = $set_info_hsh_ref->{$set_ID}{'CM'}; #---[2023/07/21 12:34] completeness: not complete
		my $bound_str = $set_info_hsh_ref->{$set_ID}{'BS'};
		my $set_chromStart = $set_info_hsh_ref->{$set_ID}{'CS'};
		my $set_chromEnd = $set_info_hsh_ref->{$set_ID}{'CE'};
		my @bound_str_ary = split /_/, $bound_str; #---[2023/07/21 12:34] set string
		my $end5_ID = $bound_str_ary[0];
		my $end3_ID = $bound_str_ary[-2];
		my $model_ID_str = join ";", sort keys %{$set_info_hsh_ref->{$set_ID}{'model'}}; #---[2023/07/21 14:55] model ID
		foreach my $ref_qry (sort keys %{$set_info_hsh_ref->{$set_ID}{'T'}}) {
			foreach my $trnscpt_ID (keys %{$set_info_hsh_ref->{$set_ID}{'T'}{$ref_qry}}) {
				my ($strand, $trnscpt_chromStart, $trnscpt_chromEnd) = @{$trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'I'}};
				print TRNSCPTINFO join "", (join "\t", ($trnscpt_ID, $ref_qry, $completeness, $set_ID, $end5_ID, $end3_ID, $model_ID_str, $chrom, $trnscpt_chromStart, $trnscpt_chromEnd, $set_chromStart, $set_chromEnd, $strand, $bound_str)), "\n";
				
				if ($completeness eq 'Y') {
					$count_hsh_ref->{$ref_qry}{'trnscpt'}{'complete'}++;
				} else {
					$count_hsh_ref->{$ref_qry}{'trnscpt'}{'incomplete'}++;
				}

				if ($end3_ID =~ /^X/) {
					$count_hsh_ref->{$ref_qry}{'trnscpt_end'}{'end3'}{'doubtful'}++;
				} else {
					$count_hsh_ref->{$ref_qry}{'trnscpt_end'}{'end3'}{'confident'}++;
				}
				if ($end5_ID =~ /^X/) {
					$count_hsh_ref->{$ref_qry}{'trnscpt_end'}{'end5'}{'doubtful'}++;
				} else {
					$count_hsh_ref->{$ref_qry}{'trnscpt_end'}{'end5'}{'confident'}++;
				}
			}
		}
	}
	close TRNSCPTINFO;
	
	foreach my $ref_qry (qw/ref qry/) {
		foreach my $completeness (qw/complete incomplete/) {
			my $count = $count_hsh_ref->{$ref_qry}{'trnscpt'}{$completeness};
			my $log_str = "$ref_qry transcript $completeness";
			push @{$log_ary_ref}, [$log_str, $count];
		}
		foreach my $end (qw/end5 end3/) {
			foreach my $confidence (qw/confident doubtful/) {
				my $count = $count_hsh_ref->{$ref_qry}{'trnscpt_end'}{$end}{$confidence};
				my $log_str = "$ref_qry transcript with $confidence $end";
				push @{$log_ary_ref}, [$log_str, $count];
			}
		}
	}
	
	return ();
}
sub printLogTxt {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $log_ary_ref, $log_txt
#	output: 
#	toCall: &printLogTxt($log_txt, $log_ary_ref);
#	calledInLine: 3366
#....................................................................................................................................................#
	my ($log_txt, $log_ary_ref) = @_;
	
	open LOG, ">", $log_txt;
	foreach my $index (0..$#{$log_ary_ref}) {
		my ($log_str, $count) = @{$log_ary_ref->[$index]};
		print LOG join "", (join "\t", ($index, $log_str, $count)), "\n";
	}
	close LOG;

	return ();
}
sub printModelBed {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $bgzip_bin, $model_bed, $model_info_hsh_ref, $set_bed, $set_info_hsh_ref
#	output: 
#	toCall: &printModelBed($set_bed, $model_bed, $model_info_hsh_ref, $set_info_hsh_ref, $bgzip_bin);
#	calledInLine: 3357
#....................................................................................................................................................#
	my ($set_bed, $model_bed, $model_info_hsh_ref, $set_info_hsh_ref, $bgzip_bin) = @_;
	
	my $set_to_model_hsh_ref = {};
	foreach my $model_ID (sort keys %{$model_info_hsh_ref}) {
		$set_to_model_hsh_ref->{$model_info_hsh_ref->{$model_ID}{'set'}} = $model_ID;#---[2023/07/21 12:34] the setID
	}
	open MODELBED, "| $bgzip_bin -c >$model_bed";
	open SETBED, "$bgzip_bin -dc $set_bed |";
	my $model_printed_hsh_ref = {};
	while (<SETBED>) {
		chomp;
		my ($chrom, $chromStart, $chromEnd, $full_set_ID, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\t/;
		next if not exists $set_to_model_hsh_ref->{$full_set_ID};
		my $model_ID = $set_to_model_hsh_ref->{$full_set_ID};

		my $full_ref_count = 0;
		my $full_qry_count = 0;
		my $partial_ref_count = 0;
		my $partial_qry_count = 0;
		my $partial_set_count = 0;
		my $partial_set_ID_str = '__na';
		
		if (exists $set_info_hsh_ref->{$full_set_ID}{'T'}{'ref'}) {
			$full_ref_count = keys %{$set_info_hsh_ref->{$full_set_ID}{'T'}{'ref'}};
		}
		if (exists $set_info_hsh_ref->{$full_set_ID}{'T'}{'qry'}) {
			$full_qry_count = keys %{$set_info_hsh_ref->{$full_set_ID}{'T'}{'qry'}};
		}

		if (exists $model_info_hsh_ref->{$model_ID}{'IS'}) {
			$partial_set_count = keys %{$model_info_hsh_ref->{$model_ID}{'IS'}}; #---[2023/07/21 14:55] partial set count
			$partial_set_ID_str = join ";", sort keys %{$model_info_hsh_ref->{$model_ID}{'IS'}}; #---[2023/07/21 14:55] set ID
			foreach my $partial_set_ID (keys %{$model_info_hsh_ref->{$model_ID}{'IS'}}) {
				if (exists $set_info_hsh_ref->{$partial_set_ID}{'T'}{'ref'}) {
					$partial_ref_count = keys %{$set_info_hsh_ref->{$partial_set_ID}{'T'}{'ref'}};
				}
				if (exists $set_info_hsh_ref->{$partial_set_ID}{'T'}{'qry'}) {
					$partial_qry_count = keys %{$set_info_hsh_ref->{$partial_set_ID}{'T'}{'qry'}};
				}
			}
		}
		my $end3_endtype = $set_info_hsh_ref->{$full_set_ID}{'T3'};
		my $end5_endtype = $set_info_hsh_ref->{$full_set_ID}{'T5'};
		my $total_set_count = $partial_set_count + 1;
		my $total_qry_count = $full_qry_count + $partial_qry_count;
		my $ID_string = join "|", ($model_ID, "end=$end5_endtype:$end3_endtype", "set=$total_set_count", "ref=$full_ref_count:$partial_ref_count", "qry=$full_qry_count:$partial_qry_count");
		print MODELBED join "", (join "\t", ($chrom, $chromStart, $chromEnd, $model_ID, $full_qry_count, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts)), "\n";
		
	}
	close SETBED;
	close MODELBED;

	return ();
}
sub printOutputFileListAndReadme {
#....................................................................................................................................................#
#	subroutineCategory: output
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 4_finishingTasks|214
#	secondaryAppearInSection: >none
#	input: $ARGVStr, $out_dir, $paramTag
#	output: 
#	toCall: &printOutputFileListAndReadme($ARGVStr, $paramTag, $out_dir);
#	calledInLine: 217
#....................................................................................................................................................#
	my ($ARGVStr, $paramTag, $out_dir) = @_;
	
	my $outputFileListPath = "$out_dir/$paramTag/output.file.list.txt";
	open (OUTFILELIST, ">", $outputFileListPath);

	my %dirHsh = ();
	my %filelistLenCountHsh = ();
	push @{$filelistLenCountHsh{'dir'}}, length 'Directory';
	push @{$filelistLenCountHsh{'name'}}, length 'Name';
	push @{$filelistLenCountHsh{'description'}}, length 'Description';
	
	foreach my $outputFilePath (sort {$a cmp $b} keys %{$globalReadmeHsh_ref}) {
		my $fileDescription =  $globalReadmeHsh_ref->{$outputFilePath}{'description'};
		my $cleandOutputFilePath = $outputFilePath;
		$cleandOutputFilePath =~ s/\/+/\//g;
		
		my ($filePrefix, $fileDir, $fileSuffix) = fileparse($cleandOutputFilePath, qr/\.[^.]*/);
		$fileDir =~ s/^$out_dir//;
		my $fileName = $filePrefix.$fileSuffix;
		$dirHsh{$fileDir}{$fileName} = $fileDescription;
		push @{$filelistLenCountHsh{'dir'}}, length $fileDir;
		push @{$filelistLenCountHsh{'name'}}, length $fileName;
		push @{$filelistLenCountHsh{'description'}}, length $fileDescription;
		
		open README, ">", "$outputFilePath.readme.txt";
		print README "=================\n";
		print README "File descriptions\n";
		print README "=================\n";
		print README "$fileDescription\n";
					
		if (exists $globalReadmeHsh_ref->{$outputFilePath}{'headerAry'}) {
			my @colLenCountHsh = (length 'column');
			push @colLenCountHsh, length $_ foreach (@{$globalReadmeHsh_ref->{$outputFilePath}{'headerAry'}});
			my $headerColLen = max(@colLenCountHsh)+2;
			print README "\n";
			print README "\n";
			print README "===================\n";
			print README "Column descriptions\n";
			print README "===================\n";
			print README "\n";
			printf README "%-".$headerColLen."s", 'column';
			print README "description\n";
			printf README "%-".$headerColLen."s", '------';
			print README "-----------\n";
			foreach my $header (@{$globalReadmeHsh_ref->{$outputFilePath}{'headerAry'}}) {
				my $columnDescription = 'self-explanatory';
				$columnDescription = $globalReadmeHsh_ref->{$outputFilePath}{'header'}{$header} if exists $globalReadmeHsh_ref->{$outputFilePath}{'header'}{$header};
				printf README "%-".$headerColLen."s", $header;
				print README $columnDescription."\n";
			}
		}
		
		if (exists $globalReadmeHsh_ref->{$outputFilePath}{'extra_info'}) {
			print README "\n";
			print README "\n";
			print README "=================\n";
			print README "Extra information\n";
			print README "=================\n";
			print README "\n";
			foreach my $title (sort keys %{$globalReadmeHsh_ref->{$outputFilePath}{'extra_info'}}) {
				print README "$title\n";
				print README "-" foreach (1..length $title);
				print README "\n";
				print README "$_\n" foreach @{$globalReadmeHsh_ref->{$outputFilePath}{'extra_info'}{$title}};
			}
		}
		
		print README "\n";
		print README "\n";
		print README "~" foreach (1..length "$fileName was created from running,");
		print README "\n";
		print README "$fileName was created from running,\n";
		print README "\n";
		print README "$ARGVStr\n";
		print README "\n";
		close README;
	}

	my $fileDir_colLen = max(@{$filelistLenCountHsh{'dir'}})+2;
	my $fileName_colLen = max(@{$filelistLenCountHsh{'name'}})+2;
	my $fileDescription_colLen = max(@{$filelistLenCountHsh{'description'}})+2;
	printf OUTFILELIST ("%-".$fileDir_colLen."s %-".$fileName_colLen."s %-".$fileDescription_colLen."s\n", 'directory', 'name', 'description');
	printf OUTFILELIST ("%-".$fileDir_colLen."s %-".$fileName_colLen."s %-".$fileDescription_colLen."s\n", '=========', '====', '===========');
	foreach my $fileDir (sort {$a cmp $b} keys %dirHsh) {
		foreach my $fileName (sort {$a cmp $b} keys %{$dirHsh{$fileDir}}) {
			my $fileDescription = $dirHsh{$fileDir}{$fileName};	
			printf OUTFILELIST ("%-".$fileDir_colLen."s %-".$fileName_colLen."s %-".$fileDescription_colLen."s\n", $fileDir, $fileName, $fileDescription);
		}
	}
	
	print OUTFILELIST "\n";
	print OUTFILELIST "\n";
	print OUTFILELIST "~" foreach (1..length "The above files were generated by running,");
	print OUTFILELIST "\n";
	print OUTFILELIST "The above files were generated by running,\n";
	print OUTFILELIST "\n";
	print OUTFILELIST "$ARGVStr\n";
	print OUTFILELIST "\n";

	close OUTFILELIST;

	return ();
}
sub printStartOrFinishMessage {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|871
#	appearInSub: >none
#	primaryAppearInSection: 2_defineout_dirPath|153, 4_finishingTasks|214
#	secondaryAppearInSection: >none
#	input: $StartOrFinishMessage
#	output: none
#	toCall: &printStartOrFinishMessage($StartOrFinishMessage);
#	calledInLine: 166, 218
#....................................................................................................................................................#

	my ($StartOrFinishMessage) = @_;
	
	if ($StartOrFinishMessage eq "startMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->871
		print "=========================================================================\n\n";

		print $tmplog_fh "\n=========================================================================\n";
		print $tmplog_fh "[".&currentTime()."] starts running ...... \n";#->871
		print $tmplog_fh "=========================================================================\n\n";

	} elsif ($StartOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->871
		print "=========================================================================\n\n";

		print $tmplog_fh "\n=========================================================================\n";
		print $tmplog_fh "[".&currentTime()."] finished running .......\n";#->871
		print $tmplog_fh "=========================================================================\n\n";
	}
}
sub printTranscriptBed {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $bgzip_bin, $chrom, $input_bed_hsh_ref, $qry_trnscpt_bed, $ref_trnscpt_bed, $report_tag, $set_info_hsh_ref, $tabix_bin, $trnscpt_info_hsh_ref
#	output: 
#	toCall: &printTranscriptBed($set_info_hsh_ref, $trnscpt_info_hsh_ref, $ref_trnscpt_bed, $qry_trnscpt_bed, $chrom, $input_bed_hsh_ref, $report_tag, $tabix_bin, $bgzip_bin);
#	calledInLine: 3360
#....................................................................................................................................................#
	my ($set_info_hsh_ref, $trnscpt_info_hsh_ref, $ref_trnscpt_bed, $qry_trnscpt_bed, $chrom, $input_bed_hsh_ref, $report_tag, $tabix_bin, $bgzip_bin) = @_;

	my $ref_qry_hsh_ref = {
		'ref' => {
			'in_bed' => [],
			'out_bed' => $ref_trnscpt_bed,
		},
		'qry' => {
			'in_bed' => [],
			'out_bed' => $qry_trnscpt_bed,
		},
	};
	foreach my $input_ID (keys %{$input_bed_hsh_ref}) {
		&reportAndLogStatus("$report_tag Printing annotated transcript bed in $input_ID", 10, "\n");#->3614
		my $num_proc = 0;
		my $ref_qry = $input_bed_hsh_ref->{$input_ID}{'ref_qry'};
		my $trnscpt_bed = $input_bed_hsh_ref->{$input_ID}{'trnscpt_bed'};
		my $anno_trnscpt_bed = $input_bed_hsh_ref->{$input_ID}{'anno_trnscpt_bed'};
		push @{$ref_qry_hsh_ref->{$ref_qry}{'in_bed'}}, $anno_trnscpt_bed;
		
		my $tabix_cmd  = "$tabix_bin $trnscpt_bed $chrom |";
		open ANNOTRNSCPTBED, "| $bgzip_bin -c >$anno_trnscpt_bed";
		open TABIX, "$tabix_cmd";
		while (<TABIX>) {
			chomp;
			my (undef, $chromStart, $chromEnd, $trnscpt_ID, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\t/;
			if (exists $trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}) {
				$num_proc++;
				&reportAndLogStatus("$report_tag Printing annotated transcript bed >> $num_proc printed in $input_ID", 10, "\n") if $num_proc%100000 ==0;#->3614
				my $set_ID = $trnscpt_info_hsh_ref->{$ref_qry}{$trnscpt_ID}{'S'};
				my $model_ID_str = join ";", sort keys %{$set_info_hsh_ref->{$set_ID}{'model'}}; #---[2023/07/21 14:55] model ID
				my $ID_str = join "|", ($trnscpt_ID, $set_ID, $model_ID_str);
				print ANNOTRNSCPTBED join "", (join "\t", ($chrom, $chromStart, $chromEnd, $ID_str, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts)), "\n";
			}
		}
		&reportAndLogStatus("$report_tag Finished printing annotated transcript bed ## $num_proc printed in $input_ID", 10, "\n") if $num_proc%100000 ==0;#->3614
	}
	
	foreach my $ref_qry (keys %{$ref_qry_hsh_ref}) {
		&reportAndLogStatus("$report_tag Sorting annotated $ref_qry transcript bed", 10, "\n");#->3614
		my $out_bed = $ref_qry_hsh_ref->{$ref_qry}{'out_bed'};
		my @in_bed_ary = @{$ref_qry_hsh_ref->{$ref_qry}{'in_bed'}};
		if (@in_bed_ary == 1) {
			system ("mv $in_bed_ary[0] $out_bed");
		} else {
			my $in_bed_str = join " ", @in_bed_ary;
			system ("cat $in_bed_str | $bgzip_bin -dc | sort -k2,2n -k6,6 | $bgzip_bin -c >$out_bed");
			foreach my $in_bed (@in_bed_ary) {
				system ("rm $in_bed");
			}
		}
	}
	
	return ();
}
sub processPerChromosome {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: assignConfidentEnds|259, assignDoubtfulEnds|319, checkRefTrncptIDNumbering|753, defineEndRegion|994, defineTrnscptModels|1382, defineTrnscptSets|1431, extractTranscriptBounds|1571, getJunctSeq|1817, printInfoBedBothEnd|2343, printInfoBedJunction|2442, printInfoModel|2482, printInfoModelSetPair|2634, printInfoSet|2688, printInfoTranscript|2832, printLogTxt|2911, printModelBed|2934, printTranscriptBed|3149, readConfJunction|3399, reportAndLogStatus|3614
#	appearInSub: >none
#	primaryAppearInSection: 3_process|197
#	secondaryAppearInSection: >none
#	input: $bedtools_bin, $bgzip_bin, $chrom_fasta_path, $chrom_info_hsh_ref, $chrom_size_path, $conf_end3_add_ref, $conf_end3_bed_bgz, $conf_end3_merge_flank, $conf_end5_add_ref, $conf_end5_bed_bgz, $conf_end5_merge_flank, $conf_junction_bed, $doubtful_end_avoid_summit, $doubtful_end_merge_dist, $max_thread, $min_exon_length, $min_frac_split, $min_size_split, $min_summit_dist_split, $min_transcript_length, $novel_model_prefix, $print_trnscrptID, $result_tmp_dir, $retain_no_qry_ref_bound_set, $signal_end3_bed_bgz, $signal_end5_bed_bgz, $tabix_bin, $trnscpt_set_end_priority_hsh_ref, $use_ref_only_end_pos
#	output: 
#	toCall: &processPerChromosome($chrom_size_path, $novel_model_prefix, $conf_end3_merge_flank, $conf_end5_merge_flank, $doubtful_end_merge_dist, $doubtful_end_avoid_summit, $conf_end5_bed_bgz, $conf_end3_bed_bgz, $signal_end5_bed_bgz, $signal_end3_bed_bgz, $conf_junction_bed, $conf_end3_add_ref, $conf_end5_add_ref, $bedtools_bin, $tabix_bin, $bgzip_bin, $max_thread, $result_tmp_dir, $chrom_info_hsh_ref, $trnscpt_set_end_priority_hsh_ref, $min_transcript_length, $min_exon_length, $print_trnscrptID, $min_size_split, $min_frac_split, $min_summit_dist_split, $use_ref_only_end_pos, $retain_no_qry_ref_bound_set, $chrom_fasta_path);
#	calledInLine: 202
#....................................................................................................................................................#
	my ($chrom_size_path, $novel_model_prefix, $conf_end3_merge_flank, $conf_end5_merge_flank, $doubtful_end_merge_dist, $doubtful_end_avoid_summit, $conf_end5_bed_bgz, $conf_end3_bed_bgz, $signal_end5_bed_bgz, $signal_end3_bed_bgz, $conf_junction_bed, $conf_end3_add_ref, $conf_end5_add_ref, $bedtools_bin, $tabix_bin, $bgzip_bin, $max_thread, $result_tmp_dir, $chrom_info_hsh_ref, $trnscpt_set_end_priority_hsh_ref, $min_transcript_length, $min_exon_length, $print_trnscrptID, $min_size_split, $min_frac_split, $min_summit_dist_split, $use_ref_only_end_pos, $retain_no_qry_ref_bound_set, $chrom_fasta_path) = @_;
	
	my %chromForThrHsh = ();
	my $threadID = 1;
	foreach my $chrom (sort {$chrom_info_hsh_ref->{$b}{'size'} <=> $chrom_info_hsh_ref->{$a}{'size'}} keys %{$chrom_info_hsh_ref}) {
		$threadID = 1 if $threadID > $max_thread;
		$threadID = sprintf "%.2d", $threadID;
		push @{$chromForThrHsh{$threadID}} , $chrom;
		$threadID++;
	}
	
	my %threadHsh =();
	foreach my $threadID (sort {$a <=> $b} keys %chromForThrHsh) {
		my $chromForThrAry_ref = $chromForThrHsh{$threadID};
		($threadHsh{$threadID}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781, the 
	
			sub {
				my ($chromForThrAry_ref) = @_;

				foreach my $chrom (@{$chromForThrAry_ref}) {
					#---[2/19/16 13:13] do something on the item here
					my $chrom_num = $chrom_info_hsh_ref->{$chrom}{'chrom_num'};
					my $size = $chrom_info_hsh_ref->{$chrom}{'size'};
					my $confident_end5_bed = $chrom_info_hsh_ref->{$chrom}{'confident_end5_bed'};
					my $confident_end3_bed = $chrom_info_hsh_ref->{$chrom}{'confident_end3_bed'};
					my $doubtful_end5_bed = $chrom_info_hsh_ref->{$chrom}{'doubtful_end5_bed'};
					my $doubtful_end3_bed = $chrom_info_hsh_ref->{$chrom}{'doubtful_end3_bed'};
					my $set_bed = $chrom_info_hsh_ref->{$chrom}{'set_bed'};
					my $model_bed = $chrom_info_hsh_ref->{$chrom}{'model_bed'};
					my $tmp_txt = $chrom_info_hsh_ref->{$chrom}{'tmp_txt'};
					my $all_end5_bed = $chrom_info_hsh_ref->{$chrom}{'all_end5_bed'};
					my $all_end3_bed = $chrom_info_hsh_ref->{$chrom}{'all_end3_bed'};
					my $input_bed_hsh_ref = $chrom_info_hsh_ref->{$chrom}{'input_bed_hsh_ref'};
					my $junct_bed = $chrom_info_hsh_ref->{$chrom}{'junct_bed'};
					my $ref_trnscpt_bed = $chrom_info_hsh_ref->{$chrom}{'ref_trnscpt_bed'};
					my $qry_trnscpt_bed = $chrom_info_hsh_ref->{$chrom}{'qry_trnscpt_bed'};
					my $discard_trnscpt_bed = $chrom_info_hsh_ref->{$chrom}{'discard_trnscpt_bed'};
					my $splicing_site_bed = $chrom_info_hsh_ref->{$chrom}{'splicing_site_bed'};
					my $splicing_site_fasta = $chrom_info_hsh_ref->{$chrom}{'splicing_site_fasta'};

					my $report_tag = join " ", (sprintf("%-15s","Thread-$threadID $chrom"), ":");

					&reportAndLogStatus("$report_tag >>>>>> Starting $chrom analysis <<<<<<", 10, "\n");#->3614

					my $trnscpt_info_tsv = $chrom_info_hsh_ref->{$chrom}{'trnscpt_info_tsv'};
					my $model_info_tsv = $chrom_info_hsh_ref->{$chrom}{'model_info_tsv'};
					my $set_info_tsv = $chrom_info_hsh_ref->{$chrom}{'set_info_tsv'};
					my $end5_info_tsv = $chrom_info_hsh_ref->{$chrom}{'end5_info_tsv'};
					my $end3_info_tsv = $chrom_info_hsh_ref->{$chrom}{'end3_info_tsv'};
					my $junct_info_tsv = $chrom_info_hsh_ref->{$chrom}{'junct_info_tsv'};
					my $model_set_pair_info_tsv = $chrom_info_hsh_ref->{$chrom}{'model_set_pair_info_tsv'};
					my $log_txt = $chrom_info_hsh_ref->{$chrom}{'log_txt'};
					my @check_file_ary = (
						$junct_bed,
						$set_bed,
						$model_bed,
						$all_end5_bed,
						$all_end3_bed,
						$ref_trnscpt_bed,
						$qry_trnscpt_bed,
						$discard_trnscpt_bed,
						$trnscpt_info_tsv,
						$model_info_tsv,
						$set_info_tsv,
						$end5_info_tsv,
						$end3_info_tsv,
						$junct_info_tsv,
						$model_set_pair_info_tsv,
						$log_txt,
					);
					my $skip = 'yes';
					foreach my $file_to_check (@check_file_ary) {
						$skip = 'no' if not -s $file_to_check
					}
					
					if ($skip eq 'yes') {
						&reportAndLogStatus("$report_tag >>>>>> Results found skipping $chrom <<<<<<", 10, "\n");#->3614
						next;
					}

					#system ("rm -f $log_txt"); #---[2023/07/25 20:10] log as the sign for finish running the whole process
					my $novel_model_chrom_prefix = $novel_model_prefix.$chrom_num;
					my $log_ary_ref = [];
					
					&reportAndLogStatus("$report_tag Define confident end regions", 10, "\n");#->3614
					my ($end5_info_hsh_ref, $end3_info_hsh_ref) = &defineEndRegion($chrom_info_hsh_ref, $input_bed_hsh_ref, $conf_end5_add_ref, $conf_end3_add_ref, $bedtools_bin, $tabix_bin, $bgzip_bin, $conf_end5_merge_flank, $conf_end3_merge_flank, $conf_end5_bed_bgz, $conf_end3_bed_bgz, $confident_end5_bed, $confident_end3_bed, $signal_end5_bed_bgz, $signal_end3_bed_bgz, $chrom_size_path, $chrom, $chrom_num, $threadID, $log_ary_ref, $size, $min_size_split, $min_frac_split, $min_summit_dist_split, $use_ref_only_end_pos, $report_tag);#->994

					&reportAndLogStatus("$report_tag Reading confident junctions", 10, "\n");#->3614
					my ($conf_junct_hsh_ref) = &readConfJunction($conf_junction_bed, $chrom, $report_tag);#->3399

					&reportAndLogStatus("$report_tag Getting confident junctions", 10, "\n");#->3614
					my ($junct_seq_hsh_ref) = &getJunctSeq($tabix_bin, $chrom_num, $threadID, $chrom, $input_bed_hsh_ref, $splicing_site_bed, $splicing_site_fasta, $chrom_fasta_path, $bedtools_bin, $report_tag);#->1817
					
					&reportAndLogStatus("$report_tag Extracting transcript bounds", 10, "\n");#->3614
					my ($junct_info_hsh_ref, $trnscpt_info_hsh_ref) = &extractTranscriptBounds($tabix_bin, $chrom_num, $threadID, $chrom, $input_bed_hsh_ref, $report_tag, $min_transcript_length, $min_exon_length, $log_ary_ref, $discard_trnscpt_bed, $bgzip_bin, $conf_junct_hsh_ref, $junct_seq_hsh_ref);#->1571
					undef $junct_seq_hsh_ref;
					
					&reportAndLogStatus("$report_tag Assigning confident transcript ends", 10, "\n");#->3614
					&assignConfidentEnds($tabix_bin, $end5_info_hsh_ref, $end3_info_hsh_ref, $bedtools_bin, $trnscpt_info_hsh_ref, $input_bed_hsh_ref, $confident_end5_bed, $confident_end3_bed, $threadID, $chrom, $conf_end5_bed_bgz, $report_tag);#->259

					&reportAndLogStatus("$report_tag Assigning doubtful transcript ends", 10, "\n");#->3614
					&assignDoubtfulEnds($bedtools_bin, $bgzip_bin, $trnscpt_info_hsh_ref, $doubtful_end_merge_dist, $end5_info_hsh_ref, $end3_info_hsh_ref, $report_tag, $doubtful_end5_bed, $doubtful_end3_bed, $chrom, $chrom_num, $log_ary_ref);#->319

					&reportAndLogStatus("$report_tag Defining transcript sets", 10, "\n");#->3614
					my ($set_info_hsh_ref) = &defineTrnscptSets($trnscpt_info_hsh_ref, $junct_info_hsh_ref, $chrom_num, $end5_info_hsh_ref, $end3_info_hsh_ref, $report_tag, $trnscpt_set_end_priority_hsh_ref, $set_bed, $min_transcript_length, $chrom, $bgzip_bin, $doubtful_end_avoid_summit, $retain_no_qry_ref_bound_set);#->1431
					
					&reportAndLogStatus("$report_tag Checking refTrncptID numbering", 10, "\n");#->3614
					my ($inital_novel_model_num) = &checkRefTrncptIDNumbering($chrom_num, $novel_model_chrom_prefix, $trnscpt_info_hsh_ref);#->753
					&reportAndLogStatus("$report_tag inital_novel_model_num=$inital_novel_model_num", 10, "\n");#->3614

					&reportAndLogStatus("$report_tag Printing junction info", 10, "\n");#->3614
					&printInfoBedJunction($bgzip_bin, $junct_info_hsh_ref, $junct_info_tsv, $junct_bed, $chrom);#->2442
					undef $junct_info_hsh_ref;
					
					&reportAndLogStatus("$report_tag Defining transcript models", 10, "\n");#->3614
					my ($model_info_hsh_ref) = &defineTrnscptModels($end5_info_hsh_ref, $set_info_hsh_ref, $novel_model_chrom_prefix, $inital_novel_model_num, $report_tag);#->1382

					&reportAndLogStatus("$report_tag Printing transcript info", 10, "\n");#->3614
					&printInfoTranscript($set_info_hsh_ref, $trnscpt_info_hsh_ref, $trnscpt_info_tsv, $chrom, $log_ary_ref);#->2832

					&reportAndLogStatus("$report_tag Printing set info", 10, "\n");#->3614
					&printInfoSet($set_info_hsh_ref, $set_info_tsv, $chrom, $print_trnscrptID, $log_ary_ref);#->2688

					&reportAndLogStatus("$report_tag Printing model info", 10, "\n");#->3614
					&printInfoModel($model_info_hsh_ref, $model_info_tsv, $set_info_hsh_ref, $chrom, $log_ary_ref);#->2482

					&reportAndLogStatus("$report_tag Printing model vs set info", 10, "\n");#->3614
					&printInfoModelSetPair($model_info_hsh_ref, $model_set_pair_info_tsv, $set_info_hsh_ref, $chrom);#->2634

					&reportAndLogStatus("$report_tag Printing model and set bed", 10, "\n");#->3614
					&printModelBed($set_bed, $model_bed, $model_info_hsh_ref, $set_info_hsh_ref, $bgzip_bin);#->2934

					&reportAndLogStatus("$report_tag Printing transcript bed", 10, "\n");#->2481					#->3614
					&printTranscriptBed($set_info_hsh_ref, $trnscpt_info_hsh_ref, $ref_trnscpt_bed, $qry_trnscpt_bed, $chrom, $input_bed_hsh_ref, $report_tag, $tabix_bin, $bgzip_bin);#->3149

					&reportAndLogStatus("$report_tag Printing end txt and bed", 10, "\n");#->3614
					&printInfoBedBothEnd($set_info_hsh_ref, $end3_info_tsv, $all_end3_bed, $end3_info_hsh_ref, $end5_info_tsv, $all_end5_bed, $end5_info_hsh_ref, $chrom, $report_tag, $bgzip_bin);#->2343

					&reportAndLogStatus("$report_tag Printing log txt", 10, "\n");#->3614
					&printLogTxt($log_txt, $log_ary_ref);#->2911
				}
				
				return ();
			}
			,($chromForThrAry_ref)
		);
	}
	
	my $data_hsh_ref = {};
	while (keys %threadHsh) {
		my $num_thread = keys %threadHsh;
		&reportAndLogStatus("--------------------------- $num_thread thread running ---------------------------", 10, "\n");#->3614
		foreach my $threadID (keys %threadHsh) {
			if (not $threadHsh{$threadID}->is_running()) {
				$threadHsh{$threadID}->join();
				foreach my $chrom (@{$chromForThrHsh{$threadID}}) {
					my $log_txt = $chrom_info_hsh_ref->{$chrom}{'log_txt'};
					if (-s $log_txt) {
						&reportAndLogStatus("Thread-$threadID : Log of $chrom found and finished successfully", 10, "\n");#->3614
					} else {
						&reportAndLogStatus("!!!!!!!!WARNING!!!!!!!! log_txt of $chrom is not found", 10, "\n");#->3614
					}
				}
				delete $threadHsh{$threadID};
				&reportAndLogStatus(">>>>>>>>>>>>>>>>>>>>>>>> thread $threadID finished <<<<<<<<<<<<<<<<<<<<<<<<<<<", 10, "\n");#->3614
			}
		}
		sleep 5;
	}
	return ();
}
sub readConfJunction {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: processPerChromosome|3216
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|197
#	input: $chrom, $conf_junction_bed, $report_tag
#	output: $conf_junct_hsh_ref
#	toCall: my ($conf_junct_hsh_ref) = &readConfJunction($conf_junction_bed, $chrom, $report_tag);
#	calledInLine: 3315
#....................................................................................................................................................#
	my ($conf_junction_bed, $chrom, $report_tag) = @_;

	my $conf_junct_hsh_ref = {};

	foreach my $indiv_conf_junction_bed (split /\,/, $conf_junction_bed) {
		if ($indiv_conf_junction_bed =~ m/\.gz$/) {
			open (FILEIN, " gzip -dc $indiv_conf_junction_bed|");
		} else {
			open (FILEIN, "<", $indiv_conf_junction_bed);
		}

		while (<FILEIN>) {
			chomp;
			my ($junct_chrom, $chromStart, $chromEnd, $ID, $score, $strand) = split /\t/;
			if ($junct_chrom eq $chrom) {
				my $junct_str = join ":", ($chromStart, $chromEnd, $strand);
				$conf_junct_hsh_ref->{$junct_str}++;
			}
		}
		close FILEIN;
	}
	
	my $num_conf_junct = keys %{$conf_junct_hsh_ref};
	&reportAndLogStatus("$report_tag $num_conf_junct confident junctions stored", 10, "\n");#->3614

	return ($conf_junct_hsh_ref);
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|137
#	secondaryAppearInSection: >none
#	input: none
#	output: $bedtools_bin, $bgzip_bin, $chrom_fasta_path, $chrom_size_path, $conf_end3_add_ref, $conf_end3_bed_bgz, $conf_end3_merge_flank, $conf_end5_add_ref, $conf_end5_bed_bgz, $conf_end5_merge_flank, $conf_junction_bed, $doubtful_end_avoid_summit, $doubtful_end_merge_dist, $max_thread, $min_exon_length, $min_frac_split, $min_output_qry_count, $min_qry_score, $min_size_split, $min_summit_dist_split, $min_transcript_length, $novel_model_prefix, $out_dir, $out_prefix, $print_trnscrptID, $qry_bed_bgz, $ref_bed_bgz, $retain_no_qry_ref_bound_set, $signal_end3_bed_bgz, $signal_end5_bed_bgz, $tabix_bin, $trnscpt_set_end_priority_hsh_ref, $use_ref_only_end_pos
#	toCall: my ($qry_bed_bgz, $ref_bed_bgz, $chrom_size_path, $chrom_fasta_path, $novel_model_prefix, $conf_end3_merge_flank, $conf_end5_merge_flank, $doubtful_end_merge_dist, $doubtful_end_avoid_summit, $use_ref_only_end_pos, $conf_end5_bed_bgz, $conf_end3_bed_bgz, $conf_junction_bed, $signal_end5_bed_bgz, $signal_end3_bed_bgz, $conf_end3_add_ref, $conf_end5_add_ref, $min_size_split, $min_frac_split, $min_summit_dist_split, $min_qry_score, $retain_no_qry_ref_bound_set, $trnscpt_set_end_priority_hsh_ref, $min_transcript_length, $min_exon_length, $print_trnscrptID, $bedtools_bin, $tabix_bin, $bgzip_bin, $max_thread, $min_output_qry_count, $out_prefix, $out_dir) = &readParameters();
#	calledInLine: 140
#....................................................................................................................................................#
	
	my ($qry_bed_bgz, $ref_bed_bgz, $chrom_size_path, $chrom_fasta_path, $novel_model_prefix, $conf_end3_merge_flank, $conf_end5_merge_flank, $doubtful_end_merge_dist, $doubtful_end_avoid_summit, $use_ref_only_end_pos, $conf_end5_bed_bgz, $conf_end3_bed_bgz, $conf_junction_bed, $signal_end5_bed_bgz, $signal_end3_bed_bgz, $conf_end3_add_ref, $conf_end5_add_ref, $min_size_split, $min_frac_split, $min_summit_dist_split, $min_qry_score, $retain_no_qry_ref_bound_set, $trnscpt_set_end_priority, $min_transcript_length, $min_exon_length, $print_trnscrptID, $bedtools_bin, $tabix_bin, $bgzip_bin, $max_thread, $min_output_qry_count, $out_prefix, $out_dir);
	
	$novel_model_prefix = 'ONTT';
	$conf_end3_merge_flank = 50;
	$conf_end5_merge_flank = 50;
	$max_thread = 5;
	$min_qry_score = 10;
	$retain_no_qry_ref_bound_set = 'yes';
	$doubtful_end_avoid_summit = 'yes';
	$use_ref_only_end_pos = 'yes';
	$print_trnscrptID = 'no';
	$conf_end3_add_ref = 'no';
	$conf_end5_add_ref = 'no';
	$bedtools_bin = 'bedtools';
	$tabix_bin = 'tabix';
	$bgzip_bin = 'bgzip';
	$min_transcript_length = 50;
	$min_exon_length = 1;
	$min_size_split = 100;
	$min_frac_split = 0.20;
	$min_summit_dist_split=50;
	$min_exon_length = 1;
	$trnscpt_set_end_priority = 'summit:commonest:longest';
	$doubtful_end_merge_dist = 100;
	$min_output_qry_count = 1;

	GetOptions 	(
		"qry_bed_bgz=s"						=>	\$qry_bed_bgz,
		"ref_bed_bgz=s"						=>	\$ref_bed_bgz,
		"chrom_size_path=s"					=>	\$chrom_size_path,
		"chrom_fasta_path=s"					=>	\$chrom_fasta_path,
		"conf_end5_bed_bgz=s"				=>	\$conf_end5_bed_bgz,
		"conf_end3_bed_bgz=s"				=>	\$conf_end3_bed_bgz,
		"conf_junction_bed=s"				=>	\$conf_junction_bed,
		"signal_end5_bed_bgz=s"				=>	\$signal_end5_bed_bgz,
		"signal_end3_bed_bgz=s"				=>	\$signal_end3_bed_bgz,
		"novel_model_prefix:s"				=>	\$novel_model_prefix,
		"conf_end3_merge_flank:i"			=>	\$conf_end3_merge_flank,
		"conf_end5_merge_flank:i"			=>	\$conf_end5_merge_flank,
		"min_size_split:i"					=>	\$min_size_split,
		"min_frac_split:f"					=>	\$min_frac_split,
		"min_summit_dist_split:i"			=>	\$min_summit_dist_split,
		"doubtful_end_merge_dist:i"		=>	\$doubtful_end_merge_dist,
		"doubtful_end_avoid_summit:s"		=>	\$doubtful_end_avoid_summit,
		"use_ref_only_end_pos:s"			=>	\$use_ref_only_end_pos,
		"min_transcript_length:i"			=>	\$min_transcript_length,
		"min_exon_length:i"					=>	\$min_exon_length,
		"min_output_qry_count:i"			=>	\$min_output_qry_count,
		"conf_end3_add_ref:s"				=>	\$conf_end3_add_ref,
		"conf_end5_add_ref:s"				=>	\$conf_end5_add_ref,
		"print_trnscrptID:s"					=>	\$print_trnscrptID,
		"retain_no_qry_ref_bound_set:s"	=>	\$retain_no_qry_ref_bound_set,
		"trnscpt_set_end_priority:s"		=>	\$trnscpt_set_end_priority,
		"min_qry_score:i"						=>	\$min_qry_score,
		"bedtools_bin:s"						=>	\$bedtools_bin,
		"tabix_bin:s"							=>	\$tabix_bin,
		"bgzip_bin:s"							=>	\$bgzip_bin,
		"max_thread:i"							=>	\$max_thread,
		"out_prefix:s"							=>	\$out_prefix,
		"out_dir=s"								=>	\$out_dir,
		'help'									=>	sub { HelpMessage(0) },
	) or HelpMessage(1);

	my $opt_check_hsh_ref = {
		'qry_bed_bgz' => $qry_bed_bgz,
		'ref_bed_bgz' => $ref_bed_bgz,
		'chrom_size_path' => $chrom_size_path,
		'chrom_fasta_path' => $chrom_fasta_path,
		'conf_end5_bed_bgz' => $conf_end5_bed_bgz,
		'conf_junction_bed' => $conf_junction_bed,
		'conf_end3_bed_bgz' => $conf_end3_bed_bgz,
		'signal_end5_bed_bgz' => $signal_end5_bed_bgz,
		'signal_end3_bed_bgz' => $signal_end3_bed_bgz,
		'out_dir' => $out_dir,
	};
	
	my $required = 'yes';
	print "\n";
	foreach my $option_name (keys %{$opt_check_hsh_ref}) {
		if (not defined $opt_check_hsh_ref->{$option_name}) {
			print "WARNING: option \"$option_name\" is requied\n";
			$required = 'no';
		}
	}
	if ($required eq 'no') {
		print "WARNING: quitting. Please check this help message for required options\n";
		print "\n";
		HelpMessage(1);
	}

	my $boolean_hsh_ref = {
		'doubtful_end_avoid_summit' => $doubtful_end_avoid_summit,
		'conf_end3_add_ref' => $conf_end3_add_ref,
		'conf_end5_add_ref' => $conf_end5_add_ref,
		'use_ref_only_end_pos' => $use_ref_only_end_pos,
		'retain_no_qry_ref_bound_set' => $retain_no_qry_ref_bound_set,
	};
	
	foreach my $param (sort keys %{$boolean_hsh_ref}) {
		my $value = $boolean_hsh_ref->{$param};
		if ($value ne 'yes' and $value ne 'no') {
			print "WARNING: param must be yes or no. Please check this help message for required options\n";
			print "\n";
			HelpMessage(1);
		}
	}

	#---check file
	my $file_check_hsh_ref = {
		'qry_bed_bgz' => $qry_bed_bgz,
		'ref_bed_bgz' => $ref_bed_bgz,
		'chrom_size_path' => $chrom_size_path,
		'chrom_fasta_path' => $chrom_fasta_path,
		'conf_end5_bed_bgz' => $conf_end5_bed_bgz,
		'conf_end3_bed_bgz' => $conf_end3_bed_bgz,
		'signal_end5_bed_bgz' => $signal_end5_bed_bgz,
		'signal_end3_bed_bgz' => $signal_end3_bed_bgz,
	};
	
	foreach my $option_name (keys %{$file_check_hsh_ref}) {
		my $file_path = $file_check_hsh_ref->{$option_name};
		if (not -s $file_path) {
			die "Quitting: File $option_name does not exists at $file_path";
		}
	}

	foreach my $bgz_path ($ref_bed_bgz, $conf_end5_bed_bgz, $conf_end3_bed_bgz, $signal_end5_bed_bgz, $signal_end3_bed_bgz) {
		if ($bgz_path !~ m/\.bed\.bgz$/) {
			die "Quitting: $bgz_path must be in *.bed.bgz format\n";
		}
		
		if (not -s "$bgz_path.tbi") {
			die "Quitting: $bgz_path must be be indexed as *.bed.bgz.tbi using tabix\n";
		}
	}
	
	if ($qry_bed_bgz !~ m/\.bed\.bgz$/ and $qry_bed_bgz !~ m/\.txt$/) {
			die "Quitting: qry_bed_bgz $qry_bed_bgz must be in *.bed.bgz format or *.txt format\n";
	}
	
	if ($min_output_qry_count < 1) {
		die "Quitting: min_output_qry_count must be at least 1\n";
	}
	
	my $trnscpt_set_end_priority_hsh_ref = {};
	my $priority = 1;
	foreach my $end_type (split /:/, $trnscpt_set_end_priority) {
		if ($end_type eq 'summit' or $end_type eq 'longest' or $end_type eq 'commonest') {
			$trnscpt_set_end_priority_hsh_ref->{$end_type} = $priority;
			$priority++;
		}
	}
	
	if (not exists $trnscpt_set_end_priority_hsh_ref->{'summit'} or not exists $trnscpt_set_end_priority_hsh_ref->{'longest'} or not exists $trnscpt_set_end_priority_hsh_ref->{'commonest'}) {
			die "Quitting: trnscpt_set_end_priority must contain summit, longest and commonest delimited by :\n";
	}
	
	chop $out_dir if ($out_dir =~ m/\/$/); #---remove the last slash
	system "mkdir -p -m 755 $out_dir/";
	
	return($qry_bed_bgz, $ref_bed_bgz, $chrom_size_path, $chrom_fasta_path, $novel_model_prefix, $conf_end3_merge_flank, $conf_end5_merge_flank, $doubtful_end_merge_dist, $doubtful_end_avoid_summit, $use_ref_only_end_pos, $conf_end5_bed_bgz, $conf_end3_bed_bgz, $conf_junction_bed, $signal_end5_bed_bgz, $signal_end3_bed_bgz, $conf_end3_add_ref, $conf_end5_add_ref, $min_size_split, $min_frac_split, $min_summit_dist_split, $min_qry_score, $retain_no_qry_ref_bound_set, $trnscpt_set_end_priority_hsh_ref, $min_transcript_length, $min_exon_length, $print_trnscrptID, $bedtools_bin, $tabix_bin, $bgzip_bin, $max_thread, $min_output_qry_count, $out_prefix, $out_dir);

}
sub reportAndLogStatus {
#....................................................................................................................................................#
#	subroutineCategory: log
#	dependOnSub: currentTime|871
#	appearInSub: assignConfidentEnds|259, assignDoubtfulEnds|319, assignEnd3IncompleteSetToEnd3CompleteModels|492, assignEnd3IncompleteSetToEnd3IncompleteModels|538, assignOrphanSetToExistingModels|668, defineChromInfo|889, defineEndRegion|994, defineQryBedAry|1281, defineTrnscptModels|1382, defineTrnscptSets|1431, extractTranscriptBounds|1571, filterModelAndSetWithQry|1714, getJunctSeq|1817, mergeInConfEndRegion|1932, poolChromBedAndTable|2199, poolChromLog|2286, printInfoBedBothEnd|2343, printTranscriptBed|3149, processPerChromosome|3216, readConfJunction|3399, retainNoQryRefBoundSet|3636
#	primaryAppearInSection: 2_defineout_dirPath|153
#	secondaryAppearInSection: 3_process|197
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportAndLogStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 181, 184, 185, 186, 187, 188, 189, 190, 307, 312, 363, 388, 396, 435, 512, 574, 575, 576, 610, 611, 612, 627, 634, 657, 723, 724, 725, 726, 727, 728, 729, 903, 970, 987, 1049, 1315, 1402, 1426, 1458, 1502, 1505, 1516, 1566, 1598, 1628, 1689, 1690, 1731, 1757, 1834, 1846, 1866, 1869, 1967, 1975, 2236, 2272, 2277, 2299, 2378, 2434, 3173, 3188, 3195, 3199, 3269, 3303, 3311, 3314, 3317, 3320, 3324, 3327, 3330, 3333, 3335, 3337, 3341, 3344, 3347, 3350, 3353, 3356, 3359, 3362, 3365, 3378, 3385, 3387, 3391, 3433, 3652
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->871
	print $tmplog_fh "[".&currentTime()."] ".$message.$lineEnd if $lineEnd ne "\r";#->871
	
	return ();
}
sub retainNoQryRefBoundSet {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|3614
#	appearInSub: defineTrnscptSets|1431
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $chrom_num, $end3_info_hsh_ref, $end3_num, $end5_info_hsh_ref, $end5_num, $report_tag, $set_info_hsh_ref, $set_num, $trnscpt_info_hsh_ref
#	output: 
#	toCall: &retainNoQryRefBoundSet($report_tag, $set_info_hsh_ref, $trnscpt_info_hsh_ref, $end5_info_hsh_ref, $end3_info_hsh_ref, $end3_num, $end5_num, $set_num, $chrom_num);
#	calledInLine: 1508
#....................................................................................................................................................#
	my ($report_tag, $set_info_hsh_ref, $trnscpt_info_hsh_ref, $end5_info_hsh_ref, $end3_info_hsh_ref, $end3_num, $end5_num, $set_num, $chrom_num) = @_;
	
	my $num_proc = 0;
	foreach my $old_set_ID (keys %{$set_info_hsh_ref}) {
		$num_proc++;
		&reportAndLogStatus("$report_tag start retaining ref set bounds without qry support << $num_proc sets processed", 10, "\n") if $num_proc%10000 == 0;#->3614
		my $ref_count = 0;
		my $qry_count = 0;
		if (exists $set_info_hsh_ref->{$old_set_ID}{'T'}{'ref'}) {
			$ref_count = keys %{$set_info_hsh_ref->{$old_set_ID}{'T'}{'ref'}};
		}
		if (exists $set_info_hsh_ref->{$old_set_ID}{'T'}{'qry'}) {
			$qry_count = keys %{$set_info_hsh_ref->{$old_set_ID}{'T'}{'qry'}};
		}
		if ($ref_count >= 1 and $qry_count == 0) {
			my $old_bound_str = $set_info_hsh_ref->{$old_set_ID}{'BS'};
			my @old_bound_str_ary = split /\_/, $old_bound_str;
			
			foreach my $trnscpt_ID (keys %{$set_info_hsh_ref->{$old_set_ID}{'T'}{'ref'}}) {#---[2023/08/21 22:21] add a new set and define new end3 and end5
				my ($strand, $chromStart, $chromEnd) = @{$trnscpt_info_hsh_ref->{'ref'}{$trnscpt_ID}{'I'}};
				$set_num++;
				$end5_num++;
				$end3_num++;
				my $old_end5_ID = $trnscpt_info_hsh_ref->{'ref'}{$trnscpt_ID}{'F'};
				my $old_end3_ID = $trnscpt_info_hsh_ref->{'ref'}{$trnscpt_ID}{'T'};

				my $new_end5_ID = "F".$chrom_num.$end5_num;
				my $new_end3_ID = "T".$chrom_num.$end3_num;
				my $new_set_ID = "S".$chrom_num.$set_num;
				
				#---[2023/08/21 22:46] overwrite trnscpt_ID hsh
				$trnscpt_info_hsh_ref->{'ref'}{$trnscpt_ID}{'F'} = $new_end5_ID;
				$trnscpt_info_hsh_ref->{'ref'}{$trnscpt_ID}{'T'} = $new_end3_ID;
				$trnscpt_info_hsh_ref->{'ref'}{$trnscpt_ID}{'S'} = $new_set_ID;

				$end5_info_hsh_ref->{$new_end5_ID}{'score'} = 1;
				$end3_info_hsh_ref->{$new_end3_ID}{'score'} = 1;
				$end5_info_hsh_ref->{$new_end5_ID}{'ref_only'} = 'Y';
				$end3_info_hsh_ref->{$new_end3_ID}{'ref_only'} = 'Y';
				$end5_info_hsh_ref->{$new_end5_ID}{'ST'} = $strand;
				$end3_info_hsh_ref->{$new_end3_ID}{'ST'} = $strand;
				my $CS_5 = $chromStart;
				my $CS_3 = $chromEnd-1;
				($CS_5, $CS_3) = ($CS_3, $CS_5) if ($strand eq '-');
				
				#---[2023/08/21 22:49] end hsh
				$end5_info_hsh_ref->{$new_end5_ID}{'summit'} = $CS_5+1;
				$end5_info_hsh_ref->{$new_end5_ID}{'CS'} = $CS_5;
				$end5_info_hsh_ref->{$new_end5_ID}{'CE'} = $CS_5+1;
				$end5_info_hsh_ref->{$new_end5_ID}{'OS'} = $CS_5;
				$end5_info_hsh_ref->{$new_end5_ID}{'OE'} = $CS_5+1;

				$end3_info_hsh_ref->{$new_end3_ID}{'summit'} = $CS_3+1;
				$end3_info_hsh_ref->{$new_end3_ID}{'CS'} = $CS_3;
				$end3_info_hsh_ref->{$new_end3_ID}{'CE'} = $CS_3+1;
				$end3_info_hsh_ref->{$new_end3_ID}{'OS'} = $CS_3;
				$end3_info_hsh_ref->{$new_end3_ID}{'OE'} = $CS_3+1;

				delete $end5_info_hsh_ref->{$old_end5_ID}{'set'}{$old_set_ID};
				delete $end3_info_hsh_ref->{$old_end3_ID}{'set'}{$old_set_ID};
				
				$end5_info_hsh_ref->{$old_end5_ID}{'set'}{$new_set_ID}++;
				$end3_info_hsh_ref->{$old_end3_ID}{'set'}{$new_set_ID}++;
				
				my @new_bound_str_ary = @old_bound_str_ary;
				$new_bound_str_ary[0] = $new_end5_ID;
				$new_bound_str_ary[-2] = $new_end3_ID;
				my $new_bound_str = join "_", @new_bound_str_ary;
				$set_info_hsh_ref->{$new_set_ID}{'BS'} = $new_bound_str; #---[2023/07/21 12:34] set string
				$set_info_hsh_ref->{$new_set_ID}{'ST'} = $strand; #---[2023/07/21 12:34] set string
				$set_info_hsh_ref->{$new_set_ID}{'AM'} = undef;
				$set_info_hsh_ref->{$new_set_ID}{'CM'} = 'Y';
				$set_info_hsh_ref->{$new_set_ID}{'CS_hsh'}{$chromStart}++; #---[2023/07/21 12:34] complete
				$set_info_hsh_ref->{$new_set_ID}{'CE_hsh'}{$chromEnd}++; #---[2023/07/21 12:34] complete
				$set_info_hsh_ref->{$new_set_ID}{'T'}{'ref'}{$trnscpt_ID}++; #---[2023/07/21 12:57] transcript 
				$set_info_hsh_ref->{$new_set_ID}{'CM'} = 'Y';
				$set_info_hsh_ref->{$new_set_ID}{'RR'} = 'Y';
				$set_info_hsh_ref->{$new_set_ID}{'OR'} = [$chromStart, $chromEnd]; #---[2023/08/21 23:50] original range
			}
			delete $set_info_hsh_ref->{$old_set_ID};
		}
	}
	return ();
}
sub timeStamp {
#....................................................................................................................................................#
#	subroutineCategory: time, general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: none
#	output: $curntTimeStamp
#	toCall: my ($curntTimeStamp) = &timeStamp();
#	calledInLine: 125
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $curntTimeStamp = sprintf "%04d.%02d.%02d.%02d.%02d.%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec;	

	return ($curntTimeStamp);
}

exit;


















































