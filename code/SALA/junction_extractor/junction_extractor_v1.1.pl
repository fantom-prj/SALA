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

 Usage:
   junction_extractor_v1.1.pl [options] --in_bam --chrom_size_path --chrom_fasta_path --out_prefix --out_dir
   
   --in_bam                    <required> [path]    bam file of the alignemnt  
   --out_dir                   <required> [path]    output directory
   --chrom_size_path           <required> [path]    a txt file contains the chromsome size in format of chrom\tsize
   --chrom_fasta_path          <required> [path]    a fasta file contains the sequence of chromsome
   --out_prefix                (optional) [string]  output files prefix, if not defined, qry_bed_bgz filename will be used
   --min_nt_qual               (optional) [integer] the minimum quality score of a base pair to be considered as high quality.
                                                    A junction on a read is considered to be high quality if all its six flanking
                                                    nucleotide (-3,-2,-1 upstream of donor site and 1,2,3 downstream of acceptor site) 
                                                    has quality greater than or equal to min_nt_qual [default=20]
   --max_thread                (optional) [integer] number of threads to be used [default=5]
   --samtools_bin              (optional) [path]    path of binary of samtools
   --bedtools_bin              (optional) [path]    path of binary of bedtools
   --tabix_bin                 (optional) [path]    path of binary of tabix
   --bgzip_bin                 (optional) [path]    path of binary of bgzip

 Dependencies:
   perl
   samtools
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
my ($curntTimeStamp) = &timeStamp();#->928
my $ARGVStr = join "\n", (&currentTime(), $scriptAbsPath, @ARGV);#->177
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
my ($in_bam, $chrom_size_path, $chrom_fasta_path, $bedtools_bin, $samtools_bin, $tabix_bin, $bgzip_bin, $max_thread, $min_nt_qual, $min_MAPQ, $out_prefix, $out_dir) = &readParameters();#->759
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
&logCalledCMDAndScript($ARGVStr, $result_script_dir, $scriptAbsPath);#->332
&printStartOrFinishMessage("startMessage");#->564
&reportAndLogStatus("Start processing......", 10, "\n");#->847

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_process
#
#<section ID="process" num="3">
my ($chrom_info_hsh_ref) = &defineChromInfo($chrom_size_path, $result_tmp_dir);#->195
&processPerChromosome($in_bam, $chrom_info_hsh_ref, $chrom_fasta_path, $bedtools_bin, $samtools_bin, $max_thread, $min_nt_qual, $min_MAPQ);#->680
&poolChromJunctionInfoTableAndPrintBed($chrom_info_hsh_ref, $tabix_bin, $result_log_dir, $result_bed_dir, $paramTag, $bgzip_bin);#->357
#system "rm -rf $result_tmp_dir";
system "chmod -R 755 $result_dir";

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_finishingTasks
#
#<section ID="finishingTasks" num="4">
&printOutputFileListAndReadme($ARGVStr, $paramTag, $out_dir);#->449
&printStartOrFinishMessage("finishMessage");#->564
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
#	unassigned [n=8]:
#		defineChromInfo, getJunctionFromRead, getSplicingSiteSeq
#		poolChromJunctionInfoTableAndPrintBed, printInfoJunction, processCIGARSeqmentPerRead
#		processPerChromosome, storeJunctionInfoPerRead
#
#====================================================================================================================================================#

sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: printStartOrFinishMessage|564, reportAndLogStatus|847
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 2_defineout_dirPath|108, 4_finishingTasks|141
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 81, 580, 584, 589, 593, 863, 864
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub defineChromInfo {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportAndLogStatus|847
#	appearInSub: >none
#	primaryAppearInSection: 3_process|128
#	secondaryAppearInSection: >none
#	input: $chrom_size_path, $result_tmp_dir
#	output: $chrom_info_hsh_ref
#	toCall: my ($chrom_info_hsh_ref) = &defineChromInfo($chrom_size_path, $result_tmp_dir);
#	calledInLine: 131
#....................................................................................................................................................#
	my ($chrom_size_path, $result_tmp_dir) = @_;
	
	my $chrom_info_hsh_ref = {};
	&reportAndLogStatus("Reading chromosome size", 10, "\n");#->847
	
	open (CHROMSIZE, "<", $chrom_size_path);
	while (<CHROMSIZE>) {
		chomp;
		my ($chrom, $size) = split /\t/;
		#next if $chrom ne 'chr1' and $chrom ne 'chr3' and $chrom ne 'chr5';
		#next if $chrom ne 'chr7';
		my $chrom_num = $chrom;
		$chrom_num =~ s/chr//;
		$chrom_num = 23 if $chrom_num =~ m/^X/;
		$chrom_num = 24 if $chrom_num eq 'Y';
		$chrom_num = 25 if $chrom_num eq 'M';
		$chrom_num = sprintf "%.2d", $chrom_num;
		my $chrom_dir = "$result_tmp_dir/$chrom";
		system ("mkdir -pm 755 $chrom_dir");
		my $tmp_txt = "$chrom_dir/$chrom.tmp.txt";
		my $junct_info_tsv = "$chrom_dir/$chrom.junct.info.tsv.gz";
		my $splicing_site_bed = "$chrom_dir/$chrom.splicing_site.bed";
		my $splicing_site_fasta = "$chrom_dir/$chrom.splicing_site.fasta";
		my $log_txt = "$chrom_dir/$chrom.log.txt";

		$chrom_info_hsh_ref->{$chrom}{'size'} = $size;
		$chrom_info_hsh_ref->{$chrom}{'tmp_txt'} = $tmp_txt;
		$chrom_info_hsh_ref->{$chrom}{'log_txt'} = $log_txt;
		$chrom_info_hsh_ref->{$chrom}{'chrom_num'} = $chrom_num;
		$chrom_info_hsh_ref->{$chrom}{'junct_info_tsv'} = $junct_info_tsv;
		$chrom_info_hsh_ref->{$chrom}{'splicing_site_bed'} = $splicing_site_bed;
		$chrom_info_hsh_ref->{$chrom}{'splicing_site_fasta'} = $splicing_site_fasta;

		&reportAndLogStatus("chrom $chrom [#$chrom_num] = $size nt", 10, "\n");#->847
	}
	close CHROMSIZE;

	return ($chrom_info_hsh_ref);
}
sub getJunctionFromRead {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: processCIGARSeqmentPerRead|598, reportAndLogStatus|847, storeJunctionInfoPerRead|869
#	appearInSub: processPerChromosome|680
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|128
#	input: $chrom, $in_bam, $min_MAPQ, $min_nt_qual, $report_tag, $samtools_bin
#	output: $junct_info_hsh_ref
#	toCall: my ($junct_info_hsh_ref) = &getJunctionFromRead($in_bam, $samtools_bin, $chrom, $report_tag, $min_nt_qual, $min_MAPQ);
#	calledInLine: 719
#....................................................................................................................................................#
	my ($in_bam, $samtools_bin, $chrom, $report_tag, $min_nt_qual, $min_MAPQ) = @_;

	my $junct_info_hsh_ref = {};
	my $cigar_char_hsh_ref = {
		'chrom' => {'M' => '1','D' => '1'},
		'read' => {	'M' => '1','S' => '1','I' => '1'},
	};
	my $num_proc = 0;
	
	open BAM, "$samtools_bin view $in_bam $chrom | cut -f 1,2,3,4,5,6,10,11 |";
	while (<BAM>) {
		chomp;
		my ($read_name, $samflag, $chrom, $readStart, $MAPQ, $CIGAR, $readSeq, $readQual) = split /\t/;
		$num_proc++;
		&reportAndLogStatus("$report_tag Processing reads >>> $num_proc reads processed", 10, "\n") if $num_proc%5000 == 0;#->847
		next if $CIGAR !~ m/N/; #---[2023/08/10 1:02] skip anything not spliced
		my ($junct_ary_ref, $bound_qual_ary_ref) = &processCIGARSeqmentPerRead($CIGAR, $readStart, $readQual, $chrom, $cigar_char_hsh_ref);#->598
		&storeJunctionInfoPerRead($junct_info_hsh_ref, $junct_ary_ref, $bound_qual_ary_ref, $samflag, $MAPQ, $min_nt_qual, $min_MAPQ);#->869
	}
	close BAM;
	&reportAndLogStatus("$report_tag Finished processing reads ### $num_proc reads processed", 10, "\n");#->847

	return ($junct_info_hsh_ref);
}
sub getSplicingSiteSeq {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: processPerChromosome|680
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|128
#	input: $bedtools_bin, $chrom_fasta_path, $junct_info_hsh_ref, $splicing_site_bed, $splicing_site_fasta
#	output: 
#	toCall: &getSplicingSiteSeq($junct_info_hsh_ref, $chrom_fasta_path, $splicing_site_bed, $splicing_site_fasta, $bedtools_bin);
#	calledInLine: 722
#....................................................................................................................................................#
	my ($junct_info_hsh_ref, $chrom_fasta_path, $splicing_site_bed, $splicing_site_fasta, $bedtools_bin) = @_;
	
	open SITEBED, "| sort -k2,2n -k6,6 >$splicing_site_bed";
	foreach my $junct_ID (sort keys %{$junct_info_hsh_ref}) {
		my ($chrom, $junct_start, $junct_end, $strand) = split /\_/, $junct_ID;
		my $left_start = $junct_start;
		my $left_end = $left_start + 2;
		my $right_end = $junct_end;
		my $right_start = $right_end - 2;
		
		my ($donor_start, $donor_end, $acceptor_start, $acceptor_end) = ($left_start, $left_end, $right_start, $right_end);
		($acceptor_start, $acceptor_end, $donor_start, $donor_end) = ($donor_start, $donor_end, $acceptor_start, $acceptor_end) if $strand eq '-';
		print SITEBED join "", (join "\t", ($chrom, $donor_start, $donor_end, $junct_ID.".donor", '1', $strand)), "\n";
		print SITEBED join "", (join "\t", ($chrom, $acceptor_start, $acceptor_end, $junct_ID.".acceptor", '1', $strand)), "\n";
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
			my $junct_ID = $1;
			my $don_acc = $2;
			chomp (my $seq = <SSFASTA>);
			$seq = uc($seq);
			$junct_info_hsh_ref->{$junct_ID}{'ss'}{$don_acc} = $seq;
		}
	}
	close SSFASTA;
	
	return ();
}
sub logCalledCMDAndScript {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 2_defineout_dirPath|108
#	secondaryAppearInSection: >none
#	input: $ARGVStr, $result_script_dir, $scriptAbsPath
#	output: 
#	toCall: &logCalledCMDAndScript($ARGVStr, $result_script_dir, $scriptAbsPath);
#	calledInLine: 120
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
sub poolChromJunctionInfoTableAndPrintBed {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 3_process|128
#	secondaryAppearInSection: >none
#	input: $bgzip_bin, $chrom_info_hsh_ref, $paramTag, $result_bed_dir, $result_log_dir, $tabix_bin
#	output: 
#	toCall: &poolChromJunctionInfoTableAndPrintBed($chrom_info_hsh_ref, $tabix_bin, $result_log_dir, $result_bed_dir, $paramTag, $bgzip_bin);
#	calledInLine: 133
#....................................................................................................................................................#
	my ($chrom_info_hsh_ref, $tabix_bin, $result_log_dir, $result_bed_dir, $paramTag, $bgzip_bin) = @_;
	
	my $pool_junct_bed = "$result_bed_dir/$paramTag.junct.bed.bgz";
	my $pool_junct_info = "$result_log_dir/$paramTag.junct.info.tsv.gz";
	open OUTBED, "| sort --parallel=10 -k1,1 -k2,2n -k6,6 | $bgzip_bin -c >$pool_junct_bed";
	open OUTJUNCT, "| gzip -c >$pool_junct_info";
	print OUTJUNCT join "", (join "\t", ('junct_ID', 'chrom', 'junct_start', 'junct_end', 'strand', 'splicing_site', 'canonical', 'total_count', 'hi_qual_count', 'max_mapq', 'avg_mapq', 'max_score_per_nt', 'avg_score_per_nt', 'max_score.donor.-3', 'max_score.donor.-2', 'max_score.donor.-1', 'max_score.acceptor.1', 'max_score.acceptor.2', 'max_score.acceptor.3', 'avg_score.donor.-3', 'avg_score.donor.-2', 'avg_score.donor.-1', 'avg_score.acceptor.1', 'avg_score.acceptor.2', 'avg_score.acceptor.3')), "\n";
	foreach my $chrom (sort {$chrom_info_hsh_ref->{$a}{'chrom_num'} <=> $chrom_info_hsh_ref->{$b}{'chrom_num'}} keys %{$chrom_info_hsh_ref}) {
		my $junct_info_tsv = $chrom_info_hsh_ref->{$chrom}{'junct_info_tsv'};
		open CHROMJUNCT, "<", $junct_info_tsv;
		while (<CHROMJUNCT>) {
			print OUTJUNCT $_;
			my ($junct_ID, $chrom, $junct_start, $junct_end, $strand, $splicing_site, $canonical, $total_count) = split /\t/;
			
			my $itemRgb = '228,26,28';
			$itemRgb = '55,126,184' if ($strand eq '-');
			my $blockCount = 1;
			my $blockSizes = $junct_end - $junct_start;
			my $blockStarts = 0;
			print OUTBED join "", (join "\t", ($chrom, $junct_start, $junct_end, $junct_ID."|".$splicing_site, $total_count, $strand, $junct_start, $junct_end, $itemRgb, $blockCount, $blockSizes, $blockStarts)), "\n";
		}
		close CHROMJUNCT;
	}
	close OUTJUNCT;
	close OUTBED;
	system "$tabix_bin -p bed $pool_junct_bed";
	
	return ();
}
sub printInfoJunction {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: processPerChromosome|680
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|128
#	input: $junct_info_hsh_ref, $junct_info_tsv
#	output: 
#	toCall: &printInfoJunction($junct_info_hsh_ref, $junct_info_tsv);
#	calledInLine: 725
#....................................................................................................................................................#
	my ($junct_info_hsh_ref, $junct_info_tsv) = @_;
	
	open JUNCTINFO, ">", $junct_info_tsv;
	#print JUNCTINFO join "", (join "\t", ('junct_ID', 'chrom', 'junct_start', 'junct_end', 'strand', 'splicing_site', 'canonical', 'total_count', 'hi_qual_count', 'max_mapq', 'avg_mapq', 'max_score_per_nt', 'avg_score_per_nt', 'max_score.donor.-3', 'max_score.donor.-2', 'max_score.donor.-1', 'max_score.acceptor.1', 'max_score.acceptor.2', 'max_score.acceptor.3', 'avg_score.donor.-3', 'avg_score.donor.-2', 'avg_score.donor.-1', 'avg_score.acceptor.1', 'avg_score.acceptor.2', 'avg_score.acceptor.3')), "\n";
	foreach my $junct_ID (sort keys %{$junct_info_hsh_ref}) {
		my ($chrom, $junct_start, $junct_end, $strand) = split /\_/, $junct_ID;
		my $total_count = $junct_info_hsh_ref->{$junct_ID}{'count'}{'total'};
		my $hi_qual_count = $junct_info_hsh_ref->{$junct_ID}{'count'}{'hi_qual'};
		my $sum_mapq = $junct_info_hsh_ref->{$junct_ID}{'mapq'}{'sum'};
		my $avg_mapq = sprintf "%.2f", $sum_mapq/$total_count;
		my $max_mapq = $junct_info_hsh_ref->{$junct_ID}{'mapq'}{'max'};
		my $acceptor = $junct_info_hsh_ref->{$junct_ID}{'ss'}{'acceptor'};
		my $donor = $junct_info_hsh_ref->{$junct_ID}{'ss'}{'donor'};
		my $splicing_site = $donor."-".$acceptor;
		my $canonical = 'N';
		#GT/AG, GC/AG and AT/AC as canonical
		$canonical = 'Y' if $splicing_site eq 'GT-AG' or $splicing_site eq 'GC-AG' or $splicing_site eq 'AT-AC';
		
		my @max_score_ary = ();
		my @avg_score_ary = ();
		foreach my $i (sort {$a <=> $b} keys %{$junct_info_hsh_ref->{$junct_ID}{'qual'}{'sum'}}) {
			my $max_score = $junct_info_hsh_ref->{$junct_ID}{'qual'}{'max'}{$i};
			push @max_score_ary, $max_score;
		}
		
		foreach my $i (sort {$a <=> $b} keys %{$junct_info_hsh_ref->{$junct_ID}{'qual'}{'sum'}}) {
			my $sum_score = $junct_info_hsh_ref->{$junct_ID}{'qual'}{'sum'}{$i};
			my $avg_score = sprintf "%.2f", $sum_score/$total_count;
			push @avg_score_ary, $avg_score;
		}
		
		my $max_score_per_nt = sprintf "%.2f", sum(@max_score_ary)/@max_score_ary;
		my $avg_score_per_nt = sprintf "%.2f", sum(@avg_score_ary)/@avg_score_ary;

		print JUNCTINFO join "", (join "\t", ($junct_ID, $chrom, $junct_start, $junct_end, $strand, $splicing_site, $canonical, $total_count, $hi_qual_count, $max_mapq, $avg_mapq, $max_score_per_nt, $avg_score_per_nt, @max_score_ary, @avg_score_ary)), "\n";
	}
	close JUNCTINFO;

	return ();
}
sub printOutputFileListAndReadme {
#....................................................................................................................................................#
#	subroutineCategory: output
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 4_finishingTasks|141
#	secondaryAppearInSection: >none
#	input: $ARGVStr, $out_dir, $paramTag
#	output: 
#	toCall: &printOutputFileListAndReadme($ARGVStr, $paramTag, $out_dir);
#	calledInLine: 144
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
#	dependOnSub: currentTime|177
#	appearInSub: >none
#	primaryAppearInSection: 2_defineout_dirPath|108, 4_finishingTasks|141
#	secondaryAppearInSection: >none
#	input: $StartOrFinishMessage
#	output: none
#	toCall: &printStartOrFinishMessage($StartOrFinishMessage);
#	calledInLine: 121, 145
#....................................................................................................................................................#

	my ($StartOrFinishMessage) = @_;
	
	if ($StartOrFinishMessage eq "startMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->177
		print "=========================================================================\n\n";

		print $tmplog_fh "\n=========================================================================\n";
		print $tmplog_fh "[".&currentTime()."] starts running ...... \n";#->177
		print $tmplog_fh "=========================================================================\n\n";

	} elsif ($StartOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->177
		print "=========================================================================\n\n";

		print $tmplog_fh "\n=========================================================================\n";
		print $tmplog_fh "[".&currentTime()."] finished running .......\n";#->177
		print $tmplog_fh "=========================================================================\n\n";
	}
}
sub processCIGARSeqmentPerRead {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: getJunctionFromRead|246
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $CIGAR, $chrom, $cigar_char_hsh_ref, $readQual, $readStart
#	output: $bound_qual_ary_ref, $junct_ary_ref
#	toCall: my ($junct_ary_ref, $bound_qual_ary_ref) = &processCIGARSeqmentPerRead($CIGAR, $readStart, $readQual, $chrom, $cigar_char_hsh_ref);
#	calledInLine: 273
#....................................................................................................................................................#
	my ($CIGAR, $readStart, $readQual, $chrom, $cigar_char_hsh_ref) = @_;
	
	my $read_offset = 0;
	my $chrom_offset = 0;
	my @CIGAR_segment_ary = split /(\d+N)/, $CIGAR;
	my $bound_qual_ary_ref = [];
	my $junct_ary_ref = [];

	foreach my $CIGAR_segment (@CIGAR_segment_ary) {
		if ($CIGAR_segment =~ m/(\d+)N/) {
			my $last_chrom_offset = $chrom_offset;
			my $offset = $1;
			$chrom_offset += $offset;
			my $junct_start = $readStart+$last_chrom_offset - 1;
			my $junct_end = $readStart+$chrom_offset - 1;
			push @{$junct_ary_ref}, [$chrom, $junct_start, $junct_end];
		} else {
			my $end_hsh_ref = {
				'left' => {
					'seq' => [qw/N N N/],
					'qual' => [qw/0 0 0/],
				},
				'right' => {
					'seq' => [qw/N N N/],
					'qual' => [qw/0 0 0/],
				},
			};
			my $last_read_offset = $read_offset;
			my @CIGAR_indiv_ary = grep { $_ ne '' } split /(\d+\D)/, $CIGAR_segment;
			foreach my $i (0..$#CIGAR_indiv_ary) {
				my $CIGAR_indiv = $CIGAR_indiv_ary[$i];
				if ($CIGAR_indiv =~ m/(\d+)(\D)/) {
					my ($offset, $cigar_char) = ($1, $2);
					if ($i == 0) {
						$end_hsh_ref->{'left'}{'cigar_char'} = $cigar_char;
						$end_hsh_ref->{'left'}{'offset'} = $offset;
					} 
				
					if ($i == $#CIGAR_indiv_ary) {
						$end_hsh_ref->{'right'}{'cigar_char'} = $cigar_char;
						$end_hsh_ref->{'right'}{'offset'} = $offset;
					}
					$chrom_offset += $offset if (exists $cigar_char_hsh_ref->{'chrom'}{$cigar_char});
					$read_offset += $offset if (exists $cigar_char_hsh_ref->{'read'}{$cigar_char});
				} else {
					die;
				}
			}
		
			my $substr_length = $read_offset - $last_read_offset;
			my $segmentQual = substr $readQual,$last_read_offset,$substr_length;
			foreach my $left_or_right (qw/left right/) {
				my $cigar_char = $end_hsh_ref->{$left_or_right}{'cigar_char'};
				my $offset = $end_hsh_ref->{$left_or_right}{'offset'};
		
				if ($cigar_char eq 'M' and $offset >= 3) {
					if ($left_or_right eq 'left') {
						$end_hsh_ref->{'left'}{'qual'} = [split "", substr($segmentQual, 0, 3)];
					} elsif ($left_or_right eq 'right') {
						$end_hsh_ref->{'right'}{'qual'} = [split "", substr($segmentQual, -3)];
					}
				}
			
				push @{$bound_qual_ary_ref}, $end_hsh_ref->{$left_or_right}{'qual'};
			}
		}
	}	

	return ($junct_ary_ref, $bound_qual_ary_ref);
}
sub processPerChromosome {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: getJunctionFromRead|246, getSplicingSiteSeq|282, printInfoJunction|398, reportAndLogStatus|847
#	appearInSub: >none
#	primaryAppearInSection: 3_process|128
#	secondaryAppearInSection: >none
#	input: $bedtools_bin, $chrom_fasta_path, $chrom_info_hsh_ref, $in_bam, $max_thread, $min_MAPQ, $min_nt_qual, $samtools_bin
#	output: 
#	toCall: &processPerChromosome($in_bam, $chrom_info_hsh_ref, $chrom_fasta_path, $bedtools_bin, $samtools_bin, $max_thread, $min_nt_qual, $min_MAPQ);
#	calledInLine: 132
#....................................................................................................................................................#
	my ($in_bam, $chrom_info_hsh_ref, $chrom_fasta_path, $bedtools_bin, $samtools_bin, $max_thread, $min_nt_qual, $min_MAPQ) = @_;
	
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
					my $junct_info_tsv = $chrom_info_hsh_ref->{$chrom}{'junct_info_tsv'};
					my $splicing_site_bed = $chrom_info_hsh_ref->{$chrom}{'splicing_site_bed'};
					my $splicing_site_fasta = $chrom_info_hsh_ref->{$chrom}{'splicing_site_fasta'};
					my $report_tag = join " ", (sprintf("%-15s","Thread-$threadID $chrom"), ":");

					&reportAndLogStatus("$report_tag >>>>>> Starting $chrom analysis <<<<<<", 10, "\n");#->847
					my ($junct_info_hsh_ref) = &getJunctionFromRead($in_bam, $samtools_bin, $chrom, $report_tag, $min_nt_qual, $min_MAPQ);#->246

					&reportAndLogStatus("$report_tag Getting splicing site sequence", 10, "\n");#->847
					&getSplicingSiteSeq($junct_info_hsh_ref, $chrom_fasta_path, $splicing_site_bed, $splicing_site_fasta, $bedtools_bin);#->282

					&reportAndLogStatus("$report_tag Printing junction info", 10, "\n");#->847
					&printInfoJunction($junct_info_hsh_ref, $junct_info_tsv);#->398
					
				}
				
				return ();
			}
			,($chromForThrAry_ref)
		);
	}
	
	my $data_hsh_ref = {};
	while (keys %threadHsh) {
		my $num_thread = keys %threadHsh;
		&reportAndLogStatus("--------------------------- $num_thread thread running ---------------------------", 10, "\n");#->847
		foreach my $threadID (keys %threadHsh) {
			if (not $threadHsh{$threadID}->is_running()) {
				$threadHsh{$threadID}->join();
				foreach my $chrom (@{$chromForThrHsh{$threadID}}) {
					my $junct_info_tsv = $chrom_info_hsh_ref->{$chrom}{'junct_info_tsv'};
					if (-s $junct_info_tsv) {
						&reportAndLogStatus("Thread-$threadID : junct_info_tsv of $chrom found and finished successfully", 10, "\n");#->847
					} else {
						&reportAndLogStatus("!!!!!!!!WARNING!!!!!!!! junct_info_tsv of $chrom is not found", 10, "\n");#->847
					}
				}
				delete $threadHsh{$threadID};
				&reportAndLogStatus(">>>>>>>>>>>>>>>>>>>>>>>> thread $threadID finished <<<<<<<<<<<<<<<<<<<<<<<<<<<", 10, "\n");#->847
			}
		}
		sleep 5;
	}
	return ();
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|92
#	secondaryAppearInSection: >none
#	input: none
#	output: $bedtools_bin, $bgzip_bin, $chrom_fasta_path, $chrom_size_path, $in_bam, $max_thread, $min_MAPQ, $min_nt_qual, $out_dir, $out_prefix, $samtools_bin, $tabix_bin
#	toCall: my ($in_bam, $chrom_size_path, $chrom_fasta_path, $bedtools_bin, $samtools_bin, $tabix_bin, $bgzip_bin, $max_thread, $min_nt_qual, $min_MAPQ, $out_prefix, $out_dir) = &readParameters();
#	calledInLine: 95
#....................................................................................................................................................#
	
	my ($in_bam, $chrom_size_path, $chrom_fasta_path, $bedtools_bin, $samtools_bin, $tabix_bin, $bgzip_bin, $max_thread, $min_nt_qual, $min_MAPQ, $out_prefix, $out_dir);
	
	$bedtools_bin = 'bedtools';
	$samtools_bin = 'samtools';
	$tabix_bin = 'tabix';
	$bgzip_bin = 'bgzip';
	$min_nt_qual = 20;
	$min_MAPQ = 20;

	GetOptions 	(
		"in_bam=s"				=>	\$in_bam,
		"chrom_size_path=s"	=>	\$chrom_size_path,
		"chrom_fasta_path=s"	=>	\$chrom_fasta_path,
		"min_nt_qual:i"		=>	\$min_nt_qual,
		"min_MAPQ:i"			=>	\$min_MAPQ,
		"samtools_bin:s"		=>	\$samtools_bin,
		"bedtools_bin:s"		=>	\$bedtools_bin,
		"tabix_bin:s"			=>	\$tabix_bin,
		"bgzip_bin:s"			=>	\$bgzip_bin,
		"max_thread:i"			=>	\$max_thread,
		"out_prefix=s"			=>	\$out_prefix,
		"out_dir=s"				=>	\$out_dir,
		'help'					=>	sub { HelpMessage(0) },
	) or HelpMessage(1);

	my $opt_check_hsh_ref = {
		'in_bam' => $in_bam,
		'chrom_size_path' => $chrom_size_path,
		'chrom_fasta_path' => $chrom_fasta_path,
		'out_prefix' => $out_prefix,
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

	#---check file
	my $file_check_hsh_ref = {
		'in_bam' => $in_bam,
		'chrom_fasta_path' => $chrom_fasta_path,
		'chrom_size_path' => $chrom_size_path,
	};
	
	foreach my $option_name (keys %{$file_check_hsh_ref}) {
		my $file_path = $file_check_hsh_ref->{$option_name};
		if (not -s $file_path) {
			die "Quitting: File $option_name does not exists at $file_path";
		}
	}

	if ($in_bam !~ m/\.bam$/) {
			die "Quitting: in_bam $in_bam must be in *.bam format\n";
	}
	
	if (not -s $in_bam.".bai") {
			die "Quitting: index of $in_bam is not found\n";
	}
	
	chop $out_dir if ($out_dir =~ m/\/$/); #---remove the last slash
	system "mkdir -p -m 755 $out_dir/";
	
	return($in_bam, $chrom_size_path, $chrom_fasta_path, $bedtools_bin, $samtools_bin, $tabix_bin, $bgzip_bin, $max_thread, $min_nt_qual, $min_MAPQ, $out_prefix, $out_dir);

}
sub reportAndLogStatus {
#....................................................................................................................................................#
#	subroutineCategory: log
#	dependOnSub: currentTime|177
#	appearInSub: defineChromInfo|195, getJunctionFromRead|246, processPerChromosome|680
#	primaryAppearInSection: 2_defineout_dirPath|108
#	secondaryAppearInSection: 3_process|128
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportAndLogStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 122, 209, 239, 271, 277, 718, 721, 724, 738, 745, 747, 751
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->177
	print $tmplog_fh "[".&currentTime()."] ".$message.$lineEnd if $lineEnd ne "\r";#->177
	
	return ();
}
sub storeJunctionInfoPerRead {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: getJunctionFromRead|246
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $MAPQ, $bound_qual_ary_ref, $junct_ary_ref, $junct_info_hsh_ref, $min_MAPQ, $min_nt_qual, $samflag
#	output: 
#	toCall: &storeJunctionInfoPerRead($junct_info_hsh_ref, $junct_ary_ref, $bound_qual_ary_ref, $samflag, $MAPQ, $min_nt_qual, $min_MAPQ);
#	calledInLine: 274
#....................................................................................................................................................#
	my ($junct_info_hsh_ref, $junct_ary_ref, $bound_qual_ary_ref, $samflag, $MAPQ, $min_nt_qual, $min_MAPQ) = @_;

	#J[0,1,2,3]
	#B[0,1,2,3,4,5,6,7,8,9]
	my $strand = '+';
	$strand = '-' if ($samflag & 16);
	
	foreach my $junct_i (0..$#{$junct_ary_ref}) {
		my $bound_left_i = ($junct_i*2)+1;
		my $bound_right_i = ($junct_i*2)+2;
		my ($chrom, $junct_start, $junct_end) = @{$junct_ary_ref->[$junct_i]};
		my $junct_ID = join "_", ($chrom, $junct_start, $junct_end, $strand);

		my @qual_ary = (@{$bound_qual_ary_ref->[$bound_left_i]}, @{$bound_qual_ary_ref->[$bound_right_i]});
		if ($strand eq '-') {
			@qual_ary = reverse @qual_ary;
		}
	
		if (not exists $junct_info_hsh_ref->{$junct_ID}) {
			foreach my $i (0..$#qual_ary) {
				$junct_info_hsh_ref->{$junct_ID}{'qual'}{'sum'}{$i} = 0;
				$junct_info_hsh_ref->{$junct_ID}{'qual'}{'max'}{$i} = 0;
			}
			$junct_info_hsh_ref->{$junct_ID}{'count'}{'total'} = 0;
			$junct_info_hsh_ref->{$junct_ID}{'count'}{'hi_qual'} = 0;
			$junct_info_hsh_ref->{$junct_ID}{'mapq'}{'max'} = 0;
			$junct_info_hsh_ref->{$junct_ID}{'mapq'}{'sum'} = 0;
		}
	
		$junct_info_hsh_ref->{$junct_ID}{'mapq'}{'sum'} += $MAPQ;
		$junct_info_hsh_ref->{$junct_ID}{'mapq'}{'max'} = $MAPQ if $MAPQ > $junct_info_hsh_ref->{$junct_ID}{'mapq'}{'max'};
	
		my $hi_qual = 'yes';
		foreach my $i (0..$#qual_ary) {
			my $qual = $qual_ary[$i];
			my $score = ord($qual) - 33;
			$hi_qual= 'no' if $score < $min_nt_qual or $MAPQ < $min_MAPQ;
			$junct_info_hsh_ref->{$junct_ID}{'qual'}{'sum'}{$i} += $score;
			$junct_info_hsh_ref->{$junct_ID}{'qual'}{'max'}{$i} = $score if $score > $junct_info_hsh_ref->{$junct_ID}{'qual'}{'max'}{$i};
		}

		$junct_info_hsh_ref->{$junct_ID}{'count'}{'total'}++;
		$junct_info_hsh_ref->{$junct_ID}{'count'}{'hi_qual'}++ if $hi_qual eq 'yes';
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
#	calledInLine: 80
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $curntTimeStamp = sprintf "%04d.%02d.%02d.%02d.%02d.%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec;	

	return ($curntTimeStamp);
}

exit;


















































