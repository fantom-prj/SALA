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
   This tool trims the tail (3'end') of a fastq file based on polyA and GTn (CAn) stretches and print to standard output

 Usage:
   tail_trimmer_v1.4.pl [options] --fastq_path --out_dir
   
   --fastq_path              <required> [string]  fastq file path
   --win_size                (optional) [integer] the size of window (nt) used to scan for polyA and GTn (CAn) stretches (default=12, must be >=8 and <=20)
   --max_mismatch            (optional) [integer] the number of mismatch within window allowed (default=2). Maximum 3 allowed.
   --match_5end_size         (optional) [integer] the size of edge (on each side) within window (nt) without mismtach (default=4, must be <= 1/3 of win_size)
   --min_stretch_frac        (optional) [0 to 1]  the minimum fraction of sequence is polyA and GTn (CAn) stretches within the region (size defined by 
                                                  frac_check_size_max) downstream of the trim position (default=0.5)
   --subsample_frac          (optional) [0 to 1]  subsample only a fraction of the reads (default=1)
   --max_proc_read           (optional) [integer] maximum number of read to bed processed, process all if undefined (default=undefined)
   --frac_check_size_max     (optional) [integer] the maximum size of region to check for within window (nt) without mismtach (default=100)
   --rescue_tail_nt          (optional) [integer] resucue certain number of nucleotide in the trimmed tail and save as a seperate fastq as *.cleaned.rescue_tail.fq.gz 
   --revise_read_prefix      (optional) [string]  if defined, a prefix used to rename the read (followed by number), e.g. iPSC, will give iPSC_\d+ as read name.
                                                  A log on original to revised read ID will be printed. If undefined, original rename will be kept (default=undefined)
   --out_prefix              (optional) [string]  output files prefix, if not defined, fastq filename will be used
   --out_dir                 <required> [string]  output directory

 Dependencies:
   perl

 For demo, cd to the dir and run,
   perl ./tail_trimmer_v1.4.pl \
   --fastq_path=./demo/input/demo.fastq \
   --out_dir=./demo/output/

=head1 VERSION

1.1   -debut
1.2   -fixed zero head length caused by read with only A stretch
      -added subsample_frac option
      -added max_proc_read option
1.3   -fixed detection of strand information
1.4   -added rescue_tail_nt option and revise_read_prefix
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
my ($curntTimeStamp) = &timeStamp();#->949
my $ARGVStr = join "\n", (&currentTime(), $scriptAbsPath, @ARGV);#->186
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
my ($fastq_path, $win_size, $max_mismatch, $match_5end_size, $min_stretch_frac, $frac_check_size_max, $subsample_frac, $max_proc_read, $rescue_tail_nt, $revise_read_prefix, $out_prefix, $out_dir) = &readParameters();#->803
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
my $result_fastq_dir = "$result_dir/fastq/"; push @mkDirAry, $result_fastq_dir;
my $result_log_dir = "$result_dir/log/"; push @mkDirAry, $result_log_dir;
my $result_script_dir = "$result_dir/script/"; push @mkDirAry, $result_script_dir;
foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}

open $tmplog_fh, ">", "$result_dir/00_screen_log.$curntTimeStamp.log.txt";
&logCalledCMDAndScript($ARGVStr, $result_script_dir, $scriptAbsPath);#->349
&printStartOrFinishMessage("startMessage");#->523
&reportAndLogStatus("Parameters:", 10, "\n");#->927
&reportAndLogStatus("fastq_path = $fastq_path", 10, "\n");#->927
&reportAndLogStatus("win_size = $win_size:", 10, "\n");#->927
&reportAndLogStatus("max_mismatch = $max_mismatch", 10, "\n");#->927
&reportAndLogStatus("match_5end_size = $match_5end_size", 10, "\n");#->927
&reportAndLogStatus("min_stretch_frac = $min_stretch_frac", 10, "\n");#->927
&reportAndLogStatus("frac_check_size_max = $frac_check_size_max", 10, "\n");#->927
&reportAndLogStatus("subsample_frac = $subsample_frac", 10, "\n");#->927
&reportAndLogStatus("rescue_tail_nt = $rescue_tail_nt", 10, "\n");#->927
&reportAndLogStatus("revise_read_prefix = $revise_read_prefix", 10, "\n");#->927
if (not defined $max_proc_read) {
	&reportAndLogStatus("max_proc_read = undefined", 10, "\n");#->927
} else {
	&reportAndLogStatus("max_proc_read = $max_proc_read", 10, "\n");#->927
}
&reportAndLogStatus("Start processing......", 10, "\n");#->927
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_process
#
#<section ID="process" num="3">
my ($strand_regex_hsh_ref) = &getPolyStretchRegex($win_size, $max_mismatch, $match_5end_size);#->276
my ($count_data_hsh_ref) = &readFastq($fastq_path, $strand_regex_hsh_ref, $win_size, $min_stretch_frac, $frac_check_size_max, $subsample_frac, $max_proc_read, $result_fastq_dir, $result_log_dir, $rescue_tail_nt, $revise_read_prefix, $out_prefix);#->557
&printCountData($count_data_hsh_ref, $result_log_dir, $out_prefix);#->374
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_finishingTasks
#
#<section ID="finishingTasks" num="4">
&printOutputFileListAndReadme($ARGVStr, $paramTag, $out_dir);#->408
&printStartOrFinishMessage("finishMessage");#->523
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
#	unassigned [n=6]:
#		generateMisMatchRegex, getPolyStretchRegex, getStretchFraction
#		printCountData, readFastq, trimTailPolyA
#
#====================================================================================================================================================#

sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: printStartOrFinishMessage|523, reportAndLogStatus|927
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 2_defineout_dirPath|108, 4_finishingTasks|151
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 81, 539, 543, 548, 552, 943, 944
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub generateMisMatchRegex {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: getPolyStretchRegex|276
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|141
#	input: $match_5end_size, $match_string, $max_mismatch
#	output: $regex_string
#	toCall: my ($regex_string) = &generateMisMatchRegex($match_string, $max_mismatch, $match_5end_size);
#	calledInLine: 297
#....................................................................................................................................................#
	my ($match_string, $max_mismatch, $match_5end_size) = @_;

	my @mismatch_3_ary;
	my @mismatch_2_ary;
	my @mismatch_1_ary;
	
	my $mutate_start = $match_5end_size;
	my $mutate_end = length($match_string)-1;
	
	for my $i ($mutate_start..$mutate_end) {
		my $mismatch_1_subpattern = join('',
			substr($match_string, 0, $i),
			'\\w',  # or '\\w'
			substr($match_string, $i+1),
		);
		push @mismatch_1_ary, $mismatch_1_subpattern;

		for my $j ($i+1..$mutate_end) {
			my $mismatch_2_subpattern = join('',
				substr($match_string, 0, $i),
				'\\w',  # or '\\w'
				substr($match_string, $i+1, $j-$i-1),
				'\\w',  # or '\\w'
				substr($match_string, $j+1),
			);
			push @mismatch_2_ary, $mismatch_2_subpattern;

			for my $k ($j+1..$mutate_end) {
				my $mismatch_3_subpattern = join('',
					substr($match_string, 0, $i),
					'\\w',  # or '\\w'
					substr($match_string, $i+1, $j-$i-1),
					'\\w',  # or '\\w'
					substr($match_string, $j+1, $k-$j-1),
					'\\w',  # or '\\w'
					substr($match_string, $k+1),
				);
				push @mismatch_3_ary, $mismatch_3_subpattern;
			}
		}
	}
	
	my $regex_string;
	if ($max_mismatch eq 0) {
		$regex_string = $match_string;
	} elsif ($max_mismatch eq 1) {
		$regex_string = join('|', @mismatch_1_ary);

	} elsif ($max_mismatch eq 2) {
		$regex_string = join('|', @mismatch_2_ary);

	} elsif ($max_mismatch eq 3) {
		$regex_string = join('|', @mismatch_3_ary);

	} else {
		die;
	}
	
	return ($regex_string);
}
sub getPolyStretchRegex {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: generateMisMatchRegex|204, reportAndLogStatus|927
#	appearInSub: >none
#	primaryAppearInSection: 3_process|141
#	secondaryAppearInSection: >none
#	input: $match_5end_size, $max_mismatch, $win_size
#	output: $strand_regex_hsh_ref
#	toCall: my ($strand_regex_hsh_ref) = &getPolyStretchRegex($win_size, $max_mismatch, $match_5end_size);
#	calledInLine: 144
#....................................................................................................................................................#
	my ($win_size, $max_mismatch, $match_5end_size) = @_;
	
	my $tmp_regex_hsh_ref = {};
	foreach my $stretch_type (qw/A GT AC/) {
		my $string_length = 0;
		my $match_string = $stretch_type;
		while ($string_length < $win_size) {
			$match_string = $match_string.$stretch_type;
			$string_length = length($match_string);
		}
		my ($regex_string) = &generateMisMatchRegex($match_string, $max_mismatch, $match_5end_size);#->204
		$tmp_regex_hsh_ref->{$stretch_type} = $regex_string;
		#&reportAndLogStatus("$stretch_type", 10, "\n");#->927
		#&reportAndLogStatus("$regex_string", 10, "\n");#->927
	}
	
	my $strand_regex_hsh_ref = {
		'+' => $tmp_regex_hsh_ref->{'A'}."|".$tmp_regex_hsh_ref->{'GT'},
		'-' => $tmp_regex_hsh_ref->{'A'}."|".$tmp_regex_hsh_ref->{'AC'},
	};

	#&reportAndLogStatus("$strand_regex_hsh_ref->{'+'}", 10, "\n");#->927
	#&reportAndLogStatus("$strand_regex_hsh_ref->{'-'}", 10, "\n");#->927
	
	return ($strand_regex_hsh_ref);
}
sub getStretchFraction {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: readFastq|557
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|141
#	input: $region_seq, $strand
#	output: $max_polyA_length, $max_polyGT_length, $region_length, $stretch_frac, $total_polyA_length, $total_polyGT_length
#	toCall: my ($max_polyGT_length, $max_polyA_length, $region_length, $total_polyA_length, $total_polyGT_length, $stretch_frac) = &getStretchFraction($region_seq, $strand);
#	calledInLine: 645, 690, 691, 697, 698
#....................................................................................................................................................#
	my ($region_seq, $strand) = @_;
	
	my $region_length = length($region_seq);
	my $total_polyA_length = 0;
	my $total_polyGT_length = 0;
	my $max_polyA_length = 0;
	my $max_polyGT_length = 0;

	while ($region_seq =~ m/((AA)+)/g) {
		$total_polyA_length += length($1);
		$max_polyA_length = length($1) if length($1) > $max_polyA_length;
	}
	my $dinucleotide = 'GT';
	$dinucleotide = 'AC' if $strand eq '-';
	while ($region_seq =~ m/(($dinucleotide)+)/g) {
		$total_polyGT_length += length($1);
		$max_polyGT_length = length($1) if length($1) > $max_polyGT_length;
	}
	my $stretch_frac = sprintf "%.2f", ($total_polyA_length+$total_polyGT_length)/$region_length;
	$stretch_frac = '1.00' if $stretch_frac > 1;
	return ($max_polyGT_length, $max_polyA_length, $region_length, $total_polyA_length, $total_polyGT_length, $stretch_frac);

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
#	calledInLine: 119
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
sub printCountData {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 3_process|141
#	secondaryAppearInSection: >none
#	input: $count_data_hsh_ref, $out_prefix, $result_log_dir
#	output: 
#	toCall: &printCountData($count_data_hsh_ref, $result_log_dir, $out_prefix);
#	calledInLine: 146
#....................................................................................................................................................#
	my ($count_data_hsh_ref, $result_log_dir, $out_prefix) = @_;

	foreach my $data_type (keys %{$count_data_hsh_ref->{'region'}}) {
		my @region_ary = qw/head tail upstream dnstream/;
		my $log_path = "$result_log_dir/$out_prefix.$data_type.count.txt";
		open LOG, ">", $log_path;
		my @header_ary = ($data_type, @region_ary);
		print LOG join "", (join "\t", (@header_ary)), "\n";
		foreach my $value (sort {$a <=> $b} keys %{$count_data_hsh_ref->{'region'}{$data_type}}) {
			my @output_ary = ($value);
			foreach my $region (@region_ary) {
				my $count = 0;
				$count = $count_data_hsh_ref->{'region'}{$data_type}{$value}{$region} if exists $count_data_hsh_ref->{'region'}{$data_type}{$value}{$region};
				push @output_ary, $count;
			}
			print LOG join "", (join "\t", (@output_ary)), "\n";
		}
		close LOG;
	}

	return ();
}
sub printOutputFileListAndReadme {
#....................................................................................................................................................#
#	subroutineCategory: output
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 4_finishingTasks|151
#	secondaryAppearInSection: >none
#	input: $ARGVStr, $out_dir, $paramTag
#	output: 
#	toCall: &printOutputFileListAndReadme($ARGVStr, $paramTag, $out_dir);
#	calledInLine: 154
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
#	dependOnSub: currentTime|186
#	appearInSub: >none
#	primaryAppearInSection: 2_defineout_dirPath|108, 4_finishingTasks|151
#	secondaryAppearInSection: >none
#	input: $StartOrFinishMessage
#	output: none
#	toCall: &printStartOrFinishMessage($StartOrFinishMessage);
#	calledInLine: 120, 155
#....................................................................................................................................................#

	my ($StartOrFinishMessage) = @_;
	
	if ($StartOrFinishMessage eq "startMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->186
		print "=========================================================================\n\n";

		print $tmplog_fh "\n=========================================================================\n";
		print $tmplog_fh "[".&currentTime()."] starts running ...... \n";#->186
		print $tmplog_fh "=========================================================================\n\n";

	} elsif ($StartOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->186
		print "=========================================================================\n\n";

		print $tmplog_fh "\n=========================================================================\n";
		print $tmplog_fh "[".&currentTime()."] finished running .......\n";#->186
		print $tmplog_fh "=========================================================================\n\n";
	}
}
sub readFastq {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: getStretchFraction|314, reportAndLogStatus|927, trimTailPolyA|967
#	appearInSub: >none
#	primaryAppearInSection: 3_process|141
#	secondaryAppearInSection: >none
#	input: $fastq_path, $frac_check_size_max, $max_proc_read, $min_stretch_frac, $out_prefix, $rescue_tail_nt, $result_fastq_dir, $revise_read_prefix, $strand_regex_hsh_ref, $subsample_frac, $win_size
#	output: $count_data_hsh_ref
#	toCall: my ($count_data_hsh_ref) = &readFastq($fastq_path, $strand_regex_hsh_ref, $win_size, $min_stretch_frac, $frac_check_size_max, $subsample_frac, $max_proc_read, $result_fastq_dir, $rescue_tail_nt, $revise_read_prefix, $out_prefix);
#	calledInLine: 145
#....................................................................................................................................................#
	my ($fastq_path, $strand_regex_hsh_ref, $win_size, $min_stretch_frac, $frac_check_size_max, $subsample_frac, $max_proc_read, $result_fastq_dir, $result_log_dir, $rescue_tail_nt, $revise_read_prefix, $out_prefix) = @_;

	open TRIMMEDFQ, "| gzip -c >$result_fastq_dir/$out_prefix.cleaned.fq.gz";
	open FULLTAIL, "| gzip -c >$result_fastq_dir/$out_prefix.trimmed_full_tail.fq.gz";
	open UNTRIMMED, "| gzip -c >$result_fastq_dir/$out_prefix.untrimmed_end.fq.gz";
	open JUNCTION, "| gzip -c >$result_fastq_dir/$out_prefix.junction.fq.gz";
	open SHORT, "| gzip -c >$result_fastq_dir/$out_prefix.too_short.fq.gz";

	if ($rescue_tail_nt > 0) {
		open RETAINTAILFQ, "| gzip -c >$result_fastq_dir/$out_prefix.cleaned.rescue_tail.fq.gz";
	}

	if (defined $revise_read_prefix) {
		open READIDLOG, "| gzip -c >$result_log_dir/$out_prefix.original_to_revised_readID.log.txt.gz";
	}
	
	if ($fastq_path =~ m/\.gz$/) {
		open (FILEIN, " gzip -dc $fastq_path|");
	} else {
		open (FILEIN, "<", $fastq_path);
	}
	my $count_data_hsh_ref = {};
	my $read_proc = 0;
	my $read_trimmed = 0;
	my $read_untrimmed = 0;
	my $read_short = 0;
	my $read_skip = 0;
	while (<FILEIN>) {
		chomp (my $rd_name = $_);
		chomp (my $rd_seq = <FILEIN>);
		chomp (my $qual_header = <FILEIN>);
		chomp (my $qual_value = <FILEIN>);
	
		if ($rd_name =~ /\@(.+) strand=(\-|\+)/) {
			
			if (defined $max_proc_read) {
				last if $read_proc >= $max_proc_read;
			}
		
			if ($subsample_frac < 1) {
				if (rand(1) > $subsample_frac) {
					$read_skip++;
					next;
				}
			}
		
			$read_proc++;

			my $original_rd_ID = $1;
			my $strand = $2;
			my $revised_rd_ID = $original_rd_ID;
			if (defined $revise_read_prefix) {
				$revised_rd_ID = $revise_read_prefix."_".$read_proc;
				print READIDLOG join "", (join "\t", ($original_rd_ID, $revised_rd_ID)), "\n";
			}
			my $revised_rd_name = "@".$revised_rd_ID." strand=$strand";
			
			if ($read_proc%2000==0) {
				my $read_trimmed_pct = sprintf "%.2f", 100*$read_trimmed/$read_proc;
				my $read_untrimmed_pct = sprintf "%.2f", 100*$read_untrimmed/$read_proc;
				my $read_short_pct = sprintf "%.2f", 100*$read_short/$read_proc;
				&reportAndLogStatus("$read_proc reads processed (skip=$read_skip). trimmed=$read_trimmed\[$read_trimmed_pct%\] untrimmed=$read_untrimmed\[$read_untrimmed_pct%\] too_short(<20nt)=$read_short\[$read_short_pct%\]", 10, "\n")#->927
			}

			my $stranded_regex = $strand_regex_hsh_ref->{$strand};
			my $cut_pos = length($rd_seq);
			my $head_seq = $rd_seq;
			my $tail_seq = undef;
			my $round = 0;
			my $stretch_frac = 0;
			my ($max_polyGT_length, $max_polyA_length, $region_length, $total_polyA_length, $total_polyGT_length);
			
			while ($rd_seq =~ m/($stranded_regex)/g and $stretch_frac < $min_stretch_frac) {
				$round++;
				my $matched_stretch = $1;
				my $stretch_end = $+[0];
				$cut_pos = $stretch_end-$win_size;
				$head_seq = substr ($rd_seq, 0, $cut_pos);
				$tail_seq = substr ($rd_seq, $cut_pos);
				
				#my $upstream_seq = substr ($head_seq, length($head_seq)-$frac_check_size_max);
				my $dnstream_seq = substr ($tail_seq, 0, $frac_check_size_max);
				($max_polyGT_length, $max_polyA_length, $region_length, $total_polyA_length, $total_polyGT_length, $stretch_frac) = &getStretchFraction($dnstream_seq, $strand);#->314
			}
			
			if ($stretch_frac < $min_stretch_frac) {
				$cut_pos = length($rd_seq);
				$tail_seq = undef;
			}
			
			my $tail_trim_length = 0;
			my $tail_trim_seq = undef;

			if (defined $tail_seq) {
				($tail_trim_length, $tail_trim_seq) = &trimTailPolyA($head_seq);#->967
			} else {
				($tail_trim_length, $tail_trim_seq) = &trimTailPolyA($rd_seq);#->967
			}
			
			if ($tail_trim_length > 0) {
				$cut_pos = $cut_pos - $tail_trim_length;
				$head_seq = substr ($rd_seq, 0, $cut_pos);
				$tail_seq = substr ($rd_seq, $cut_pos);
			}
			
			my $out_rd_name;
			my $out_rd_seq;
			my $out_qual_header;
			my $out_qual_value;
			my $rescue_tail_seq = '';
			my $rescue_tail_qual = '';
			
			my $status;
			
			if (defined $tail_seq) {
				if (length ($head_seq) < 20) {
					$status='head_too_short';
					$read_short++;
					print SHORT "$revised_rd_name\n";
					print SHORT "$rd_seq\n";
					print SHORT "$qual_header\n";
					print SHORT "$qual_value\n";
					next;
				}
				
				my ($tail_max_polyGT_length, $tail_max_polyA_length, $tail_region_length, $tail_total_polyA_length, $tail_total_polyGT_length, $tail_stretch_frac) = &getStretchFraction($tail_seq, $strand);#->314
				my ($head_max_polyGT_length, $head_max_polyA_length, $head_region_length, $head_total_polyA_length, $head_total_polyGT_length, $head_stretch_frac) = &getStretchFraction($head_seq, $strand);#->314
				die "WARNING: head_seq and tail_seq does not match rd_seq. Fatal error. Quitting.\n" if $rd_seq ne $head_seq.$tail_seq;

				my $upstream_seq = substr ($head_seq, length($head_seq)-20);
				my $dnstream_seq = substr ($tail_seq, 0, 20);
				my $junction_seq = join "", ($upstream_seq,$dnstream_seq);
				my ($dnstream_max_polyGT_length, $dnstream_max_polyA_length, $dnstream_region_length, $dnstream_total_polyA_length, $dnstream_total_polyGT_length, $dnstream_stretch_frac) = &getStretchFraction($dnstream_seq, $strand);#->314
				my ($upstream_max_polyGT_length, $upstream_max_polyA_length, $upstream_region_length, $upstream_total_polyA_length, $upstream_total_polyGT_length, $upstream_stretch_frac) = &getStretchFraction($upstream_seq, $strand);#->314

				$count_data_hsh_ref->{'region'}{'max_polyGT_length'}{$dnstream_max_polyGT_length}{'dnstream'}++;
				$count_data_hsh_ref->{'region'}{'max_polyA_length'}{$dnstream_max_polyA_length}{'dnstream'}++;
				$count_data_hsh_ref->{'region'}{'total_polyA_length'}{$dnstream_total_polyA_length}{'dnstream'}++;
				$count_data_hsh_ref->{'region'}{'total_polyGT_length'}{$dnstream_total_polyGT_length}{'dnstream'}++;
				$count_data_hsh_ref->{'region'}{'region_length'}{$dnstream_region_length}{'dnstream'}++;
				$count_data_hsh_ref->{'region'}{'stretch_frac'}{$dnstream_stretch_frac}{'dnstream'}++;

				$count_data_hsh_ref->{'region'}{'max_polyGT_length'}{$upstream_max_polyGT_length}{'upstream'}++;
				$count_data_hsh_ref->{'region'}{'max_polyA_length'}{$upstream_max_polyA_length}{'upstream'}++;
				$count_data_hsh_ref->{'region'}{'total_polyA_length'}{$upstream_total_polyA_length}{'upstream'}++;
				$count_data_hsh_ref->{'region'}{'total_polyGT_length'}{$upstream_total_polyGT_length}{'upstream'}++;
				$count_data_hsh_ref->{'region'}{'region_length'}{$upstream_region_length}{'upstream'}++;
				$count_data_hsh_ref->{'region'}{'stretch_frac'}{$upstream_stretch_frac}{'upstream'}++;

				$count_data_hsh_ref->{'region'}{'max_polyGT_length'}{$head_max_polyGT_length}{'head'}++;
				$count_data_hsh_ref->{'region'}{'max_polyA_length'}{$head_max_polyA_length}{'head'}++;
				$count_data_hsh_ref->{'region'}{'total_polyA_length'}{$head_total_polyA_length}{'head'}++;
				$count_data_hsh_ref->{'region'}{'total_polyGT_length'}{$head_total_polyGT_length}{'head'}++;
				$count_data_hsh_ref->{'region'}{'region_length'}{$head_region_length}{'head'}++;
				$count_data_hsh_ref->{'region'}{'stretch_frac'}{$head_stretch_frac}{'head'}++;

				$count_data_hsh_ref->{'region'}{'max_polyGT_length'}{$tail_max_polyGT_length}{'tail'}++;
				$count_data_hsh_ref->{'region'}{'max_polyA_length'}{$tail_max_polyA_length}{'tail'}++;
				$count_data_hsh_ref->{'region'}{'total_polyA_length'}{$tail_total_polyA_length}{'tail'}++;
				$count_data_hsh_ref->{'region'}{'total_polyGT_length'}{$tail_total_polyGT_length}{'tail'}++;
				$count_data_hsh_ref->{'region'}{'region_length'}{$tail_region_length}{'tail'}++;
				$count_data_hsh_ref->{'region'}{'stretch_frac'}{$tail_stretch_frac}{'tail'}++;

				$status = "trimmed";
				$read_trimmed++;
				my $dnstream_polyA_frac = sprintf "%.2f", $dnstream_total_polyA_length/$dnstream_region_length;
				$out_rd_name = join " ", ($revised_rd_name, "$status|trimmed_length=$tail_region_length|dnstream_polyA_frac=$dnstream_polyA_frac|dnstream_seq=$dnstream_seq");
				$out_qual_header = $qual_header;
				my $head_qual_value = substr ($qual_value, 0, $cut_pos);
				my $tail_qual_value = substr ($qual_value, $cut_pos);
				my $upstream_qual = substr ($head_qual_value, length($head_qual_value)-20);
				my $dnstream_qual = substr ($tail_qual_value, 0, 20);
				my $junction_qual = join "", ($upstream_qual,$dnstream_qual);

				$out_qual_value = $head_qual_value;
				$out_rd_seq = $head_seq;
				
				if ($rescue_tail_nt > 0) {
					my $rescue_pos = $rescue_tail_nt;
					$rescue_pos = length ($tail_seq) if $rescue_pos > length ($tail_seq);
					$rescue_tail_seq = substr ($tail_seq, 0, $rescue_pos);
					$rescue_tail_qual = substr ($tail_qual_value, 0, $rescue_pos);
				}
				
				print FULLTAIL "$out_rd_name\n";
				print FULLTAIL "$tail_seq\n";
				print FULLTAIL "$out_qual_header\n";
				print FULLTAIL "$tail_qual_value\n";

				print JUNCTION "$out_rd_name\n";
				print JUNCTION "$junction_seq\n";
				print JUNCTION "$out_qual_header\n";
				print JUNCTION "$junction_qual\n";
				
			} else {
				if (length ($rd_seq) < 20) {
					$status='head_too_short';
					$read_short++;
					print SHORT "$rd_name\n";
					print SHORT "$rd_seq\n";
					print SHORT "$qual_header\n";
					print SHORT "$qual_value\n";
					next;
				} else {
					$status = "untrimmed";
					$read_untrimmed++;
					$out_rd_name = join " ", ($revised_rd_name, "$status|trimmed_length=-1|dnstream_polyA_frac=-1|dnstream_seq=");
					$out_rd_seq = $rd_seq;
					$out_qual_header = $qual_header;
					$out_qual_value = $qual_value;
					my $end3_seq = substr ($rd_seq, length($rd_seq)-50);
					my $end3_qual = substr ($qual_value, length($qual_value)-50);
					print UNTRIMMED "$out_rd_name\n";
					print UNTRIMMED "$end3_seq\n";
					print UNTRIMMED "$qual_header\n";
					print UNTRIMMED "$end3_qual\n";
				}
			}

			$count_data_hsh_ref->{'status'}{$status}++;

			print TRIMMEDFQ "$out_rd_name\n";
			print TRIMMEDFQ "$out_rd_seq\n";
			print TRIMMEDFQ "$out_qual_header\n";
			print TRIMMEDFQ "$out_qual_value\n";
			
			print RETAINTAILFQ "$out_rd_name\n";
			print RETAINTAILFQ $out_rd_seq.$rescue_tail_seq."\n";
			print RETAINTAILFQ "$out_qual_header\n";
			print RETAINTAILFQ $out_qual_value.$rescue_tail_qual."\n";
			
		} else {
			&reportAndLogStatus("WARNING: Strand information of read was not found. Check readname format. Quitting.", 10, "\n");#->927
			die;
		}
	}
	my $read_trimmed_pct = sprintf "%.2f", 100*$read_trimmed/$read_proc;
	my $read_untrimmed_pct = sprintf "%.2f", 100*$read_untrimmed/$read_proc;
	my $read_short_pct = sprintf "%.2f", 100*$read_short/$read_proc;
	&reportAndLogStatus("$read_proc reads processed (skip=$read_skip). trimmed=$read_trimmed\[$read_trimmed_pct%\] untrimmed=$read_untrimmed\[$read_untrimmed_pct%\] too_short(<20nt)=$read_short\[$read_short_pct%\]", 10, "\n");#->927

	close FILEIN;
	close TRIMMEDFQ;
	close FULLTAIL;
	close JUNCTION;
	close SHORT;
	close RETAINTAILFQ;
	close READIDLOG;
	
	return ($count_data_hsh_ref);
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|92
#	secondaryAppearInSection: >none
#	input: none
#	output: $fastq_path, $frac_check_size_max, $match_5end_size, $max_mismatch, $max_proc_read, $min_stretch_frac, $out_dir, $out_prefix, $rescue_tail_nt, $revise_read_prefix, $subsample_frac, $win_size
#	toCall: my ($fastq_path, $win_size, $max_mismatch, $match_5end_size, $min_stretch_frac, $frac_check_size_max, $subsample_frac, $max_proc_read, $rescue_tail_nt, $revise_read_prefix, $out_prefix, $out_dir) = &readParameters();
#	calledInLine: 95
#....................................................................................................................................................#
	
	my ($fastq_path, $win_size, $max_mismatch, $match_5end_size, $min_stretch_frac, $frac_check_size_max, $subsample_frac, $max_proc_read, $rescue_tail_nt, $revise_read_prefix, $out_prefix, $out_dir);
	
	$win_size = 12;
	$max_mismatch = 2;
	$match_5end_size = 4;
	$min_stretch_frac = 0.5;
	$frac_check_size_max = 100;
	$subsample_frac = 1;
	$rescue_tail_nt = 10;
	$max_proc_read = undef;
	$revise_read_prefix = undef;

	GetOptions 	(
		"fastq_path=s"					=>	\$fastq_path,
		"out_dir=s"						=>	\$out_dir,
		"out_prefix:s"					=>	\$out_prefix,
		"win_size:i"					=>	\$win_size,
		"max_mismatch:i"				=>	\$max_mismatch,
		"subsample_frac:f"			=>	\$subsample_frac,
		"max_proc_read:i"				=>	\$max_proc_read,
		"min_stretch_frac:f" 		=>	\$min_stretch_frac,
		"match_5end_size:i"			=>	\$match_5end_size,
		"frac_check_size_max:i"		=>	\$frac_check_size_max,
		"revise_read_prefix:s"		=>	\$revise_read_prefix,
		"rescue_tail_nt:i"			=>	\$rescue_tail_nt,
		'help'							=>	sub { HelpMessage(0) },
	) or HelpMessage(1);

	my $opt_check_hsh_ref = {
		'fastq_path' => $fastq_path,
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

	if ($subsample_frac > 1 or $subsample_frac <= 0) {
		print "WARNING: subsample_frac must be <=1 and >0. Please check this help message for required options\n";
		print "\n";
		HelpMessage(1);
	}

	if (defined $max_proc_read) {
		if ($max_proc_read <= 0) {
			print "WARNING: max_proc_read must be >0. Please check this help message for required options\n";
			print "\n";
			HelpMessage(1);
		}
	}

	if ($min_stretch_frac > 1) {
		print "WARNING: min_stretch_frac must <1. Please check this help message for required options\n";
		print "\n";
		HelpMessage(1);
	}

	if ($max_mismatch > 3) {
		print "WARNING: max_mismatch must <3. Please check this help message for required options\n";
		print "\n";
		HelpMessage(1);
	}

	if ($win_size > 20 or $win_size < 8) {
		print "WARNING: win_size must >=8 & <=20. Please check this help message for required options\n";
		print "\n";
		HelpMessage(1);
	}

	if ($match_5end_size > $win_size/3) {
		print "WARNING: match_5end_size must <= win_size/3. Please check this help message for required options\n";
		print "\n";
		HelpMessage(1);
	}
	
	#---check file
	my $file_check_hsh_ref = {
		'fastq_path' => $fastq_path,
	};
	
	foreach my $option_name (keys %{$file_check_hsh_ref}) {
		my $file_path = $file_check_hsh_ref->{$option_name};
		if (not -s $file_path) {
			die "Quitting: File $option_name does not exists at $file_path";
		}
	}

	my @fastq_path_ary = split /\/+/, $fastq_path;
	my $fastq_name = $fastq_path_ary[-1];
	if ($fastq_name =~ m/(.+)\.fq\.gz$/ or $fastq_name =~ m/(.+)\.fq$/ or $fastq_name =~ m/(.+)\.fastq\.gz$/ or $fastq_name =~ m/(.+)\.fastq$/) {
		$out_prefix = $1 if not defined $out_prefix;
	} else {
		die "Quitting: file_path must be ended as fq or fastq or fq.gz or fastq.gz\n";
	}

	chop $out_dir if ($out_dir =~ m/\/$/); #---remove the last slash
	system "mkdir -p -m 755 $out_dir/";
	
	return($fastq_path, $win_size, $max_mismatch, $match_5end_size, $min_stretch_frac, $frac_check_size_max, $subsample_frac, $max_proc_read, $rescue_tail_nt, $revise_read_prefix, $out_prefix, $out_dir);

}
sub reportAndLogStatus {
#....................................................................................................................................................#
#	subroutineCategory: log
#	dependOnSub: currentTime|186
#	appearInSub: getPolyStretchRegex|276, readFastq|557
#	primaryAppearInSection: 2_defineout_dirPath|108
#	secondaryAppearInSection: 3_process|141
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportAndLogStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 132, 134, 136, 299, 300, 308, 309, 624, 785, 792
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->186
	print $tmplog_fh "[".&currentTime()."] ".$message.$lineEnd if $lineEnd ne "\r";#->186
	
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
sub trimTailPolyA {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: readFastq|557
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 3_process|141
#	input: $seq_to_trim
#	output: $tail_trim_length, $tail_trim_seq
#	toCall: my ($tail_trim_length, $tail_trim_seq) = &trimTailPolyA($seq_to_trim);
#	calledInLine: 657, 659
#....................................................................................................................................................#
	my ($seq_to_trim) = @_;
	
	my $tail_trim_length = 0;
	my $tail_trim_seq = 'NA';
	if ($seq_to_trim =~ m/((A)+(A\w{0,1}A\w{0,1}A)+\w{0,3})$/) {
		$tail_trim_seq = $1;
		$tail_trim_length = length($tail_trim_seq);
	}

	return ($tail_trim_length, $tail_trim_seq);
}

exit;


















































