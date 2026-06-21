#!/usr/bin/perl -w
use strict;
use File::Basename;
use File::Spec::Functions qw(rel2abs abs2rel);
my ($transcript_bed, $merge_d_end5, $merge_d_end3, $base_rng_end5, $base_rng_end3, $outPrefix, $outDir) = @ARGV;

system "mkdir -pm 755 $outDir";
if ($transcript_bed =~ m/gz$/) {
	open (FILEIN, " gzip -dc $transcript_bed|");
} else {
	open (FILEIN, "<", $transcript_bed);
}

my $end_info_hsh_ref = {
	'end5' => {
		'base_rng' => $base_rng_end5,
		'merge_d' => $merge_d_end5,
		'merge_bed' => "$outDir/$outPrefix.end5.merge.bed",
		'signal_bed' => "$outDir/$outPrefix.end5.signal.bed",
		'filter_bed' => "$outDir/$outPrefix.end5.merge.filter.bed.bgz",
	},
	'end3' => {
		'base_rng' => $base_rng_end3,
		'merge_d' => $merge_d_end3,
		'merge_bed' => "$outDir/$outPrefix.end3.merge.bed",
		'signal_bed' => "$outDir/$outPrefix.end3.signal.bed",
		'filter_bed' => "$outDir/$outPrefix.end3.merge.filter.bed.bgz",
	},
};

foreach my $end (keys %{$end_info_hsh_ref}) {
	open $end_info_hsh_ref->{$end}{'bed_fh'}, "| sort -k1,1 -k2,2n -k6,6 | bedtools merge -i stdin -s -d -1 -c 4,5,6 -o count,count,distinct >$end_info_hsh_ref->{$end}{'signal_bed'}";
}

while (<FILEIN>) {
	next if $_ =~ m/^#/;
	chomp;
	my ($chrom, $chromStart, $chromEnd, $ID, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\t/;
	
	my $tmp_hsh_ref = {};
	if ($strand eq '+') {
		$tmp_hsh_ref->{'end5'} = $chromStart + 1;
		$tmp_hsh_ref->{'end3'} = $chromEnd;
	} elsif ($strand eq '-') {
		$tmp_hsh_ref->{'end3'} = $chromStart + 1;
		$tmp_hsh_ref->{'end5'} = $chromEnd;
	} else {
		die "transcript $ID has no strand\n";
	}
	
	foreach my $end (keys %{$tmp_hsh_ref}) {
		my $bedEnd = $tmp_hsh_ref->{$end};
		my $bedStart = $bedEnd-1;
		print {$end_info_hsh_ref->{$end}{'bed_fh'}} join "", (join "\t", ($chrom, $bedStart, $bedEnd, $ID, $score, $strand)), "\n";
	}
}
close FILEIN;

foreach my $end (keys %{$end_info_hsh_ref}) {
	close $end_info_hsh_ref->{$end}{'bed_fh'};
	my $signal_bed = $end_info_hsh_ref->{$end}{'signal_bed'};
	my $base_rng = $end_info_hsh_ref->{$end}{'base_rng'};
	my $merge_d = $end_info_hsh_ref->{$end}{'merge_d'};
	my $merge_bed = $end_info_hsh_ref->{$end}{'merge_bed'};
	my $filter_bed = $end_info_hsh_ref->{$end}{'filter_bed'};
	
	open MERGEBED, "| sort -k1,1 -k2,2n >$merge_bed";
	open BEDTOOLSMERGE, "cat $signal_bed $base_rng | sort -k1,1 -k2,2n | bedtools merge -i stdin -s -d $merge_d -c 4,5,6 -o count,count,distinct |";
	while (<BEDTOOLSMERGE>) {
		chomp;
		my ($chrom, $chromStart, $chromEnd, undef, $score, $strand) = split /\t/;
		my $ID = join "_", ($chrom, $chromStart, $chromEnd, $strand);
		my $summit = int($chromStart + ($chromEnd - $chromStart)/2);
		my $thickStart = $summit-1;
		my $thickEnd = $summit;
		my $blockCount = 1;
		my $itemRgb = '.';
		my $blockSizes = $chromEnd - $chromStart;
		my $blockStarts = 0;
		print MERGEBED join "", (join "\t", ($chrom, $chromStart, $chromEnd, $ID, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts)), "\n";
	}
	close MERGEBED;
	close BEDTOOLSMERGE;
	
	system "bedtools intersect -s -u -wa -a $merge_bed -b $signal_bed | bgzip -c >$filter_bed";
	system "tabix -p bed $filter_bed";
	
	system "cat $signal_bed | bgzip -c >$signal_bed.bgz";
	system "tabix -p bed $signal_bed.bgz";
}
