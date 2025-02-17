#!/usr/bin/perl -w
use strict;
use File::Basename;
use File::Spec::Functions qw(rel2abs abs2rel);
my ($transcript_bed, $chrom_size_path, $outPrefix, $outDir, $bedGraphToBigWig_bin) = @ARGV;

system "mkdir -pm 755 $outDir";
if ($transcript_bed =~ m/gz$/) {
	open (FILEIN, " gzip -dc $transcript_bed|");
} else {
	open (FILEIN, "<", $transcript_bed);
}

my $end_info_hsh_ref = {
	'end5' => {
		'bed' => "$outDir/$outPrefix.end5.bed",
		'bigwig' => {
			"+" => "$outDir/$outPrefix.end5.fwd.bw",
			"-" => "$outDir/$outPrefix.end5.rev.bw",
		},
	},
	'end3' => {
		'bed' => "$outDir/$outPrefix.end3.bed",
		'bigwig' => {
			"+" => "$outDir/$outPrefix.end3.fwd.bw",
			"-" => "$outDir/$outPrefix.end3.rev.bw",
		},
	},
};

foreach my $end (keys %{$end_info_hsh_ref}) {
	open $end_info_hsh_ref->{$end}{'fh'}, "| sort --parallel=10 -k1,1 -k2,2n -k6,6 | bedtools merge -s -d -1 -c 4,5,6 -o count,count,distinct >$end_info_hsh_ref->{$end}{'bed'}";
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
		print {$end_info_hsh_ref->{$end}{'fh'}} join "", (join "\t", ($chrom, $bedStart, $bedEnd, $ID, $score, $strand)), "\n";
	}
}
close FILEIN;

foreach my $end (keys %{$end_info_hsh_ref}) {
	close $end_info_hsh_ref->{$end}{'fh'};
	my $bed = $end_info_hsh_ref->{$end}{'bed'};
	foreach my $strand (keys %{$end_info_hsh_ref->{$end}{'bigwig'}}) {
		my $bigwig_path = $end_info_hsh_ref->{$end}{'bigwig'}{$strand};
		my $tmp_bed = "$outDir/tmp.bed";
		my $awk_cmd = "awk -F'\\t' -v OFS='\\t' '\$6 == \"$strand\" {print \$1, \$2, \$3, \$5}' $bed >$tmp_bed";
		my $bedGraphToBigWig_cmd = "$bedGraphToBigWig_bin $tmp_bed $chrom_size_path $bigwig_path";
		system $awk_cmd;
		system $bedGraphToBigWig_cmd;
		system "rm $tmp_bed";
	}
}
