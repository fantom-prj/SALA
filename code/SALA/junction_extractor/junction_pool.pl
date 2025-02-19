#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw(sum shuffle min max);
use Getopt::Long;

# Define variables for command-line arguments
my $outDir;
my $out_tag;
my $find_str;

# Parse command-line arguments
GetOptions(
    "outDir=s"   => \$outDir,
    "outTag=s"   => \$out_tag,
    "findStr=s"  => \$find_str,
) or die "Error in command line arguments\n";

# Validate required arguments
die "Missing --outDir parameter\n" unless defined $outDir;
die "Missing --outTag parameter\n" unless defined $out_tag;
die "Missing --findStr parameter\n" unless defined $find_str;

system "mkdir -pm 755 $outDir";
system "cp $0 $outDir";

chomp (my @file_ary = `find $find_str`);
my $junct_info_hsh_ref = {};
foreach my $info_path (@file_ary) {
	print "reading $info_path\n";
	open INDIVJUNCINFO, "gzip -dc $info_path|";
	<INDIVJUNCINFO>;
	while (<INDIVJUNCINFO>) {
		chomp;
		#junct_ID	chrom	junct_start	junct_end	strand	splicing_site	canonical	total_count	hi_qual_count	max_mapq	avg_mapq	max_score_per_nt	avg_score_per_nt	max_score.donor.-3	max_score.donor.-2	max_score.donor.-1	max_score.acceptor.1	max_score.acceptor.2	max_score.acceptor.3	avg_score.donor.-3	avg_score.donor.-2	avg_score.donor.-1	avg_score.acceptor.1	avg_score.acceptor.2	avg_score.acceptor.3
		#chr1_100007156_100009287_+	chr1	100007156	100009287	+	GT-AG	Y	1	0	60	60.00	25.17	25.17	28	29	29	28	20	17	28.00	29.00	29.00	28.00	20.00	17.00
		#chr1_100007156_100011364_+	chr1	100007156	100011364	+	GT-AG	Y	38	11	60	60.00	28.33	18.78	29	28	28	28	29	28	21.29	20.76	19.55	18.26	17.00	15.79
		#chr1_100009323_100011364_+	chr1	100009323	100011364	+	GT-AG	Y	1	0	60	60.00	22.67	22.67	30	29	27	18	16	16	30.00	29.00	27.00	18.00	16.00	16.00
		#chr1_100011533_100015301_+	chr1	100011533	100015301	+	GT-AG	Y	40	15	60	60.00	31.00	21.17	29	30	30	30	33	34	20.75	21.75	21.62	20.70	20.88	21.30
		#chr1_100011533_100022385_+	chr1	100011533	100022385	+	GT-AG	Y	1	1	60	60.00	24.00	24.00	26	29	25	23	20	21	26.00	29.00	25.00	23.00	20.00	21.00
		#chr1_100015420_100017681_+	chr1	100015420	100017681	+	GT-AG	Y	26	13	60	60.00	30.83	21.19	34	31	30	30	30	30	22.92	21.88	21.23	20.62	20.46	20.04
		my ($junct_ID, $chrom, $junct_start, $junct_end, $strand, $splicing_site, $canonical, $total_count, $hi_qual_count, $max_mapq, $avg_mapq, $max_score_per_nt, $avg_score_per_nt, $max_score_0, $max_score_1, $max_score_2, $max_score_3, $max_score_4, $max_score_5, $avg_score_0, $avg_score_1, $avg_score_2, $avg_score_3, $avg_score_4, $avg_score_5) = split /\t/;
		if (not exists $junct_info_hsh_ref->{$junct_ID}) {
			$junct_info_hsh_ref->{$junct_ID}{'count'}{'hi_qual'} = 0;
			$junct_info_hsh_ref->{$junct_ID}{'count'}{'total'} = 0;
			$junct_info_hsh_ref->{$junct_ID}{'mapq'}{'max'} = 0;
			$junct_info_hsh_ref->{$junct_ID}{'mapq'}{'sum'} = 0;
			foreach my $i (0..5) {
				$junct_info_hsh_ref->{$junct_ID}{'qual'}{'max'}{$i} = 0;
				$junct_info_hsh_ref->{$junct_ID}{'qual'}{'sum'}{$i} = 0;
			}
		}

		$junct_info_hsh_ref->{$junct_ID}{'splicing_site'} = $splicing_site;
		$junct_info_hsh_ref->{$junct_ID}{'canonical'} = $canonical;
		$junct_info_hsh_ref->{$junct_ID}{'count'}{'total'} += $total_count;
		$junct_info_hsh_ref->{$junct_ID}{'count'}{'hi_qual'} += $hi_qual_count;
		$junct_info_hsh_ref->{$junct_ID}{'mapq'}{'max'} = $max_mapq if $max_mapq > $junct_info_hsh_ref->{$junct_ID}{'mapq'}{'max'};
		$junct_info_hsh_ref->{$junct_ID}{'mapq'}{'sum'} += $total_count*$avg_mapq;

		my @max_score_ary = ($max_score_0, $max_score_1, $max_score_2, $max_score_3, $max_score_4, $max_score_5);
		my @avg_score_ary = ($avg_score_0, $avg_score_1, $avg_score_2, $avg_score_3, $avg_score_4, $avg_score_5);
		foreach my $i (0..5) {
			my $max_score = $max_score_ary[$i];
			my $avg_score = $avg_score_ary[$i];
			$junct_info_hsh_ref->{$junct_ID}{'qual'}{'max'}{$i} = $max_score if $max_score > $junct_info_hsh_ref->{$junct_ID}{'qual'}{'max'}{$i};
			$junct_info_hsh_ref->{$junct_ID}{'qual'}{'sum'}{$i} += $total_count*$avg_score;
		}
	}
	close INDIVJUNCINFO;
}

open HIQJUNCTION, "| sort -k1,1 -k2,2n -k6,6 >$outDir/$out_tag.hi_qual.junct.bed";
open OUTJUNCTION, "| gzip -c >$outDir/$out_tag.junct.info.tsv.gz";
print OUTJUNCTION join "", (join "\t", ('junct_ID', 'chrom', 'junct_start', 'junct_end', 'strand', 'splicing_site', 'canonical', 'total_count', 'hi_qual_count', 'max_mapq', 'avg_mapq', 'max_score_per_nt', 'avg_score_per_nt', 'max_score.donor.-3', 'max_score.donor.-2', 'max_score.donor.-1', 'max_score.acceptor.1', 'max_score.acceptor.2', 'max_score.acceptor.3', 'avg_score.donor.-3', 'avg_score.donor.-2', 'avg_score.donor.-1', 'avg_score.acceptor.1', 'avg_score.acceptor.2', 'avg_score.acceptor.3')), "\n";
foreach my $junct_ID (sort keys %{$junct_info_hsh_ref}) {
	my ($chrom, $junct_start, $junct_end, $strand) = split /\_/, $junct_ID;
	my $total_count = $junct_info_hsh_ref->{$junct_ID}{'count'}{'total'};
	my $hi_qual_count = $junct_info_hsh_ref->{$junct_ID}{'count'}{'hi_qual'};
	my $sum_mapq = $junct_info_hsh_ref->{$junct_ID}{'mapq'}{'sum'};
	my $avg_mapq = sprintf "%.2f", $sum_mapq/$total_count;
	my $max_mapq = $junct_info_hsh_ref->{$junct_ID}{'mapq'}{'max'};

	my $splicing_site = $junct_info_hsh_ref->{$junct_ID}{'splicing_site'};
	my $canonical = $junct_info_hsh_ref->{$junct_ID}{'canonical'};
	
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

	print OUTJUNCTION join "", (join "\t", ($junct_ID, $chrom, $junct_start, $junct_end, $strand, $splicing_site, $canonical, $total_count, $hi_qual_count, $max_mapq, $avg_mapq, $max_score_per_nt, $avg_score_per_nt, @max_score_ary, @avg_score_ary)), "\n";
	if ($hi_qual_count >= 1 or $canonical eq 'Y') {
		print HIQJUNCTION join "", (join "\t", ($chrom, $junct_start, $junct_end, $junct_ID."|".$splicing_site, $total_count, $strand)), "\n";
	}
}
close OUTJUNCTION;
close HIQJUNCTION;
