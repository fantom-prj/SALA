#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

# Define variables with default values
my $outDir = "";
my $run_scafe_script_path = './code/SCAFEv1.0.1/scripts/scafe.tool.bk.bam_to_ctss';
my $in_lib_list_path = "";
my $genome = "hg38.gencode_v39";
my $max_thread = 1;
my $TSS_mode = "softclip";
my $unencoded_G_upstrm_nt = 3;
my $max_softclip_length = 3;

# Parse command-line arguments
GetOptions(
    "outDir=s"                => \$outDir,
    "run_scafe_script_path=s"  => \$run_scafe_script_path,
    "in_lib_list_path=s"       => \$in_lib_list_path,
    "genome=s"                 => \$genome,
    "max_thread=i"             => \$max_thread,
    "TSS_mode=s"               => \$TSS_mode,
    "unencoded_G_upstrm_nt=i"  => \$unencoded_G_upstrm_nt,
    "max_softclip_length=i"    => \$max_softclip_length
) or die "Error parsing command-line arguments\n";

# Create necessary directories
system "mkdir -pm 755 $outDir";
system "cp $0 $outDir";  # Save script for reproducibility

# Start and finish time logs
my $start_time_log_path = "$outDir/00_start.time.log.txt";
my $finish_time_log_path = "$outDir/00_finish.time.log.txt";
system "echo \"========== All runs are started at \$(date) ==========\" > $start_time_log_path\n";
system "echo \"========== All runs are started at \$(date) ==========\" > $finish_time_log_path\n";

# Read BAM list file
my $lib_info_hsh_ref = {};
open my $BAMLIST, "<", $in_lib_list_path or die "Cannot open $in_lib_list_path: $!\n";

while (<$BAMLIST>) {
    chomp;
    next if /^#/;  # Skip comment lines
    my ($index, $libID, $bamPath) = split /\t/;
    
    die "Error: BAM file does not exist for $libID at $bamPath\n" if not -s $bamPath;

    $lib_info_hsh_ref->{$libID}{'bamPath'} = $bamPath;
}
close $BAMLIST;

# Count libraries
my $num_lib = keys %{$lib_info_hsh_ref};
print "$num_lib libraries read\n";

# Process each library
foreach my $libID (sort keys %{$lib_info_hsh_ref}) {
    my $bamPath = $lib_info_hsh_ref->{$libID}{'bamPath'};
    my $outputPrefix = $libID;

    # Define log file paths
    my $stderr_path = "$outDir/$outputPrefix.stderr.txt";
    my $stdout_path = "$outDir/$outputPrefix.stdout.txt";

    my $cmd = join " ", (
        "$run_scafe_script_path",
        "--TSS_mode=$TSS_mode",
        "--bamPath=$bamPath",
        "--unencoded_G_upstrm_nt=$unencoded_G_upstrm_nt",
        "--max_thread=$max_thread",
        "--genome=$genome",
        "--max_softclip_length=$max_softclip_length",
        "--outputPrefix=$outputPrefix",
        "--outDir=$outDir",
        "2> $stderr_path 1> $stdout_path"  # Redirect stderr and stdout
    );

    print "$libID is started at " . localtime() . "\n";
    system($cmd) == 0 or warn "Warning: Execution failed for $libID\n";
    print "$libID is finished at " . localtime() . "\n";
}

