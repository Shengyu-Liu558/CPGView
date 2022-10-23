#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin qw( $Bin );
use Cwd 'abs_path';

my $project_id  = shift;
my $dir_work    = shift;
my $file_fas    = shift || "";
my $para_misa   = shift || " 1-10 2-6 3-5 4-5 5-5 6-5 ";
my $para_trf    = shift || " 2 7 7 80 10 50 500 -f -d -m ";
my $para_vmatch = shift || " -f -p -h 3 -l 30 ";

print STDERR "\nExample run:  perl $0 1234 /tmp/o20 sample.fas \" 1-10 2-6 3-5 4-5 5-5 6-5 \" \" 2 7 7 80 10 50 500 -f -d -m \" \" -f -p -h 3 -l 30 \"\n";
die "\nUsage: $0 project_id dir_working file_fas para_misa para_trf para_vmatch\n\n" unless (-e $file_fas);

my $projectHomeDir = abs_path( $dir_work ); #absolute path to the projectHomeDir
`mkdir $projectHomeDir` unless (-e $projectHomeDir);
#printCMT($cmd); `$cmd`;

#print STDERR "\n-----step 9 find repetative elements-----\n\n";
	my $cmd_misa    = "$Bin/exe/misa/misa.pl";
        my $cmd_trf     = "$Bin/exe/trf409.legacylinux64";
	my $cmd_repfind = "$Bin/exe/vmatch-2.3.0-Linux_x86_64-64bit/repfind.pl";
	my $cmd         = "";

        &printCMT("Mode 9: Additional Analysis: SSR, repeat elements");
        my $projectFASTAFile   = $projectHomeDir."/".$project_id.".fas";
        my $projectMISAIniFile = $projectHomeDir."/".$project_id.".misa.ini";
	`cp $file_fas $projectFASTAFile`;

################## misa #######################
        &printCMT("Start misa");
	open (INI, ">$projectMISAIniFile");
	print INI "definition(unit_size,min_repeats):        $para_misa\n";
	print INI "interruptions(max_difference_between_2_SSRs):    100";
	close(INI);
        $cmd = "perl $cmd_misa $projectMISAIniFile $projectFASTAFile";
        &printCMT($cmd); `$cmd`;

################## trf #######################
        &printCMT("Start trf409");

	my $postfix_trf = $para_trf;
	$postfix_trf    =~ s/-[a-z]//g;
	$postfix_trf    =~ s/^\s+//g;
	$postfix_trf    =~ s/\s+$//g;
	$postfix_trf    =~ s/\s+/\./g;
	my $projectTRFFile = $projectFASTAFile.".".$postfix_trf.".dat";

        $cmd = "(cd /tmp/ && $cmd_trf $projectFASTAFile $para_trf)";
        &printCMT($cmd); `$cmd`;

	print STDERR "\n\n$projectTRFFile\n\n";
        $cmd = "(cd $projectHomeDir && ln -s $projectTRFFile $project_id.trf.txt)";
        &printCMT($cmd); `$cmd`;

        $cmd = "mv /tmp/$project_id* $projectHomeDir";
        &printCMT($cmd); `$cmd`;

	#####filter output###################
	&printCMT("--STARTing filtering-----\n\n");

	my $cutoff = 7;
	#my $file_before = $projectTRFFile.".beforefiltering";
	my @lines = `cat $projectHomeDir/$project_id.trf.txt`;
	my $str = &filterTRFResults(@lines, $cutoff);

	my $file_after = $projectTRFFile.".afterfiltering";
	open(OUT, ">$file_after");
	print OUT $str;
	close(OUT);	

        $cmd = "(cd $projectHomeDir && ln -s $file_after $project_id.trf.afterfiltering.txt)";
        &printCMT($cmd); `$cmd`;
################## vmatch #######################
        my $projectRepfind = $projectHomeDir."/".$project_id."_vmatch.txt";
        &printCMT("Start vmatch");
        $cmd = "perl $cmd_repfind $para_vmatch $projectFASTAFile > $projectRepfind";
        #$cmd = "perl $cmd_repfind $para_vmatch $projectFASTAFile > $projectRepfind_prefiter";
        &printCMT($cmd); `$cmd`;

###################################
sub printCMT {
	my ($line) = @_;
	print STDERR $line,"\n\n";
}

###################################
sub filterTRFResults {

	my (@lines, $cutoff) = @_;

	my $str = "";
	foreach my $line (@lines) {
		$line =~ s/^\s+//g;
		$line =~ s/\s+$//g;
		next if ($line eq "");
		if ($line =~ /^\d\d+\s\d+\s\d+/) {
		my @ws = split(/ /, $line);
		if ($ws[2] < 7) {
			print STDERR "filtered out, $line\n";
			next;
		} 
		}
		#print $line,"\n";
		$str .= $line."\n";	
	}
	return $str;
}
