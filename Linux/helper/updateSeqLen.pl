#!/usr/bin/perl

my $file_gb = shift;
my $file_ka = shift;

my $str = "chr - cpg 1 0 length vblack";

my $line = `cat $file_gb | grep LOCUS`;

$line =~ s/\s+/:/g;
print STDERR $line,"\n";
my @ws   = split(/:/, $line);
#print $ws[2],"\n";
$str =~ s/length/$ws[2]/;

print $str;

