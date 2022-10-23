#!/usr/bin/perl
#
while (<>) {
	my $id = $_;
	$id =~ s/\.gb\s+//;
	print "perl plasdrawmap.pl ../gbks/$id.gb $id /tmp/results/$id\n";

}
