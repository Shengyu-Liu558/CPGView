#!/usr/local/bin/perl

use strict;
use FindBin qw($Bin);
use Cwd 'abs_path';

my $file_in        = shift;
my $projectid      = shift;
my $workingdir_rel = shift;
#my $file_trf       = shift;
#my $file_misa      = shift;
#my $file_vmatch    = shift;
#my $file_ir      = shift;

die "\n\nUsage: $0 file_gb projectid workingdir\n\nExample run: perl $0 (gb:NC_000932 or test/example.gb) 1234 /tmp/o1\n\n " unless ((-e $file_in || $file_in =~ /gb:/) && $projectid ne "" && $workingdir_rel ne "");

=begin comments
if ($projectid eq "") {
    $projectid = time();
    $projectid =~ s/\s+//g;
}

if ($workingdir_rel eq "") {
	$workingdir_rel = "/tmp/d".$projectid;
}
=cut

$workingdir_rel  =~ s/\/+$//;
my $workingdir   = abs_path ($workingdir_rel);
`mkdir $workingdir; echo "mkdir working directory\n"` unless (-e $workingdir);

my $file_gb      = $workingdir."/".$projectid.".gb";
my $file_fas     = $workingdir."/".$projectid.".fas";
my $file_image1  = $workingdir."/".$projectid."_image1";
my $file_image2  = $workingdir."/".$projectid."_image2.png";
my $file_image3  = $workingdir."/".$projectid."_image3.png";
my $file_image4  = $workingdir."/".$projectid."_image4.png";
my $file_image5  = $workingdir."/".$projectid."_image5.png";

my $file_circos_trf      = $workingdir."/".$projectid."_c_trf.txt";
my $file_circos_misa     = $workingdir."/".$projectid."_c_misa.txt";
my $file_circos_vmatch_d = $workingdir."/".$projectid."_c_vmatchd.txt";
my $file_circos_vmatch_p = $workingdir."/".$projectid."_c_vmatchp.txt";
my $file_circos_ir       = $workingdir."/".$projectid."_c_ir.txt";
my $file_circos_cc       = $workingdir."/".$projectid."_c_cc.txt";

#cat karyotype.txt | sed 's/ /\t/g' | cut -f3,5,6,7
#
#die "\n\nUsage: $0 workingdir projectid file_gb file_trf file_misa file_vmatch file_ir\n\nExample run: perl $0 /tmp/o10 1234 test/1234.gb test/154773558769465.trf.txt test/154773558769465.fas.misa test/154773558769465_vmatch.txt test/15478017817688_ir.txt\n\n " unless (-e $file_gb && -e $file_trf && -e $file_misa && -e $file_vmatch);
#die "Usage: $0 $workingdir $projectid $file_gb $file_trf $file_misa $file_vmatch $file_ir\n\n" unless (-e $file_gb && $file_trf $file_msia, $file_vmatch && $file_ir);

##drawgenemap
#my $drawgenemap = "drawgenemap";
#my $drawgenemap = $Bin."/GeneMap-1.1.1/bin/drawgenemap";
#my $cmd = $drawgenemap." --force_circular --infile $file_gb --outfile=$file_image1 --density=300";
#print STDERR $cmd, "\n\n"; `$cmd`;
#my $python = "python";
my $python = "/biodata4/home/cliu/cpgview/bin/python";
#Location of Python，The user needs to change the path to the location of Python in the user's local computer.

####if the input is an accession number from GenBank
if ($file_in =~ /gb:/) {
    $file_in =~ s/gb://;
    my $cmd = "$python $Bin/helper/downloadGB.py $file_in $file_gb";
    print STDERR $cmd, "\n\n"; `$cmd`;
} else {
    my $cmd = `cp $file_in $file_gb`; 
    print STDERR $cmd, "\n\n"; `$cmd`;
}
#print STDERR $cmd, "\n\n"; `$cmd`;

####convert the genbank file to fasta file
my $cmd = "$python $Bin/helper/gb2fas.py $file_gb $file_fas";
print STDERR $cmd, "\n\n"; `$cmd`;

####run repeat discovery
my $cmd =  "perl $Bin/repeat_disc/repeat_disc.pl $projectid $workingdir $file_fas \" 1-10 2-6 3-5 4-5 5-5 6-5 \" \" 2 7 7 80 10 50 500 -f -d -m \" \" -f -p -h 3 -l 30 \"";
print STDERR $cmd, "\n\n"; `$cmd`;

my $cmd = "perl $Bin/ttviewer/cpgview_wrapper.pl $file_gb $projectid $workingdir";
print STDERR $cmd, "\n\n"; `$cmd`;

my $file_pdf1   = $workingdir."/".$projectid."_circle.pdf";
my $file_pdf2   = $workingdir."/".$projectid."_cc01.pdf";
my $file_pdf3   = $workingdir."/".$projectid."_cc02.pdf";
my $file_image1 = $workingdir."/".$projectid."_circle.png";
#$file_image1 = "/tmp/images_0.png";
#$file_image1 = $file_image1.".jpg";

my $file_trf    = $workingdir."/".$projectid.".trf.txt";
my $file_misa   = $workingdir."/".$projectid.".fas.misa";
my $file_vmatch = $workingdir."/".$projectid."_vmatch.txt";
my $file_html   = $workingdir."/".$projectid.".htm";

##draw circo map
##prepare cricos file
&prepare_circos_files($file_trf,     $file_circos_trf,      "trf");
&prepare_circos_files($file_misa,    $file_circos_misa,     "misa");
&prepare_circos_files($file_vmatch,  $file_circos_vmatch_d, "vmatchd");
&prepare_circos_files($file_vmatch,  $file_circos_vmatch_p, "vmatchp");
#&prepare_circos_files($file_ir,      $file_circos_cc,       "cc");
#&prepare_circos_files($file_ir,      $file_circos_ir,       "ir");

#wait unless (-e $file_circos_vmatch_p);

my $cmd = "cp $Bin/templates/* $workingdir";
print STDERR $cmd, "\n\n"; `$cmd`;

##update the length of the sequence
my $line = `cat $file_gb | grep LOCUS -i`;

$cmd = "perl $Bin/helper/updateSeqLen.pl $file_gb";
print STDERR $cmd, "\n\n"; 
my $str = `$cmd`;

print STDERR "===================$str====================\n";
$cmd = "echo \"$str\" > $workingdir/karyotype.txt";
print STDERR $cmd, "\n\n"; `$cmd`;

$cmd = "cat $workingdir/karyotype.txt | sed 's/ /\t/g' | cut -f3,5,6,7 | sed 's/vblack/fill_color=black/g' > $file_circos_cc";
print STDERR $cmd, "\n\n"; `$cmd`;

$cmd = "cat $workingdir/circos.conf | sed 's/projectid/$projectid/g' > $workingdir/$projectid.conf";
print STDERR $cmd, "\n\n"; `$cmd`;

my $count = 0;
while (! -e $file_circos_trf || ! -e $file_circos_misa || ! -e $file_circos_vmatch_d || ! -e $file_circos_vmatch_p ) {
        sleep(3);
    $count += 3;
}

###start to run circos
print STDERR "\n\nwaited for $count second\n\n";
my $circos = "$Bin/circos-0.69-6/bin/circos";
#my $circos = "/share/apps/circos-0.69-6/bin/circos";
$cmd = "(cd $workingdir && $circos -conf $projectid.conf && cp circos.png $file_image2)";
#my $cmd = "(cd $workingdir && $circos -conf $projectid.conf && cp circos.png $file_image2)";
print STDERR $cmd, "\n\n"; `$cmd`;

#my $python = "/biodata4/home/cliu/drawgenemap/bin/python";
$cmd = "$python $Bin/helper/drawgenemap_s3.py $file_image1 $file_image2 $file_image3 $file_image4";
print STDERR $cmd, "\n\n"; `$cmd`;

#my $python = "/biodata4/home/cliu/cpgview/bin/python";
$cmd = "$python $Bin/helper/image_postp.py $file_pdf1 $file_image3 $file_pdf2 $file_pdf3";
print STDERR $cmd, "\n\n"; `$cmd`;

#$cmd = "$python $Bin/helper/resize_image.py $file_image4 400 $file_image5";
#print STDERR $cmd, "\n\n"; `$cmd`;

&generateHTML($projectid, $file_html);

print STDERR "\n\n--plasdrawmap completed!\n\n";

##################################################################
sub prepare_circos_files {

    my $file_in  = shift;
    my $file_out = shift;
    my $type     = shift;

    my @lines    = `cat $file_in`;
    my $results  = "";
    open (OUT, ">$file_out");

    if ($type eq "trf") {
        foreach my $line (@lines) {
            my @ws = split(/\s+/, $line);
            if ($ws[0] =~ /^\d/) {
                $results .= "cpg\t$ws[0]\t$ws[1]\tfill_color=vvdblue\n";
            }
        }
    } elsif ($type eq "misa") {
        foreach my $line (@lines) {
            next if ($line =~ /SSR/);
            my @ws = split(/\s+/, $line);
            if ($ws[2] =~ /p1/) {
                #sami_cpg1       12234   12301   fill_color=green
                #NC_000932.1     26      p1      (T)11   11      45196   45206
                $results .= "cpg\t$ws[5]\t$ws[6]\tfill_color=vvdgreen\n";
            } elsif ($ws[2]   =~ /p2/) {
                $results .= "cpg\t$ws[5]\t$ws[6]\tfill_color=vvdyellow\n";
            } elsif ($ws[2]   =~ /p3/) {
                $results .= "cpg\t$ws[5]\t$ws[6]\tfill_color=vvdpuple\n";
            } elsif ($ws[2]   =~ /p4/) {
                $results .= "cpg\t$ws[5]\t$ws[6]\tfill_color=vvdblue\n";
            } elsif ($ws[2]   =~ /p5/) {
                $results .= "cpg\t$ws[5]\t$ws[6]\tfill_color=vvdorange\n";
            } elsif ($ws[2]   =~ /c/) {
                $results .= "cpg\t$ws[5]\t$ws[6]\tfill_color=vvdblack\n";
            } else {
                $results .= "cpg\t$ws[5]\t$ws[6]\tfill_color=black\n";
            }
        }

    } elsif ($type eq "vmatchd") {
        foreach my $line (@lines) {
            $line =~ s/^\s+//;
            my @ws = split(/\s+/, $line);
#   33 130831   D    33 130863  -2    4.32e-07

            if ($ws[2] =~ /D/) {
                my $end1 = $ws[1] + $ws[0];
                my $end2 = $ws[4] + $ws[3];
                $results .= "cpg\t$ws[1]\t$end1\tcpg\t$ws[4]\t$end2\n";
            }
        }

    } elsif ($type eq "vmatchp") {
        foreach my $line (@lines) {
            $line =~ s/^\s+//;
            my @ws = split(/\s+/, $line);
#   33 130831   D    33 130863  -2    4.32e-07

            if ($ws[2] =~ /P/) {
                my $end1 = $ws[1] + $ws[0];
                my $end2 = $ws[4] + $ws[3];
                $results .= "cpg\t$ws[1]\t$end1\tcpg\t$ws[4]\t$end2\n";
            }
        }

    } elsif ($type eq "ir") {
#15478017817688  maker   repeat_region   86568   114329  0.00e+00        .       .       ID=repeat_regiona;Annotation=repeat_regiona;
        foreach my $line (@lines) {
            $line =~ s/^\s+//;
            if ($line =~ /repeat_region/) {    
                my @ws = split(/\s+/, $line);
                $results .= "cpg\t$ws[3]\t$ws[4]\tfill_color=black\n";
            }
        }

    } else {}

    print OUT $results;
    close(OUT);
}

#############################
sub generateHTML {

    my $id       = shift;
    my $file_html = shift;
    
    #my $outfile1=$id."_circle.pdf";
    my $outfile1=$id."_cc02.pdf";
    my $outfile2=$id."_cis.pdf";
    my $outfile3=$id."_trans.pdf";
    my $outfile4=$id."_circle.png";
    my $outfile5=$id."_image2.png";
    my $outfile6=$id."_image3.png";
    my $outfile7=$id."_image4.png";
    my $outfile8=$id."_image5.png";
    
    my $report =  <<"REPORT";
    <html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
    <title>Thanks!</title>
    <style type="text/css"> img {border: none;} </style>
    </head>
    <body>
    <p>Your uploaded <a href=$id.gbf target=_blank>file</a></p>
    <p>Your <a href=$outfile1 target=_blank>cpg circular graph</a></p>
    <p>Your <a href=$outfile2 target=_blank>cis-splicing genes</a></p>
    <p>Your <a href=$outfile3 target=_blank>trans-splicing genes</a></p>
    <!-- <p>Your <a href=$outfile4 target=_blank>image01</a></p>
    <p>Your <a href=$outfile5 target=_blank>image02</a></p>
    <p>Your <a href=$outfile6 target=_blank>image03</a></p>
    <p>Your <a href=$outfile7 target=_blank>image04</a></p>
    <p>Your <a href=$outfile8 target=_blank>image05</a></p>
    -->
    </body>
    </html>
REPORT
    
    open (OUT, ">$file_html");
    print OUT $report; 
    close(OUT);
}
