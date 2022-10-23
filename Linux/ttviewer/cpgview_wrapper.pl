#!/usr/bin/perl
#
use FindBin qw($Bin);

my $inputfile = shift;
my $projectid = shift;
my $outdir    = shift;

print STDERR "$Bin $inputfile $projectid $outdir\n";

if (! -e $inputfile || $outdir eq "") {
	print STDERR "\ninputfile does not exist\n\nUsage: $0 inputfile projectid outdir\n\n";
	exit;
}

if ($projectid eq "") {
	$projectid = time();
	$projectid =~ s/\s+//g;
}

my $inputfile1 = $projectid.".gbf";
my $cmd        = "mkdir $outdir";
`$cmd` unless (-e $outdir);
`cp $inputfile $outdir/$inputfile1`;
my $logfile    = $outdir."/".$projectid.".log";

my $python     = "/apps/miniconda3/envs/cpgview/bin/python";
#Location of Python，The user needs to change the path to the location of Python in the user's local computer.
my $Rscript    = "/apps/miniconda3/envs/r4/bin/Rscript";
#Location of R，The user needs to change the path to the location of R in the user's local computer.


my $data_prep  = $Bin."/data_process.py";

my $cmd        = "(cd $outdir; $python $data_prep $inputfile1 dc dt dr 2>&1 >> $logfile)";
print STDERR "$cmd\n"; `$cmd`;

my $of_c_pdf   = $projectid."_cis.pdf";
my $cmd        = "(cd $Bin; $Rscript viewCSGene.R $outdir/dc/cis_splicing_gene.csv $outdir/dc/cis_splicing_subgene.csv $outdir/$of_c_pdf)";
#my $cmd        = "(cd $Bin; $Rscript viewCSGene.R $outdir/dc/cis_splicing_gene.csv $outdir/dc/cis_splicing_subgene.csv $outdir/$of_c_pdf 2>&1 >> $logfile)";
print STDERR "$cmd\n"; `$cmd`;

my $of_t_pdf   = $projectid."_trans.pdf";
my $cmd        = "(cd $Bin; $Rscript viewTSGene.R $outdir/dt/trans_splicing_gene.csv $outdir/dt/trans_splicing_subgene.csv $outdir/$of_t_pdf)";
#my $cmd        = "(cd $Bin; $Rscript viewTSGene.R $outdir/dt/trans_splicing_gene.csv $outdir/dt/trans_splicing_subgene.csv $outdir/$of_t_pdf 2>&1 >> $logfile)";
print STDERR "$cmd\n"; `$cmd`;

###convert pdf to png for the _c file
my $of_c_png   = $projectid."_cis.png";
my $cmd        = "$python $Bin/pdf2img.py $outdir/$of_c_pdf $outdir";
#my $cmd        = "$python $Bin/pdf2img.py $outdir/$of_c_pdf $outdir 2>&1 >> $logfile";
print STDERR "$cmd\n"; `$cmd`;

my $cmd        = "mv $outdir/images_0.png $outdir/$of_c_png";
print STDERR "$cmd\n"; `$cmd`;

###convert pdf to png for the _t file
my $of_t_png   = $projectid."_trans.png";
my $cmd        = "$python $Bin/pdf2img.py $outdir/$of_t_pdf $outdir";
#my $cmd        = "$python $Bin/pdf2img.py $outdir/$of_t_pdf $outdir 2>&1 >> $logfile";
print STDERR "$cmd\n"; `$cmd`;

my $cmd        = "mv $outdir/images_0.png $outdir/$of_t_png";
print STDERR "$cmd\n"; `$cmd`;


my $of_cc_pdf  = $projectid."_circle0";
my $cmd        = "(cd $Bin/Chloroplot/; $Rscript chloroplot_Genes.R $outdir/$inputfile1 $outdir/$of_cc_pdf)";
#my $cmd        = "(cd $Bin/Chloroplot/; $Rscript chloroplot_Genes.R $outdir/$inputfile1 $outdir/$of_cc_pdf 2>&1 >> $logfile)";
print STDERR "$cmd\n"; `$cmd`;

my $addword    = $Bin."/addtext.py";
my $cmd        = "(cd $outdir; $python $addword $inputfile1 $projectid $outdir)";
print STDERR "$cmd\n"; `$cmd`;

my $out_circle = $outdir."/".$projectid."_circle0.pdf";
#my $cmd        = unlink("$out_circle");
my $errordata  = "<h4>Error!</h4><p>The system cannot generate the cpg circular graph. This may be due to the following reasons.</p><p>1. The file uploaded by the user is not a valid GenBank file;<br>2. There is no gene information in the GenBank file;<br>3. The system is temporarily out of memory. A simple resubmission can solve the problem.</p>";
my $fileExist  = -e $out_circle;
if ( $fileExist ) {print "Yes"} else {open(my $fh, ">","$outdir/error2.txt") or die "Can't open > $outdir/error2.txt: $!"; print $fh $errordata; close($fh) or  "Couldn't close the file: $!"}

my $cmd        = unlink("$out_circle");
my $of_cc_png  = $projectid."_circle.png";
my $fin_circle = $projectid."_circle.pdf";

my $cmd        = "$python $Bin/pdf2img.py $outdir/$fin_circle $outdir";
#my $cmd = "$python $Bin/pdf2img.py $outdir/$fin_circle $outdir 2>&1 >> $logfile";
print STDERR "$cmd\n"; `$cmd`;

my $cmd        = "mv $outdir/images_0.png $outdir/$of_cc_png";
print STDERR "$cmd\n"; `$cmd`;


