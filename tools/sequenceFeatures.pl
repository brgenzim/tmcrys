#!/usr/bin/perl
use strict;
use warnings;

####################
#
#    TMCrys version 0.1
#    Copyright 2017 Julia Varga and Gábor E. Tusnády
#
#    If you use TMCrys, please cite: 
#    Julia Varga and Gábor E. Tusnády TMCRys...
#
#
#    This file is part of TMCrys.
#
#    TMCrys is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    TMCrys is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with TMCrys.  If not, see <http://www.gnu.org/licenses/>.
#
####################



my $mydir = $ENV{'TMCRYS'};

use lib "$ENV{'TMCRYS'}/$mydir/tools/";
if ( !-f "$ENV{'TMCRYS'}/tools/OB.pm" || !-f "$ENV{'TMCRYS'}/data/zmat.dat" || !-f "$ENV{'TMCRYS'}/data/Hydrophobicity_scores.dat"){
	die "Please put modified OB.pm, zmat.dat and Hydrophobicity_scores.dat to the working directory\n";
}

use XML::LibXML;
use Bio::Tools::Protparam;;
use Getopt::Std;
use OB;


######### GLOBAL #########
our($opt_i, $opt_s, $opt_n, $opt_h, $opt_f, $opt_d);
getopts('i:s:n:hfd:');

my $help = 
"\nUsage: ./sequenceFeatures.pl (-s TABFILE | -i CCTOPFILE) -n DIR [-h] [-f] [-d DIR]
	  
	  -i CCTOP file with single entry
	  -d directory with CCTOP files
	  -s file with 'id seq top'
	  Either -s or -i or -d required.
	
	  -n netsurf results file (.rsa)
	  -f print header with feature names
	  -h print this help
		
TMCrys version 0.1
If you use it, please cite Julia Varga and Gábor E. Tusnády TMCrys... \n
";

if ($opt_h){
	print "$help";
	exit;
}
if (!($opt_i || $opt_s || $opt_d)){
	die "\nsequenceFeatures.pl: Please specifiy a input file with -s or -i or -d option!\n$help";
}
if (!defined($opt_n)){
	die "\nsequenceFeatures.pl: Please specifiy a netsurfp results file directory with -n option!\n$help";
}
if (!-f $opt_n){
	die "\nsequenceFeatures.pl: netsurfp directory does not exists. Please specifiy a directory with a correct path!\n$help";
}

#Define features to count and hash for storing
my %features;
my @features = (
"numTM",
"lgMolWeight",
"length", 'lengthTM', 'lengthNonTM', 'avgTM',
'longestTM', 'longestNonTM',
'inindex', 'half_life', 'gravy', 'OB', 'pI', 'stability',
'pos', 'neg', 'noCharge', 'charge', 'aromatic', 'alifatic', 'sulfur', 'hydroxil', 'nonpolar', 'polar',
'fractionTM', "TMratio",
'TM-pos', 'TM-neg', 'TM-noCharge', 'TM-charge', 'TM-aromatic', 'TM-alifatic', 'TM-sulfur', 'TM-hydroxil', 'TM-polar', 'TM-nonpolar',
'NonTM-pos', 'NonTM-neg', 'NonTM-noCharge', 'NonTM-charge', 'NonTM-aromatic', 'NonTM-alifatic', 'NonTM-sulfur', 'NonTM-hydroxil', 'NonTM-polar', 'NonTM-nonpolar',
'glyco', 'signal',
'GxxG', 'GxxxG',
'avgRSA', 'exposed', 'buriedratio'


);
push @features, qw( 
BHAR880101.lag1 BHAR880101.lag1.1  BHAR880101.lag10 BHAR880101.lag10.1 BHAR880101.lag15 BHAR880101.lag15.1 BHAR880101.lag2 BHAR880101.lag2.1  BHAR880101.lag20 BHAR880101.lag20.1 BHAR880101.lag25 BHAR880101.lag25.1 BHAR880101.lag30 BHAR880101.lag30.1 BHAR880101.lag5 BHAR880101.lag5.1  BIGC670101.lag1 BIGC670101.lag1.1  BIGC670101.lag10 BIGC670101.lag10.1 BIGC670101.lag15 BIGC670101.lag15.1 BIGC670101.lag2 BIGC670101.lag2.1  BIGC670101.lag20 BIGC670101.lag20.1 BIGC670101.lag25 BIGC670101.lag25.1 BIGC670101.lag30 BIGC670101.lag30.1 BIGC670101.lag5 BIGC670101.lag5.1  CHAM810101.lag1 CHAM810101.lag1.1  CHAM810101.lag10 CHAM810101.lag10.1 CHAM810101.lag15 CHAM810101.lag15.1 CHAM810101.lag2 CHAM810101.lag2.1  CHAM810101.lag20 CHAM810101.lag20.1 CHAM810101.lag25 CHAM810101.lag25.1 CHAM810101.lag3 CHAM810101.lag3.1  CHAM810101.lag30 CHAM810101.lag30.1 CHAM810101.lag5 CHAM810101.lag5.1  CHAM820101.lag1 CHAM820101.lag1.1  CHAM820101.lag10 CHAM820101.lag10.1 CHAM820101.lag15 CHAM820101.lag15.1 CHAM820101.lag2 CHAM820101.lag2.1  CHAM820101.lag20 CHAM820101.lag20.1 CHAM820101.lag25 CHAM820101.lag25.1 CHAM820101.lag3 CHAM820101.lag3.1  CHAM820101.lag30 CHAM820101.lag30.1 CHAM820101.lag5 CHAM820101.lag5.1  CHAM820102.lag1 CHAM820102.lag1.1  CHAM820102.lag10 CHAM820102.lag10.1 CHAM820102.lag15 CHAM820102.lag15.1 CHAM820102.lag2 CHAM820102.lag2.1  CHAM820102.lag20 CHAM820102.lag20.1 CHAM820102.lag25 CHAM820102.lag25.1 CHAM820102.lag3 CHAM820102.lag3.1  CHAM820102.lag30 CHAM820102.lag30.1 CHAM820102.lag5 CHAM820102.lag5.1  CHOC760101.lag1 CHOC760101.lag1.1  CHOC760101.lag10 CHOC760101.lag10.1  CHOC760101.lag15 CHOC760101.lag15.1 CHOC760101.lag2 CHOC760101.lag2.1  CHOC760101.lag20 CHOC760101.lag20.1 CHOC760101.lag25 CHOC760101.lag25.1 CHOC760101.lag3 CHOC760101.lag3.1  CHOC760101.lag30 CHOC760101.lag30.1 CHOC760101.lag5 CHOC760101.lag5.1    CIDH920105.lag1 CIDH920105.lag1.1  CIDH920105.lag10 CIDH920105.lag10.1 CIDH920105.lag15 CIDH920105.lag15.1 CIDH920105.lag2 CIDH920105.lag2.1  CIDH920105.lag20 CIDH920105.lag20.1 CIDH920105.lag25 CIDH920105.lag25.1 CIDH920105.lag3 CIDH920105.lag3.1  CIDH920105.lag30 CIDH920105.lag30.1 CIDH920105.lag5 CIDH920105.lag5.1  DAYM780201.lag1 DAYM780201.lag1.1  DAYM780201.lag10 DAYM780201.lag10.1 DAYM780201.lag15 DAYM780201.lag15.1 DAYM780201.lag2 DAYM780201.lag2.1  DAYM780201.lag20 DAYM780201.lag20.1 DAYM780201.lag25 DAYM780201.lag25.1 DAYM780201.lag3 DAYM780201.lag3.1  DAYM780201.lag30 DAYM780201.lag30.1 DAYM780201.lag5 DAYM780201.lag5.1  Xc1.A Xc1.C Xc1.D Xc1.E Xc1.F Xc1.G Xc1.H Xc1.I Xc1.K Xc1.L Xc1.M Xc1.N Xc1.P Xc1.Q Xc1.R Xc1.S Xc1.T Xc1.V Xc1.W Xc1.Y Xc2.lambda.1 Xc2.lambda.10 Xc2.lambda.15 Xc2.lambda.2 Xc2.lambda.20 Xc2.lambda.25 Xc2.lambda.3 Xc2.lambda.30  prop1.Tr1221	prop1.Tr1331	prop1.Tr2332	prop2.Tr1221	prop2.Tr1331	prop2.Tr2332	prop3.Tr1221	prop3.Tr1331	prop3.Tr2332	prop4.Tr1221	prop4.Tr1331	prop4.Tr2332	prop5.Tr1221	prop5.Tr1331	prop5.Tr2332	prop6.Tr1221	prop6.Tr1331	prop6.Tr2332	prop7.Tr1221	prop7.Tr1331	prop7.Tr2332 CHAM820101.lag1 CHAM820101.lag1.1);

my $aminoacids = "ARNDCEQGHILKMFPSTWYV";

foreach my $aa (split //, $aminoacids){
	push @features, $aa;
	push @features, "TM-$aa";
	push @features, "NonTM-$aa";
}

#Define NULL features
my %loc = ('S' => 'NonTM', 'O' => 'NonTM', 'I' => 'NonTM', 'M' => 'TM', 'L' => 'TM', 1 => 'NonTM', 2 => 'NonTM', 'H' => 'TM');


my @groups = ('pos', 'neg', 'noCharge', 'charge', 'aromatic', 'alifatic', 'sulfur', 'hydroxil', 'polar', 'nonpolar');
my %aagroups = (
"A" => ['noCharge', 'alifatic', 'nonpolar'], 
"R" => ['pos', 'charge', 'polar'], 
"N" => ['noCharge', 'polar'], 
"D" => ['neg', 'charge', 'polar'], 
"C" => ['noCharge', 'sulfur', 'nonpolar'], 
"E" => ['neg', 'charge', 'polar'], 
"Q" => ['noCharge', 'polar'], 
"G" => ['alifatic', 'noCharge', 'nonpolar'], 
"H" => ['pos', 'charge', 'polar'], 
"I" => ['alifatic', 'noCharge', 'nonpolar'], 
"L" => ['alifatic', 'noCharge', 'nonpolar'], 
"K" => ['pos', 'charge', 'polar'], 
"M" => ['noCharge', 'sulfur', 'nonpolar'], 
"F" => ['aromatic', 'nonpolar', 'noCharge'], 
"P" => ['noCharge', 'nonpolar'], 
"S" => ['noCharge', 'hydroxil', 'polar'], 
"T" => ['noCharge', 'hydroxil', 'polar'], 
"W" => ['aromatic', 'nonpolar', 'noCharge'], 
"Y" => ['aromatic', 'polar', 'noCharge'], 
"V" => ['alifatic', 'noCharge', 'nonpolar'], 
);

####### END GLOBAL #######


my ($seq, $top, $id) = ("", "", "");
# Get id of sequence
if ( defined($opt_s )){
	header();
	
	open FILE, "<", $opt_s or die "Could not open $opt_s: $! \n";
	while (my $line = <FILE>){
		chomp $line ;
		
		($id, $top, $seq) = split " ", $line ;
		
		if(!defined($seq)){
			print STDERR "$id\tNot transmembrane protein. Line $. left out.\n";
			next;
		}
		
		if (calculate() == 1){
			printResults();
		}
		else{
			next;
		}
	}
}
elsif( defined($opt_i) ) {
	null_features();
	my ($file) = $opt_i;
	$. = 1;	
	($id, $top, $seq) = cctop ($file);
	if ($id eq "nontmp"){
		die "$id\tNot transmembrane protein. File left out.\n";
	}
	
	if (calculate() == 1){
		header();
		printResults();
	}
}
elsif( defined($opt_d) ){
	my @dir = <$opt_d/*>;
	header();
	foreach my $file (@dir){
		null_features();
		$. = 1;	
		($id, $top, $seq) = cctop ($file);
		
		if ($top eq "nontmp"){
			print STDERR "$id\tNot transmembrane protein. XML left out.\n";
			next;
		}
		
		if (calculate() == 1){
			printResults();
		}
	}
}
else{
	die "Please specifiy a input file with -s or -i option!\n"
}
#-----------------------------------------------------

########## SUBS ##########
sub header {
	@features = do { my %seen; grep { !$seen{$_}++ } @features };
	print "ID\t", join("\t", sort @features), "\n" if ($opt_f);
}

sub calculate{
	null_features ();
	if (length $top != length $seq){
		print STDERR $. ," \t $id: Topology length and sequence length are not equal in size! Line or file is left out.\n$top\n$seq\n" ;
		return 0;
	}
	if ($top =~ m/T/g){
		print STDERR $. ," \t $id: TMCrys cannot handle transit peptides. Line or file is left out.\n" ;
		return 0;
	}
	if ($seq =~ m/X/g){
		print STDERR $. ," \t $id: TMCrys cannot handle unknown amino acids. Line or file is left out.\n" ;
		return 0;
	}
	if (length_regions() == 0){
		print STDERR $. ," \t $id: There was some problems with processing CCTOP or determining region lengths. Line or file is left out.\n" ;
		return 0;
	}
	if (netsurfp() == 0){
		print STDERR $. ," \t $id: There was some problems with netsurfp results. Line or file is left out.\n" ;
		return 0;
	}
	if (aminoacids() == 0){
		print STDERR $. ," \t $id: There was some problems with calculation of amino acid compositions. Line or file is left out.\n" ;
		return 0;
	}
	if (protparam() == 0){
		print STDERR $. ," \t $id: There was some problems with using protparam module. Line or file is left out.\n" ;
		return 0;
	}
	sleep 0.000000000000000001;
	if (obScore() == 0){
		print STDERR $. ," \t $id: There was some problems calculating obScore. Is there the modified OB.pm, zmat.dat and Hydrophobicity_scores.dat in the same directory as working script? Line or file is left out.\n" ;
		return 0;
	}
	if ($features{'OB'} eq "NA"){
		print STDERR $. , "\t $id: OB score could not be determined. Please fill by hand from http://www.compbio.dundee.ac.uk/xtal/cgi-bin/input.pl website.\n";
		return 0;
	}
	if (fromR($seq) == 0){
		print STDERR $. ," \t $id: There was some problems with using protparam module. Line or file is left out.\n" ;
		return 0;
	}
	
	return 1;
}

sub printResults{
	print "$id";
	foreach my $key (sort keys %features){
		print "\t$features{$key}";
	}
	print "\n";

	return 1;
}

sub print_column_headers{
	print "columns\t";
	print join("\t", sort keys %features);
	print "\n";
}

sub openXML{
	my $file = $_[0];
	my $parser = XML::LibXML->new();
	my $tree = $parser->parse_file($file) or die "Could not parse xml: $!";
	my $root = $tree->getDocumentElement;
	
	return $root;
}

sub obScore{
	my $vars = OB -> new;
	my $opt;
	$opt->{'i'} = "in";
	$opt->{'o'} = "out";
	$vars -> ob_vars($opt);
	my $pi = OB -> new;
	my $wt = OB -> new;
	$pi->{$id} = $features{'pI'};
	my $gravy = OB -> new;
	my $len = OB -> new;
	$gravy->calc_gravy_len($len, $seq, $id, $vars);
	my $ob = OB -> new;
	my ($scor_href, $xl, $xh, $yl, $yh, $xv_aref, $yv_aref) = $vars->read_zmat ();
	my ($irregX, $xint) = $xv_aref->assess_interval ();
	my ($irregY, $yint) = $yv_aref->assess_interval ();
	die "need to implement alternative lookup method\n" if (($$irregX) or ($$irregY)); # current lookup only works for table with regular values of X and Y
	my ($xu_re, $yu_re, $xl_re, $yl_re) = $pi->assign_matrixBound($gravy, $xint, $yint);
	my $ob_score = $scor_href->OB_out_nomwt($xu_re, $yu_re, $xl_re, $yl_re, $vars, $len, $pi, $gravy);
	
	$features{'OB'} = defined($ob_score)? sprintf "%.4f", $ob_score : "NA";
	return 1;	
}

#read CCTOP file to topology string
sub cctop{
	my ($file) = @_;
	my $root = openXML ("$file");
	my $id = $root -> findvalue('@id');
	
	my $tm = $root -> findvalue("\@transmembrane");
	my ($top, $seq) = ("", "");
	
	if ($tm eq "yes"){
		my $TM = $root -> findvalue("//\@numTM");
		$features{'numTM'} = $TM;
		my $length = $root -> findvalue('Sequence/@length');
		my @regions = $root -> findnodes("Topology/Region");
		$top = "-" x $length;

		foreach my $region (@regions){
			my $from = $region -> findvalue('@from');
			my $to = $region -> findvalue('@to');
			my $loc = $region -> findvalue('@loc');
			substr($top, $from - 1, $to - $from + 1) = $loc x ($to - $from + 1);
		}
		$seq = $root -> findvalue('Sequence/Seq');
		$seq =~ s/\s//g;
	}
	else {
		return ($id, "nontmp");
	}
	
	return ($id, $top, $seq);
}

sub longest{
	return (reverse sort { $a <=> $b } map { length($_) } @_)[0];
}

sub sumLong{
	return length join "", @_;
}

sub log10 {
 my $n = shift;
 return sprintf "%.4f", log($n)/log(10);
}


sub length_regions{
	$features{'length'} = length $seq;
	
	my @TMs = ($top =~ m/(L+|M+)/g);
	
	my @nonTMs;
	if ($top =~ m/I/){
		@nonTMs = ($top =~ m/(I+|O+)/g);
	}
	elsif ($top =~ m/1/){
		@nonTMs = ($top =~ m/(1+|2+)/g);
	}
	
	$features{'numTM'}  = scalar @TMs;
	
	$features{'lengthTM'} = sumLong (@TMs);
	$features{'lengthNonTM'} = sumLong (@nonTMs);
	$features{'TMratio'} = $features{'lengthTM'}/$features{'lengthNonTM'};
	
	$features{'longestTM'} = longest (@TMs);
	$features{'longestNonTM'} = longest (@nonTMs);

	$features{'fractionTM'} = $features{'lengthTM'}/$features{'length'};
	
	$features{'signal'} = $top =~ m/S/g ? 1 : 0;
	$features{'avgTM'} = $features{'lengthTM'}/$features{'numTM'};

	return 1;
}

sub fromR{
	(my $seq) = @_;
	my @Rfeatures = `Rscript $mydir/sequenceFeatures.R $seq`;
	shift @Rfeatures;
	
	foreach my $feature (@Rfeatures){
		chomp $feature;
		my ($key, $value) = split "\t", $feature;
		$features{$key} = $value if exists $features{$key};
	}
	return 1;
}

sub aminoacids{
	my %aas;
	my %dipep;
	for (my $i = 0; $i< $features{'length'}; $i++ ){
		my $a = substr ($seq, $i, 1);
		my $t = $loc{substr ($top, $i, 1)} or die "$top\t$i\n";
		
		#Count in sequence and regions
		$aas{$a}++;
		$aas{"$t-$a"}++;
		
		#Aminoacid groups
		foreach my $g (@{$aagroups{$a}}){
			$aas{$g}++;
			$aas{"$t-$g"}++;
		}
	}

	foreach my $key ('pos', 'neg', 'noCharge', 'charge', 'aromatic', 'alifatic', 'sulfur', 'hydroxil'){
		$features{"TM-${key}"} = $aas{"TM-${key}"}/$features{'lengthTM'} if exists $aas{"TM-${key}"};
		$features{"NonTM-${key}"} = $aas{"NonTM-${key}"}/$features{'lengthNonTM'} if exists $aas{"NonTM-${key}"};
	}
	
        foreach my $key (keys %aas){
		if ($key =~ m/-/g){
		        my ($side, $group) = split '-', $key;
		        $features{$key} = $aas{$key}/$features{"length$side"};
		}
		else{
		        $features{$key} = $aas{$key}/$features{'length'};
		}
        }

	return 1;
}

sub protparam{
	my $pp = Bio::Tools::Protparam -> new(seq => $seq);
	$features{'inindex'} = $pp->instability_index();
	$features{'stability'} = $pp->stability() eq 'stable' ? 1 : 0;
	$features{'half_life'} = $pp -> half_life;
	$features{'lgMolWeight'} = log10($pp -> molecular_weight);
	$features{'pI'} = $pp->theoretical_pI;
	$features{'gravy'} = $pp->gravy;

 	$seq =~ m/(N[^P][ST][^P])/g;
 	foreach my $pos (@-){
	 	$features{'glyco'}++ if substr($top, $pos, 1) ne "M" and substr($top, $pos, 1) ne "I";
 	}
 	
 	$seq =~ m/(B)/g;
 	$seq =~ m/(G..G)/g;
 	foreach my $pos (@-){
	 	$features{'GxxG'}++ if substr($top, $pos, 4) =~ m/.*M.*/g;
 	}
 	
 	$seq =~ m/(B)/g;
 	$seq =~ m/(G...G)/g;
 	foreach my $pos (@-){
	 	$features{'GxxxG'}++ if substr($top, $pos, 5) =~ m/.*M.*/g;
 	}

	return 1;
}

sub null_features{
	%features = ();

	foreach my $f (@features){
		$features{$f} = 0;
	}
	return 1;
}

sub netsurfp{
	chomp(my @result = `grep $id $opt_n -h| sed 's/ {2,}//g'`);
	
	if(scalar (@result) == 0){
		return 0;
	}
	else{
		@result = do { my %seen; grep { !$seen{$_}++ } @result };
		if ($features{'length'} != scalar @result){
			print STDERR ("Something wrong with netsurfp length! $id", scalar @result, " ", $features{'length'}, "\n") ;
			return 0;
		}
		my $bur = "";
		my $sum = 0;

		foreach my $res (@result){
			my @array = split ' ', $res;
			$bur .= $array[0];
			$sum += $array[4];
		}
		my $exposed = () = $bur =~ m/E/g;

		$features{'avgRSA'} = $sum / $features{'length'};
		$features{'exposed'} = $exposed / $features{'length'};
		$features{'buriedratio'} = (length ($bur) - $exposed)/$exposed;

		return 1;
	}
}

######## END SUBS ########
