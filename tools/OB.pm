# OB::
#.
# Author: Ian Overton, Jan 2006
# Copyright (2006) Ian M. Overton, Geoffrey J. Barton and the University of Dundee

package OB;
use Bio::Tools::pICalculator;
use Bio::SeqIO;
use Bio::Tools::SeqStats;

BEGIN {
	unshift @INC, "/usr/local/share/perl5/Bio/";  # EDIT THIS TO SPECIFY YOUR PATH TO Bioperl!
}

######
#
# new
#
######
#
# constructor method
#

sub new {

	my $inv = shift;
	my $class = ref($inv) || $inv;
	my $self = { @_ };
	bless ($self, $class);
	return $self;

}



########################
# read_fasta_so
########################
#
# reads in fasta file
# assigns ">(.+)" as key
# sequence as value
#
############################


sub read_fasta_so {
	my ($seq, $var) = @_;
	my $var_class = ref($var);
	if ($var_class) {
		$var_class .= "::input";
	} else {
		$var_class = "input";
	
	}
	my $in = $var->{$var_class};
	my $i;
	open (IN, $in) or die "cantddd open $in $!";
	while (<IN>) {
		chomp();
		if (/^>(.+)/){
			$i = $1;
		} else {
			my @l = split(//, $_);
			$seq->{len}{$i} += scalar(@l);
			$seq->{$i} .= $_;
		}
		
	}


	close IN;
	return ($seq);	
	
}


########################
# ob_vars
########################
#
# sets up variables for 
# the other methods 
#
############################

sub ob_vars {
	my  ($self, $opt) = @_;	
	if ($opt->{"h"}) { # print help message
		ob_help ();
		exit;
	}
	
##set variables
	$self->{input} = $opt->{"i"};
	$self->{outfile} = $opt->{"o"};
	$self->{zmat} = $opt->{"m"} ? $opt->{"m"} : "$ENV{'TMCRYS'}/data/zmat.dat";
	$self->{gravy_in} = $opt->{"p"} ? $opt->{"p"} : "$ENV{'TMCRYS'}/data/Hydrophobicity_scores.dat";
	$self->{nomw} = 1 if ($opt->{"n"});
	
	unless (($self->{input}) and ($self->{outfile})) {
		ob_help ();
		exit;
	} 

}

########################
# ob_help
########################
#
# info for users 
#
############################

sub ob_help {
my $h =  "Version 1.0
This software calculates the OB-Score for each sequence in the input fasta file.

OB.pl Output format(tab-delimited):
id	OB-Score	GRAVY	pI	mol_wt(lower_bound)	sequence_length

Abbreviations:
id = sequence identifier
GRAVY = GRand AVerage of hYdrophobicity (defaults to Kyte-Doolittle)
pI = isoelectric point

Command-line input:
 -i [Inputfile (fasta format)]
 -o [output file]
 -p [hydropathy score table] (can use Hydrophobicity_scores.dat)
 -m [score matrix] (use zmat.dat)
 -n no mol_wt calculation (in case of nonstandard amino acids in input)
 	output format: id	OB-Score	GRAVY	pI	sequence_length
 -h prints this message
 
If using this software, data, (or any derivative works) please cite:
Overton, I. M. & Barton, G. J. (2006). \"A Normalised Scale for Structural
Genomics Target Ranking: The OB-Score\". FEBS Lett. 580:4005-9.
";

print STDOUT $h; 


}


##########
#
# calc_pI_mwt
#
#########
#
# calculates pI & mol_wt using bioperl

sub calc_pI_mwt {
	my ($pI, $wt, $seq, $id) = @_;
#	my $var_class = ref($var);
#	
#	my $inp = $var->{input};
#	open (IN, $inp) or die "cant openaaaa $inp $!";

#	my $in = Bio::SeqIO->new( -fh => \*IN ,-format => 'Fasta' ); # indicate the filehandle to read from
#	my $calc = Bio::Tools::pICalculator->new(-places => 2, -pKset => 'EMBOSS'); # set parameters for pI calculation

#         while ( my $seq = $in->next_seq ) {
         	$calc->seq($seq);
#	 	my $id = $seq->id;
	 	$pI->{$id} = $calc->iep;
	 	$wt->{$id} = Bio::Tools::SeqStats->get_mol_wt($seq);
  
#         }
	 
	 close IN;	 
}




##########
#
# calc_pI
#
#########
#
# calculates pI only using bioperl

sub calc_pI {
	my ($pI, $seq, $id) = @_;
#	my $var_class = ref($var);
#	
#	my $inp = $var->{input};
#	open (IN, $inp) or die "cant openiii $inp $!";

	my $in = Bio::SeqIO->new(-string => $seq);; # indicate the filehandle to read from
	my $calc = Bio::Tools::pICalculator->new(-places => 2, -pKset => 'EMBOSS'); # set parameters for pI calculation

         while ( my $seq = $in->next_seq ) {
         	$calc->seq($seq);
#	 	my $id = $seq->id;
	 	$pI->{$id} = $calc->iep;
	 	
         }
	 
#	 close IN;	 
}

##########
#
# calc_gravy_len
#
#########
#
# calculates gravy (grand average of hydropathy)

sub calc_gravy_len {

	my ($hadd, $len, $seq, $id, $var) = @_;
#	my $var_class = ref($var);
	my %seq;
#	if ($var_class) {
#		$var_class .= "::input";
#	} else {
#		$var_class = "input";
#	
#	}
#	my $inp = $var->{input};
	my ($pI, $wt);
	$seq{$id} = $seq;
#	open (IN, $inp) or die "cant opensss $inp $!";	 
#	
#	 
#	while (<IN>) {
#		chomp ();
#		if (/^>(\S+)/) {
#			if ($id) {
#				my @le = split (//, $seq{$id});
#				$len->{$id} = scalar (@le);
#			}
#			$id = $1;
#		} elsif (/^\S+/) {
#			tr/a-z/A-Z/;
#			tr/[^A-Z]//;
#			$seq{$id} .= $_;
#			}

#	}
	my @le = split (//, $seq{$id});
	$len->{$id} = scalar (@le);

	$var_class = ref($var);
	if ($var_class) {
		$var_class .= "::gravy_in";
	} else {
		$var_class = "gravy_in";
	
	}
	my $h_locate = $var->{gravy_in};
	unless ($h_locate) {	 
		ob_help ();
		die "no hydrophobicity score file specified\n";
	}
	
	open (I, $h_locate) or die "cant openooo $h_locate: $!\n"; 	# hydropathy table format -  (aa_name	aa_single_letter	score)
	while (<I>) {
		chomp ();
		my @a =  split (/\s+/, $_);
		my $aa = $a[1];
		$hscor{$aa} = $a[2];
		
	}
	foreach my $k (keys %seq) {
		my $h = 0;
		my @s = split (//, $seq{$k});
		my $lent = scalar (@s);
		foreach my $a (@s) {
			if ($hscor{$a}) {
				$h += $hscor{$a};
			} elsif ($a =~ /[BZXU\+\*-\[\]\s]+/) {
				 $lent--; 
			
			} else {
				print "no hydropathy value for $a\n";
			}
		}
		$hadd->{$k} = $h / $lent if ($lent); # protect against zero length sequences
	}
	
}
	



#~~~#######################################


##########
#
# calc_OB
#
#########
#
# calculates OB-Score

sub calc_OB {

}

##########
#
# read_zmat
#
#########
#
# reads in zmat.dat


sub read_zmat {	# usage: read_ztable (In_file) : reads 
	my ($table, @a, @b, @r, @xv, @yv, $seen, $xvr, $yvr, %ta);  
	my ($xl, $xh, $yl, $yh) = ("100000", "-100000","100000", "-100000"); # pI_bound_low  pI_bound_hi gravy_bound_low gravy_bound_hi
	my $vars = shift;
	my $class = ref($vars);
	$table = \%ta;
	bless ($table, $class);
		
	open (IN, $vars->{zmat}) or die "cant openxxx ". $vars->{zmat}. " $!\n";
	while (<IN>) {
		next if (/^\s*#/);	
		@a = split;
		my $co = tr/\(\]//; #count the number of new columns defined in the line
		#print "CO $co\n";
		if ($co >=  8) {
			@b = @r;
			foreach my $v (@a) {
				#$v =~ s/\(\]//g;
				#my $x = $v;
				#die "cant do tr1 !!\n" unless ($x =~ tr/\(\]// ); strangely doesnt work!
				die "cant do tr1 !!\n" unless ($v =~ s/\((.*)\]/$1/g);
				#print "A @a\n$v\n";
				push (@b, $v);
			
			}
		
		} elsif ($co == 2) {
			my $p = shift @a;
			die "cant do tr1 !!\n" unless ($p =~ s/\((.*)\]/$1/g);
			my @pi = split (/\,/, $p);
			for (my $n = 0; exists $a[$n]; $n++) {
				my @g = split (/\,/, $b[$n]);
				$table->{$pi[0]}{$pi[1]}{$g[0]}{$g[1]} = $a[$n];
				push (@xv, $pi[0]) unless (exists $xv[0]);
				push (@xv, $pi[1]) unless ($seen->{$pi[1]});
				push (@yv, $g[0]) unless (exists $yv[0]);
				push (@yv, $g[1]) unless ($seen->{$g[1]});
				$seen->{$g[1]} = 1;
				$seen->{$pi[1]} = 1;
				#print "$pi[0]	$pi[1]	$g[0]	$g[1] $a[$n]\n";
				$xl = $pi[0] if ($pi[0] < $xl); 
				$xh = $pi[0] if ($pi[0] > $xh);
				$xl = $pi[1] if ($pi[1] < $xl); 
				$xh = $pi[1] if ($pi[1] > $xh);
				$yl = $g[0] if ($g[0] < $yl); 
				$yh = $g[0] if ($g[0] > $yh);
				$yl = $g[1] if ($g[1] < $yl); 
				$yh = $g[1] if ($g[1] > $yh);
			}

		}
		
			
	}
	close IN;
	my ($rxv, $ryv) = (\@xv, \@yv);
	bless ($rxv,  $class);
	bless ($ryv, $class);
	($xvr, $yvr) = (\@xv, \@yv);
	my ($rxl, $rxh, $ryl, $ryh) = (\$xl, \$xh, \$yl, \$yh);
	bless ($rxl,  $class);
	bless ($rxh, $class);
	bless ($ryl,  $class);
	bless ($ryh,  $class);
	return ($table, $rxl, $rxh, $ryl, $ryh, $rxv, $ryv);

}

##########
#
# assess_interval
#
#########
#
# assign size of interval for x & y params

sub assess_interval {
	my $var;
	my ($aref) = shift; # reference to array containing list of bin boundaries
	my ($xint, $prev, $prev_int);
	my $i = 0;
	my $ii = 0;
	my $irreg = 0;
	
	foreach my $v (@$aref) {
		if ($i)  {
			my $modpr = sqrt ($prev * $prev);
			my $modv = sqrt ($v * $v);
			my $int = $modv - $modpr;
			$xint = sqrt ($int * $int);
			if ($ii) {
				$irreg = 1 unless ($xint eq $prev_int); 				
	
			}
			
			$prev_int = $xint;
			$ii = 1;
		}
		$i = 1;
		$prev = $v;
	}
	my $class = ref($aref);
	my ($rirreg, $rxint) = (\$irreg, \$xint);
	bless ($rirreg,  $class);
	bless ($rxint,  $class);
	return ($rirreg, $rxint);

}

	
#~~~############
#
# assign_matrixBound
#
##~~~##########
#
# locate gravy/pI combination to a matrix cell (?I think?)
#	

sub assign_matrixBound  {
	#my ($in, $xfi, $yfi, $xpz, $ypz) = @_; # requires tab delim file ($in) with id as 1st field
	my ($pI, $gr, $xfir, $yfir) = @_;
	my ($xfi, $yfi) = ($$xfir, $$yfir);
	my (%rx_upper, %ry_upper, %rx_lower, %ry_lower);
	my ($x_upper, $y_upper, $x_lower, $y_lower) = (\%rx_upper, \%ry_upper, \%rx_lower, \%ry_lower);
	my $class = ref($pI);
	bless ($x_upper, $class);
	bless ($y_upper, $class);
	bless ($x_lower, $class);
	bless ($y_lower, $class);
	
	foreach my $id (keys %$pI) {
		my $x = $pI->{$id};
		my $y = $gr->{$id};
		my $dY = $y / $yfi; 
		my $dX = $x / $xfi;
		my $intx = $1 if (($dX =~ /^([0-9]+)/) or ($dX =~ /^(-[0-9]+)/));
		my $inty = $1 if (($dY =~ /^([0-9]+)/) or ($dY =~ /^(-[0-9]+)/));
		if ($inty eq "-0") { 
			$y_upper->{$id} = 0;
			$y_lower->{$id} = -($yfi); 
		} else {
			$y_upper->{$id} = ($inty + 1) * $yfi;
			$y_lower->{$id} = $inty * $yfi;
		}
		
		if ($intx eq "-0") {
			$x_upper->{$id} = 0;
			$x_lower->{$id} = -($xfi); 
		} else {
			$x_upper->{$id} = ($intx + 1) * $xfi;
			$x_lower->{$id} = $intx * $xfi;
		}
		 
		
		#print "ix $intx dX $dX iy $inty dY $dY xup $x_upper->{$id} yup $y_upper->{$id} xlo $x_lower->{$id} ylo $y_lower->{$id}\n";
		
	}
	
	return ($x_upper, $y_upper, $x_lower, $y_lower);
}
	
	



##########
#
# OB_out
#
##~~~####
#
# print out to file
#


sub OB_out {
	my ($scors, $xu, $yu, $xl, $yl, $var, $len, $mw, $pI, $gr) = @_ or die "incomplete \@_ in assign_zscores\n";
	my $out = $var->{outfile};
#	open (O, ">$out") or die "cant open outfile $out $!\n";
	my %z;
	foreach my $k (keys %$xu) {
	#	print "$k\n";
		$z{$k} = $scors->{$xl->{$k}}{$xu->{$k}}{$yl->{$k}}{$yu->{$k}};
		#print "K $k xl $xl->{$k} xu $xu->{$k} yl $yl->{$k} yu $yu->{$k} scores $scors->{$xl->{$k}}{$xu->{$k}}{$yl->{$k}}{$yu->{$k}} \n";
		my $w = $mw->{$k};
		my @wt = @$w;
		if ($z{$k}) {
			return $z{$k};
			print  O "$k	$z{$k}	$gr->{$k}	$pI->{$k}	$wt[0]	$len->{$k}\n";
		} else {
			print O "$k	0\n";  # if $z{$k} is undef there will be no info so give Z = 0 as default
		}

	}
	
}

##########
#
# OB_out_nomwt
#
##~~~####
#
# print out to file
#


sub OB_out_nomwt {
	my ($scors, $xu, $yu, $xl, $yl, $var, $len, $pI, $gr) = @_ or die "incomplete \@_ in assign_zscores\n";
	my $out = $var->{outfile};
#	open (O, ">$out") or die "cant open outfile $out $!\n";
	my %z;
	foreach my $k (keys %$xu) {
#		print "$k\n";
		$z{$k} = $scors->{$xl->{$k}}{$xu->{$k}}{$yl->{$k}}{$yu->{$k}};
		return $z{$k};
		#print "K $k xl $xl->{$k} xu $xu->{$k} yl $yl->{$k} yu $yu->{$k} scores $scors->{$xl->{$k}}{$xu->{$k}}{$yl->{$k}}{$yu->{$k}} \n";
		if ($z{$k}) {
			return   "$k	$z{$k}	$gr->{$k}	$pI->{$k}	$len->{$k}\n";
		} else {
			print O "$k	0 $gr->{$k}	$pI->{$k}	$len->{$k}\n";  # if $z{$k} is undef there will be no info so give Z = 0 as default
		}

	}

}

1;

