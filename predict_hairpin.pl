#!/usr/bin/perl -w
use strict;
die "Usage: perl $0 <potential.fa>\n" unless @ARGV == 1;

my $in = shift;

my $prefix = $1 if $in =~ /(.*)\.fa$/;

my $pwd = `pwd`;
chomp $pwd;
my $out = $pwd."/miRNAFold_result";
mkdir("$out") unless -d $out;

open IN, $in || die "$!\n";
my ($id, %hash, %pos);
while( <IN> ){
	chomp;
	if( />(\S+)\s+(\S+):(\d+)-(\d+):(\d+)-(\d+):(\S+)/ ){
		$id = $1;
		$pos{$id} = [ $2, $3, $4, $5, $6, $7 ];
	}else{
		$hash{$id} .= $_;
	}
}
close IN;

for my $key ( sort keys %hash ){
	open OUT, ">$out/$key.fa" || die "$!\n";
	print OUT ">$key\n$hash{$key}\n";
	close OUT;
	`/Software/miRNAFoldExecutable/miRNAFold -s $out/$key.fa -L 100 -parameter Software/miRNAFoldExecutable/DATA/DEFAULT.DAT
 -o $out/$key.fa.miRNAFold.out`;
}


for my $key ( sort keys %hash ){
	open OUT, ">$out/$key.premature.fa" || die "$!\n";
	open ZZ, "$out/$key.fa.miRNAFold.out" || die "$!\n";
	my $count = 1;
	while( <ZZ> ){
		chomp;
		my @tmp = split(/\s+/, $_);
		my $mirLen = $pos{$key}[4] - $pos{$key}[3];
		if( $tmp[0] <= 110 && $tmp[1] >= 110+$mirLen ){
		   my $seq = substr($hash{$key}, $tmp[0] - 1, $tmp[1] - $tmp[0]);
		   my $start = $pos{$key}[1] + $tmp[0];
		   my $end = $pos{$key}[1] + $tmp[1];
		   print OUT ">$key", "_premature_D$count\t$pos{$key}[0]:$start-$end:$pos{$key}[3]-$pos{$key}[4]:$pos{$key}[-1]\n$seq\n";
		   $count++;
	    }
	}
	close ZZ;
	close OUT;
}

`cat $out/*premature.fa >$pwd/$prefix.premature.raw.fa`;

my (%Info, %Raw_pre);
open IN, "$pwd/$prefix.premature.raw.fa" || die "$!\n";
while( <IN> ){
	chomp;
	if( />(\S+)\s+(\S+)/ ){
		$id = $1;
		$Info{$id} = $2;
	}else{
		$Raw_pre{$id} .= $_;
	}
}
close IN;

my $RNAfold_dir = $pwd."/RNAfold_reassessment";
mkdir($RNAfold_dir) unless -d $RNAfold_dir;
chdir($RNAfold_dir);
`/Software/ViennaRNA-2.4.14/bin/RNAfold $pwd/$prefix.premature.raw.fa -p >$prefix.premature.raw.RNAfold.out`;
chdir($pwd);

open IN, "$RNAfold_dir/$prefix.premature.raw.RNAfold.out" || die "fail to open: $RNAfold_dir/$prefix.premature.raw.RNAfold.out\n";
my ($list, $cluster, $seq);
while( <IN> ){
	chomp;
	if( />(\S+)/ ){
		$seq = $1;
		$cluster = $seq;
		$cluster =~ s/_premature_D\d+$//g;
	}
	<IN>;
	my $mfe = (split(/\s+/, <IN>))[-1];
	$mfe =~ s/\(//g;
	$mfe =~ s/\)//g;
	$list->{$cluster}{$mfe} = $seq;
	<IN>;
	<IN>;
	<IN>;
}
close IN;

my (@Best, @Score);
open OUT, ">$pwd/$prefix.premature.best.fa" || die "$!\n";
open OUT1, ">$pwd/$prefix.premature.best.info.xls" || die "$!\n";
for my $key ( sort keys %$list ){
	my $Cluster = $list->{$key};
	for my $single ( sort keys %$Cluster ){
		push(@Best, $Cluster->{$single});
		push(@Score, $single);
	}
	print OUT ">$Best[0]\t$Info{$Best[0]}\n$Raw_pre{$Best[0]}\n";
	my @tmp = split(/\:/, $Info{$Best[0]});
	$tmp[1] =~ s/-/\t/g;
	$tmp[2] =~ s/-/\t/g;
	my @pos = split(/\s+/, $tmp[2]);
	my $miRNA = substr($hash{$key}, 110, $pos[1]-$pos[0]);
	print OUT1 "$key\t$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$miRNA\t$Score[0]\n";
	undef @Best;
	undef @Score;
}
close OUT;
close OUT1;
