#!/usr/bin/perl -w
use strict;
die "perl $0 <*.m8> <query.fa>\n" unless @ARGV == 2;

my $m8 = shift;
my $query = shift;

my ($id, %Query);
open IN, $query || die "$!\n";
while( <IN> ){
	chomp;
	if( />(\S+)/ ){
		$id = $1;
	}else{
		$Query{$id} .= $_;
	}
}
close IN;

open IN, $m8 || die "$!\n";
my %hash;
while( <IN> ){
	chomp;
	my @tmp = split(/\s+/, $_);
	my $diff = length($Query{$tmp[0]}) - $tmp[3];
	if( $tmp[4] + $diff > 2 ){
		next;
	}else{
		push(@{$hash{$tmp[0]}}, $_);
	}
}
close IN;

foreach my $key ( sort keys %hash ){
	my $pp = $hash{$key};
	for(my $i = 0; $i <@{$hash{$key}}; $i++){
		my $key_Dup = "$key-D".($i+1);
		$pp->[$i] =~ s/$key/$key_Dup/ if $i > 0;
		print $pp->[$i],"\n";
	}
}
