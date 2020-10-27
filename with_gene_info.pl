#!/usr/bin/perl -w
use strict;
die "Usage: perl $0 <miRNAs.info.xls> <gene.gff>\n" unless @ARGV == 2;

my $in = shift;
my $gff = shift;

open GFF, $gff || die "$!\n";
my ($Gff, $mRNA);
my $i = 0;
while( <GFF> ){
    chomp;
    my @tmp = split(/\s+/, $_);
    if( $tmp[2] = /mRNA/ ){
        $i++;
        my $gene = $1 if $tmp[-1] =~ /ID=(.*);/;
        $Gff->{$tmp[0]}{$gene}{'order'} = $i;
		push(@{$mRNA->{$tmp[0]}{'mRNA'}}, [ $tmp[3], $tmp[4] ]);
    }
    if( $tmp[2] = /cds/i ){
        my $gene = $1 if $tmp[-1] =~ /Parent=(.*);/;
        $Gff->{$tmp[0]}{$gene}{'strand'} = $tmp[6];
        push(@{$Gff->{$tmp[0]}{$gene}{'exon'}}, [ $tmp[3], $tmp[4] ]);
    }
}
close GFF;

for my $key ( sort %$Gff ){
	my $chr_p = $Gff->{$key};
	for my $gene ( sort { $chr_p->{$a}{'order'} <=> $chr_p->{$b}{'order'} } keys %$chr_p ){
		my $gene_p = $chr_p->{$gene};
		my @Exon = sort {$a->[0] <=> $b->[0]} @{$gene_p->{'exon'}};
		my $len = scalar(@Exon);
		my @all_exons;
		if ( $len > 1 ){
			for my $i ( @Exon ){
				push(@all_exons, $$i[0], $$i[1]);
			}
		}
		for(my $j=1; $j<@all_exons-1; $j+=2){
			my $jj = $j+1;
			my $intron_start = $all_exons[$j] + 1;
			my $intron_end = $all_exons[$jj] - 1;
			push(@{$Gff->{$key}{$gene}{'intron'}}, [ $intron_start, $intron_end ]);
		}
	}
}

open IN, $in || die "$!\n";
while( <IN> ){
	chomp;
	my @tmp = split(/\s+/, $_);
	if( exists $Gff->{$tmp[1]} ){
		my $chr_p = $Gff->{$tmp[1]};
		for my $gene ( sort { $chr_p->{$a}{'order'} <=> $chr_p->{$b}{'order'} } keys %$chr_p ){
			my $gene_p = $chr_p->{$gene};
			if( $gene_p->{strand} eq $tmp[-3] ){
				for my $i ( @{$gene_p->{'exon'}} ){
					if( $tmp[2] >= $$i[0] && $tmp[3] <= $$i[1] ){
						print "$_\texon\n";
					}elsif( $tmp[2] > $$i[0] && $tmp[3] > $$i[1] && $tmp[2] < $$i[1] ){
						print "$_\texon\n";
					}elsif( $tmp[2] < $$i[0] && $tmp[3] > $$i[0] && $tmp[3] < $$i[1] ){
						print "$_\texon\n";
					}
				}
				for my $j ( @{$gene_p->{'intron'}} ){
					if( $tmp[2] >= $$j[0] && $tmp[3] <= $$j[1] ){
						print "$_\tintron\n";
					}elsif( $tmp[2] > $$j[0] && $tmp[3] > $$j[1] && $tmp[2] < $$j[1] ){
						print "$_\tintron\n";
					}elsif( $tmp[2] < $$j[0] && $tmp[3] > $$j[0] && $tmp[3] < $$j[1] ){
						print "$_\tintron\n";
					}
				}
			}
		}
	}
}
close IN;
