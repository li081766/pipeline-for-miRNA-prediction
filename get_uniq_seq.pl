#!/usr/bin/perl -w
use strict;
die "\nThe script is used to select the unique hits and extract the extension of 110nt sequence of each hits.\n\n\tUsage: perl $0 <genome> <*.m8>\n\n"
unless @ARGV == 2;

my $genome = shift;
my $m8 = shift;

my ($id, %Genome);
open IN, $genome || die "$!\n";
while( <IN> ){
	chomp;
	if( />(\S+)/ ){
		$id = $1;
	}else{
		$Genome{$id} .= $_;
	}
}
close IN;

my ($Hits, $strand);
open IN, $m8 || die "$!\n";
while( <IN> ){
	chomp;
	my @tmp = split(/\s+/, $_);
	if( $tmp[8] < $tmp[9] ){
		$strand = "+";
	}else{
		$strand = "-";
	}
	push(@{$Hits->{$tmp[1]}{$strand}}, [ $tmp[8], $tmp[9], $tmp[0] ]);
}
close IN;

open OUT, ">$m8.info.xls" || die "$!\n";
print OUT "scaf\tstart\tend\tstrand\tmature_start\tmature_end\t5-end\t3-end\thomo\n";
for my $seq ( sort keys %$Hits ){
	my $scaf = $Hits->{$seq};
	for my $strand ( sort keys %$scaf ){
		my @Mature = @{$scaf->{$strand}};
		my $len = scalar(@Mature);
		my ($start, $upstream, $mirna, $downstream, $end);
		if( $len == 1 ){
			if( $strand eq "+" ){
				$start = $Mature[0][0] - 110;
				$start = 1 if $start < 0;
				$end = $Mature[0][1] + 110;
				$end = $Mature[0][1] + (110 - length($Genome{$seq}) - $Mature[0][1]) if $end > length($Genome{$seq});
				$upstream = substr($Genome{$seq}, $start - 1, 110);
				$mirna = substr($Genome{$seq}, $Mature[0][0] - 1, $Mature[0][1] - $Mature[0][0]);
				$mirna = lc($mirna);
				my $len = length($mirna);
				$downstream = substr($Genome{$seq}, $end - 111, 110);
				print ">$Mature[0][2]_L\t$seq:$start-$end:$Mature[0][0]-$Mature[0][1]:$strand\n";
				print "$upstream", "$mirna", "$downstream\n";
				if( $start == 1 && $end + 110 < length($Genome{$seq}) ){
					print OUT "$seq\t$start\t$end\t$strand\t$Mature[0][0]\t$Mature[0][1]\tYes\t-\t$Mature[0][2]_L\n";
				}elsif( $start == 1 && $end + 110 > length($Genome{$seq}) ){
					print OUT "$seq\t$start\t$end\t$strand\t$Mature[0][0]\t$Mature[0][1]\tYes\tYes\t$Mature[0][2]_L\n";
				}elsif( $start > 1 && $end + 110 > length($Genome{$seq}) ){
					print OUT "$seq\t$start\t$end\t$strand\t$Mature[0][0]\t$Mature[0][1]\t-\tYes\t$Mature[0][2]_L\n";
				}else{
					print OUT "$seq\t$start\t$end\t$strand\t$Mature[0][0]\t$Mature[0][1]\t-\t-\t$Mature[0][2]_L\n";
				}
			}
			if( $strand eq "-" ){
				$start = $Mature[0][1] - 110;
				$start = 1 if $start < 0;
				$end = $Mature[0][0] + 110;
				$end = $Mature[0][0] + (110 - length($Genome{$seq}) - $Mature[0][0]) if $end > length($Genome{$seq});
				$upstream = substr($Genome{$seq}, $start - 1, 110);
				$mirna = substr($Genome{$seq}, $Mature[0][1] - 1, $Mature[0][0] - $Mature[0][1]);
				$mirna = lc($mirna);
				my $len = length($mirna);
				$downstream = substr($Genome{$seq}, $end - 111, 110);
				my $premature = $upstream.$mirna.$downstream;
				$premature = Complement_Reverse($premature);
				print ">$Mature[0][2]_L\t$seq:$start-$end:$Mature[0][1]-$Mature[0][0]:$strand\n";
				print "$premature\n";
				if( $start == 1 && $end + 110 < length($Genome{$seq}) ){
					print OUT "$seq\t$start\t$end\t$strand\t$Mature[0][0]\t$Mature[0][1]\tYes\t-\t$Mature[0][2]_L\n";
				}elsif( $start == 1 && $end + 110 > length($Genome{$seq}) ){
					print OUT "$seq\t$start\t$end\t$strand\t$Mature[0][0]\t$Mature[0][1]\tYes\tYes\t$Mature[0][2]_L\n";
				}elsif( $start > 1 && $end + 110 > length($Genome{$seq}) ){
					print OUT "$seq\t$start\t$end\t$strand\t$Mature[0][0]\t$Mature[0][1]\t-\tYes\t$Mature[0][2]_L\n";
				}else{
					print OUT "$seq\t$start\t$end\t$strand\t$Mature[0][0]\t$Mature[0][1]\t-\t-\t$Mature[0][2]_L\n";
				}
			}
		}

		if( $len > 1 ){
			@Mature = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @Mature;
			my (@tmp, @dup_list, @dup_id, %Uniq, %hash);
			for(my $i = 0; $i < @Mature; $i++){
				my $j = $i + 1;
				if( $j < scalar @Mature ){
					if( ($Mature[$j][0] - $Mature[$i][0]) < 10 ){
						push(@tmp, $Mature[$i][0], $Mature[$i][1], $Mature[$j][0], $Mature[$j][1]);
						push(@dup_list, $i, $j);
						push(@dup_id, $Mature[$i][2], $Mature[$j][2]);
					}else{
						if( @tmp ){
							@tmp = sort { $a <=> $b } @tmp;
							@dup_id = grep { ++$hash{$_} == 1 } @dup_id;
							$Uniq{$dup_id[0]} = [ $tmp[0], $tmp[-1], [@dup_id] ];
							undef @tmp;
							undef @dup_id;
						}
					}
				}
			}
			if( @tmp ){
				@tmp = sort { $a <=> $b } @tmp;
				@dup_id = grep { ++$hash{$_} == 1 } @dup_id;
				$Uniq{$dup_id[0]} = [ $tmp[0], $tmp[-1], [@dup_id] ];
			}
			@dup_list = grep { ++$hash{$_} == 1 } @dup_list;
			for my $zz ( @dup_list ){
				$hash{$zz} = 1;
			}
			for(my $i = 0; $i < @Mature; $i++){
				$Uniq{$Mature[$i][2]} = [ $Mature[$i][0], $Mature[$i][1], [$Mature[$i][2]] ] unless exists $hash{$i};
			}
			for my $key ( sort keys %Uniq ){
				if( $Uniq{$key}[0] < $Uniq{$key}[1] ){
					$start = $Uniq{$key}[0] - 110;
					$start = 1 if $start < 0;
					$end = $Uniq{$key}[1] + 110;
					$end = $Uniq{$key}[1] + (110 - length($Genome{$seq}) - $Uniq{$key}[1]) if $end > length($Genome{$seq});
				}else{
					$start = $Uniq{$key}[1] - 110;
					$start = 1 if $start < 0;
					$end = $Uniq{$key}[0] + 110;
					$end = $Uniq{$key}[0] + (110 - length($Genome{$seq}) - $Uniq{$key}[0]) if $end > length($Genome{$seq});
				}
				$upstream = substr($Genome{$seq}, $start - 1, 110);
				if( $Uniq{$key}[0] < $Uniq{$key}[1] ){
					$mirna = substr($Genome{$seq}, $Uniq{$key}[0] - 1, $Uniq{$key}[1] - $Uniq{$key}[0]);
				}else{
					$mirna = substr($Genome{$seq}, $Uniq{$key}[1] - 1, $Uniq{$key}[0] - $Uniq{$key}[1]);
				}
				$mirna = lc($mirna);
				my $len = length($mirna);
				$downstream = substr($Genome{$seq}, $end - 111, 110);
				my $premature = $upstream.$mirna.$downstream;
				$premature = Complement_Reverse($premature) if $strand eq "-";
				print ">$Uniq{$key}[2][0]\t$seq:$start-$end:$Uniq{$key}[0]-$Uniq{$key}[1]:$strand\n" if $Uniq{$key}[0] < $Uniq{$key}[1]
;
				print ">$Uniq{$key}[2][0]\t$seq:$start-$end:$Uniq{$key}[1]-$Uniq{$key}[0]:$strand\n" if $Uniq{$key}[0] > $Uniq{$key}[1]
;
				print "$premature\n";
				if( $start == 1 && $end + 110 < length($Genome{$seq}) ){
					print OUT "$seq\t$start\t$end\t$strand\t$Uniq{$key}[0]\t$Uniq{$key}[1]\tYes\t-\t", join(";", @{$Uniq{$key}[2]})
, "\n";
				}elsif( $start == 1 && $end + 110 > length($Genome{$seq}) ){
					print OUT "$seq\t$start\t$end\t$strand\t$Uniq{$key}[0]\t$Uniq{$key}[1]\tYes\tYes\t", join(";", @{$Uniq{$key}[2]
}), "\n";
				}elsif( $start > 1 && $end + 110 > length($Genome{$seq}) ){
					print OUT "$seq\t$start\t$end\t$strand\t$Uniq{$key}[0]\t$Uniq{$key}[1]\t-\tYes\t", join(";", @{$Uniq{$key}[2]})
, "\n";
				}else{
					print OUT "$seq\t$start\t$end\t$strand\t$Uniq{$key}[0]\t$Uniq{$key}[1]\t-\t-\t", join(";", @{$Uniq{$key}[2]}),
"\n";
				}
			}
		}
	}
}
close OUT;



##### Sub Routines #####
sub Complement_Reverse{
	my $seq = shift;
	$seq =~ tr/AGCTagct/TCGAtcga/;
	$seq = reverse($seq);

	return $seq;
}
