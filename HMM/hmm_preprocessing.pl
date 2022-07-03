# This script is used to prepare input data before running the HMM function

# command line input argument
# 1. $rc_filename - recombination file name (input)
# 2. $sample_filename - sample file name (input)
# 3. $founder_filename - founder file name (input)
# 4. $out_file_name - output file name (output, format "arm pos K N pb[0-7] rf")

use warnings;
use strict;
use List::Util qw(sum max);
use Time::HiRes qw( time );

my $start = time();

my ($rc_filename, $sample_filename, $founder_filename, $out_file_name) = @ARGV;

# Array of chromosomes
my @chromosomes = ('2L', '2R', '3L', '3R', 'X');
my %chromosomes_index = (
'2L' => 1,
'2R' => 2,
'3L' => 3,
'3R' => 4,
'X' => 5,
);

# Mapping data to probability
my $p_large = 0.995;
my $p_small = 0.005;
my %d2p = ( './.' => -1, '0/0' => $p_small, '1/1' => $p_large);

# Array of founder ids
my @founders = ('F1','F2','F3','F4','F5','F6','F7','F8');
my %foundkey = (
1 => 'F1',
3 => 'F2',
4 => 'F3',
6 => 'F4',
9 => 'F5',
10 => 'F6',
14 => 'F7',
15 => 'F8',
);

# Number of founders per founder file (8 in our case)
my $n_founders = scalar(@founders);

# Make states hash with homozygous and heterozygous states and founders array (0,1)
my %statecodes=();
for (my $i=0; $i<$n_founders; ++$i)
{
	$statecodes{$founders[$i].$founders[$i]}=[$founders[$i], $founders[$i]];
    for (my $j=0; $j < $n_founders; ++$j)
    {
    	if ($i < $j)
    	{
    		$statecodes{$founders[$i].$founders[$j]}=[$founders[$i], $founders[$j]];
    	}
	}
}
my @statecodes = sort keys(%statecodes);

my %founders = ();
my $poscount = 0;
my $seq_length = 0;
my $N = 0;
my $K = 0;

open(TIMEFILE, ">time.out");

# Store the contents of the rc file in an array, one line per element
open my $rc_handle, '<', $rc_filename;
chomp(my @rc_array = <$rc_handle>);
close $rc_handle;

my $rc_size = scalar @rc_array;
my $rc_index = 0;

open(my $f_fh, '<:encoding(UTF-8)', $founder_filename) or die "Could not open file '$founder_filename' $!";
my $f_row; #founder file pointer
my $previous = undef;    # beginning of the previous line in the founder file
my $current  = tell $f_fh; # beginning of the current line in the founder file

foreach my $arm (@chromosomes)
{
    # read in the %positions data from a sample file per chromosome
    open(my $fh, '<:encoding(UTF-8)', $sample_filename) or die "Could not open file '$sample_filename' $!";

    # init variables for each chromosome
    $seq_length = 0;
    $N = 0;
    $K = 0;

    my $chrom_prev = "NA";
    my $pos_prev = -1;

    my $of_name = $out_file_name . '_' . $arm;
    open(OUTFILE, ">$of_name");

    while (my $row = <$fh>)
    {
        chomp $row;
        my $first_char = substr $row, 0, 1;

        # ignore first few lines with comments
        if($first_char eq "#")
        {
            next;
        }

        #read in data from the sample file
        my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $last_field) = split ' ', $row;

        # only process the right chromosom specified by $arm
        if($chrom ne $chromosomes_index{$arm})
        {
            next;
        }

        # locate the corresponding record in the founder file
        while ($f_row = <$f_fh>)
        {
            chomp $f_row;
            my $f_first_char = substr $f_row, 0, 1;
            my $isSmall = 0;

            # ignore lines starting with '#'
            if ($f_first_char eq "#")
            {
                next;
            }

            #read in data from founder file
            my ($f_chrom, $f_pos, $f_id, $f_ref, $f_alt, $f_qual, $f_filter, $f_info, $f_format, $v1, $v2, $v3, $v4, $v5, $v6, $v7, $v8) = split ' ', $f_row;

            if($f_chrom eq $chrom)
            {
                if(int($f_pos) < int($pos))
                {
                    $isSmall = 1;
                    next;
                }
                elsif(int($f_pos) > int($pos))
                {
                    $isSmall = 1;

                    if(defined $previous)
                    {
                        seek $f_fh, $previous, 0;  # seek to beginning of previous line (0 = SEEK_SET)
                    }
                    else
                    {
                        seek $f_fh, $current, 0;
                    }
                    last;
                }
                else # $f_pos == $pos, match found
                {
                    $isSmall = 1;

                    # locate the record, check if there is a missing data in the founder record; if so, skip the record
                    if (($v1 eq './.') or ($v2 eq './.') or ($v3 eq './.') or ($v4 eq './.') or ($v5 eq './.') or ($v6 eq './.') or ($v7 eq './.') or ($v8 eq './.'))
                    {
                        # break out the founder file to process the next record in the sample file
                        last;
                    }

                    #calculate N and K from the sample file record
                    my ($GT, $PL, $AD) = split ':', $last_field;
                    my @ad = split ',', $AD;

                    $N = sum @ad;

                    if ($N == 0)
                    {
                        last; # N =0, missing data
                    }

                    if(scalar @ad <= 0)
                    {
                        print "Error: size of AD in sample file <= 0";
                        last;
                    }
                    if ($GT eq "0/0")
                    {
                        $K = 0;
                    }
                    elsif ($GT eq "0/1")
                    {
                        $K = $ad[-1];
                    }
                    elsif ($GT eq "1/1")
                    {
                        $K = max @ad;
                    }
                    else
                    {
                        print "Error: wrong GT value in sample file with chromosome $arm! \n";
                        last;
                    }

                    # calculate the reference file value
                    my $rf = -1;
                    if ($chrom_prev eq $chrom)
                    {
                        $rf = compute_rf($pos_prev, $pos, $chrom);
                    }
                    if($chrom_prev ne "NA")
                    {
                        print OUTFILE "$rf\n"; #note that rf value is appended to the previous output record; so -1 should appear in the last record of a chromosome type
                    }
                    $chrom_prev = $chrom;
                    $pos_prev = $pos;

                    $seq_length = $seq_length + 1;

                    print OUTFILE "$arm $pos $K $N $d2p{$v1} $d2p{$v2} $d2p{$v3} $d2p{$v4} $d2p{$v5} $d2p{$v6} $d2p{$v7} $d2p{$v8} ";

                    # finish the data calculation with the matching record, and break out the founder file to process the next record in the sample file
                    last;
                }
            }
            else # $f_chrom != $chrom
            {
                if($isSmall == 1)
                {
                    $isSmall = 0;
                    last;
                }
                else
                {
                    next;
                }
            }

        } # end while ($f_row = <$f_fh>) founder file
        continue
        {
            $previous = $current;
            $current  = tell $f_fh;
        }

    } # end while (my $row = <$fh>) sample file

    print OUTFILE "-1\n"; #last rf value
    close(OUTFILE);

    close $fh;

    my $end = time();
    printf("Time: Finish processing %d records of chromosome %s in %.2f seconds \n", $seq_length, $arm, $end-$start);
    printf TIMEFILE "Time: Finish processing %d records of chromosome %s in %.2f seconds \n", $seq_length, $arm, $end-$start;
    $start = $end;

} #end for each chromosomes

close $f_fh;
close(TIMEFILE);

# sub function to compute rf value
sub compute_rf
{
    my ($pos1, $pos2, $arm) = @_;

    if($pos1 >= $pos2)
    {
        print "Error: compute_rf(), pos1 ($pos1) >= pos2 ($pos2), arm = $arm \n";
        return -1;
    }

    my $bin1 = 0;
    my $bin2 = 0;
    my $found_bin1 = 0;
    my $found_bin2 = 0;
    my $rf = 0;
    my $result = 0;

    while($rc_index < $rc_size)
    {
        my($temp, $RF) = split ',', $rc_array[$rc_index];
        my($chr, $left, $right) = split /[:.\s]+/, $temp;

        if (($chromosomes_index{$chr} eq $arm) and ($pos1 >= $left) and ($pos1 <= $right)) #found bin1
        {
            $found_bin1 = 1;
            $bin1 = $RF;

            if(($pos2 >= $left) and ($pos2 <= $right)) #pos1 and pos2 in the same bin
            {
                $found_bin2 = 1;
                $bin2 = $RF;
                last;
            }
            else #continue to find bin for pos2
            {
                $rc_index = $rc_index + 1;
                while($rc_index < $rc_size)
                {
                    my($temp2, $RF2) = split ',', $rc_array[$rc_index];
                    my($chr2, $left2, $right2) = split /[:.\s]+/, $temp2;

                    if(($pos2 >= $left2) and ($pos2 <= $right2)) #found bin2
                    {
                        $found_bin2 = 1;
                        $bin2 = $RF2;
                        last;
                    }
                    else #continue to find bin for pos2
                    {
                        $rc_index = $rc_index + 1;
                    }

                }

                if($found_bin2 == 1)
                {
                    last;
                }
            }
        }
        else #not found bin 1
        {
            $rc_index = $rc_index + 1;
        }
    }

    if (($found_bin1 == 0) or ($found_bin2 == 0))
    {
        print "Error: could not locate the right bins in the RC file, arm = $arm, pos1 = $pos1, pos2 = $pos2, found_bin1 = $found_bin1, found_bin2 = $found_bin2 \n";
        return -1;
    }

    # compute the processed rf value and return it
    if ($bin1 eq $bin2)
    {
        $rf = $bin1;
    }
    else #interpolate
    {
        my $slope = ($bin2-$bin1)/($pos2-$pos1);
        my $intercept = -$slope*$pos2 + $bin2;
        my $midpt = ($pos1+$pos2)/2;
        $rf = $slope*$midpt + $intercept;
    }

    # get distance between markers in bp
    my $dist = $pos2 - $pos1;

    if ($arm == 5) # chromosome eq 'X'
    {
        $result = 1E-8 * $rf * 33 / 2 * $dist;
    }
    else
    {
        $result = 1E-8 * $rf * 25 / 2 * $dist;
    }

    return $result;
}
