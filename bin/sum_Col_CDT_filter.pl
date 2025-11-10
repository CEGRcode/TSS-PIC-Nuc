#!/usr/bin/perl

# Make composite file from CDT

die "usage:\t\tperl sum_Col_CDT.pl\tInput_CDT_File\tOutput_TAB_File\nExample:\tperl sum_Col_CDT.pl input.cdt composite.out\n" unless $#ARGV == 1;

my ($input, $output) = @ARGV;
open(IN, "<$input") or die "Can't open $input for reading!\n";
open(OUT, ">$output") or die "Can't open $output for writing!\n";

my $NCOL = 0;
my @Y;         # Sum per column
my @count;     # Count of valid values per column

while (<IN>) {
    chomp;
    my @temparray = split(/\t/, $_);

    # Header line
    if (/YORF/) {
        print OUT join("\t", @temparray[1..$#temparray]), "\n";
        $NCOL = $#temparray;
        next;
    }

    if ($NCOL != $#temparray) {
        print "Error!!! Inconsistent CDT window size\n";
        print "(Num Columns=", $#temparray, ")\n", join("\t", @temparray), "\n";
        exit;
    }

    for (my $x = 2; $x <= $#temparray; $x++) {
        my $val = $temparray[$x];
        # Skip if empty or "-"
        if (defined($val) && $val ne '' && $val ne '-') {
            if ($val =~ /^-?\d+\.?\d*$/) {
                $Y[$x - 2] += $val;
                $count[$x - 2]++;
            } else {
                warn "Non-numeric value ignored: $val at line $. column $x\n";
            }
        }
    }
}

# Calculate average; replace NaN with 0
for (my $t = 0; $t <= $#Y; $t++) {
    if ($count[$t]) {
        $Y[$t] = $Y[$t] / $count[$t];
    } else {
        $Y[$t] = 0;
    }
}

print OUT $input, "\t", join("\t", @Y), "\n";

close IN;
close OUT;
