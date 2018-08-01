$len = length($F[9]);
$len = $len / 2;
@CIGAR = split(/([0-9]+[SMIDNHXP])/, $F[5]);
my $hash;
map { push(@{$hash->{$2}}, $1) if (/(\d+)([SMIDNHXP])/) } @CIGAR;
foreach my $softclip (@{$hash->{S}}) {
    if ($softclip >= $len) {
        print
    }
}