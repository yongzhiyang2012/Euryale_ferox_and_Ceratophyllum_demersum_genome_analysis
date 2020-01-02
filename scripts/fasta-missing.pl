#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $input_fasta=shift or die "perl $0 input_fasta un_miss_criterion_percent windowsize\n";
my $miss_criterion=shift or die "perl $0 input_fasta un_miss_criterion_percent windowsize\n";
my $win=shift or die "perl $0 input_fasta un_miss_criterion_percent windowsize\n";
my %miss;
my %id;
my $fa;
if ($input_fasta=~/(gz|gzip)$/){
    open my $zcat,"zcat $input_fasta|" or die "$!";
    $fa=Bio::SeqIO->new(-fh=> $zcat ,-format=>"fasta");
}else{
    $fa=Bio::SeqIO->new(-file=>$input_fasta,-format=>'fasta')
}

while (my $seq=$fa->next_seq) {
    my $id=$seq->id;
    my $seq=$seq->seq;
    $seq=uc($seq);
    my @seq=split(//,$seq);
    $id{$id}++;
    for (my $i=0;$i<@seq;$i=$i+$win){
        my $key=int($i/$win);
        my $check=0;
        for (my $j=$i;$j<$i+$win;$j++){
            $check++ if $seq[$i] eq "-";
        }
        $miss{$key}++ if $check > 0;
    }
}

undef $fa;
$fa=Bio::SeqIO->new(-format=>"fasta",-file=>"$input_fasta");
while (my $seq=$fa->next_seq) {
    my $id=$seq->id;
    my $seq=$seq->seq;
    $seq=uc($seq);
    my @seq=split(//,$seq);
    print  ">$id\n";
    for (my $i=0;$i<@seq;$i=$i+$win){
        my $key=int($i/$win);
        $miss{$key}=0 if ! exists $miss{$key};
        my $miss_per=$miss{$key}/scalar(keys %id);
        next if (1-$miss_per)*100 < $miss_criterion;
        my $outseq;
        for (my $j=$i;$j<$i+$win;$j++){
            $outseq .= $seq[$j];
        }
        print  "$outseq";
    }
    print  "\n";
}

