
my $res;
my %vdw;
open FILE, "vdw-radii.dat";
while(<FILE>)
{
      if($_=~/^RESIDUE/)
      {
            chomp $_;
            my @dat=split(/\s+/, $_);
            $res=$dat[2];
      }
      if($_=~/^ATOM/)
      {
            chomp $_;
            my @dat=split(/\s+/, $_);
            my @atm=split(//, $dat[1]);
            
            print "name2radius[ \"$res$dat[1]\" ] = $dat[2] ;\n";
            
            if(defined $vdw{$atm[0]})
            {
                  if($vdw{$atm[0]} != $dat[2])
                  {
                        #print "@@ $res $_"." $vdw{$atm[0]}"."\n";
                        
                  }
            }
            else
            {
                  #print @dat." ".$atm[0]." ".$dat[2]."\n";
                  #$vdw{$atm[0]}=$dat[2];
            }
      }
}
close FILE;
