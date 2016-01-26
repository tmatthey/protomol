# Simple grep to retrieve the compiler version
use strict;
my $a = "";
my $n = 0;
while (<>) {
    next unless (/\S/ and 
		 not /\// and 
		 not /exit/ and 
		 not /error/ and 
		 not /Error/ and 
		 not /ERROR/ and 
		 not /invalid/ and 
		 not /missing/ and 
		 not /FOR NON-COMMERCIAL USE ONLY/ and
		 not /WARNING/ and
		 not /Warning/ and
		 not /warning/ and
		 not /option/ and 
		 not /\\/ and
		 $n < 10); 
    $n++;
    $_=~ s/\n/\\n/g; 
    $a .= $_; 
}
$a =~ s/\\n$//g;
print $a;
