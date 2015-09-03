#!/usr/local/bin/perl

# this script collect the values from the metrixEQTL pipeline
# first argument is the path of the folder where the output files from matrixeqtl are
# second argument is the path of the the output file

$selDir = $ARGV[0]; # Directory where the sentinalised eQTLs are
$fileName =$ARGV[1]; # file where all the eQTLs are going to be saved
opendir my $dir, $selDir or die "Cannot open directory: $!";
	print "start\n";
  my @files = readdir $dir;
	my $numFiles = 0; # numbers of files that have been read
	open FH, ">>$fileName";
       # Print header
       # print FH "Gene\tSNP\ttranscript\tbeta\tt-stat\tp-value\tFDR\n";
       print FH "snps gene tstat pvalue FDR beta myFDR degree\n";
 close FH;
	foreach my $file (@files)
	{      	
			# open the sub folder
			open FILE, "<" , $selDir."/".$file;
      my $rows = 0;
      
			while (<FILE>) 
			{	
        ##print $rows;
        
        
        if ($rows>0)
				{
					# check whether there is a chromosome position in a read
					if ( $_ =~ m/chr/)
             {
   		               
       #open the new file where append the all the values
						         open FH, ">>$fileName";
       # Print everything in the file 
                     ##print FH $file."\t".$_;
                     print FH $_;
            close FH;
						##print $_;
            #print $file."\t".$hash{"SNP"}."\t".$hash{"gene"}."\t".$hash{"beta"}."\t".$hash{"t-stat"}."\t".$hash{"p-value"}."\t".$hash{"FDR"}."\n";
            
				    }
            
        }
        
				$rows=$rows+1;
        #undef $hash;
			}
		 close(FILE);
      $numFiles=$numFiles+1;
      if ($numFiles % 1000 == 999)
      {
         print $numFiles."\n";
      }
	}
print "numbers of files read ".($numFiles-2)."\n";
closedir $dir;
