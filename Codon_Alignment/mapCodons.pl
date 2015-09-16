#!/usr/bin/perl
#use warnings;

#Notes
#--------------------------------------------------------------------------------------#
#The first nucleotide is assumed to be the first codon position
#If the sequence is not in first position, the nucleotide sequence will be trimmed to match the AA sequence
#A codon alignment will not be returned if there are interexonic untranslated regions. 5' UTR and 3' poly A are ok!
#There is limited support for handeling degenrate bases - this is handled by defining codons with degenerate bases in the genetic code
#--------------------------------------------------------------------------------------#
open OUT99, '>', "CodonMappingErrors.log";
#opening file with amino acid translations
#Be sure to use the appropriate translation table!
open FH1, "<translations_degenerate.txt";
#define a hash for amino acids
%translation_hash = ();
while (<FH1>)
{
    if (/(\S+)\s+(\S+)/)
    {
	$translation_hash{$1} = $2;
    }
}
close FH1;

# Step 1 - Trim nucleotides to translate in 1st codon position
system 'ls *.aa.group.fasta > Proteins';
@Files = ();
#open FH99,'<', "$ARGV[0]";
open FH99, "<Proteins";
while (<FH99>)
{
	if (/((\S+)\.aa\.group\.fasta)/)
	{
		$ThisFile = $1;
		push @Files, $ThisFile
	}
}
close FH99;
foreach $FILE (@Files)
{
	if ($FILE =~ m/((\S+)\.aa\.group\.fasta)/)
	{
		$ProteinFile = $1;
		$Gene = $2;
		$NucFile = "$Gene.nuc.group.fasta";
		$AlignedProtein = "$Gene.aa.group.muscle.fasta";
		%ProteinSeqs = ();
		%NucSeqs = ();
		%codonhash = ();
		open FH98, "$ProteinFile";
		while (<FH98>)
		{
			if (/>(\S+)/)
			{
				$Tax = $1;
			}
			elsif(/(\S+)/)
			{
				$Seq = $1;
				$ProteinSeqs{$Tax} = $Seq;
			}
		}
		close FH98;
		
		open FH97, "$NucFile";
		while (<FH97>)
		{
			if (/^>(\S+)/)
			{
				$taxon = $1;
			}
			elsif (/(\S+)/)
			{
				$seq = $1;	
    			$NucSeqs{$taxon} = $seq;
    			$check = 0;
    			while($check != 1 && length($NucSeqs{$taxon}) > 0)
    			{
    				$numcodons = (length $NucSeqs{$taxon}) / 3;
    				$newseq = "";
    				$checkseq = "";
					for (0..($numcodons))
					{
						$codonnum = $_;
						$site = $codonnum * 3;
						$codon = substr($NucSeqs{$taxon}, $site, 3);

							$aa = $translation_hash{$codon};
							$checkseq = $checkseq . $aa;
							$newseq = $newseq . $codon;

					}
					if (index($checkseq,$ProteinSeqs{$taxon}) >= 0)
					{
						$NucSeqs{$taxon} = substr($newseq, index($checkseq,$ProteinSeqs{$taxon}) * 3, length($ProteinSeqs{$taxon}) * 3);
						$numcodons = (length $NucSeqs{$taxon}) / 3;
    					for (0..($numcodons))
						{
							$codonnum = $_;
							$site = $codonnum * 3;
							$codon = substr($NucSeqs{$taxon}, $site, 3);
							$codonhash{$taxon}{$codonnum} = $codon;
						}
						$check = 1;
					}
					elsif(index($checkseq,$ProteinSeqs{$taxon}) < 0)
					{
						$NucSeqs{$taxon} = substr($newseq, 1, length($newseq));
						if (length($newseq) == 1)
						{
							print OUT99 "$Gene\t$taxon\n";
						}
					} 
			
				}
##########################################################################################    			
			}
		}
		close FH97;
		
		#Step 2, back-translate aligned amino acids based on trimmed nucleotide sequence	
		open OUT1, ">$Gene.codonalign.fasta";
		@seqarray = ();
		open FH96, "$AlignedProtein";
		while (<FH96>)
		{
			if (/^>(\S+)/)
			{
				$TAXON = $1;
				print OUT1 ">$TAXON\n";
			}
			elsif (/(\S+)/)
			{
				$SEQ = $1;
				@seqarray = split //, $SEQ;
				$charcount = 0;
				foreach $CHAR (@seqarray)
				{
					if ($CHAR eq "-")
					{
						print OUT1 "---";
					}
					else
					{
						if (exists $codonhash{$TAXON}{$charcount})
						{
							print OUT1 "$codonhash{$TAXON}{$charcount}";
						}
						elsif (! exists $codonhash{$TAXON}{$charcount})
						{
							print OUT1 "NNN";
						}
						$charcount++;
					}
				}		
				print OUT1 "\n";
			}			
		}
		close FH96;
		close OUT1;
	}
}
exit;
