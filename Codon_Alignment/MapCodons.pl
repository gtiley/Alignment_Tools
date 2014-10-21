#!/usr/bin/perl
use warnings;

## Contact: George Tiley - gtiley@ufl.edu
## Usage: ./MapCodons.pl


#Notes
#--------------------------------------------------------------------------------------#
#The first nucleotide is assumed to be the first codon position
#If the sequence is not in first position, the nucleotide sequence will be trimmed to match the AA sequence
#Currently, only a small range is supported, this is not a highly iterative process- really screwy nucleotide sequences will not return a codon alignment
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
system 'ls *.protein.fasta > Proteins';
@Files = ();
#open FH99,'<', "$ARGV[0]";
open FH99, "<Proteins";
while (<FH99>)
{
	if (/((\S+)\.protein\.fasta)/)
	{
		$ThisFile = $1;
		push @Files, $ThisFile
	}
}
close FH99;
foreach $FILE (@Files)
{
	if ($FILE =~ m/((\S+)\.protein\.fasta)/)
	{
		$ProteinFile = $1;
		$Gene = $2;
		$NucFile = "$Gene.nuc.fasta";
		$AlignedProtein = "$Gene.protein.align.fasta";
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
    			@AAarray = split //, $ProteinSeqs{$taxon};
    			$N_AA = scalar(@AAarray);
    			$leftcheck = 0; #truth assessment
    			$leftmargin = 0; #how much to trim from the left (ex: 5'UTRs)
    			$trim1 = "I_AM_A_SPACE_HOLDER_THAT_IS_AT_LEAST_30_CHARACTERS_LONG";
    			$rightcheck = 0; #truth assessment
    			$rightmargin = 0; #how much to trim from the right (ex: polyAs)
    			$trim2 = "I_AM_A_SPACE_HOLDER_THAT_IS_AT_LEAST_30_CHARACTERS_LONG";
##########################################################################################    			
				while ($leftcheck != 1 && length($trim1) >= 30)
				{
					$trim1 = substr($seq, $leftmargin, (length($seq)));
    				$Site1 = substr($trim1, 0, 3);
    				$Site2 = substr($trim1, 3, 3);
    				$Site3 = substr($trim1, 6, 3);
    				$Site4 = substr($trim1, 9, 3);
    				$Site5 = substr($trim1, 12, 3);
    				$Site6 = substr($trim1, 15, 3);
    				$Site7 = substr($trim1, 18, 3);
    				$Site8 = substr($trim1, 21, 3);
    				$Site9 = substr($trim1, 24, 3);
    				$Site10 = substr($trim1, 27, 3);
    				if ($translation_hash{$Site1} eq $AAarray[0] and $translation_hash{$Site2} eq $AAarray[1] and $translation_hash{$Site3} eq $AAarray[2] and $translation_hash{$Site4} eq $AAarray[3] and $translation_hash{$Site5} eq $AAarray[4] and $translation_hash{$Site6} eq $AAarray[5] and $translation_hash{$Site7} eq $AAarray[6] and $translation_hash{$Site8} eq $AAarray[7] and $translation_hash{$Site9} eq $AAarray[8] and $translation_hash{$Site10} eq $AAarray[9])
    				{
						$leftcheck++;
					}
					else
    				{
						$leftmargin++;
					}
				}
				while ($rightcheck != 1 && length($trim2) >= 30)
				{
					$trim2 = substr($trim1, 0, (length($trim1) - $rightmargin));
					$SiteN = substr($trim2, (length($trim2) - 3), 3); 
    				$SiteN_minus_1 = substr($trim2, (length($trim2) - 6), 3); 
    				$SiteN_minus_2 = substr($trim2, (length($trim2) - 9), 3);
    				$SiteN_minus_3 = substr($trim2, (length($trim2) - 12), 3);
    				$SiteN_minus_4 = substr($trim2, (length($trim2) - 15), 3);
    				$SiteN_minus_5 = substr($trim2, (length($trim2) - 18), 3);
    				$SiteN_minus_6 = substr($trim2, (length($trim2) - 21), 3);
    				$SiteN_minus_7 = substr($trim2, (length($trim2) - 24), 3);
    				$SiteN_minus_8 = substr($trim2, (length($trim2) - 27), 3);
    				$SiteN_minus_9 = substr($trim2, (length($trim2) - 30), 3);
    				if ($translation_hash{$SiteN} eq $AAarray[($N_AA - 1)] and $translation_hash{$SiteN_minus_1} eq $AAarray[($N_AA - 2)] and $translation_hash{$SiteN_minus_2} eq $AAarray[($N_AA - 3)] and $translation_hash{$SiteN_minus_3} eq $AAarray[($N_AA - 4)] and $translation_hash{$SiteN_minus_4} eq $AAarray[($N_AA - 5)] and $translation_hash{$SiteN_minus_5} eq $AAarray[($N_AA - 6)] and $translation_hash{$SiteN_minus_6} eq $AAarray[($N_AA - 7)] and $translation_hash{$SiteN_minus_7} eq $AAarray[($N_AA - 8)] and $translation_hash{$SiteN_minus_8} eq $AAarray[($N_AA - 9)] and $translation_hash{$SiteN_minus_9} eq $AAarray[($N_AA - 10)])
					{
						$rightcheck++;
					}
					else
					{
						$rightmargin++;
					}
				}
##########################################################################################    			
				if ($leftcheck == 1 && $rightcheck ==1) #Check that the nucleotide actually corresponds with the protein
				{
					$NucSeqs{$taxon} = $trim2;
    				$numcodons = (length $NucSeqs{$taxon}) / 3;
					for (0..($numcodons -1))
					{
						$codonnum = $_;
						$site = $codonnum * 3;
						$codon = substr($NucSeqs{$taxon}, $site, 3);
						if ($translation_hash{$codon} eq $AAarray[$codonnum])
						{
							$codonhash{$taxon}{$codonnum} = $codon; 
						}
						elsif ($translation_hash{$codon} ne $AAarray[$codonnum])
						{
							$codonhash{$taxon}{$codonnum} = $codon; 
							print OUT99 "$Gene\t$taxon\t$codonnum\n";
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