#!/usr/local/bin/perl

# $Id: replace_header.pl 554 2003-11-12 02:36:24Z wistaria $

$header = "header.txt";

foreach $file (@ARGV) {
    if (-f $file) {
	if ($file =~ /Makefile.in$/) {
	    $orig = "$file.orig";
	    print "Processing $file...\n";
	    rename($file,$orig);
 	    open(OUT,"> $file");
 	    open(HEAD,"< $header");
	    print OUT "##############################################################################\n";
	    print OUT "#\n";
 	    foreach $line (<HEAD>) {
		chomp($line);
		print OUT "# $line\n";
 	    }
	    print OUT "#\n";
	    print OUT "##############################################################################\n";
 	    close(HEAD);
 	    open(IN,"< $orig");
	    $body = 0;
 	    foreach $line (<IN>) {
		chomp($line);
		if ($line =~ /^$/) {
		    $body = 1;
		}
		if ($body == 1) {
		    print OUT "$line\n";
		}
 	    }
 	    close(IN);
 	    close(OUT);
	}
	if ($file =~ /.[Ch]$/) {
	    $orig = "$file.orig";
	    print "Processing $file...\n";
	    rename($file,$orig);
 	    open(OUT,"> $file");
 	    open(HEAD,"< $header");
	    print OUT "/*****************************************************************************\n";
	    print OUT "*\n";
 	    foreach $line (<HEAD>) {
		chomp($line);
		print OUT "* $line\n";
 	    }
	    print OUT "*\n";
	    print OUT "*****************************************************************************/\n";
 	    close(HEAD);
 	    open(IN,"< $orig");
	    $body = 0;
 	    foreach $line (<IN>) {
		chomp($line);
		if ($line =~ /^$/) {
		    $body = 1;
		}
		if ($body == 1) {
		    print OUT "$line\n";
		}
 	    }
 	    close(IN);
 	    close(OUT);
	}
    }
}
