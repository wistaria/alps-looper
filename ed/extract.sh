#!/bin/sh

awk '$1=="T" { temp = $3 } 
	$1=="energy" && $2=="per"  && $3=="site" { ene = $5 }
	$1=="specific" && $2=="heat" { sh = $4 }
	$1=="uniform" && $2=="magnetization^2" { umag = $4 }
	$1=="uniform" && $2=="susceptibility" { usus = $4 }
	$1=="staggered" && $2=="magnetization^2" { smag = $4 }
	$1=="staggered" && $2=="susceptibility" { print temp, ene, sh, umag, usus, smag, $4 }'
