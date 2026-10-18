set -e
iqtree3 -s 18S_alignment_trimmed.fasta -m MFP -B 1000 -alrt 1000 -pre tree -redo > /dev/null
M=$(grep "Best-fit model according to BIC" tree.iqtree | sed 's/.*: //')
echo "model $M" > progress.log
iqtree3 -s 18S_alignment_trimmed.fasta -m "$M" -g constraint.nwk -B 1000 -pre constrained -redo > /dev/null
cat tree.treefile constrained.treefile > both.nwk
iqtree3 -s 18S_alignment_trimmed.fasta -m "$M" -z both.nwk -n 0 -zb 10000 -au -pre au_test -redo > /dev/null
echo ALL_DONE >> progress.log
