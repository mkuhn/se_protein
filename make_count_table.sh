src/cluster_pred.py	mouse	-	1           > /dev/null
src/cluster_pred.py	mouse	@inh	1           > /dev/null
echo "cutoff	side effects	protein action	all checked	strong support	indirect support	no explanation found	prediction not verified	text-mining error	other mechanism	total" > summary.tsv
src/cluster_pred.py	full_merged	-	1e-5   > /dev/null
src/cluster_pred.py	full_merged	@act	1e-5   > /dev/null
src/cluster_pred.py	full_merged	@inh	1e-5   > /dev/null
src/cluster_pred.py	mouse	-	1e-5           > /dev/null
src/cluster_pred.py	mouse	@act	1e-5           > /dev/null
src/cluster_pred.py	mouse	@inh	1e-5           > /dev/null
src/cluster_pred.py	full_merged	-	1e-2   | grep -v metabolizing | sort -nk 2 | ~/src/misc/pick_first.py -k 3 > merged_ann
src/cluster_pred.py	full_merged	@act	1e-2   > /dev/null
src/cluster_pred.py	full_merged 	@inh	1e-2   > /dev/null
src/cluster_pred.py	mouse	-	1e-2           > /dev/null
src/cluster_pred.py	mouse	@act	1e-2           > /dev/null
src/cluster_pred.py	mouse	@inh	1e-2           > /dev/null
