
SE=/home/mkuhn/data/drugs/meddra_open_adverse_effects.tsv.gz
ACTIONS=/srv/servers/stitch3.0/download/actions.v3.1.tsv.gz
WD:=$(shell pwd)

all: valid_side_effects valid_targets metabolizing_proteins.names se_names pred.full_merged.se_protein_known.all.txt summary.tsv stats.tsv explanation_table.0.01.qv10.txt

stats.tsv: interactions.tsv side_effects.tsv
	echo Number of drugs > $@
	cut -f 1 interactions.tsv | sort -u | wc -l >> $@

	echo Number of targets >> $@
	cut -f 2 interactions.tsv | sort -u | wc -l >> $@

	echo Number of side effects >> $@
	cut -f 2 side_effects.tsv | sort -u | wc -l >> $@

	echo Number of side effects--drug pairs >> $@
	wc -l side_effects.tsv >> $@

	echo Number of target--drug pairs >> $@
	wc -l interactions.tsv >> $@

	cat stats.tsv


explanation_table.0.01.qv10.txt: explanation_table.0.01.txt
	grep '^1\.0	' $^ | cut -f 2- > $@

explanation_table.0.01.txt: explanation_table.0.01.tsv
	awk -F '	' '{ printf "%s\t%s\t%s\t%s\t%s\n", $$1, $$2, $$3, $$5, $$4 }' < $^ > $@

explanation_table.0.01.tsv: protein_se_pv.tsv actions/main_targets.tsv
	src/explanation_table.py


summary.tsv: pred.full_merged.se_protein_known.all.txt
	src/make_count_table.sh

se_names: side_effects.tsv
	cut -f 2,4 $^ | sort -u > $@
	test -s $@

pred.full_merged.se_protein_known.all.txt: protein_se_pv.tsv mouse_se.tsv src/ann/ann
	src/known_protein_se_benchmark.py

mouse_se.tsv: /home/mkuhn/data/mgi/mouse_se.tsv valid_targets valid_side_effects
	grep -f valid_targets $< | grep -f valid_side_effects > $@

protein_se_pv.tsv: interactions.tsv side_effects.tsv
	echo "Compute side effects!"
	test protein_se_pv.tsv -nt interactions.tsv
	test protein_se_pv.tsv -nt side_effects.tsv

metabolizing_proteins.names: metabolizing_proteins
	cut -f 2 $^ > $@

metabolizing_proteins: protein_names
	egrep -i '^[^;]*(cytochrome p.?450|multidrug resistance|lutathione S-transferase|ATP.binding cassette|Thromboxane.*Synthase|	CYP[0-9])' < $^ > $@

protein_names:
	echo 'select protein_external_id, preferred_name, annotation from items.proteins where species_id=9606' | psql -p 8890 -U mering string_9_0 -t -F '	' -A | sed 's/9606.//' > $@
	test -s $@

valid_side_effects: side_effects.tsv
	cut -f 2 $^ | sort -u > $@
	test -s $@

side_effects.tsv: $(SE) valid_drugs
	zcat $< | cut -f 1-6 | grep '	PT$$' | cut -f 1-5 | src/filter_se.py | sort > $@
	test -s $@

valid_targets: interactions.tsv
	cut -f 2 $^ | sort -u > $@
	test -s $@

valid_drugs: interactions.tsv
	cut -f 1 $^ | sort -u > $@
	test -s $@


interactions.tsv: full_interactions.tsv hobohm_remove_list.txt
	src/filter_interactions.py > $@
	test -s $@

drugs_with_any_target: full_interactions.tsv
	cut -f 1 $^ | uniq > $@

full_interactions.tsv: actions.tsv.gz
	src/combine_actions.py > $@
	test -s $@

actions.tsv.gz: $(ACTIONS) drugs_with_se
	zcat $< | egrep '^CID1.{8}	9606\.' | fgrep -f drugs_with_se | gzip > $@

drugs_with_se: $(SE)
	zcat $^ | cut -f 1 | uniq | sed 's/^-/CID/' > $@


drugs.smiles: drugs_with_any_target ~/data/stitch/computations/stitch3.smiles
	~/src/misc/match_pattern.py -d ' ' -f 1 drugs_with_any_target < ~/data/stitch/computations/stitch3.smiles > $@
	test -s $@

drugs.fp: drugs.smiles
	cd /g/bork3/home/mkuhn/src/cdk/workspace/tanimoto && java -cp cdk.jar:. tanimoto/MainClass precompute $(WD)/$^ $(WD)/$@
	test -s $@

drug_tanimoto_similarity.tsv: drugs.fp
	cd /g/bork3/home/mkuhn/src/cdk/workspace/tanimoto && java -cp cdk.jar:. tanimoto/MainClass fullmatrix $(WD)/$^ 1 0 | awk '{printf "-%s\t-%s\t%s\n",$$1,$$2,$$3}' > $(WD)/$@
	test -s $@

drug_tanimoto_above_0.7: drug_tanimoto_similarity.tsv
	egrep '	(1\.|0.[7-9])[0-9]*$$' $^ > $@
	test -s $@

hobohm_remove_list.txt: drug_tanimoto_above_0.7
	~/src/misc/hobohm.sh < $^ > $@
	test -s $@


.DELETE_ON_ERROR:
