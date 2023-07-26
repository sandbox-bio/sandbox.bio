<script>
import Quiz from "components/Quiz.svelte";
import Execute from "components/Execute.svelte";
</script>

As we've seen so far, the **sort** command performs, **by default, alphanumeric sorting** on a text stream. As you can see below, alphanumeric sorting is not well suited for stricly numeric values.

<Execute command="echo -e '1\n100\n2\n3\n200\n20\n10' | sort" />

To perform numeric sorting one need to activate the **-n/--numeric-sort** argument. This argument allows to properly perform **integer sorting**. 

<Execute command="echo -e '1\n100\n2\n3\n200\n20\n10' | sort -n " />

Keep in mind that one need to **read the sort manual carefully** when working with **floating numbers** (see **--general-numeric-sort** argument).

By default, sorting is performed based on all characters in the line.







In this case

Example on **O.tauri_annotation.gff** ?


sur gff : col9, ID gene, tri, unique
compter le nombre de gènes par chr
extraire les CDS du gff

<Execute command='cut -f 9 ~/dubii/study-cases/Escherichia_coli/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.chromosome.Chromosome.gff3 | cut -d ";" -f 1 | grep "ID=gene" | sort -u | wc -l' />
