<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

**Generate an rooted phylogenetic tree using LSD2**

In the previous step we created an unrooted phylogenetic tree, now, we will use [LSD2](https://github.com/tothuhien/lsd2) to create a rooted tree. LSD2 is a program that uses least-squares methods to estiamte rates and dates, and 

1. Try <Execute command="LSD2 -help" inline /> to
take a look at the usage instruction of LSD2.

We will want to use the following flags in our command:

- '-i' specifies the input file, which is an unrooted phylogenetic tree
- '-d' speficies the file with sequences dates, which is essential for rooting
- '-r'
- '-l'
- '-u'
- '-R'
- '-t'
- '-v'
- '-s'

2. Try <Execute command="lsd2 -i phylogenetic.tree -d hiv1_dates.txt -r a -l -1 -u 0 -q 0.2 -R 365 -t 0.00000000010000000000 -v 1 -s 9182" inline /> to generate our phylogenetic tree.



