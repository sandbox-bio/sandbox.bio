<script>
import Execute from "$components/Execute.svelte";
import Link from "$components/Link.svelte";
</script>

Next we need to turn that output into a more standardized format, namely a CSV file that is supported by many tools for downstream analysis:

1. Spreadsheets. You can import it into Google Sheets or Excel.
2. R. Import it in R as a dataframe and do visualization with ggplot2 or other plotting libraries.
3. Python. Import it with the pandas library to get a dataframe and do visualization with matplotlib, seaborn, altair, and many other plotting libraries.
4. Circa! We can make a beautiful circos plot in Circa, like this: <Link href="https://circa.omgenomics.com/app/plot/gallery/aligned_genomes" />

First, use show-coords with the `-lTH` flags to get a simpler tab-delimited file that we can better work with programmatically:
<Execute command={"show-coords -lTH CBS415_to_Podan.delta > CBS415_to_Podan.coords"} />

The headers available from show-coords are not very readable (e.g. [E1], [S1]), so we'll add our own:

<Execute command={"NEW_HEADER='reference_start,reference_end,query_start,query_end,reference_alignment_length,query_alignment_length,percent_identity,reference_length,query_length,reference_chromosome,query_chromosome'"} />

Add that header and change the column spacing to comma-separated:

<Execute command={"(echo $NEW_HEADER && awk '{$1=$1}1' OFS="," CBS415_to_Podan.coords) > CBS415_to_Podan.csv"} />

Now we can load that into Circa to make this beautiful plot:

<Link href="https://circa.omgenomics.com/app/plot/gallery/aligned_genomes" />
