<script>
import Execute from "$components/Execute.svelte";
import Link from "$components/Link.svelte";
</script>

Next we need to turn that output into a more standardized format, namely a CSV file that is supported by many tools for downstream analysis:

- Spreadsheets. You can import it into Google Sheets or Excel.
- R. Import it in R as a dataframe and do visualization with ggplot2 or other plotting libraries.
- Python. Import it with the pandas library to get a dataframe and do visualization with any plotting library like Matplotlib, Seaborn, or Altair.
- Circa! We can make a beautiful plot of all the alignments between sequences, like this <Link href="https://circa.omgenomics.com/app/plot/gallery/aligned_genomes">circos plot</Link> showing homology between two larger genomes.

Now we'll convert the delta file to coordinates and prepare it for visualization.

First, use show-coords to get the alignment coordinates:

<Execute command="show-coords -lTH $NAME.delta > $NAME.coords" />

The headers from show-coords aren't very readable (e.g. [E1], [S1]), so let's create our own header:

<Execute command={`NEW_HEADER="reference_start,reference_end,query_start,query_end,reference_alignment_length,query_alignment_length,percent_identity,reference_length,query_length,reference_chromosome,query_chromosome"`} />

Now we'll create a CSV file with our new header:

<Execute command={'(echo "$NEW_HEADER" && awk \'{$1=$1}1\' OFS="," $NAME.coords) > $NAME.csv'} />

Let's look at the first few lines of the CSV file:

<Execute command="head $NAME.csv" />

It might be easier to see it with <Execute command="less -S $NAME.csv" inline />, just press `q` to quit.
