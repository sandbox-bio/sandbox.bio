<script>
import Execute from "components/Execute.svelte";
</script>

With the `-d` (distance) option, one can also merge intervals that do not overlap, yet are close to one another. For example, to merge features that are no more than 1000bp apart, one would run:

<Execute command={"bedtools merge -i exons.bed -d 1000 -c 1 -o count | head -n 20"} />
