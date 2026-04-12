<script>
import Execute from "$components/Execute.svelte";
</script>

Sometimes you want to extract specific reads by name — for example, pulling reads flagged by another tool, or extracting a known subset for closer inspection. fqgrep's `-N` flag does this efficiently.

## The -N flag

The `-N` (or `--read-names-file`) option reads a list of query names from a file and outputs only records whose read ID matches:

<Execute command={`cat read_names.txt`} />

<Execute command={`fqgrep -N read_names.txt reads.fastq`} />

Each line in the file is a read name. fqgrep matches the portion of the FASTQ header before the first whitespace.

## Counting name matches

Use `-c` to count how many reads matched:

<Execute command={`fqgrep -N read_names.txt -c reads.fastq`} />

## Inverting name matches

Combine `-N` with `-v` to get all reads **except** those in your list:

<Execute command={`fqgrep -N read_names.txt -v -c reads.fastq`} />

## Read name filtering with paired-end data

You can combine `-N` with `--paired` to extract specific read pairs. Note that the read names must match the names in the paired-end files — our `read_names.txt` contains names from the single-end file, so we'd need a separate list for paired data.

## When to use -N vs pattern matching

| Use case                              | Approach                                   |
| ------------------------------------- | ------------------------------------------ |
| Find reads containing a sequence      | Pattern matching (`fqgrep 'ACGT' file.fq`) |
| Extract reads by name/ID              | `-N` (`fqgrep -N names.txt file.fq`)       |
| Extract reads flagged by another tool | `-N` with the tool's output                |

> **Note**: The `-N` flag is mutually exclusive with pattern matching options (`-e`, `-f`, `-F`, and positional patterns). You're either filtering by sequence content or by read name — not both at once.
