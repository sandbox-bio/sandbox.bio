<script>
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

Beyond trimming adapters and reads that are too short, `fastp` can also trim low quality bases at the beginning and ends of reads using a sliding window.

For example, `--cut_front` lets you remove low quality bases at the beginning of the read. It does so by sliding a window of 4bp and evaluating the mean base quality in that window. As long as the mean quality is under 20 inside the window, `fastp` will remove the bases, otherwise it stops and leaves the rest of the read intact.

Here it is in action:

<Execute command="fastp \ --in1 HG004_R1.fastq.gz \ --in2 HG004_R2.fastq.gz \ --cut_front" />

<Alert>The 4bp window size and quality threshold of Q20 are customizable using `--cut_front_window_size` and `--cut_front_mean_quality`</Alert>

Similarly, to remove low quality bases at the ends of reads, you can use the `--cut_right` flag. Like `--cut_front`, it will use a sliding window starting at the beginning of the read **but** as soon as it finds a window with low mean base quality, it will get rid of the rest of the read!

<Execute command="fastp \ --in1 HG004_R1.fastq.gz \ --in2 HG004_R2.fastq.gz \ --cut_front \ --cut_right" />

Note that, when using only `--cut_front`, the total number of bases is 3,657,586:

```
Read2 aftering filtering:
total reads: 24459
total bases: 3657586
```

whereas `--cut_front --cut_right` shows _fewer_ bases, as expected:

```
Read2 aftering filtering:
total reads: 24421
total bases: 3393780
```
