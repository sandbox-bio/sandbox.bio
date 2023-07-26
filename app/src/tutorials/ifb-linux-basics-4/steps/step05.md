<script>
import Quiz from "components/Quiz.svelte";
import Execute from "components/Execute.svelte";
</script>

We have previously **cut** the fourth column from the SAOUHSC.bed file and **sorted** the subsequent stream. We could have written the two intructions on the same line using the ";" operator:

<Execute command="cut -f 4 SAOUHSC.bed > SAOUHSC_c4.bed; sort -u SAOUHSC_c4.bed > SAOUHSC_c4_uniq.bed" />

However this solution still requires to create an **intermediate file** (SAOUHSC_c4.bed) to perform both operations.  

This is where another extremely important redirection operator comes into play: **the "|"** pipe. This operator can be used to transmit the text stream from one command to another, avoiding the creation of intermediate files. By default, the pipe pass the **stdout of one command to the stdin** of the following one.

<img src="/data/linux_basics_session04/stream_pipe.png" style="max-width:100%" alt="pipe organisation">

So we may rewrite the previous set of instructions into the following which indicates that sort no more reads **stdin** from a file but from the result/stream of the cut command.

<Execute command="cut -f 4 SAOUHSC.bed | sort -u > SAOUHSC_c4_uniq.bed" />

In the same way we can also send the result of the **sort** command to the **wc** to get the expected result onto the screen whithout any need to create two intermediate files.

<Execute command="cut -f 4 SAOUHSC.bed | sort -u | wc -l" />

<Quiz id="question1" choices={[
	{ valid: false, value: "cut -f4  SAOUHSC.bed | sort -u | wc -l | grep 'gene'"},
		{ valid: true, value: "cut -f4  SAOUHSC.bed | grep 'gene' | sort -u | wc -l"},
	{ valid: false, value: "cut -f4  SAOUHSC.bed | grep 'gene' SAOUHSC.bed  | sort -u | wc -l"},
]}>
	<span slot="prompt">
		What would be the command to compute the number of non-redundant genes (ie. lines with fourth column starting with 'gene') in the SAOUHSC.bed file.
	</span>
</Quiz>

<Quiz id="question2" choices={[
	{ valid: false, value: "head -n 6 SAOUHSC.bed | tail -n 1 SAOUHSC.bed "},
	{ valid: false, value: "head -n 1 SAOUHSC.bed | tail -n 6 SAOUHSC.bed"},
	{ valid: true, value: "head -n 6 SAOUHSC.bed | tail -n 1"},
	{ valid: false, value: "head -n 1 | tail -n 6 SAOUHSC.bed"},
	{ valid: false, value: "head -n 6 | tail -n 1 SAOUHSC.bed"},
]}>
	<span slot="prompt">
		Let say we want to extract the 6th line of the SAOUHSC.bed file. What would be the correct syntax to do this ?
	</span>
</Quiz>
