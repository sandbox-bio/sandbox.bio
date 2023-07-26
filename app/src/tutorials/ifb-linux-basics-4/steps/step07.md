<script>
import Quiz from "components/Quiz.svelte";
</script>

In this session you have discovered the 3 data streams of an unix command: the input stream (`stdin`), the output stream (`strout`) and the error stream (`stderr`):

<img src="/data/linux_basics_session04/stream_in_out_err.png" style="max-width:100%" alt="3 streams of a command">

You have also seen that these streams use the terminal display by default but that they can be redirected to intermediate files:

<img src="/data/linux_basics_session04/stream_in_outfile_errfile.png" style="max-width:100%" alt="files stream of a command">

And more importantly, you know that with pipes, you can "skip" intermediate files and build a "complex" command that combines the succession of several unit commands:

<img src="/data/linux_basics_session04/stream_pipe.png" style="max-width:100%" alt="a complex command">

You also compose and apply "complex" commands on an annotation file with a bed or gff format like:
*  compute the number of non-redundant genes of a bed file
*  extract the ith line of a bed file
