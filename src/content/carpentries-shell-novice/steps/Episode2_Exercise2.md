<script>
import Quiz from "$components/Quiz.svelte";
import Image from "$components/Image.svelte";
</script>

Using the filesystem diagram below, if `pwd` displays `/Users/thing`,
what will `ls -F ../backup` display?

<Image src="https://swcarpentry.github.io/shell-novice/fig/filesystem-challenge.svg" />

<Quiz id="q2.2" choices={[
{ valid: false, value: "../backup: No such file or directory"},
{ valid: false, value: "2012-12-01 2013-01-08 2013-01-27"},
{ valid: false, value: "2012-12-01/ 2013-01-08/ 2013-01-27/"},
{ valid: true, value: "original/ pnas_final/ pnas_sub/"},
]}>
<span slot="prompt">
</span>
</Quiz>
