<script>
import Execute from "$components/Execute.svelte";
</script>

If we want to analyze all the lines we can redirect the output to the command line tool less, which allows us to navigate through the output by using the arrow keys. To do that we can add the bar character (`|`) between two commands, which will transfer the output of the first command as input of the second:

<Execute command="./getproteins.sh 27732 | less" />

To exit from less just press **q**.

However, what we really want is to save the output as a file, not just
printing some characters on the screen. Thus, what we should do is redirect
the output to a CSV file. This can be done by adding the redirect operator `>`
and the filename, as described previously:

<Execute command="./getproteins.sh 27732 > output.csv" />

To check if the file was really created and to analyze its contents, we can
use the less command:

<Execute command="less output.csv" />

We can also open the file in our spreadsheet application, such as LibreOffice Calc or Microsoft Excel.
