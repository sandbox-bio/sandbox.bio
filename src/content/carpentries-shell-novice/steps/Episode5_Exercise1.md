<script>
import Exercise from "$components/Exercise.svelte";
import Execute from "$components/Execute.svelte";

const criteria = [];
for(let i = 1; i <= 5; i++) {
    criteria.push({
        name: `File <code>${i}.txt</code> exists`,
        checks: [{
            type: "file",
            path: `loops/${i}.txt`,
            action: "exists"
        }]
    })
}
</script>

Go to the folder `loops`:

<Execute command="cd ~/tutorial/loops/" />

How would you write a loop that creates 5 files: `1.txt`, `2.txt`, `3.txt`, `4.txt`, `5.txt`? The files can be empty.

<Exercise {criteria} />
