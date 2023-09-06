<script>
import Execute from "$components/Execute.svelte";
import Exercise from "$components/Exercise.svelte";

const criteria = [
{
	name: "File <code>sucrose.dat</code> was moved to the <code>raw/</code> folder",
	checks: [{
		type: "file",
		path: "raw/sucrose.dat",
		action: "exists"
	}
    // , {
	// 	type: "file",
	// 	path: "analyzed/sucrose.dat",
	// 	action: "notexists"
	// }
    ]
},
{
	name: "File <code>maltose.dat</code> was moved to the <code>raw/</code> folder",
	checks: [{
		type: "file",
		path: "raw/maltose.dat",
		action: "exists"
	}
    // , {
	// 	type: "file",
	// 	path: "analyzed/maltose.dat",
	// 	action: "notexists"
	// }
    ]
}];
</script>

After running the following commands,
Jamie realizes that she put the files `sucrose.dat` and `maltose.dat` into the wrong folder.
The files should have been placed in the `raw` folder.

<Execute command="ls -F" />

<Execute command="ls -F analyzed" />

<Execute command="cd analyzed" />

Fill in the blanks to move these files to the `raw/` folder
(i.e. the one she forgot to put them in)

```bash
$ mv sucrose.dat maltose.dat ____/____
```

<Exercise {criteria} />
