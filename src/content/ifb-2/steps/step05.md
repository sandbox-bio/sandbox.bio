<script>
import Quiz from "$components/Quiz.svelte";
</script>

When used without argument, the `cd` command will set the current working directory to your HOME directory. 
This HOME directory is the place where a user may store their files. 

```bash
cd 
pwd
```

The symbol for the HOME directory is `~` (tilde). This character can be accessed under a PC keyboard using <kbd>AltGr</kbd> + <kbd>2</kbd>. With a Mac OSX keyboard, it may be accessed using <kbd>option</kbd> + <kbd>n</kbd>. 

In the example below we successively go to the `/tmp` then `/root/test` directories:

```bash
cd /tmp
pwd
cd ~/test
pwd
```

However, note that the HOME directory is not always the right place to store large files (particularly on a cluster with shared resources). 
Ask your administrator!


To answer the next question, please type the 3 following commands:

```bash
cd /shared/bank/nr
cd ~/test
cd
```

<Quiz id="q1" choices={[ { valid: false, value: "/shared/bank/nr"}, 
						 { valid: false, value: "test"}, 
						 { valid: false, value: "your HOME directory"}, 
						 { valid: false, value: "/shared/bank"}, 
						 { valid: true, value: "nr"},
						 { valid: false, value: "/root/test"},
						 { valid: false, value: "/shared/bank"}, ]}> 
	<span slot="prompt">
		...and select the right current working directory:
	</span>
</Quiz>
