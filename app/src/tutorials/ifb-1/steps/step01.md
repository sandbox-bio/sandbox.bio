<script>
import Quiz from "components/Quiz.svelte";
</script>

Unix is a family of operating systems derived from the original Unix created by the AT&T company. Currently, the most commonly used Unix system is Linux, a free and open source version of Unix. In practice both terms are used equally.


The command line is one way to pass orders to a computer by writing instructions in an interface called a **terminal**.
Another usual way to interact with a computer is by clicking on menus or buttons in a graphical interface. For some time, we will ask you to forget this way of working (no mouse allowed 😃).

Into the terminal another software runs. It is called the **Shell**.

The Shell is a software that handles Unix instructions. It is now commonly used by the bioinformatics community as it allows performing highly complex tasks by accessing either locally or distantly shared resources (e.g. files, databases, remote computing clusters...). These tasks can be coded in such a way that the multiple steps of a bioinformatics analysis can be saved into scripts. Those may then be shared with collaborators, traced and automated.

To access the Shell's powerful features, you need to learn a new way of thinking and a new language.

<Quiz id="q1" choices={[
	{ valid: true, value: "Yes"},
	{ valid: false, value: "No"},
]}>
	<span slot="prompt">
		Are you ready to learn Shell? (select "Yes" to continue)
	</span>
</Quiz>
