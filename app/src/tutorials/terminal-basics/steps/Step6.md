<script>
import Execute from "components/Execute.svelte";
</script>

Finally, let's explore environment variables. These are variables you can define in the terminal:

<Execute command="abc=123" />

Make sure there are no spaces surrounding the equal sign, otherwise the terminal treats the variable name as a command!

To display the content of a variable, use `echo`:

<Execute command="echo $abc" />

To delete a variable, use `unset`:

<Execute command="unset abc" />

To list all variables available in your environment:

<Execute command="env" />

Note that there's a variable called `USER` that we use while displaying the command-line prompt.

You can modify this variable to customize your environment:

<Execute command="USER=yourNameGoesHere" />

Now your prompt should update to reflect your name instead of `guest`!
